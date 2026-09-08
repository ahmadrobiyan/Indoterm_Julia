"""
C3b — direct sparse LU at the stalled point. The control for `diag_stall_nullspace.jl`.

A large CGNR least-squares residual has TWO possible causes and they demand opposite fixes:

  (a) `J` is genuinely rank-deficient → `F ∉ range(J)`, the system is inconsistent, and no linear
      solver helps. Fix the equations.
  (b) `J` is full rank but so ill-conditioned that CGNR (a Krylov method whose convergence is
      governed by κ(J)²  on the normal equations) cannot converge in the iteration budget. Then
      the "unreachable" component is an artifact of the *iterative solver*, and a **direct**
      factorization solves it exactly.

`diag_stall_nullspace.jl` alone cannot tell these apart. A sparse LU can: at 25×6 the Ruiz-scaled
Jacobian is 83,541² with ~430k nonzeros and factorizes in ~3 s / 2.5 GB (PLAN.md's region-scaling
table), so this is cheap. Note the solver moved `lu`/`qr` → CGNR on 2026-07-28 *before* Ruiz
two-sided equilibration was added; the original reason for abandoning LU (step explosion from
near-singular columns) may no longer hold.

Reports, at the stalled point:
  * whether `lu(J_s)` succeeds, and the true linear-solve residual ‖J_s·d + b‖/‖b‖;
  * the same for CGNR at μ=0, for direct comparison;
  * whether the LU Newton step actually reduces ‖F‖∞ under a line search — i.e. whether swapping
    the linear solver back to LU is the fix.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_stall_direct.jl
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra
const MOI = JuMP.MOI
const IJ = IndotermJulia

canon(s::AbstractString) = replace(s, r"\[[0-9, ]+\]" => "[.]")
famof(s::AbstractString) = (mm = match(r"([A-Za-z_][A-Za-z0-9_]*)\[", s);
                            mm === nothing ? "(no-index)" : String(mm.captures[1]))

SHOCK_VAR = get(ENV, "STALL_SHOCK_VAR", "blabnat")
SHOCK_VAL = parse(Float64, get(ENV, "STALL_SHOCK_VAL", "0.9997"))

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

println("\n── benchmark solve ──")
r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("  status=$(r0.status)  scaled ‖F‖∞=$(r0.residual)")

println("\n── shock $SHOCK_VAR = $SHOCK_VAL, let it stall ──")
sv = vars[SHOCK_VAR]
JuMP.fix(sv isa AbstractArray ? first(sv) : sv, SHOCK_VAL; force=true)
r1 = solve_newton!(m, vars; maxit=12, verbose=false)
println("  status=$(r1.status)  iters=$(r1.iters)  stalled ‖F‖∞=$(r1.residual)")

# ── rebuild evaluator (same order as solve_newton!) ─────────────────────────
nowfree = VariableRef[]
for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr) || push!(nowfree, vr)
end
N = length(nowfree)
allv = JuMP.all_variables(m); allidx = [JuMP.index(v) for v in allv]
pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))

nlmodel = MOI.Nonlinear.Model(); rhs = Float64[]; rowsig = String[]
for (Ftype, S) in list_of_constraint_types(m)
    Ftype <: VariableRef && continue
    for con in all_constraints(m, Ftype, S)
        co = JuMP.constraint_object(con); r = IJ._rhs_of(co.set)
        MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r))
        push!(rhs, r); push!(rowsig, canon(string(co.func)))
    end
end
NC = length(rhs)
println("\n  $NC equations, $N unknowns" * (NC == N ? " (square)" : " ⚠ NOT SQUARE"))
evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
MOI.initialize(evaluator, [:Grad, :Jac])

x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
     (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]
freecols = [pos[JuMP.index(v).value] for v in nowfree]
lb = [JuMP.has_lower_bound(v) ? JuMP.lower_bound(v) : -Inf for v in nowfree]

st = MOI.jacobian_structure(evaluator)
jrows_all = getindex.(st, 1); jcols_all = getindex.(st, 2)
Jval = zeros(length(jrows_all))
col2free = zeros(Int, length(allv))
for (i, c) in enumerate(freecols); col2free[c] = i; end
keep = [col2free[c] != 0 for c in jcols_all]
jr = jrows_all[keep]; jc = [col2free[c] for c in jcols_all[keep]]

MOI.eval_constraint_jacobian(evaluator, Jval, x); vals = Jval[keep]
rowscale = ones(NC)
let mx = zeros(NC)
    for k in eachindex(jr); a = abs(vals[k]); a > mx[jr[k]] && (mx[jr[k]] = a); end
    for i in 1:NC; rowscale[i] = mx[i] > 1e-12 ? mx[i] : 1.0; end
end
gbuf = zeros(NC)
resid!(g, xv) = (MOI.eval_constraint(evaluator, g, xv); g .-= rhs; g)
sres(xv, buf) = (resid!(buf, xv); buf ./ rowscale)
gs = sres(x, gbuf)
println("  scaled ‖F‖∞ at stalled point = $(norm(gs, Inf))")

function ruiz(J0::SparseMatrixCSC{Float64,Int}; maxiter=20, tol=1e-6)
    mm, nn = size(J0); Dr = ones(mm); Dc = ones(nn); Jw = copy(J0)
    for _ in 1:maxiter
        rowmax = zeros(mm)
        for k in 1:nnz(Jw); i = Jw.rowval[k]; a = abs(Jw.nzval[k]); a > rowmax[i] && (rowmax[i] = a); end
        for i in 1:mm; Dr[i] *= rowmax[i] > 1e-12 ? sqrt(rowmax[i]) : 1.0; end
        for k in 1:nnz(Jw); i = Jw.rowval[k]; Jw.nzval[k] /= (rowmax[i] > 1e-12 ? sqrt(rowmax[i]) : 1.0); end
        colmax = zeros(nn)
        for j in 1:nn, k in Jw.colptr[j]:(Jw.colptr[j+1]-1)
            a = abs(Jw.nzval[k]); a > colmax[j] && (colmax[j] = a)
        end
        for j in 1:nn; Dc[j] *= colmax[j] > 1e-12 ? sqrt(colmax[j]) : 1.0; end
        for j in 1:nn
            s = colmax[j] > 1e-12 ? sqrt(colmax[j]) : 1.0
            for k in Jw.colptr[j]:(Jw.colptr[j+1]-1); Jw.nzval[k] /= s; end
        end
        dev = 0.0
        for i in 1:mm; s = rowmax[i] > 1e-12 ? sqrt(rowmax[i]) : 1.0; dev = max(dev, abs(s-1)); end
        for j in 1:nn; s = colmax[j] > 1e-12 ? sqrt(colmax[j]) : 1.0; dev = max(dev, abs(s-1)); end
        dev < tol && break
    end
    return Dr, Dc
end

Jbase = sparse(jr, jc, [vals[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, N)
Dr, Dc = ruiz(Jbase)
J_s = sparse(jr, jc, [vals[k] / (rowscale[jr[k]] * Dr[jr[k]] * Dc[jc[k]]) for k in eachindex(jr)], NC, N)
b = -(gs ./ Dr)
println("  J_s: $(size(J_s))  nnz=$(nnz(J_s))  ‖b‖₂=$(norm(b))")

colnrm = zeros(N)
for k in eachindex(jr); colnrm[jc[k]] += (vals[k]/(rowscale[jr[k]]*Dr[jr[k]]*Dc[jc[k]]))^2; end
colnrm .= sqrt.(colnrm)
println("  columns with norm < 1e-12: $(count(colnrm .< 1e-12))   < 1e-3: $(count(colnrm .< 1e-3))")

# ── 1. DIRECT sparse LU ─────────────────────────────────────────────────────
println("\n── direct sparse LU ──")
d_lu = nothing
try
    tlu = @elapsed F = lu(J_s)
    println("  lu() OK in $(round(tlu, digits=1))s   nnz(L)+nnz(U)=$(nnz(F.L) + nnz(F.U))")
    global d_lu = F \ b
    rr = J_s * d_lu .- b
    println("  ‖J_s·d − b‖₂ / ‖b‖₂ = $(norm(rr)/norm(b))     max|d|=$(maximum(abs, d_lu))")
catch e
    println("  lu() FAILED: $(sprint(showerror, e))")
    println("  → trying sparse QR least-squares")
    try
        global d_lu = qr(J_s) \ b
        rr = J_s * d_lu .- b
        println("  qr:  ‖J_s·d − b‖₂ / ‖b‖₂ = $(norm(rr)/norm(b))   max|d|=$(maximum(abs, d_lu))")
    catch e2
        println("  qr() FAILED too: $(sprint(showerror, e2))")
    end
end

# ── 2. CGNR at μ=0, same right-hand side, for comparison ────────────────────
println("\n── CGNR μ=0, maxit=4000 (same b) ──")
d_cg = zeros(N)
tcg = @elapsed IJ.cgnr!(d_cg, J_s, b; mu=0.0, maxit=4000, tol=1e-14)
println("  $(round(tcg, digits=1))s  ‖J_s·d − b‖₂ / ‖b‖₂ = $(norm(J_s*d_cg .- b)/norm(b))  max|d|=$(maximum(abs, d_cg))")

# ── 3. does the LU step actually reduce ‖F‖∞? ───────────────────────────────
if d_lu !== nothing
    println("\n── line search on the LU Newton step (uncapped, then capped) ──")
    for (label, capped) in (("uncapped", false), ("10% cap + fraction-to-boundary", true))
        dd = d_lu ./ Dc   # Dc is a DIVISOR (J_s = Dr⁻¹·Jbase·Dc⁻¹) ⇒ original units = d./Dc
        if capped
            for i in 1:N
                absx = abs(x[freecols[i]])
                maxs = (isfinite(lb[i]) && lb[i] >= 1e-6 - 1e-9) ? 0.10*max(absx,1e-3) : max(0.10*absx, 10.0)
                dd[i] = clamp(dd[i], -maxs, maxs)
            end
        end
        a = 1.0
        for i in 1:N
            if isfinite(lb[i]) && dd[i] < -1e-10
                gap = x[freecols[i]] - lb[i]
                if gap > 1e-8
                    cap = 0.99*gap/(-dd[i]); (cap > 0 && cap < a) && (a = cap)
                else; dd[i] = 0.0; end
            end
        end
        f0 = norm(gs, Inf); xt = copy(x); tb = zeros(NC)
        best_a = 0.0; best_g = f0
        println("  [$label] a0=$(round(a, sigdigits=3)) max|d|=$(round(maximum(abs, dd), sigdigits=4))")
        for _ in 1:24
            for i in 1:N; xt[freecols[i]] = x[freecols[i]] + a*dd[i]; end
            gt = sres(xt, tb); gi = norm(gt, Inf)
            println("     a=$(rpad(round(a, sigdigits=3), 10)) ‖F‖∞=$(round(gi, sigdigits=6))")
            isfinite(gi) && gi < best_g && (best_g = gi; best_a = a)
            a *= 0.5; a < 1e-8 && break
        end
        println("  [$label] best a=$best_a   ‖F‖∞ $(round(f0, sigdigits=6)) → $(round(best_g, sigdigits=6))" *
                (best_g < f0*0.5 ? "   ✅ LU STEP WORKS" : best_g < f0 ? "   ~ partial" : "   ❌ no descent"))
    end
end

# ── 4. if the residual really is unreachable, name the rows ─────────────────
if d_lu !== nothing
    lsr = b .- J_s * d_lu
    if norm(lsr)/norm(b) > 1e-6
        println("\n── unreachable component r = b − J·d, by equation family ──")
        byfam = Dict{String,Tuple{Int,Float64,Float64,String}}()
        for i in 1:NC
            a = abs(lsr[i]); a < 1e-14 && continue
            f = famof(rowsig[i]); n, ss, mx, ex = get(byfam, f, (0,0.0,0.0,rowsig[i]))
            byfam[f] = (n+1, ss+a^2, max(mx,a), a > mx ? rowsig[i] : ex)
        end
        tot = sum(t -> t[2], values(byfam); init=0.0)
        for (f,(n,ss,mx,ex)) in first(sort(collect(byfam); by = kv -> -kv[2][2]), 12)
            println("    $(rpad(f,16)) share=$(rpad(round(100*ss/max(tot,eps()), digits=1),6))% n=$(rpad(n,7)) max=$(round(mx,sigdigits=4))")
            println("        $(first(ex,150))")
        end
    end
end
println("\nDone.")
