"""
C1+C2+C3 — why every non-zero shock stalls at a residual floor.

See PLAN.md "Session findings (2026-07-29)". The verdict under test: the Newton solver is not
failing to converge, it is converging to the least-squares solution of an **inconsistent**
system — part of `F` lies in the left null space of `J`, so no step of any length removes it.

The decisive quantity is CGNR's own leftover. With μ = 0, CGNR minimizes `‖J·d + F‖`; its final
residual vector `r = -F - J·d` is (to convergence) exactly the projection of `F` onto
`range(J)^⊥` — i.e. the unreachable part. So:

  * `‖r‖ / ‖F‖ ≈ 0`  → J has full rank at this point, the obstruction is nonlinearity/line search.
  * `‖r‖ / ‖F‖ ≫ 0`  → J is rank-deficient here, and the large entries of `r` **name the
    equations that block the solve**. That is C2 and C3 answered by one computation, and it is
    far cheaper than a sparse-QR rank count on an 83k×83k matrix.

This also fixes the C1 diagnostic bug: residuals come from the evaluator at the solver's own `x`,
never from `JuMP.value(con)` (which reads a stale optimizer cache the Newton solver never fills —
that is why earlier reports claimed "worst |resid| = 0.0").

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_stall_nullspace.jl
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra
const MOI = JuMP.MOI
const IJ = IndotermJulia

canon(s::AbstractString) = replace(s, r"\[[0-9, ]+\]" => "[.]")
famof(s::AbstractString) = begin
    mm = match(r"([A-Za-z_][A-Za-z0-9_]*)\[", s)
    mm === nothing ? "(no-index)" : String(mm.captures[1])
end

SHOCK_VAR = get(ENV, "STALL_SHOCK_VAR", "blabnat")
SHOCK_VAL = parse(Float64, get(ENV, "STALL_SHOCK_VAL", "0.9997"))

# ── 1. model ────────────────────────────────────────────────────────────────
agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s  vars=$(num_variables(m)) cons=$(num_constraints(m; count_variable_in_set_constraints=false))")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

println("\n── benchmark solve ──")
r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("  status=$(r0.status) scaled ‖F‖∞=$(r0.residual)  ($(r0.n_eqs) eqs, $(r0.n_free) free)")
r0.residual < 1e-6 || @warn "benchmark did not converge — the stall diagnosis below is unreliable"

# ── 2. apply the shock and let the solver stall ─────────────────────────────
println("\n── shock: $SHOCK_VAR = $SHOCK_VAL ──")
sv = vars[SHOCK_VAR]
JuMP.fix(sv isa AbstractArray ? first(sv) : sv, SHOCK_VAL; force=true)
r1 = solve_newton!(m, vars; maxit=12, verbose=false)
println("  status=$(r1.status)  iters=$(r1.iters)  stalled scaled ‖F‖∞=$(r1.residual)")
# solve_newton! writes its final x back into the start values, so the model now *is* the
# stalled point and everything below reads it from there.

# ── 3. rebuild the evaluator, keeping a row → equation-signature map ────────
println("\n── rebuilding evaluator at the stalled point ──")
nowfree = VariableRef[]
for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr) || push!(nowfree, vr)
end
N = length(nowfree)

allv = JuMP.all_variables(m)
allidx = [JuMP.index(v) for v in allv]
pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))

nlmodel = MOI.Nonlinear.Model()
rhs = Float64[]; rowsig = String[]
for (Ftype, S) in list_of_constraint_types(m)
    Ftype <: VariableRef && continue
    for con in all_constraints(m, Ftype, S)
        co = JuMP.constraint_object(con)
        r = IJ._rhs_of(co.set)
        MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r))
        push!(rhs, r); push!(rowsig, canon(string(co.func)))
    end
end
NC = length(rhs)
println("  $NC equations, $N unknowns" * (NC == N ? " (square)" : " ⚠ NOT SQUARE"))

evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
MOI.initialize(evaluator, [:Grad, :Jac])   # NOTE: :Constraints is NOT a valid ReverseAD feature

x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
     (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]
freecols = [pos[JuMP.index(v).value] for v in nowfree]

st = MOI.jacobian_structure(evaluator)
jrows_all = getindex.(st, 1); jcols_all = getindex.(st, 2)
Jval = zeros(length(jrows_all))
col2free = zeros(Int, length(allv))
for (i, c) in enumerate(freecols); col2free[c] = i; end
keep = [col2free[c] != 0 for c in jcols_all]
jr = jrows_all[keep]; jc = [col2free[c] for c in jcols_all[keep]]

# ── 4. equilibrate exactly as the solver does ───────────────────────────────
MOI.eval_constraint_jacobian(evaluator, Jval, x)
vals = Jval[keep]
rowscale = ones(NC)
let mx = zeros(NC)
    for k in eachindex(jr); a = abs(vals[k]); a > mx[jr[k]] && (mx[jr[k]] = a); end
    for i in 1:NC; rowscale[i] = mx[i] > 1e-12 ? mx[i] : 1.0; end
end
g = zeros(NC); MOI.eval_constraint(evaluator, g, x); g .-= rhs
gs = g ./ rowscale
println("  scaled ‖F‖∞ at stalled point = $(norm(gs, Inf))")

Jbase = sparse(jr, jc, [vals[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, N)

# ruiz_scales is a closure inside solve_newton!; reimplement the same 20-sweep loop here.
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
Dr, Dc = ruiz(Jbase)
J_s = sparse(jr, jc, [vals[k] / (rowscale[jr[k]] * Dr[jr[k]] * Dc[jc[k]]) for k in eachindex(jr)], NC, N)
rhs_s = gs ./ Dr

# ── 5. THE TEST: unregularized least-squares leftover ───────────────────────
println("\n── CGNR with μ=0, maxit=4000 (is F in range(J)?) ──")
d = zeros(N); work_r = zeros(NC)
tcg = @elapsed IJ.cgnr!(d, J_s, -rhs_s; mu=0.0, maxit=4000, tol=1e-14, work_r=work_r)
lsr = -rhs_s .- J_s * d                       # exact leftover, recomputed (not CGNR's running r)
relres = norm(lsr) / norm(rhs_s)
println("  $(round(tcg, digits=1))s   ‖b‖₂=$(norm(rhs_s))   ‖b − J·d‖₂=$(norm(lsr))")
println("  RELATIVE LEAST-SQUARES RESIDUAL = $relres")
println(relres > 1e-3 ?
    "  ⇒ F is NOT in range(J): the system is INCONSISTENT here. J is rank-deficient at the\n     shocked point. The rows below name the blocking equations." :
    "  ⇒ F IS in range(J): J has full rank here. The stall is nonlinearity/line-search, NOT rank.\n     Verdict in PLAN.md §B is refuted — go to plan step C5 (condensation).")

# ── 6. name the blocking equations ─────────────────────────────────────────
function report(title, w, k=25)
    println("\n── $title ──")
    byfam = Dict{String,Tuple{Int,Float64,Float64,String}}()  # fam => (n, sumsq, max, example sig)
    for i in 1:NC
        a = abs(w[i]); a < 1e-14 && continue
        f = famof(rowsig[i])
        n, ss, mx, ex = get(byfam, f, (0, 0.0, 0.0, rowsig[i]))
        byfam[f] = (n+1, ss + a^2, max(mx, a), a > mx ? rowsig[i] : ex)
    end
    tot = sum(t -> t[2], values(byfam); init=0.0)
    println("  by equation family (share of squared magnitude):")
    for (f, (n, ss, mx, ex)) in first(sort(collect(byfam); by = kv -> -kv[2][2]), 12)
        println("    $(rpad(f, 16)) share=$(rpad(round(100*ss/max(tot,eps()), digits=1), 6))%  n=$(rpad(n,7)) max=$(round(mx, sigdigits=4))")
        println("        $(first(ex, 150))")
    end
    idx = sortperm(abs.(w); rev=true)[1:min(k, NC)]
    println("  top $(length(idx)) individual rows:")
    for i in idx
        println("    |w|=$(rpad(round(abs(w[i]), sigdigits=4), 12)) $(first(rowsig[i], 140))")
    end
end

report("C3: unreachable component  r = b − J·d  (the left-null projection of F)", lsr)
report("C2: raw scaled residual  F  at the stalled point", gs)

println("\nDone.")
