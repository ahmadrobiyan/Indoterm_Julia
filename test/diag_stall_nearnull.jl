"""
C3c — name the near-singular mode that wrecks the Newton step.

`diag_stall_direct.jl` settled the rank question: at the stalled point `lu(J_s)` succeeds with a
linear residual of 3.5e-10, so `J` has full rank and `F ∈ range(J)`. The system is CONSISTENT —
PLAN.md §B's inconsistency verdict is refuted, and CGNR's 0.161 leftover was non-convergence
(its rate is governed by κ(J)², and 4000 iterations were not enough), not a left null space.

What actually fails: the exact Newton step has ‖d‖∞ = 1221 for ‖b‖₂ = 0.0033 — an amplification
of ~4e5. The direction is dominated by a near-singular mode, so any α long enough to reduce `F`
linearly also pushes the nonlinear terms far outside the region where the linearization holds,
and the line search finds no descent at any α.

This script names that mode. Inverse iteration on `J_sᵀJ_s` (using the LU factor for both the
solve and the adjoint solve) converges to the smallest singular triplet:

    σ_min      — how near-singular
    v_min      — the VARIABLES that are nearly undetermined
    u_min      — the EQUATIONS that are nearly linearly dependent  (u = J·v / σ)

`v_min` and `u_min` together identify the defective block: a missing normalization, a closure
that leaves a direction unpinned, or a pair of equations that collapse onto each other away from
the benchmark. `‖d‖` is also decomposed onto the leading modes to confirm the mode found is the
one actually contaminating the step.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_stall_nearnull.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra, Random
const MOI = JuMP.MOI
const IJ = IndotermJulia

canon(s::AbstractString) = replace(s, r"\[[0-9, ]+\]" => "[.]")
famof(s::AbstractString) = (mm = match(r"([A-Za-z_][A-Za-z0-9_]*)\[", s);
                            mm === nothing ? "(no-index)" : String(mm.captures[1]))
varfam(s::AbstractString) = (mm = match(r"^([A-Za-z_][A-Za-z0-9_]*)", s);
                             mm === nothing ? s : String(mm.captures[1]))

SHOCK_VAR = get(ENV, "STALL_SHOCK_VAR", "blabnat")
SHOCK_VAL = parse(Float64, get(ENV, "STALL_SHOCK_VAL", "0.9997"))

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("benchmark: status=$(r0.status)  ‖F‖∞=$(r0.residual)")
sv = vars[SHOCK_VAR]
JuMP.fix(sv isa AbstractArray ? first(sv) : sv, SHOCK_VAL; force=true)
r1 = solve_newton!(m, vars; maxit=12, verbose=false)
println("shock $SHOCK_VAR=$SHOCK_VAL: status=$(r1.status) iters=$(r1.iters) ‖F‖∞=$(r1.residual)")

nowfree = VariableRef[]
for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr) || push!(nowfree, vr)
end
N = length(nowfree)
varname = [JuMP.name(v) for v in nowfree]
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
evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
MOI.initialize(evaluator, [:Grad, :Jac])

x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
     (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]
freecols = [pos[JuMP.index(v).value] for v in nowfree]

st = MOI.jacobian_structure(evaluator)
jrows_all = getindex.(st, 1); jcols_all = getindex.(st, 2)
Jval = zeros(length(jrows_all))
col2free = zeros(Int, length(allv)); for (i, c) in enumerate(freecols); col2free[c] = i; end
keep = [col2free[c] != 0 for c in jcols_all]
jr = jrows_all[keep]; jc = [col2free[c] for c in jcols_all[keep]]

MOI.eval_constraint_jacobian(evaluator, Jval, x); vals = Jval[keep]
rowscale = ones(NC)
let mx = zeros(NC)
    for k in eachindex(jr); a = abs(vals[k]); a > mx[jr[k]] && (mx[jr[k]] = a); end
    for i in 1:NC; rowscale[i] = mx[i] > 1e-12 ? mx[i] : 1.0; end
end
g = zeros(NC); MOI.eval_constraint(evaluator, g, x); g .-= rhs
gs = g ./ rowscale

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
println("\nJ_s $(size(J_s))  nnz=$(nnz(J_s))  ‖b‖₂=$(norm(b))")

tlu = @elapsed F = lu(J_s)
d = F \ b
println("lu $(round(tlu,digits=1))s  ‖J·d−b‖/‖b‖=$(norm(J_s*d .- b)/norm(b))  ‖d‖∞=$(maximum(abs,d))  ‖d‖₂=$(norm(d))")

# ── inverse iteration for the smallest singular triplet ─────────────────────
# v ← J⁻¹(J⁻ᵀ v) is one step of power iteration on (JᵀJ)⁻¹, whose dominant
# eigenvector is the smallest right singular vector of J. Converges linearly at
# rate (σ_min/σ_2)², which is fast precisely when the mode is well separated.
println("\n── inverse iteration on JᵀJ (smallest singular triplet) ──")
function inverse_iteration(Fac, A, n; iters=30)
    Random.seed!(20260729)
    vv = normalize!(randn(n))
    for it in 1:iters
        vv = Fac \ (Fac' \ vv)
        nv = norm(vv)
        vv ./= nv
        if it % 5 == 0 || it == iters
            println("  it $it  σ_min≈$(1/sqrt(nv))   ‖J·v‖=$(norm(A * vv))")
        end
    end
    return vv
end
v = inverse_iteration(F, J_s, N)
u = J_s * v
σ_true = norm(u)
u ./= max(σ_true, eps())
println("  converged: σ_min = $σ_true    κ(J_s) ≳ $(round(maximum(abs, J_s)*sqrt(N)/max(σ_true,eps()), sigdigits=3))")
println("  component of the Newton step along this mode: $(dot(d, v)) of ‖d‖₂=$(norm(d))" *
        "  ⇒ $(round(100*abs(dot(d,v))/norm(d), digits=1))% of the step")

# ── report ─────────────────────────────────────────────────────────────────
function topvars(title, w, k=30)
    println("\n── $title ──")
    byfam = Dict{String,Tuple{Int,Float64,Float64,String}}()
    for i in 1:N
        a = abs(w[i]); a < 1e-14 && continue
        f = varfam(varname[i]); n, ss, mx, ex = get(byfam, f, (0, 0.0, 0.0, varname[i]))
        byfam[f] = (n+1, ss + a^2, max(mx, a), a > mx ? varname[i] : ex)
    end
    tot = sum(t -> t[2], values(byfam); init=0.0)
    println("  by variable family (share of squared magnitude):")
    for (f, (n, ss, mx, ex)) in first(sort(collect(byfam); by = kv -> -kv[2][2]), 12)
        println("    $(rpad(f,16)) share=$(rpad(round(100*ss/max(tot,eps()),digits=1),6))%  n=$(rpad(n,7)) max=$(round(mx,sigdigits=4))  e.g. $ex")
    end
    idx = sortperm(abs.(w); rev=true)[1:min(k, N)]
    println("  top $(length(idx)) individual variables:")
    for i in idx
        println("    $(rpad(round(w[i], sigdigits=5), 14)) $(varname[i])")
    end
end

function toprows(title, w, k=30)
    println("\n── $title ──")
    byfam = Dict{String,Tuple{Int,Float64,Float64,String}}()
    for i in 1:NC
        a = abs(w[i]); a < 1e-14 && continue
        f = famof(rowsig[i]); n, ss, mx, ex = get(byfam, f, (0, 0.0, 0.0, rowsig[i]))
        byfam[f] = (n+1, ss + a^2, max(mx, a), a > mx ? rowsig[i] : ex)
    end
    tot = sum(t -> t[2], values(byfam); init=0.0)
    println("  by equation family (share of squared magnitude):")
    for (f, (n, ss, mx, ex)) in first(sort(collect(byfam); by = kv -> -kv[2][2]), 12)
        println("    $(rpad(f,16)) share=$(rpad(round(100*ss/max(tot,eps()),digits=1),6))%  n=$(rpad(n,7)) max=$(round(mx,sigdigits=4))")
        println("        $(first(ex,150))")
    end
    idx = sortperm(abs.(w); rev=true)[1:min(k, NC)]
    println("  top $(length(idx)) individual rows:")
    for i in idx
        println("    $(rpad(round(w[i], sigdigits=5), 14)) $(first(rowsig[i],130))")
    end
end

println("  ⇒ σ_min = $σ_true. A mode this small means the Newton step's component along v_min is")
println("    amplified by 1/σ_min; that is the term producing ‖d‖∞ = $(round(maximum(abs,d),sigdigits=5)).")
topvars("v_min — the nearly-UNDETERMINED variable direction (σ_min = $σ_true)", v)
toprows("u_min = J·v/σ — the nearly-DEPENDENT equations", u)
topvars("d — the Newton step itself, in Ruiz-scaled units", d)
topvars("d./Dc — the Newton step in ORIGINAL variable units", d ./ Dc)

println("\nDone.")
