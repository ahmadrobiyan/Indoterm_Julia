"""
Is the obstruction at t ≈ 0.3063 a FOLD, or is it the solver?

`test/run_scenario_6reg.jl` walked `blabnat` 1.0 → 0.97 adaptively and drove the
step size to `hmin` at one specific point, `t = 0.30625` (`blabnat = 0.9908125`),
after five clean steps. That is the textbook signature of a limit point — but it
is ALSO the signature of the trust-region ratchet diagnosed in PLAN.md §N.10 and
never fixed, which plateaus at a fixed *relative* accuracy no matter how small
the step. The failed-step residuals are consistent with both readings: they fall
in proportion to `h` (≈3e-3·h throughout), which says Newton reached the same
relative accuracy every time and stopped — it does not say WHY.

The two readings make opposite predictions about one measurable quantity, so
measure it. At a converged point the tangent `v = dx/dλ` solves

    J(x, λ) · v = −∂F/∂λ

where `J` is the Jacobian in the free variables and `∂F/∂λ` is the column of the
full Jacobian belonging to the (fixed) shock variable.

  * **Fold.** `J` becomes singular AT the limit point, so `‖v‖ → ∞` like
    `1/(λ* − λ)`. Walking in, the tangent norm must blow up — and `1/‖v‖` must
    fall towards zero roughly linearly in λ, which extrapolates to the fold
    location. Nothing about the solver can prevent this; it is a property of the
    equations.
  * **Solver.** `J` stays well conditioned, `‖v‖` stays bounded and boring, and
    the branch simply continues past 0.9908 — the failure is in how the step is
    accepted, not in where the branch goes.

`‖v‖` is computed by a DIRECT sparse LU, not an iterative method. PLAN.md
§N.10 records the standing lesson: never infer rank or singularity from a CGNR
leftover — at this conditioning CGNR returned a step of magnitude 0.025 against
a true 1221.

Rows are equilibrated before the factorisation because this model's rows span
1e-10 … 1e12. Row scaling leaves `v` mathematically unchanged (it scales both
sides), so it costs nothing and buys a usable LU.

Reported per point:
  * `‖v‖∞` and the variable carrying it — the DIRECTION of the singular mode,
    which is what to fix if this is a fold. Note the probe in
    `run_scenario_6reg.jl` was `pcap[2,5]`, a leftover from the refuted "cell
    [2,5]" framing; its slope DECELERATES into the wall, so it is demonstrably
    not the fold direction. Do not assume — read the argmax.
  * an estimate of `σ_min(J)` by inverse power iteration reusing the same LU.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_fold_tangent.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, LinearAlgebra, SparseArrays
const MOI = JuMP.MOI   # MathOptInterface is not a direct project dependency

const TS = let s = get(ENV, "TS", "")
    isempty(s) ? [0.05, 0.15, 0.25, 0.30, 0.30625] : parse.(Float64, split(s, ','))
end
const SPAN = -0.03          # blabnat: 1.0 → 0.97

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))
apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("benchmark: $(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= 1e-8 "benchmark regressed"

sv = vars["blabnat"]
lam = sv isa AbstractArray ? first(sv) : sv

# ── Evaluator over the SQUARED model ────────────────────────────────────────
# solve_newton! has already pinned orphans and deleted dead rows, so the row and
# column counts below should match. Build the evaluator once and reuse it: it is
# the expensive part, and the model's *structure* does not change as λ moves.
allv = JuMP.all_variables(m)
allidx = [JuMP.index(v) for v in allv]
pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))

nowfree = VariableRef[]
for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr) || push!(nowfree, vr)
end
N = length(nowfree)
freecols = [pos[JuMP.index(v).value] for v in nowfree]
col2free = zeros(Int, length(allv))
for (i, c) in enumerate(freecols); col2free[c] = i; end
lamcol = pos[JuMP.index(lam).value]
@assert col2free[lamcol] == 0 "blabnat is ENDOGENOUS — it cannot be the continuation parameter"

nlmodel = MOI.Nonlinear.Model()
# Ref, not a plain Int: a top-level `for` opens a SOFT scope, so `nrows += 1`
# inside it would create a new local and leave the global undefined.
nrows_ref = Ref(0)
for (Ftype, S) in list_of_constraint_types(m)
    Ftype <: VariableRef && continue
    for con in all_constraints(m, Ftype, S)
        co = JuMP.constraint_object(con)
        r = IndotermJulia._rhs_of(co.set)   # same helper solve_newton! uses, so
                                            # the row ORDER here matches its own
        MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r))
        nrows_ref[] += 1
    end
end
const nrows = nrows_ref[]
evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
MOI.initialize(evaluator, [:Grad, :Jac])
println("evaluator: $nrows rows × $N free columns", nrows == N ? "  (square ✓)" : "  ⚠ NOT SQUARE")

st = MOI.jacobian_structure(evaluator)
jr = getindex.(st, 1); jc_all = getindex.(st, 2)
Jbuf = zeros(length(jr))
keep = [col2free[c] != 0 for c in jc_all]
jrk = jr[keep]; jck = [col2free[c] for c in jc_all[keep]]
lamk = [c == lamcol for c in jc_all]
lamrows = jr[lamk]
println("∂F/∂λ has $(count(lamk)) structural nonzero(s)")

"""
Tangent `dx/dλ` at the current point, plus a σ_min estimate, by direct LU.
"""
function tangent_at()
    x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
         (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]
    MOI.eval_constraint_jacobian(evaluator, Jbuf, x)

    J = sparse(jrk, jck, Jbuf[keep], nrows, N)
    b = zeros(nrows)
    for (k, r) in zip(findall(lamk), lamrows); b[r] += Jbuf[k]; end

    # Equilibrate rows. Leaves v unchanged (both sides scale) but makes the LU
    # usable on rows spanning 1e-10 … 1e12.
    rs = [r > 0 ? r : 1.0 for r in vec(maximum(abs, J; dims=2))]
    Ds = sparse(Diagonal(1.0 ./ rs))
    Js = Ds * J
    bs = b ./ rs

    F = lu(Js)
    v = -(F \ bs)

    # σ_min ≈ 1/‖J⁻¹‖₂, by power iteration on (JᵀJ)⁻¹ reusing the SAME LU.
    # This is on the ROW-SCALED matrix, so compare the trend across points, not
    # the absolute number against any other script's conditioning figure.
    z = normalize!(randn(N))
    inv2 = 0.0
    for _ in 1:8
        w = F \ (F' \ z)
        inv2 = norm(w)
        z = w ./ inv2
    end
    return v, (inv2 > 0 ? 1 / sqrt(inv2) : NaN)
end

name_of = Dict(i => JuMP.name(nowfree[i]) for i in 1:N)

println("\n" * "="^78)
println("TANGENT NORM ALONG THE BRANCH")
println("="^78)
rows = Tuple{Float64,Float64,Float64,String,Float64}[]
for tt in TS
    JuMP.fix(lam, 1.0 + tt * SPAN; force=true)
    r = solve_newton!(m, vars; maxit=40, verbose=false)
    if r.residual > 1e-8
        println("\nt=$tt  ❌ did not converge (‖F‖∞=$(r.residual)) — cannot read a tangent here")
        break
    end
    v, smin = tangent_at()
    k = argmax(abs.(v))
    push!(rows, (tt, 1.0 + tt * SPAN, abs(v[k]), name_of[k], smin))
    println("\n── t=$tt  (blabnat=$(round(1.0 + tt*SPAN; sigdigits=8)))  " *
            "converged in $(r.iters) it, ‖F‖∞=$(round(r.residual; sigdigits=4))")
    println("   ‖dx/dλ‖∞ = $(round(abs(v[k]); sigdigits=6))   at  $(name_of[k])")
    println("   σ_min(J_scaled) ≈ $(round(smin; sigdigits=4))")
    ord = sortperm(abs.(v); rev=true)
    println("   top 5 tangent components:")
    for i in ord[1:min(5, N)]
        println("     $(rpad(name_of[i], 26)) $(round(v[i]; sigdigits=6))")
    end
end

println("\n" * "="^78)
println("VERDICT")
println("="^78)
println("t          blabnat        ‖dx/dλ‖∞        1/‖dx/dλ‖∞     σ_min      argmax")
for (tt, lv, nv, nmx, sm) in rows
    println(rpad(tt, 11), rpad(round(lv; sigdigits=8), 15), rpad(round(nv; sigdigits=6), 16),
            rpad(round(1 / nv; sigdigits=6), 15), rpad(round(sm; sigdigits=4), 11), nmx)
end

if length(rows) >= 2
    growth = rows[end][3] / rows[1][3]
    println("\n‖dx/dλ‖∞ grew $(round(growth; sigdigits=4))× from t=$(rows[1][1]) to t=$(rows[end][1]).")
    println("""
    HOW TO READ THIS
      * Blowing up (growth ≫ 10, and 1/‖dx/dλ‖∞ heading linearly to zero):
        genuine FOLD. Extrapolate the 1/‖dx/dλ‖∞ column to zero for its location,
        and read the argmax column for the DIRECTION — that names the block to
        look at, and it is the first honest localisation of this obstruction.
        Remedy is pseudo-arclength continuation, which parameterises by arclength
        and stays nonsingular through the turn.
      * Flat and boring: NOT a fold. The branch continues and the failure is the
        solver — specifically the trust-region ratchet at solve_newton!.jl:636-638
        (PLAN.md §N.10), which contracts the radius using a `moved` that is small
        only BECAUSE the step was damped. Fix that before touching the model.""")
end
println("\nDone.")
