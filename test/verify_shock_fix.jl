"""
Gate: does a non-zero shock now converge?

Three defects were found and fixed after `diag_stall_direct.jl` / `diag_step_quality.jl`:

 1. `solve_newton!.jl` converted the CGNR direction back to original units with
    `dy .* Dc_ruiz`; `ruiz_scales` returns Dr/Dc as DIVISORS, so the correct
    transform is `dy ./ Dc`. Every step was corrupted by Dc² per variable, giving a
    non-descent direction. (Verified: with `d = dy./Dc`, ‖Jbase·d + gs‖/‖gs‖ = 3.8e-10.)
 2. The `ppur_s` zero-flow guard used `PUR_S > 1e-10`, but the smallest strictly
    positive PUR_S in the 25x6 database is 2.7e-8 — so the guard only ever pinned
    cells that were already exactly zero. `PUR_S[5,16,6] = 9.95e-8` left a price
    index amplified by 1/X ≈ 1e7: its Newton step was -769, which alone throttled
    fraction-to-boundary to α ≈ 1.2e-3 and drove its own log argument negative.
    Threshold is now `_MIN_PRICED_FLOW = 1e5*_TINY`.
 3. The solver line-searched a direction CGNR had not actually converged to
    (200 iterations at κ ≈ 2.6e10 gave max|dy| = 0.025 against a true 1221), and
    then clipped it with `max(0.10|x|, 10.0)` — a cap unrelated to the scale of the
    value-denominated `del*` variables, which need steps of 300–1221. Now: direct
    sparse LU first, CGNR/LM as fallback, and no magnitude cap on a verified step.

Pass criteria:
  * benchmark still replicates: ‖F‖∞ <= 1e-8;
  * blabnat = 0.9997 CONVERGES (it used to stall at 2.126e-4, then 1.870e-4);
  * a 30x larger shock (blabnat = 0.99) also converges.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_shock_fix.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, LinearAlgebra

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

nfree = 0
for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr) || (global nfree += 1)
end
ncon = sum(num_constraints(m, Ft, S) for (Ft, S) in list_of_constraint_types(m)
           if !(Ft <: VariableRef); init=0)
println("free variables = $nfree   equality constraints = $ncon   square = $(nfree == ncon)")

println("\n══ benchmark ══")
r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("status=$(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= 1e-8 "benchmark replication regressed: $(r0.residual)"

sv = vars["blabnat"]
svr = sv isa AbstractArray ? first(sv) : sv

for shock in (0.9997, 0.99)
    println("\n══ shock blabnat = $shock ══")
    JuMP.fix(svr, shock; force=true)
    el = @elapsed r = solve_newton!(m, vars; maxit=60, verbose=true)
    println("status=$(r.status)  iters=$(r.iters)  ‖F‖∞=$(r.residual)  ($(round(el,digits=1))s)")
    if r.status == :converged
        println("  ✅ CONVERGED")
    else
        println("  ❌ did not converge")
    end
end

println("\nDone.")
