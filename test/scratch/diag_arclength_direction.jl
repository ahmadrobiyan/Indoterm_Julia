"""
Fast repro for the task #39 arclength-direction bug (PLAN.md §N.22, VV_PLAN.md "Open items").

`verify_arclength_fallback.jl` takes 35 minutes to reach the failure. This script skips
straight to a point near the documented fold (t≈0.30, blabnat≈0.991) with a single Newton
solve, then runs `arclength_solve!` with a SMALL `maxsteps` so the first several accepted
steps — and in particular the sign of dλ on step 1 — are visible in under a minute.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_arclength_direction.jl
"""
nothing
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, Printf

agg6, params = cached_pipeline(6)

m, vars = build_model_full!(agg6, params)
initialize_model!(m, vars; bmk_levels=benchmark_levels(params), numeraire=:gdppi)
apply_swaps!(m, vars, collect(TERM_CMF_SWAPS); verbose=false)

r0 = solve_newton!(m, vars; maxit=5, tol=1e-8, verbose=false)
@printf("benchmark: status=%s ||F||inf=%.3e\n", r0.status, r0.residual)

blabnat = vars["blabnat"]
b = 1.0
tg = pct(-3)   # 0.97
t0 = 0.30
lam0 = b + t0 * (tg - b)
@printf("stepping directly to t=%.4f  (blabnat %.6f -> %.6f)\n", t0, b, lam0)
JuMP.fix(blabnat, lam0; force=true)
r1 = solve_newton!(m, vars; maxit=30, tol=1e-8, verbose=false)
@printf("t=%.4f solve: status=%s ||F||inf=%.3e\n", t0, r1.status, r1.residual)
r1.residual <= 1e-8 || error("direct jump to t=$t0 did not converge — need a smaller t0")

println("\n" * "="^72)
println("arclength_solve! from here, target blabnat = $tg, maxsteps=12")
println("="^72)
ra = arclength_solve!(m, vars, blabnat, tg;
                      maxsteps=12, tol=1e-8, verbose=true,
                      report=("blabnat",))

@printf("\nreached=%s turned=%s lam_fold=%.6g lam_reached=%.10g\n",
        ra.reached_target, ra.turned, ra.lam_fold, ra.lam_reached)
println("path (lam, iters, resid):")
for (l, it, rr) in ra.path
    @printf("  lam=%.8f  it=%d  resid=%.3g\n", l, it, rr)
end
println("Done.")
