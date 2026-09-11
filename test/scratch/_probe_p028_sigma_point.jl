"""
One point of the σ-continuation, verbose. Diagnoses why `_probe_p028_sigma_continuation.jl`
stalled at `maxit` for every f < 0.95: the residual at the cap fell roughly in proportion to
the step (4.2e-2, 3.2e-3, 1.7e-3, 8.9e-4, 1.0e-4) — the signature of slow damped progress,
not of a fold. This runs a single f with `verbose=true` and a large `maxit`, warm-started from
the f=1.0 anchor exactly as the walk did, and prints the per-iteration trace.

Run:  julia --project=. test/scratch/_probe_p028_sigma_point.jl [f=0.9375] [maxit=200]
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, Printf

const F     = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.9375
const MAXIT = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200
const TOL   = 1e-8

function warm_start!(vars, values)
    for (nm, v) in vars
        haskey(values, nm) || continue
        src = values[nm]
        for (k, vr) in enumerate(v isa AbstractArray ? vec(v) : [v])
            JuMP.is_fixed(vr) && continue
            JuMP.set_start_value(vr, src isa AbstractArray ? vec(src)[k] : src)
        end
    end
end

agg6, params = cached_pipeline(6)
sc = COALPRICE_REFERENCE
(spec, target) = only(sc.shocks)

println("anchor f=1.0 …"); flush(stdout)
r1 = run_model!(agg6, params, sc; h0=0.25, tol=TOL, gdp=false, verbose=false)
r1.solved || error("anchor failed")
println("anchor ‖F‖∞=$(r1.residual)"); flush(stdout)

agg2 = copy(agg6); agg2["P028"] = vec(parent(agg6["P028"])) .* F
p2 = prepare_parameters!(agg2)
m, vars = build_model_full!(agg6, p2)
initialize_model!(m, vars; bmk_levels=benchmark_levels(p2), numeraire=sc.numeraire)
apply_swaps!(m, vars, collect(sc.swaps); verbose=false)
r0 = solve_newton!(m, vars; maxit=5, tol=TOL, verbose=false)
println("benchmark at f=$F: ‖F‖∞=$(r0.residual)"); flush(stdout)
_, _, vr = IndotermJulia._shock_ref(vars, spec)
JuMP.fix(vr, Float64(target); force=true)
warm_start!(vars, r1.values)

println("\n── Newton at f=$F, full shock, warm from f=1.0, maxit=$MAXIT"); flush(stdout)
t0 = time()
r = solve_newton!(m, vars; maxit=MAXIT, tol=TOL, verbose=true)
@printf("\nRESULT f=%.4f status=%s iters=%d ‖F‖∞=%.3e  %.0fs\n", F, r.status, r.iters, r.residual, time()-t0)
flush(stdout)
