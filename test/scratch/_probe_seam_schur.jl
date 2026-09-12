# Stage 2a: schur-only verbose leg to localize the stall (scratch, logs only).
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
say(msg) = (println(msg); flush(stdout))

agg6, params = cached_pipeline(6)
t = @elapsed r = run_model!(agg6, params, COALPRICE_REFERENCE;
                            h0=0.25, tol=1e-8, gdp=false, verbose=true,
                            linsolve=:schur)
say("schur leg: solved=$(r.solved) t_reached=$(r.t_reached) " *
    "steps/rejects=$(r.nsteps)/$(r.nrejects) res=$(r.residual) $(round(t/60;digits=1))min")
