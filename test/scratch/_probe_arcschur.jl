# Arclength fallback under linsolve=:schur on the INCOMPLETE fixture (scratch).
# Pass = TURNED near λ*≈0.9908083 with continuity, via bordered solves.
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
say(msg) = (println(msg); flush(stdout))

agg6, params = cached_pipeline(6)
el = @elapsed r = run_model!(agg6, params, TERM_CMF_NO_DELUNITY;
                             h0=0.25, hmin=1e-3, tol=1e-8,
                             arclength=true, gdp=false, verbose=true,
                             linsolve=:schur)
say("SCHUR-ARC RESULT ($(round(el/60;digits=1)) min): solved=$(r.solved) " *
    "t=$(round(r.t_reached;digits=5)) steps/rejects=$(r.nsteps)/$(r.nrejects) " *
    "res=$(r.residual) pathpoints=$(length(r.path))")
