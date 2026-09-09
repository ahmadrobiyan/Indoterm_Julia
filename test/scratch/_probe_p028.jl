# Diagnostic for the sole incomplete point of the V8 sweep (P028 x0.5, folded at t = 0.6426).
#
# Two hypotheses, distinguished by this script:
#   (A) continuation limit — the path exists but the stepper cannot track it. A finer hmin
#       and a larger Newton budget should then reach t = 1.
#   (B) genuine fold — no setting reaches t = 1, and the largest factor that still solves
#       locates the boundary. Halving P028 takes the primary-factor CES to ~0.1 (near
#       Leontief), so a fold is economically plausible, not a bug.
#
# Reports t_reached for every attempt. It never claims success from a small residual: a
# stalled branch reports the residual at its last ACCEPTED point, which is always small.

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using Printf

function attempt(agg6, params, factor, label; kw...)
    agg2 = copy(agg6)
    agg2["P028"] = vec(parent(agg6["P028"])) .* factor
    p2 = prepare_parameters!(agg2)
    @assert vec(parent(p2["P028"])) != vec(parent(params["P028"]))
    println("── $label")
    # Redirected stdout is block-buffered and this probe prints too little to overflow the
    # buffer, so without explicit flushes the log stays empty for the entire run and progress
    # is unobservable — which is exactly what happened on the 2026-09-08 launch.
    flush(stdout)
    t0 = time()
    res = run_model!(agg2, p2, COALPRICE_REFERENCE;
                     tol = 1e-8, gdp = false, verbose = false, kw...)
    dt = time() - t0
    if res.solved
        @printf("   SOLVED in %.0fs, residual %.2e\n", dt, res.residual)
    else
        @printf("   INCOMPLETE after %.0fs — t = %.4f, %d steps, %d rejections\n",
                dt, res.t_reached, res.nsteps, res.nrejects)
    end
    flush(stdout)
    return res.solved, res.t_reached
end

function main()
    agg6, params = cached_pipeline(6)
    println("\nP028 x0.5 diagnostic — baseline run reached t = 0.6426 (h0=0.25, hmin=1e-4)\n")

    ok, t = attempt(agg6, params, 0.5, "P028 x0.5, hmin=1e-6, h0=0.05, maxit=60";
                    h0 = 0.05, hmin = 1e-6, maxit = 60)
    if ok
        println("\nHypothesis A: continuation limit. The point is recoverable with a finer")
        println("floor; the sweep can be certified once sensitivity.jl uses these settings")
        println("for this point.")
        return
    end
    @printf("\nRetry also folded (t = %.4f). Bracketing the boundary.\n\n", t)

    for f in (0.6, 0.75)
        ok2, t2 = attempt(agg6, params, f, "P028 x$f, hmin=1e-6, h0=0.05, maxit=60";
                          h0 = 0.05, hmin = 1e-6, maxit = 60)
        if ok2
            @printf("\nHypothesis B: genuine fold. P028 x%.2f solves; x0.5 does not.\n", f)
            println("The lower half of the P028 range is not attainable under this scenario.")
            println("Report the sweep's P028 coverage as [x$f, x1.5], not [x0.5, x1.5].")
            return
        end
    end
    println("\nNeither 0.6 nor 0.75 completed — the fold is not narrowly located.")
end

main()
