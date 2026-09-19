# Stage 2b.4: 34-region small-shock homotopy probe (scratch, logs only)
# Benchmark is converged (8.6e-9); 1% blabnat shock stalled the line search (memledger_34_b.log).
# This tries 0.1% first, then 1% with small h0, on the :schur path.
# Run: julia --project=. test/scratch/_solve_34_smallshock.jl
# Log:  logs/solve_34_smallshock_<date>.log
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using IndotermJulia

say(msg) = (println(msg); flush(stdout))
say("--- Stage 2b.4: 34-region small-shock probe ---")

t_pipe = @elapsed begin
    agg34, params34 = cached_pipeline(34)
end
say("pipeline: $(round(t_pipe; digits=1))s (cache HIT expected)")

for (tag, val, h0) in (("0p1pct", 0.999, 0.05), ("1pct", 0.99, 0.05))
    say("\n=== shock $tag: blabnat 1.0 -> $val (h0=$h0) ===")
    t = @elapsed begin
        r = run_model!(agg34, params34;
            shocks = [Pair("blabnat", val)],
            swaps = (),
            name = "34reg_$tag",
            h0 = h0, hmin = 1e-4, hmax = 0.25,
            tol = 1e-8, maxit = 30,
            linsolve = :schur, verbose = true, gdp = false)
    end
    say("RESULT $tag: solved=$(r.solved) residual=$(r.residual) t_reached=$(r.t_reached) time=$(round(t; digits=1))s")
    flush(stdout)
end
say("\n--- PROBE DONE ---")
