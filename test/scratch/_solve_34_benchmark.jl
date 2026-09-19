# 34-Region Benchmark Solve
# Goal: Prove that the Schur solver allows the full model to solve without OOM.
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using IndotermJulia, JuMP

say(msg) = (println(msg); flush(stdout))

say("--- Stage 2b.3: 34-Region Benchmark Solve ---")

# 1. Load 34-region data
t_pipe = @elapsed begin
    agg34, params34 = cached_pipeline(34)
end
say("Pipeline load: $(round(t_pipe, digits=1))s")

# 2. Define a simple benchmark scenario
# Using the Scenario constructor properly: name and source are required.
sc = Scenario(
    name = "34reg_benchmark",
    source = "benchmark",
    shocks = [Pair("blabnat", 0.99)], # -1% labor supply
    swaps = [],
    pct_shocks = []
)

# 3. Execute Solve with Schur
say("\nLaunching solve with linsolve = :schur...")
t_solve = @elapsed begin
    r = run_model!(agg34, params34, sc; 
                   linsolve = :schur, 
                   verbose = true)
end

say("\n--- RESULT ---")
say("Solved: $(r.solved)")
say("Residual: $(r.residual)")
say("T reached: $(r.t_reached)")
say("Solve time: $(round(t_solve, digits=1))s")

if r.solved
    say("✅ SUCCESS: 34-region model solved using Schur complement.")
else
    say("❌ FAILURE: Solve did not converge.")
end
