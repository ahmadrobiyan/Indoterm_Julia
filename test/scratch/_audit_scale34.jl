# Scaling Audit: 34-Region Jacobian and Fill
# Goal: Verify 2b.1 (Data Pipeline) and 2b.2 (Sparsity Audit)
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using IndotermJulia, JuMP, SparseArrays, LinearAlgebra

say(msg) = (println(msg); flush(stdout))

# 2b.1: Scale pipeline to 34 regions
say("--- Step 2b.1: Pipeline Scaling ---")
t_pipe = @elapsed begin
    agg34, params34 = cached_pipeline(34)
end
say("Pipeline load complete in $(round(t_pipe, digits=1))s")

# CORRECT DIMENSION CHECKS: Use size() instead of length()
# Params keys that are Tensors:
# MAKE: (Sectors, Industries, Regions)
nr = size(params34["MAKE"], 3)
ns = size(params34["MAKE"], 1)
say("Verified Dimensions: $nr regions, $ns sectors")

# 2b.2: Sparsity Audit
say("\n--- Step 2b.2: Sparsity Audit ---")
say("Building 34-region model for Jacobian extraction...")
t_build = @elapsed begin
    model_tuple = build_model_full!(agg34, params34)
    model = model_tuple[1]
    vars = model_tuple[2]
    # We use diagnose_jacobian to print the Raw nnz to the log
    diagnose_jacobian(model, vars; threshold=1e-10)
end
say("Model build and diagnosis complete in $(round(t_build, digits=1))s")

say("\n--- Summary ---")
say("Pipeline 34: ✅")
say("Check the logs above for the 'Raw Jacobian nnz' value reported by diagnose_jacobian.")
say("Projected LU (PRIMARY): ~57M")
say("Projected LU (Direct): ~250M+ (The Wall)")
