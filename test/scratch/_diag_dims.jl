# Dimension Diagnostic: 34-Region Params Structure
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using IndotermJulia

say(msg) = (println(msg); flush(stdout))

say("--- 34-Region Structure Diagnostic ---")
agg34, params34 = cached_pipeline(34)

# Check a few key arrays
keys_to_check = ["MAKE", "TRADE", "1LAB", "SLAB", "P028", "LAB"]

for k in keys_to_check
    if haskey(params34, k)
        val = params34[k]
        say("Key: $k")
        say("  Type: $(typeof(val))")
        say("  Size: $(size(val))")
        say("  Length: $(length(val))")
    else
        say("Key: $k -> NOT FOUND")
    end
end

# Check if they are flattened
is_flattened = false
if haskey(params34, "MAKE") && ndims(params34["MAKE"]) == 1
    is_flattened = true
end

say("\nFlattening Detected: $is_flattened")
