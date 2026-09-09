include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
agg6, params = cached_pipeline(6)

function same_value(a, b)
    try
        a === b && return true
        if a isa AbstractArray && b isa AbstractArray
            return size(a) == size(b) && all(a .== b)
        end
        return isequal(a, b)
    catch
        return false
    end
end

KEYS = ["P015", "P028", "SLAB", "SCET", "SMAR", "P018"]
println()
for key in KEYS, f in (0.5, 1.5)
    agg2 = copy(agg6)
    agg2[key] = vec(parent(agg6[key])) .* f
    p2 = prepare_parameters!(agg2)
    changed = sort([k for k in keys(p2)
                    if k != key && !same_value(get(params, k, nothing), get(p2, k, nothing))])
    println(rpad("$key x$f", 12), " -> ", length(changed), "  ", join(changed, ", "))
end
