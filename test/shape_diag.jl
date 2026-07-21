using IndotermJulia, NamedArrays

println("Data shapes diagnostic")
println("="^60)

nat = read_national_data()
regsup = read_regsupp_data()
distgone = read_distgone_data()

# Print shape of each national header
ks = collect(keys(nat)); sort!(ks)
for k in ks; v = nat[k]
    if v isa NamedArray
        println("nat[\"$k\"]: ndims=$(ndims(v)), size=$(size(v))")
    elseif v isa Float64
        println("nat[\"$k\"]: scalar=$v")
    else
        println("nat[\"$k\"]: $(typeof(v))")
    end
end
println()
ks2 = collect(keys(regsup)); sort!(ks2)
for k in ks2; v = regsup[k]
    if v isa NamedArray
        println("regsup[\"$k\"]: ndims=$(ndims(v)), size=$(size(v))")
        println("  dimnames: $(v.dimnames)")
    elseif v isa Float64
        println("regsup[\"$k\"]: scalar=$v")
    else
        println("regsup[\"$k\"]: $(typeof(v))")
    end
end
println()
ks3 = collect(keys(distgone)); sort!(ks3)
for k in ks3; v = distgone[k]
    if v isa NamedArray
        println("distgone[\"$k\"]: ndims=$(ndims(v)), size=$(size(v))")
    end
end
