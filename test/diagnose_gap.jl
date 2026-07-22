push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP

nat = read_national_data()
regsup = read_regsupp_data()
distgone = read_distgone_data()
reg0_r = build_reg0!(nat, regsup)
reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
ras_r  = ras_balance!(reg2_r.reg2)
pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone)
prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast)
agg_r  = aggregate_model!(prem_r.premod)
params = prepare_parameters!(agg_r.agg)

m, vars = build_model_full!(agg_r.agg, params)
initialize_model!(m, vars)

println("Collecting variables touched by any affine constraint...")
touched = Set{VariableRef}()
for (F, S) in list_of_constraint_types(m)
    F <: GenericAffExpr || continue
    for c in all_constraints(m, F, S)
        for v in keys(constraint_object(c).func.terms)
            push!(touched, v)
        end
    end
end
println("Touched: $(length(touched))")

free_vars = filter(v -> !is_fixed(v), all_variables(m))
println("Free (non-fixed): $(length(free_vars))")

untouched = [v for v in free_vars if !(v in touched)]
println("Untouched free vars (no defining/using equation at all): $(length(untouched))")

# breakdown by base name (strip trailing [indices])
counts = Dict{String,Int}()
for v in untouched
    nm = name(v)
    base = split(nm, "[")[1]
    counts[base] = get(counts, base, 0) + 1
end
for (k, v) in sort(collect(counts), by = x -> -x[2])
    println("  $k => $v")
end
