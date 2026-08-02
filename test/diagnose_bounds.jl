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
agg6 = aggregate_regions!(agg_r.agg, REG_MAP_34_to_6, 6)
params = prepare_parameters!(agg6)
m, vars = build_model_full!(agg6, params)
bmk = benchmark_levels(params)
shock = Dict{String,Any}("blabnat" => exp(-3.0 / 100.0))
initialize_model!(m, vars; bmk_levels=bmk, shocks=shock)

# Find variables that are free, have a finite lower bound, and are close to it
println("Free variables near lower bound:")
for (nm, v) in vars
    for idx in eachindex(v isa AbstractArray ? v : [v])
        vr = v isa AbstractArray ? v[idx] : v
        JuMP.is_fixed(vr) && continue
        lb = JuMP.has_lower_bound(vr) ? JuMP.lower_bound(vr) : -Inf
        isfinite(lb) || continue
        sv = JuMP.start_value(vr)
        gap = sv - lb
        if gap < 0.1
            idx_str = v isa AbstractArray ? "[$idx]" : ""
            println("  $nm$idx_str: lb=$lb  start=$sv  gap=$gap")
        end
    end
end
println("Done.")
