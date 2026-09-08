push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP, Statistics

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

println("Building model...")
m, vars = build_model_full!(agg6, params)
bmk = benchmark_levels(params)
initialize_model!(m, vars; bmk_levels=bmk)

println("Checking closure...")
for nm in ["alab_o", "blabnat", "blab", "blab_d", "xlab_o"]
    v = vars[nm]
    if v isa AbstractArray
        fixed = [JuMP.is_fixed(v[idx]) for idx in eachindex(v)]
        n_fixed = sum(fixed)
        n_total = length(fixed)
        println("  $nm: $n_fixed/$n_total fixed")
        if n_fixed > 0
            println("    first fixed value = ", JuMP.fix_value(v[1]))
        end
    else
        println("  $nm: fixed=$(JuMP.is_fixed(v)), value=$(JuMP.is_fixed(v) ? JuMP.fix_value(v) : JuMP.start_value(v))")
    end
end

# Now run the continuation and dump alab_o
blabnat = vars["blabnat"]
bmk_blabnat = JuMP.fix_value(blabnat)
println("\nStarting from blabnat = $bmk_blabnat")

# Quick test: move blabnat just 1 step and check
res = euler_continuation!(m, vars, blabnat, 0.97, 1; tol=1e-7, verbose=true)
vals = res.values

println("\nAfter 1 step:")
println("  blabnat = ", JuMP.fix_value(blabnat))
println("  alab_o range: ", round(minimum(vals["alab_o"]);digits=6), " to ", round(maximum(vals["alab_o"]);digits=6))
println("  alab_o first 5: ", round.(vals["alab_o"][1:5,1];digits=6))
println("  xlab_o first 5: ", round.(vals["xlab_o"][1:5,1];digits=6))
