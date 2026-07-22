push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia

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

LAB_O = params["LAB_O"]
println("LAB_O size: ", size(LAB_O))
println("nonzero (>1e-10) count: ", count(x -> x > 1e-10, LAB_O), " / ", length(LAB_O))
println("sum(LAB_O) = ", sum(LAB_O))
println("max(LAB_O) = ", maximum(LAB_O))
println("sample LAB_O[1:3,1:3] = ", LAB_O[1:3,1:3])

SUPPMAR_D_arr = get(params, "SUPPMAR_D_idx", nothing)
println("\n--- margin data checks ---")
println("has TMAR: ", haskey(agg_r.agg, "TMAR"))
if haskey(agg_r.agg, "TMAR")
    TMAR = agg_r.agg["TMAR"]
    println("TMAR size: ", size(TMAR), " sum: ", sum(TMAR))
end
