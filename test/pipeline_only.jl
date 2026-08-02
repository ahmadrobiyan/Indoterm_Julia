push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia

println("Step 1: read_national_data")
nat = read_national_data(); println("  OK")

println("Step 2: read_regsupp_data")
regsup = read_regsupp_data(); println("  OK")

println("Step 3: read_distgone_data")
distgone = read_distgone_data(); println("  OK")

println("Step 4: build_reg0!")
reg0_r = build_reg0!(nat, regsup); println("  OK")

println("Step 5: build_reg1!")
reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone); println("  OK")

println("Step 6: build_reg2!")
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag); println("  OK")

println("Step 7: ras_balance!")
ras_r  = ras_balance!(reg2_r.reg2); println("  OK")

println("Step 8: build_pstras!")
pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone); println("  OK")

println("Step 9: build_premod!")
prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast); println("  OK")

println("Step 10: aggregate_model!")
agg_r  = aggregate_model!(prem_r.premod); println("  OK")

println("Step 11: aggregate_regions!")
agg6 = aggregate_regions!(agg_r.agg, REG_MAP_34_to_6, 6); println("  OK")

println("Step 12: prepare_parameters!")
params = prepare_parameters!(agg6); println("  OK")

println("ALL PIPELINE OK")
