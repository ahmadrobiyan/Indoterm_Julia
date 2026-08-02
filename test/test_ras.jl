push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia

nat = read_national_data(); regsup = read_regsupp_data(); distgone = read_distgone_data()
reg0_r = build_reg0!(nat, regsup); reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
ras_r = ras_balance!(reg2_r.reg2)
