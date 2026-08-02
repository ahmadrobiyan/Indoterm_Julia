push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia

nat = read_national_data()
regsup = read_regsupp_data()
distgone = read_distgone_data()
reg0_r = build_reg0!(nat, regsup)
reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)

# Warm up by calling once
println("First call...")
t0 = time()
ras_r1 = ras_balance!(reg2_r.reg2)
println("First: ", time()-t0, "s")

# Re-build the input (ras_balance! mutates reg2)
reg2_r2 = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
println("Second call...")
t0 = time()
ras_r2 = ras_balance!(reg2_r2.reg2)
println("Second: ", time()-t0, "s")
