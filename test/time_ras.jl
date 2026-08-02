push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia

nat = read_national_data()
regsup = read_regsupp_data()
distgone = read_distgone_data()

println("\n--- build_reg0! ---")
t0 = time(); reg0_r = build_reg0!(nat, regsup); println("build_reg0!: $(round(time()-t0,digits=2))s")

println("\n--- build_reg1! ---")
t0 = time(); reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone); println("build_reg1!: $(round(time()-t0,digits=2))s")

println("\n--- build_reg2! ---")
t0 = time(); reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag); println("build_reg2!: $(round(time()-t0,digits=2))s")

println("\n--- ras_balance! ---")
t0 = time(); ras_r = ras_balance!(reg2_r.reg2); println("ras_balance!: $(round(time()-t0,digits=2))s")

println("\nDone.")
