using IndotermJulia, Printf

println("Loading data...")
nat = read_national_data()
regsup = read_regsupp_data()
distgone = read_distgone_data()

println("Running reg0...")
@time reg0_r = build_reg0!(nat, regsup)

println("Running reg1...")
@time reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
println("reg1 done!")
