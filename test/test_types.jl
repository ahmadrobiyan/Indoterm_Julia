push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia

nat = read_national_data()
regsup = read_regsupp_data()
distgone = read_distgone_data()
reg0_r = build_reg0!(nat, regsup)
reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)

td = reg2_r.reg2["TRAD"]
println("TRADE_data type: ", typeof(td))
println("TD parent type: ", typeof(parent(td)))
println("copy(parent) type: ", typeof(copy(parent(td))))

# Time the copy
t0 = time()
for _ in 1:1000
    c = copy(parent(td))
end
println("1000 copies: ", time()-t0, "s")

# Time a simple sum
t0 = time()
for _ in 1:1000
    s = sum(td)
end
println("1000 sums: ", time()-t0, "s")
