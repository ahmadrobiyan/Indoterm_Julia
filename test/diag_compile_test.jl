push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP

println("pipeline..."); flush(stdout)
nat = read_national_data(); regsup = read_regsupp_data(); distgone = read_distgone_data()
reg0_r = build_reg0!(nat, regsup)
reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
ras_r  = ras_balance!(reg2_r.reg2)
pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone)
prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast)
agg6 = aggregate_regions!(aggregate_model!(prem_r.premod).agg, REG_MAP_34_to_6, 6)
params = prepare_parameters!(agg6)

println("build..."); flush(stdout)
t0 = time()
m, vars = build_model_full!(agg6, params)
println("Build: $(round(time()-t0;digits=1))s"); flush(stdout)

println("init..."); flush(stdout)
bmk = benchmark_levels(params)
initialize_model!(m, vars; bmk_levels=bmk)

println("solve..."); flush(stdout)
t0 = time()
res = solve_newton!(m, vars; maxit=1, tol=1e-8, verbose=true)
println("Solve: $(round(time()-t0;digits=1))s status=$(res.status) ||F||=$(res.residual)")
println("Done!")
