push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP

println("="^60)
println("INDOTERM FULL MODEL BUILD (Step 5a) — full 25×34 scale")
println("="^60)

println("\n--- Running data pipeline (Steps 0-4) ---")
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
println("Pipeline OK — agg keys: $(length(agg_r.agg)), params: $(length(params))")

println("\n--- build_model_full! (25×34 scale) ---")
t0 = time()
m, vars = build_model_full!(agg_r.agg, params)
elapsed = time() - t0
println("Build completed in $(round(elapsed, digits=1))s")

nvars = num_variables(m)
ncons = num_constraints(m; count_variable_in_set_constraints=false)
println("Variables:   $nvars")
println("Constraints: $ncons")
println("vars dict entries: $(length(vars))")

println("\n"^2)
println("="^60)
println("FULL MODEL BUILD COMPLETED — vars=$nvars, cons=$ncons")
println("="^60)
