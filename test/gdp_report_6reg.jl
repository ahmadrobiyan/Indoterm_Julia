push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP

println("="^72)
println("STEP 8 — GDP REPORTING + BOTH-SIDES CHECK @ 25×6")
println("="^72)

nat = read_national_data(); regsup = read_regsupp_data(); distgone = read_distgone_data()
reg0_r = build_reg0!(nat, regsup)
reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
ras_r  = ras_balance!(reg2_r.reg2)
pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone)
prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast)
agg6 = aggregate_regions!(aggregate_model!(prem_r.premod).agg, REG_MAP_34_to_6, 6)
params = prepare_parameters!(agg6)

m, vars = build_model_full!(agg6, params)
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))
res = solve_newton!(m, vars; verbose=false)
println("solve: status=$(res.status)  scaled ||F||_inf=$(res.residual)")
@assert res.solved "benchmark solve failed — cannot report GDP"

report = calculate_gdp(res.values, params)
print_gdp_report(report; regions=REG6)

# ── Benchmark sanity: at a zero-shock solve every (change) component is ~0, so
#    the reported levels must equal the benchmark GDPINCSUM/GDPEXPSUM exactly.
GDPINCSUM = parent(params["GDPINCSUM"]); GDPEXPSUM = parent(params["GDPEXPSUM"])
max_inc_dev = maximum(abs.(report.inc .- GDPINCSUM))
max_exp_dev = maximum(abs.(report.exp_ .- GDPEXPSUM))
println("\nBenchmark deviation — income components: $max_inc_dev, expenditure: $max_exp_dev")

scale = max(maximum(abs.(GDPINCSUM)), 1.0)
@assert max_inc_dev / scale < 1e-6 "income components moved at a zero shock"
@assert max_exp_dev / scale < 1e-6 "expenditure components moved at a zero shock"

# ── Slot-order regression guard ─────────────────────────────────────────────
# TERM.TAB indexes delGDPINC by CATEGORY NAME, so the equation letters do not run
# parallel to GDPINCCAT: E_delGDPINCb writes "Capital" (slot 3) and E_delGDPINCc
# writes "Labour" (slot 2). The port originally mapped a/b/c positionally, which
# silently swapped Labour and Capital. Totals were unaffected, so benchmark
# replication could not detect it — this component-level check can.
LAB_IO = parent(params["LAB_IO"]); CAP_I = parent(params["CAP_I"])
lab_slot = findfirst(==("Labour"), GDPINCCAT)
cap_slot = findfirst(==("Capital"), GDPINCCAT)
for d in 1:length(LAB_IO)
    @assert isapprox(report.inc[d, lab_slot], LAB_IO[d]; rtol=1e-6) """
        GDPINCCAT slot mismatch: "Labour" (slot $lab_slot) reads \
        $(report.inc[d, lab_slot]) but LAB_IO[$d] = $(LAB_IO[d]). \
        Check E_delGDPINCc writes slot $lab_slot."""
    @assert isapprox(report.inc[d, cap_slot], CAP_I[d]; rtol=1e-6) """
        GDPINCCAT slot mismatch: "Capital" (slot $cap_slot) reads \
        $(report.inc[d, cap_slot]) but CAP_I[$d] = $(CAP_I[d]). \
        Check E_delGDPINCb writes slot $cap_slot."""
end
println("GDPINCCAT slot order verified (Labour=$lab_slot, Capital=$cap_slot) ✓")

# ── The two GDP checks measure different things — keep them apart ───────────
# 1) LEVELS (income vs expenditure totals): a property of the calibrated
#    BENCHMARK DATABASE. INDOTERM's database does not satisfy the regional GDP
#    identity exactly (same imbalance that forces the `K` constant in
#    E_pdomA_sum! and shows as DIFFCOM_sc in the pipeline). Reported, not
#    asserted — asserting it would be asserting the input data is perfect.
# 2) CHANGES (TERM's own wgdpdiff = wgdpinc − wgdpexp): the MODEL-side check.
#    Both are bmk=1 ratio indices, so this must be 0 at a zero shock — that IS
#    assertable, and it is the meaningful test under a real shock.
lev_ok, lev_d, lev_worst = gdp_both_sides_check(report)
chg_ok, chg_d, chg_worst = gdp_change_consistency(res.values)

println("\nLevels  income-vs-expenditure : worst rel.diff $(round(lev_worst, sigdigits=6))" *
        (lev_ok ? " ✓" : " at $(REG6[lev_d])  [benchmark-data property, not asserted]"))
println("Changes wgdpdiff (model check) : worst |wgdpdiff| $(round(chg_worst, sigdigits=6))" *
        (chg_ok ? " ✓" : " at $(REG6[chg_d])"))

@assert chg_ok "wgdpdiff = $(chg_worst) at a ZERO shock — the income and expenditure sides of the MODEL disagree"

println("\n" * "="^72)
println("STEP 8 PASSED: GDP reporting works; model-side GDP check consistent ✓")
if !lev_ok
    println("NOTE: the benchmark database's regional GDP identity is off by up to")
    println("      $(round(100*lev_worst, sigdigits=3))% (worst: $(REG6[lev_d])). Weigh this when reading regional results.")
end
println("="^72)
