push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP

println("="^64)
println("INDOTERM BENCHMARK SOLVE @ 25×6 (Phase 1 validation gate)")
println("="^64)

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

# ── Phase 1: collapse 34 provinces → 6 island groups ─────────────────────────
println("--- aggregate_regions! (34→6) ---")
agg6 = aggregate_regions!(agg_r.agg, REG_MAP_34_to_6, 6)
params = prepare_parameters!(agg6)
@assert size(agg6["MAKE"], 3) == 6 "region collapse failed"

println("--- build_model_full! @ 25×6 ---")
t0 = time()
m, vars = build_model_full!(agg6, params)
println("Build: $(round(time()-t0, digits=1))s, vars=$(num_variables(m)), " *
        "cons=$(num_constraints(m; count_variable_in_set_constraints=false))")

println("--- initialize_model! (benchmark, all shocks = 0.0) ---")
bmk = benchmark_levels(params)
initialize_model!(m, vars; bmk_levels=bmk)

# Snapshot the benchmark seed BEFORE solving — benchmark replication means the
# solution returns to these values. (Comparing every `>= 1e-6` variable to 1.0
# would be wrong: bmk_levels deliberately seeds genuine value FLOWS at their
# benchmark flow, not at 1.)
seed = Dict{String,Any}()
for (nm, v) in vars
    seed[nm] = v isa AbstractArray ?
        [JuMP.start_value(v[i]) for i in eachindex(v)] : JuMP.start_value(v)
end

println("--- solve_newton! (square levels system F(x)=0) ---")
t0 = time()
res = solve_newton!(m, vars)
println("Solve: $(round(time()-t0, digits=1))s  status=$(res.status)  " *
        "iters=$(res.iters)  scaled ||F||_inf=$(res.residual)")

@assert res.solved "Benchmark solve did NOT converge at 25×6 — status=$(res.status)"
@assert res.n_free == res.n_eqs "system not square: $(res.n_free) vars vs $(res.n_eqs) eqs"

# ── Benchmark replication: the solution must return to the benchmark seed ────
max_abs, worst_abs = 0.0, ""
max_rel, worst_rel = 0.0, ""
for (nm, v) in vars
    v isa AbstractArray || continue
    sv = seed[nm]; solved = res.values[nm]
    for (k, i) in enumerate(eachindex(v))
        JuMP.is_fixed(v[i]) && continue
        s0 = sv[k]; s0 === nothing && continue
        dev = abs(solved[i] - s0)
        if dev > max_abs; global max_abs, worst_abs = dev, "$nm[$i]"; end
        rel = dev / max(abs(s0), 1.0)
        if rel > max_rel; global max_rel, worst_rel = rel, "$nm[$i]"; end
    end
end

println()
println("Max |Δ| vs benchmark (absolute): $max_abs  at $worst_abs")
println("Max |Δ| vs benchmark (relative): $max_rel  at $worst_rel")
println("Raw ||F||_inf (unscaled, expected ~1e-4 on 1e6-magnitude value")
println("              identities — judge on the scaled residual): $(res.raw_residual)")

@assert res.residual < 1e-7 "scaled residual $(res.residual) too large — not converged"
@assert max_rel < 1e-4 "solution deviates from benchmark by $max_rel (relative) — not replicated"

println()
println("="^64)
println("PHASE 1 GATE PASSED: 25×6 levels model replicates its benchmark ✓")
println("="^64)
