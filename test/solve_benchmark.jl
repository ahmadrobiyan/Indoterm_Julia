push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP

println("--- Running data pipeline (Steps 0-4) ---")
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

println("--- build_model_full! ---")
t0 = time()
m, vars = build_model_full!(agg_r.agg, params)
println("Build: $(round(time()-t0, digits=1))s, vars=$(num_variables(m)), cons=$(num_constraints(m; count_variable_in_set_constraints=false))")

println("--- initialize_model! (benchmark, all shocks = 0.0) ---")
initialize_model!(m, vars)

println("--- solve_model! (Ipopt feasibility solve) ---")
t0 = time()
result = solve_model!(m)
println("Solve: $(round(time()-t0, digits=1))s")
println("Status: $(result.status), solved=$(result.solved), primal_status=$(result.primal_status)")

if result.solved
    # Levels-model benchmark check (see PLAN.md): unlike the old %-change
    # formulation, "close to 0" is only the right benchmark value for the
    # shifter/`(change)`-type free variables in ZERO_START_NAMES
    # (initialize_model!.jl). Ratio-type variables (declared `>= 1e-6`, e.g.
    # prices, a-parameters, w-shares) should replicate at ~1. Genuine
    # quantity/value variables (declared `>= 0`, or unbounded ones like
    # realwage/gret/ggro/pfexp/xstocks) have no generic ground truth wired in
    # here to compare against, so they're only summarized, not asserted.
    zero_start = Set(ZERO_START_NAMES)

    max_zero_dev, worst_zero = 0.0, ""
    max_ratio_dev, worst_ratio = 0.0, ""
    unchecked_max, unchecked_worst = 0.0, ""

    for (nm, v) in vars
        for vr in (v isa AbstractArray ? v : (v,))
            is_fixed(vr) && continue
            val = value(vr)
            if nm in zero_start
                dev = abs(val)
                if dev > max_zero_dev
                    global max_zero_dev, worst_zero = dev, name(vr)
                end
            elseif has_lower_bound(vr) && isapprox(lower_bound(vr), 1e-6; atol=1e-9)
                dev = abs(val - 1.0)
                if dev > max_ratio_dev
                    global max_ratio_dev, worst_ratio = dev, name(vr)
                end
            else
                dev = abs(val)
                if dev > unchecked_max
                    global unchecked_max, unchecked_worst = dev, name(vr)
                end
            end
        end
    end

    println("Max |value| among shifter/(change) free variables (should be ~0): $max_zero_dev (at $worst_zero)")
    println("Max |value - 1| among ratio-type free variables (should be ~0): $max_ratio_dev (at $worst_ratio)")
    println("Max |value| among unchecked (genuine quantity/value, no ground truth wired in) free variables: $unchecked_max (at $unchecked_worst)")
end
