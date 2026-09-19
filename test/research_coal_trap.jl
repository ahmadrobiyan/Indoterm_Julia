include(joinpath(@__DIR__, "pipeline_cache.jl"))
using Printf

# Research: The Coal Trap — asymmetric regional pain from a world coal price drop.
#
# Vehicle: 6-region locked validation model, COALPRICE closure + numeraire
# (COALPRICE_SWAPS, :exrate) — the only closure with an external reference
# result (V9). Two negative shocks mirror COALPRICE_REFERENCE (+50%) in the
# opposite direction: fpexp_d[coal] = logpct(-20) and logpct(-50).
#
# fpexp_d is a LOG-space additive shifter (benchmark 0), so the shock target is
# log(1+p/100), NOT p itself. Coal = sector 5 (aggregation_data.jl).
#
# Outputs per region:
# - RealGDP  via MainMacro (bmk=1 index)  -> 100*(v-1)
# - AggEmploy via MainMacro (bmk=1 index)  -> 100*(v-1)
# - Coal employment via xlab_o[coal,d] (flow, bmk=LAB_O) -> 100*(v/b-1)
#
# National row via NatMacro. Regional signs screened with
# print_regional_confidence_report (V4 limitation).
# Run: julia --project=. test/research_coal_trap.jl


agg6, params = cached_pipeline(6)

COAL = findfirst(==("Coal"), AGGCOM)
@assert COAL == 5 "AGGCOM order changed — coal is no longer sector 5"
NR = length(REG6)

kGDP = findfirst(==("RealGDP"), MAINMACROS)
kEMP = findfirst(==("AggEmploy"), MAINMACROS)

scenarios = [
    ("Coal_Drop_20", -20.0),
    ("Coal_Drop_50", -50.0),
]

# results[scenario][region] = (gdp, emp_total, emp_coal)
results = Dict{String,Any}()

for (name, p) in scenarios
    println("\n" * "="^72)
    println("Running scenario: $name  (world coal price $(p)%)")
    println("="^72)
    flush(stdout)
    sc = Scenario(
        name   = "research — world coal export price $(Int(p))%",
        source = "constructed; mirrors COALPRICE_REFERENCE (src/scenarios.jl) with logpct($(Int(p)))",
        swaps  = COALPRICE_SWAPS,
        shocks = [("fpexp_d", COAL) => logpct(p)],
        numeraire = :exrate,
        notes = "Comparative static like COALPRICE_REFERENCE (no blabnat, no delUnity). "
              * "fpexp_d is LOG-additive — shock is log(1+p/100).",
    )
    try
        r = run_model!(agg6, params, sc; h0=0.25, tol=1e-8, gdp=false, verbose=true)
        if !r.solved
            println("❌ $name did NOT reach t=1 (t=$(round(r.t_reached; sigdigits=4)), "
                  * "‖F‖∞=$(r.residual)). No numbers reported — partial homotopy "
                  * "is a different experiment, not a scaled-down shock.")
            results[name] = nothing
            continue
        end
        @printf("solved: %d steps, %d rejects, ‖F‖∞ = %.3g\n", r.nsteps, r.nrejects, r.residual)

        # V4 screen before any regional sign is read
        print_regional_confidence_report(r.values, benchmark_levels(params), NR)

        MM = r.values["MainMacro"]   # nmac x (nr+1)
        NM = r.values["NatMacro"]    # nmac
        XL = r.values["xlab_o"]      # na x nr (flow)
        bXL = benchmark_levels(params)["xlab_o"]

        rows = Dict{String,Any}()
        for d in 1:NR
            gdp = 100 * (MM[kGDP, d] - 1)
            emptot = 100 * (MM[kEMP, d] - 1)
            base = bXL[COAL, d]
            empcoal = base > 1e-10 ? 100 * (XL[COAL, d] / base - 1) : NaN
            rows[REG6[d]] = (gdp=gdp, emptot=emptot, empcoal=empcoal)
        end
        rows["NATIONAL"] = (gdp=100*(NM[kGDP]-1), emptot=100*(NM[kEMP]-1),
            empcoal=let tot=sum(XL[COAL, d] for d in 1:NR),
                        b=sum(bXL[COAL, d] for d in 1:NR);
                        b > 1e-10 ? 100*(tot/b-1) : NaN end)
        results[name] = rows
        println("✓ $name converged.")
    catch e
        println("❌ $name failed: $e")
        showerror(stdout, e, catch_backtrace())
        println()
        results[name] = nothing
    end
    flush(stdout)
end

println("\n" * "="^80)
println("RESEARCH REPORT: THE COAL TRAP (6-REGION IMPACT, % change vs benchmark)")
println("="^80)
@printf("%-12s | %-12s | %12s | %12s | %12s\n", "Region", "Scenario", "Δ RealGDP %", "Δ TotEmp %", "Δ CoalEmp %")
println("-"^80)
for reg in vcat(REG6, ["NATIONAL"])
    for (name, _) in scenarios
        rows = get(results, name, nothing)
        if rows === nothing || !haskey(rows, reg)
            @printf("%-12s | %-12s | %12s | %12s | %12s\n", reg, name, "n/a", "n/a", "n/a")
        else
            v = rows[reg]
            ec = isnan(v.empcoal) ? "zero base" : @sprintf("%12.3f", v.empcoal)
            @printf("%-12s | %-12s | %12.3f | %12.3f | %s\n", reg, name, v.gdp, v.emptot, ec)
        end
    end
    println("-"^80)
end
println("Done.")
