"""
Gate 7 — end-to-end scenario run at 25 sectors × 6 regions.

This is the first script that runs the model the way a user would: database in,
policy scenario out, with no diagnostic scaffolding. It runs TERM.CMF's own
scenario, `Shock blabnat = -3`, under TERM.CMF's own closure, and prints the GDP
decomposition at the end.

⚠️ `blabnat` is labour-AUGMENTING TECHNICAL CHANGE, not labour supply
(`TERM.TAB:460` `alab_o # Labor-augmenting technical change #`; `:468` `blabnat
# Driver for alab_o #`). `TERM.CMF:112` states the reading outright:
`Shock blabnat = -3; ! 3% increase labour productivity`. So the shock is a 3%
productivity IMPROVEMENT and it is EXPANSIONARY. This file previously described
it as a labour-supply cut and asserted real GDP must FALL — an assertion that
would now fail on a correct model. The sign check below is inverted from what it
was, and the reason is written into it so it cannot be flipped back by accident.

Pass criteria, in order — each one is a separate claim and a failure of an
earlier one makes the later ones meaningless:

  1. The swapped closure reproduces the benchmark (`run_model!` enforces this
     itself and raises if it does not).
  2. The homotopy reaches `t = 1`, i.e. the FULL −3% shock, not a fraction.
  3. `wgdpdiff` — TERM's own income-vs-expenditure check on CHANGES — stays at
     zero. This is the model-side identity and it must hold under a shock, not
     only at the benchmark.
  4. The result is economically sane: a 3% labour-productivity GAIN must raise
     real GDP, and must not move it by an implausible multiple of the shock.

Criterion 4 is deliberately loose. It is a smell test for a sign error or a
runaway, not a validation of the elasticities — the point is that a wrong sign
here would otherwise sail through every residual-based gate, all of which are
satisfied by any solution of the system regardless of which way it moved.

The LEVELS income-vs-expenditure gap is reported but NOT asserted: INDOTERM's
benchmark database does not satisfy the regional GDP identity exactly (up to
~36% in the worst region), so asserting it would be asserting the input data is
perfect. See `test/gdp_report_6reg.jl` for the same distinction at the benchmark.

Run:  julia --project=IndotermJulia IndotermJulia/test/run_scenario_6reg.jl
      SHOCK=-1 julia ...          # smaller shock
      SWAPS=0  julia ...          # bare automatic closure (expect a fold ~0.99975)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP

const SHOCK_PCT = parse(Float64, get(ENV, "SHOCK", "-3"))
const USE_SWAPS = get(ENV, "SWAPS", "1") != "0"

agg6, params = cached_pipeline(6)

r = run_model!(agg6, params;
               name    = "blabnat = $(SHOCK_PCT)%",
               shocks  = ["blabnat" => pct(SHOCK_PCT)],
               swaps   = USE_SWAPS ? TERM_CMF_SWAPS : (),
               # Start small. The branch is linear in the swapped closure out to
               # at least 1e-3, but h doubles after every easy step, so a small
               # h0 costs a handful of extra solves and buys a clean first step.
               h0      = 0.05,
               report  = [("pcap", 2, 5)])

println("\n" * "="^72)
println("SCENARIO SUMMARY")
println("="^72)
println("closure         : ", USE_SWAPS ? "TERM.CMF all six swaps" : "bare automatic")
println("shock           : blabnat = $(SHOCK_PCT)%  (levels $(pct(SHOCK_PCT)))")
println("benchmark ‖F‖∞  : ", r.bench_residual)
println("reached t       : ", r.t_reached, r.solved ? "  ✅ full shock applied" : "  ❌ partial")
println("steps/rejects   : $(r.nsteps) / $(r.nrejects)")
println("final ‖F‖∞      : ", r.residual)

println("\n── homotopy path ──")
println("t              iters   ‖F‖∞")
for (t_, it, res) in r.path
    println(rpad(round(t_; sigdigits=6), 15), rpad(it, 8), round(res; sigdigits=6))
end

@assert r.solved "homotopy stopped at t = $(r.t_reached) — the full shock was never applied"

# ── 3. Model-side GDP identity, under the shock ─────────────────────────────
chg_ok, chg_d, chg_worst = gdp_change_consistency(r.values)
println("\nChanges wgdpdiff (model check): worst |wgdpdiff| = $(round(chg_worst; sigdigits=6))" *
        (chg_ok ? " ✓" : " at $(REG6[chg_d])"))
@assert chg_ok """
    wgdpdiff = $chg_worst at a converged shock: the income and expenditure sides of
    the MODEL disagree. This is not a data-quality issue — it held at the benchmark
    (test/gdp_report_6reg.jl asserts it), so something in the shocked block breaks it."""

print_gdp_report(r.gdp; regions=REG6)

lev_ok, lev_d, lev_worst = gdp_both_sides_check(r.gdp)
println("\nLevels income-vs-expenditure: worst rel.diff $(round(lev_worst; sigdigits=6))" *
        (lev_ok ? " ✓" : " at $(REG6[lev_d])  [benchmark-data property, not asserted]"))

# ── 4. Direction and magnitude smell test ───────────────────────────────────
# A residual-based gate is satisfied by ANY solution of the system, including one
# that moved the wrong way. Real GDP is the one number whose sign we can predict
# from theory alone: more output per worker, more output. A NEGATIVE blabnat is a
# POSITIVE productivity shock (TERM.CMF:112), so the expected signs are opposite.
if haskey(r.values, "NatMacro")
    nm = r.values["NatMacro"]
    k = findfirst(==("RealGDP"), MAINMACROS)
    if k !== nothing
        realgdp = nm[k]                        # bmk = 1 index
        chg = 100 * (realgdp - 1)
        println("\nNational real GDP: $(round(realgdp; sigdigits=8))  ($(round(chg; digits=4))%)")
        println("elasticity to the productivity gain: $(round(chg / -SHOCK_PCT; digits=4))")
        @assert sign(chg) == sign(-SHOCK_PCT) """
            real GDP moved $(chg)% for a $(-SHOCK_PCT)% labour-PRODUCTIVITY gain —
            WRONG SIGN. blabnat is labour-augmenting technical change, so a negative
            shock makes labour more productive and cannot reduce real GDP."""
        @assert abs(chg) <= 3 * abs(SHOCK_PCT) """
            real GDP moved $(chg)% for a $(-SHOCK_PCT)% productivity gain (ratio
            $(chg / -SHOCK_PCT)). Labour's cost share bounds this below 1 with capital
            fixed; a ratio above 3 means something is amplifying, not an elasticity."""
        println("sign and magnitude plausible ✓")
    else
        println("\n(RealGDP not in MAINMACROS — skipping the direction check)")
    end
else
    println("\n(NatMacro not solved back — skipping the direction check)")
end

println("\n" * "="^72)
println("GATE 7 PASSED — full scenario run end to end ✓")
println("="^72)
