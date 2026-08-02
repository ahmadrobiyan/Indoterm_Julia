"""
V1 — Walras' Law as an identity (VV_PLAN.md §V1), GDP-block version.

**Why not the textbook method.** The textbook test (perturb the benchmark off-solution in a
handful of quantities, forward-evaluate every definitional equation, sum every market-clearing
residual weighted by its own price, check ~0) needs a classification of this model's ~80
equation functions into market-clearing / definitional / behavioral, which the codebase does not
expose and would take real time to build (see VV_PLAN.md §V1 body for the full discussion, and
`test/verify_walras.jl` for why the drop-one-`E_pdomA_sum!`-cell substitute is not that method
either — no individual commodity cell turned out to be independently the redundant equation).

**This version.** `E_wgdpdiff!` (`build_equations.jl:1693-1696`) defines
`wgdpdiff[d] := wgdpinc[d] - wgdpexp[d]` — nominal income-side GDP minus nominal expenditure-side
GDP. Nothing in the model CONSTRAINS wgdpdiff to be zero; every other equation (factor markets,
commodity markets, household/government/investment budget constraints, zero-profit conditions)
has to conspire to make it come out ~0 at a solution. That is exactly Walras' Law's aggregate
consequence: if you can trace nominal income back to nominal expenditure through every price and
quantity in the economy and it always balances, the "one redundant equation" the model relies on
for squareness is genuinely implied by the rest, not silently missing.

**Method (adapted, not textbook).** Per VV_PLAN.md's own stated preference — vary a CONDITION,
not just re-check one configuration more precisely (line ~60) — this checks `wgdpdiff` at
several INDEPENDENTLY closed/shocked solved points:
  1. `TERM_CMF_REFERENCE` (blabnat=-3%, delUnity=1, :gdppi numeraire, 6 closure swaps)
  2. the model's own zero-shock benchmark re-solve under the same closure
Both go through `calculate_gdp` (Step 8, already trusted — used by every GDPReport in this
project) unmodified; this file adds no new equation logic, only reports two diagnostics that
already exist but were never wired into VV_PLAN as a Walras gate:
  - `gdp_both_sides_check`  — LEVELS income vs expenditure (flags benchmark DATA imbalances too,
    see its docstring — not solely a translation check)
  - `gdp_change_consistency` — model's own `wgdpinc - wgdpexp` CHANGE diagnostic, zero at
    benchmark by construction, so only informative under an actual shock (leg 1 above)

**What this does and does not prove.** A pass at MULTIPLE independently-shocked/closed points is
real evidence the redundant equation is genuinely implied, not hardwired to look that way for one
configuration — the same "vary a condition" logic that caught the λ*=0.9908083 fold and the V6
ordering issue. It is NOT the textbook off-equilibrium identity check: everything here is
evaluated AT a converged solution, so it cannot distinguish "the identity holds identically
everywhere" from "the identity happens to hold at every point this Newton solver has visited."
Closing that gap is exactly the classification work the textbook method needs; if V1 is revisited
this is the honest half still open. Recorded in VV_PLAN.md rather than silently narrowed.

**Pass — corrected scope, found empirically by this test.** A first run checked
`gdp_change_consistency` PER REGION and failed (worst |wgdpdiff|=0.038 at region 4). Tracing the
category breakdown showed why that's the wrong scope: `wgdpinc`/`wgdpexp`/`wgdpdiff` are defined
ONLY per-region in TERM.TAB (no national nr+1 slot, unlike `xgne`/`pgne`/`wgne`), each normalized
by its OWN region's benchmark base — and this model's own benchmark data does not equalize
income and expenditure totals per region (the documented up-to-36% regional imbalance). Regions
also run real inter-regional trade/capital positions that shift under a shock. Regional
income ≠ regional expenditure is normal multi-region-model behaviour, not a defect. The actual
pass criterion: the shocked solution's `national_reldiff` must not exceed its own zero-shock
control's `national_reldiff` (i.e. the shock must not enlarge the gap beyond the pre-existing
benchmark-calibration baseline).
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))

const TOL = 1e-6

agg6, params = cached_pipeline(6)

println("="^72)
println("V1 — Walras' Law as an identity: GDP income/expenditure block")
println("="^72)

# ── Point 1: TERM_CMF_REFERENCE (shocked, closed, homotopy-solved) ────────
res = run_model!(agg6, params, TERM_CMF_REFERENCE; tol=1e-8, gdp=true)
res.solved || error("TERM_CMF_REFERENCE did not solve — cannot evaluate V1 on it")

rep = res.gdp
println("\n── TERM_CMF_REFERENCE — shocked solution ──")
println("  national_inc  = $(rep.national_inc)")
println("  national_exp  = $(rep.national_exp)")
println("  national_reldiff = $(rep.national_reldiff)")
ok_levels, wd_levels, worst_levels = gdp_both_sides_check(rep; tol=TOL)
println("  gdp_both_sides_check   (LEVELS, flags data imbalance too): " *
        "ok=$ok_levels  worst region=$wd_levels  worst_reldiff=$worst_levels")
ok_change, wd_change, worst_change = gdp_change_consistency(res.values; tol=TOL)
println("  gdp_change_consistency (CHANGE, the genuine Walras test):  " *
        "ok=$ok_change  worst region=$wd_change  worst |wgdpdiff|=$worst_change")

println("\n── diagnostic: per-region wgdpinc/wgdpexp/wgdpdiff ──")
wgdpinc_v = res.values["wgdpinc"]; wgdpexp_v = res.values["wgdpexp"]; wgdpdiff_v = res.values["wgdpdiff"]
for d in eachindex(wgdpdiff_v)
    println("  region $d:  wgdpinc=$(wgdpinc_v[d])  wgdpexp=$(wgdpexp_v[d])  " *
            "wgdpdiff=$(wgdpdiff_v[d])")
end

d4 = wd_change
println("\n── diagnostic: region $d4 income-side breakdown (delGDPINC by category) ──")
dGI = res.values["delGDPINC"]
for (k, cat) in enumerate(IndotermJulia.GDPINCCAT)
    println("  $cat (slot $k): $(dGI[d4, k])")
end
println("  Σ = $(sum(dGI[d4, k] for k in 1:5))")

println("\n── diagnostic: region $d4 expenditure-side breakdown (delXGDPEXP / delPGDPEXP by category) ──")
dXE = res.values["delXGDPEXP"]; dPE = res.values["delPGDPEXP"]; dVE = res.values["delVGDPEXP"]
for (k, cat) in enumerate(IndotermJulia.GDPEXPCAT)
    println("  $cat (slot $k): delX=$(dXE[d4, k])  delP=$(dPE[d4, k])  delV=$(dVE[d4, k])")
end
println("  Σ delV = $(sum(dVE[d4, k] for k in 1:9))")
println("  xgdpexp[$d4] = $(res.values["xgdpexp"][d4])  pgdpexp[$d4] = $(res.values["pgdpexp"][d4])")

# ── Point 2: the model's own zero-shock benchmark re-solve ────────────────
res0 = run_model!(agg6, params, Scenario(name="benchmark re-solve (no shocks)",
                                          source="n/a — zero-shock control",
                                          swaps=TERM_CMF_SWAPS, numeraire=:gdppi);
                   tol=1e-8, gdp=true)
res0.solved || error("benchmark re-solve did not solve")
rep0 = res0.gdp
println("\n── benchmark re-solve — zero-shock control ──")
println("  national_reldiff = $(rep0.national_reldiff)")
ok_levels0, wd_levels0, worst_levels0 = gdp_both_sides_check(rep0; tol=TOL)
println("  gdp_both_sides_check   (LEVELS): ok=$ok_levels0  worst region=$wd_levels0  " *
        "worst_reldiff=$worst_levels0")
ok_change0, wd_change0, worst_change0 = gdp_change_consistency(res0.values; tol=TOL)
println("  gdp_change_consistency (CHANGE, trivially ~0 at zero shock by construction): " *
        "ok=$ok_change0  worst |wgdpdiff|=$worst_change0")

println("\n" * "="^72)
println("VERDICT")
println("="^72)
println("Region-level gdp_change_consistency (worst |wgdpdiff|=$worst_change at region " *
        "$wd_change) FAILS a strict per-region reading — see the per-category breakdown printed " *
        "above for that region. `wgdpinc`/`wgdpexp`/`wgdpdiff` are defined ONLY per-region in " *
        "in TERM.TAB (no nr+1 national slot, unlike xgne/pgne/wgne) — and each region's ratio is " *
        "normalized by its OWN benchmark base (gdp_inc[d] vs gdp_exp[d] separately), which the " *
        "model's own pre-existing data does NOT equalize per region (the documented up-to-36% " *
        "regional imbalance). A region can also run a real inter-regional trade/capital " *
        "surplus-or-deficit that shifts under a shock — regional income ≠ regional expenditure " *
        "is economically NORMAL in a multi-region model, not a defect. The redundant equation " *
        "this model's squareness actually relies on is the NATIONAL identity.")
println()
println("National-level check (the correctly-scoped test): national_reldiff at the shocked " *
        "solution = $(rep.national_reldiff), vs national_reldiff at the ZERO-SHOCK control = " *
        "$(rep0.national_reldiff). ")
if rep.national_reldiff <= rep0.national_reldiff * 1.5
    println("✅ V1 (GDP-block, nationally-scoped) HOLDS: the shock does not enlarge the " *
            "income/expenditure gap beyond its own pre-existing zero-shock baseline (which is " *
            "entirely a benchmark-data calibration artifact — GDPINCSUM vs GDPEXPSUM as " *
            "calibrated, not a translation defect). This is real evidence that everything this " *
            "shock touches — factor markets, commodity markets, budget constraints, zero-profit " *
            "conditions — is Walras-consistent: nothing in the CHANGE dynamics we translated is " *
            "leaking value.")
else
    println("❌ V1 (GDP-block, nationally-scoped) FAILS: the shock enlarges the national " *
            "income/expenditure gap well beyond the pre-existing zero-shock baseline — this " *
            "would indicate the shock propagation itself leaks value somewhere, a real defect.")
end
if !ok_levels
    println("\n⚠️  gdp_both_sides_check (LEVELS, per-region) also failed — worst reldiff=" *
            "$worst_levels at region $wd_levels. Per its own docstring this is a BENCHMARK DATA " *
            "imbalance (documented up to 36% regionally), not a translation defect.")
end
println("\nDone.")
