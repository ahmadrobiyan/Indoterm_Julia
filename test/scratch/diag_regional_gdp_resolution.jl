"""
Diagnostic (not a gate) — 2026-08-02.

V4's Leg 2 (region aggregates, option (b)) tested *flow components* (xinvi/xinv/xinv_s,
xsuppmar, xtrad, ...) for 6-vs-12-region sign/magnitude agreement under the coalprice.CMF
shock, and found 29/534 region-level aggregate sign disagreements, concentrated in BaliNusa.
GDP itself was never directly tested that way — only inferred to be "probably fine nationally,
possibly noisy regionally" because some of its components are noisy.

This script closes that gap directly: run the SAME coalprice.CMF shock at 6 regions and at
12-regions-then-aggregated-to-6, compute regional GDP (income side, `calculate_gdp`) at both,
and compare %deviation-from-benchmark by sign and magnitude per region — the same comparison
V4 makes, but on GDP itself rather than a proxy flow variable. Motivated by a user research
question: is a claim like "coal export ban moves Sumatra/Kalimantan GDP by X pp" supported?

Not a gate: does not fail the build, prints a comparison table and a verdict per region.
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))

const REG_MAP_34_to_12 = [
    1, 1, 1, 1, 1,  2, 2, 2, 2, 2,
    3, 3, 3,  4, 4, 4,
    5, 5, 5,  6, 6,
    7, 7, 7,  8, 8, 8,
    9,  10, 10,
    11, 11,  12, 12,
]
const MAP_12_to_6 = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]
const ISLAND_NAME6 = Dict(1 => "Sumatra", 2 => "Java", 3 => "Kalimantan",
                           4 => "Sulawesi", 5 => "BaliNusa", 6 => "MalukuPapua")

# Nesting sanity check (cheap, must hold before anything else is meaningful).
@assert [MAP_12_to_6[REG_MAP_34_to_12[p]] for p in 1:34] == REG_MAP_34_to_6

println("="^76)
println("Diagnostic — regional GDP, 6-region direct vs 12-region→6 aggregated")
println("(same coalprice.CMF instrument as V4 Leg 2 / round-5 dispersion diagnostic)")
println("="^76)

const LEG2_SHOCK_PCT = 12.0
mkscen() = Scenario(
    name   = "diag: coalprice.CMF +$(Int(LEG2_SHOCK_PCT))% (regional GDP check)",
    source = "test/diag_regional_gdp_resolution.jl",
    swaps  = COALPRICE_SWAPS,
    shocks = [("fpexp_d", 5) => logpct(LEG2_SHOCK_PCT)],
    numeraire = :exrate,
)

println("\n── 6-region solve ──")
agg6, params6 = cached_pipeline(6)
r6 = run_model!(agg6, params6, mkscen(); tol = 1e-8, gdp = true)
r6.solved || error("6-region solve failed")
gdp6 = calculate_gdp(r6.values, params6)

println("\n── 12-region solve ──")
agg12, params12 = cached_pipeline(12; rmap = REG_MAP_34_to_12)
r12 = run_model!(agg12, params12, mkscen(); tol = 1e-8, gdp = true)
r12.solved || error("12-region solve failed")
gdp12 = calculate_gdp(r12.values, params12)

# Benchmark GDP totals (income side) at 6 regions — the denominator for %dev.
GDPINCSUM6 = parent(params6["GDPINCSUM"])
bmk_inc6 = [sum(GDPINCSUM6[d, :]) for d in 1:6]

# Aggregate the 12-region solved GDP (income side) down to 6 by summing merged pairs —
# GDP totals are additive levels, so this is exact, unlike a ratio/index variable.
agg12to6_inc = zeros(6)
for r12i in 1:12
    agg12to6_inc[MAP_12_to_6[r12i]] += gdp12.inc_total[r12i]
end

dev6  = (gdp6.inc_total .- bmk_inc6) ./ bmk_inc6
dev12 = (agg12to6_inc   .- bmk_inc6) ./ bmk_inc6

println("\n" * "="^76)
println("Regional GDP (income side) — %dev from benchmark, 6-direct vs 12→6-aggregated")
println("="^76)
println(rpad("region", 14), lpad("6-region dev", 16), lpad("12→6 dev", 16),
        lpad("Δ (pp)", 12), lpad("sign agree?", 14))
signs_agree_vec = [sign(dev6[d]) == sign(dev12[d]) || abs(dev6[d]) < 1e-4 || abs(dev12[d]) < 1e-4
                   for d in 1:6]
for d in 1:6
    marker = signs_agree_vec[d] ? "✅" : "❌"
    println(rpad(ISLAND_NAME6[d], 14),
            lpad(string(round(dev6[d]*100, digits=3), "%"), 16),
            lpad(string(round(dev12[d]*100, digits=3), "%"), 16),
            lpad(string(round((dev12[d]-dev6[d])*100, digits=3)), 12),
            lpad(marker, 14))
end
nflip = count(!, signs_agree_vec)

println("\nNational GDP: 6-region dev = $(round(sum(gdp6.inc_total .- bmk_inc6)/sum(bmk_inc6)*100, digits=4))%" *
        "   12→6 dev = $(round(sum(agg12to6_inc .- bmk_inc6)/sum(bmk_inc6)*100, digits=4))%")

println("\n" * "="^76)
println("GDP both-sides check (income vs expenditure), for context, both resolutions")
println("="^76)
ok6, wd6, worst6   = gdp_both_sides_check(gdp6)
ok12, wd12, worst12 = gdp_both_sides_check(gdp12)
println("6-region:  worst rel.diff = $(round(worst6, sigdigits=4)) at region $wd6" *
        (ok6 ? " ✅" : " ⚠"))
println("12-region: worst rel.diff = $(round(worst12, sigdigits=4)) at region $wd12" *
        (ok12 ? " ✅" : " ⚠"))

println("\n" * "="^76)
println("VERDICT")
println("="^76)
if nflip == 0
    println("✅ Regional GDP sign agrees at both resolutions for all 6 regions.")
    println("   Directional GDP claims (which way a region moves) are NOT undermined by")
    println("   the resolution-sensitivity V4 found in investment/trade-margin components —")
    println("   that noise evidently washes out once summed into GDP.")
else
    println("❌ $nflip region(s) disagree in GDP sign between 6-region and 12→6-aggregated.")
    println("   Directional GDP claims for the flagged region(s) are NOT resolution-robust —")
    println("   same caveat V4 raised for their component flows now applies to GDP itself.")
end
println("\nSumatra (region 1):   dev6=$(round(dev6[1]*100,digits=3))%  dev12=$(round(dev12[1]*100,digits=3))%  " *
        (sign(dev6[1]) == sign(dev12[1]) ? "sign-robust" : "SIGN-SENSITIVE"))
println("Kalimantan (region 3): dev6=$(round(dev6[3]*100,digits=3))%  dev12=$(round(dev12[3]*100,digits=3))%  " *
        (sign(dev6[3]) == sign(dev12[3]) ? "sign-robust" : "SIGN-SENSITIVE"))
println("\nDone.")
