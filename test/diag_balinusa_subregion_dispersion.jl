"""
Diagnostic (not a gate) — 2026-08-02.

V4 Leg 2, even re-scoped to region-level aggregates (option (b)), still fails: 29 of 534
region-level totals disagree in sign, with `xinvi`/`xinv`/`xinv_s` and `xsuppmar` flipping for
region 5 (BaliNusa) specifically. This is a DIFFERENT mechanism from the earlier one.

Earlier diagnosis (VV_PLAN §V4, rounds 1-4): BaliNusa's 34→12 split carried the worst economic
WEIGHT imbalance of any island (3.25x by `xprim`). Rebalancing to the best possible split
(1.23x, in line with — actually better than — every other island) cut material cell-level
flips 299→165 but did not eliminate them, and did not stop the region-level aggregate flips
seen in run12. So weight imbalance cannot be the whole story anymore: BaliNusa's current split
is now the BEST-BALANCED of any island by weight, yet it's still the top offender.

**This script tests a different hypothesis: response DISPERSION, not weight imbalance.**
CES/CET aggregation bias in a nonlinear model depends on how DIFFERENT the sub-units' behavior
is under the shock, not on how equally they're weighted — two equally-sized sub-regions that
respond in opposite directions produce large aggregation bias even at a perfectly balanced
1:1 split; two very unequally-sized sub-regions that respond almost identically produce almost
none. If BaliNusa's two 12-region sub-splits (province 28 alone → sub-region 9; provinces
29+30 → sub-region 10) diverge more sharply in their %deviation-from-benchmark than other
islands' sub-region pairs do, that dispersion — not the weight ratio already ruled out — is
the mechanism, and it would explain why rebalancing by weight alone plateaued at 165/29.

Solves the 12-region model once (shocked, same LEG2 instrument as V4) and reports, per island,
the spread between its two sub-regions' %dev for exactly the variables that flipped:
`xinvi`, `xinv`, `xinv_s`, `xsuppmar`, `xsuppmar_p`. Ranks islands by spread so BaliNusa's
standing (top offender, or not) is a direct readout, not an assumption.
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))

const REG_MAP_34_to_12 = [
    1, 1, 1, 1, 1,  2, 2, 2, 2, 2,
    3, 3, 3,  4, 4, 4,
    5, 5, 5,  6, 6,
    7, 7, 7,  8, 8, 8,
    9,  10, 10,
    11, 11,  12, 12,
]
const ISLAND_NAME = Dict(1=>"Sumatra", 2=>"Java", 3=>"Kalimantan", 4=>"Sulawesi",
                          5=>"BaliNusa", 6=>"MalukuPapua")
const SUBPAIR = Dict(1=>(1,2), 2=>(3,4), 3=>(5,6), 4=>(7,8), 5=>(9,10), 6=>(11,12))

println("="^76)
println("Diagnostic — BaliNusa sub-region response DISPERSION vs other islands")
println("="^76)

agg12, params12 = cached_pipeline(12; rmap = REG_MAP_34_to_12)
bmk12 = benchmark_levels(params12)

const LEG2_SHOCK_PCT = 12.0
scen = Scenario(
    name   = "diag: coalprice.CMF +$(Int(LEG2_SHOCK_PCT))% (same instrument as V4 Leg 2)",
    source = "test/diag_balinusa_subregion_dispersion.jl",
    swaps  = COALPRICE_SWAPS,
    shocks = [("fpexp_d", 5) => logpct(LEG2_SHOCK_PCT)],
    numeraire = :exrate,
)

println("\nSolving 12-region model under the shock (this is the expensive step)...")
r12 = run_model!(agg12, params12, scen; tol = 1e-8, gdp = true)
r12.solved || error("12-region solve failed — cannot diagnose")
println("solved: ‖F‖∞ = $(r12.residual)")

const TARGET_VARS = ["xinvi", "xinv", "xinv_s", "xsuppmar", "xsuppmar_p", "xtrad", "xtradmar"]

"""Per-region %dev-from-benchmark, summed over every non-region axis (same rule as V4's
`regional_totals`), for one variable at 12 regions."""
function region_devs(nm::String, vals::Dict, bmk::Dict, nr::Int)
    haskey(vals, nm) && haskey(bmk, nm) || return nothing
    v, b = vals[nm], bmk[nm]
    (v isa AbstractArray && b isa AbstractArray && size(v) == size(b)) || return nothing
    rax = [d for d in 1:ndims(v) if size(v, d) == nr]
    isempty(rax) && return nothing
    otherax = [d for d in 1:ndims(v) if !(d in rax)]
    vr = isempty(otherax) ? v : dropdims(sum(v; dims = otherax); dims = Tuple(otherax))
    br = isempty(otherax) ? b : dropdims(sum(b; dims = otherax); dims = Tuple(otherax))
    # For 2-region-axis vars (xtrad[r,d] etc.) keep only the diagonal-ish per-origin total,
    # i.e. sum over the destination axis too, so we get one scalar-per-region like the others.
    while ndims(vr) > 1
        vr = dropdims(sum(vr; dims = 2); dims = 2)
        br = dropdims(sum(br; dims = 2); dims = 2)
    end
    dev = similar(vr, Float64)
    for i in eachindex(vr)
        dev[i] = br[i] > 1e-8 ? vr[i] / br[i] - 1.0 : NaN
    end
    dev
end

println("\nPer-island sub-region %dev spread (12-region model, shocked), by variable:")
println("(spread = |dev(subA) - dev(subB)|; large spread ⇒ the two sub-regions respond very")
println(" differently ⇒ merging them under a shock produces large aggregation bias regardless")
println(" of how evenly they're weighted)\n")

results = Vector{NamedTuple}()
for nm in TARGET_VARS
    dev = region_devs(nm, r12.values, bmk12, 12)
    dev === nothing && continue
    println("  $nm:")
    for isl in 1:6
        a, b = SUBPAIR[isl]
        (isfinite(dev[a]) && isfinite(dev[b])) || continue
        spread = abs(dev[a] - dev[b])
        push!(results, (var=nm, island=isl, name=ISLAND_NAME[isl], devA=dev[a], devB=dev[b], spread=spread))
        marker = isl == 5 ? "  ← BaliNusa" : ""
        println("     $(rpad(ISLAND_NAME[isl], 12)) sub9/A=$(round(dev[a]*100,digits=3))%  " *
                "sub10/B=$(round(dev[b]*100,digits=3))%  spread=$(round(spread*100,digits=3))pp$marker")
    end
end

println("\n" * "="^76)
println("RANKED — top 20 largest sub-region spreads across all target variables/islands")
println("="^76)
sort!(results; by = r -> -r.spread)
for r in first(results, min(20, length(results)))
    marker = r.island == 5 ? "  ← BaliNusa" : ""
    println("  $(rpad(r.var,12)) $(rpad(r.name,12)) spread=$(round(r.spread*100,digits=3))pp" *
            "  (subA=$(round(r.devA*100,digits=3))%, subB=$(round(r.devB*100,digits=3))%)$marker")
end

n_bali_top20 = count(r -> r.island == 5, first(results, min(20, length(results))))
println("\nBaliNusa entries in top 20 by spread: $n_bali_top20 / 20")
by_island_avg = Dict(isl => begin
    rs = filter(r -> r.island == isl, results)
    isempty(rs) ? NaN : sum(r.spread for r in rs) / length(rs)
end for isl in 1:6)
println("\nMean spread per island (across all target-variable entries):")
for isl in sort(collect(1:6); by = i -> -by_island_avg[i])
    marker = isl == 5 ? "  ← BaliNusa" : ""
    println("  $(rpad(ISLAND_NAME[isl],12)) $(round(by_island_avg[isl]*100,digits=3))pp$marker")
end
println("\nDone.")
