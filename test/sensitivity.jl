"""
V8 — elasticity sensitivity analysis (VV_PLAN.md §V8; design in V8_ELASTICITY_SENSITIVITY.md).

Every result this project reports is a point estimate from ONE set of elasticities baked into
the 2016 HAR data. Benchmark replication is structurally blind to them: a mistranslated
substitution elasticity changes nothing at the benchmark and everything off it. This gate asks
which conclusions survive that uncertainty.

METHOD (Phase A of the design — GTAP-style, group-level, one-at-a-time).
Six elasticity groups x {x0.5, x1.5}, plus the unperturbed centre point = 13 solves against
`COALPRICE_REFERENCE`, the project's external-validation scenario, so every number here is
directly comparable to the V9 baseline. Each point re-runs `run_model!` with a shallow-copied
`params` in which exactly one elasticity vector is scaled; no `src/` change and no data-pipeline
rebuild is needed (`cached_pipeline(6)` is reused across all points).

THERE IS NO PASS/FAIL ON THE ECONOMICS. Per VV_PLAN.md the deliverable is a table: centre value,
range across the sweep, and whether the SIGN is stable. A conclusion whose sign flips inside the
sweep is reported UNRESOLVED — that is the finding, not a defect. The gate's own ✅/❌ reports
only whether the sweep itself is sound: every key present, every perturbation material, every
point solved.

COST. ~4-6 min per point at 6 regions => ~55-85 min for the full 13. Set V8_SMOKE=1 for a
2-point smoke run (one group, one factor) to validate the plumbing first — do that before
spending the full hour.

Run:  julia --project=. test/sensitivity.jl
      V8_SMOKE=1 julia --project=. test/sensitivity.jl
"""

using Printf

include(joinpath(@__DIR__, "pipeline_cache.jl"))

# ── sweep definition ────────────────────────────────────────────────────────────────────────
# params key => human label. All six are national, commodity-indexed vectors; see
# V8_ELASTICITY_SENSITIVITY.md "Where the elasticities actually live".
const GROUPS = [
    "P015" => "Armington domestic/import CES",
    "P028" => "Primary-factor CES (labour/capital/land)",
    "SLAB" => "Labour-type CES (occupation mix)",
    "SCET" => "CET output mix (multi-product firms)",
    "SMAR" => "Margin-sourcing CES",
    "P018" => "Export demand elasticity",
]

const FACTORS = [0.5, 1.5]
const SMOKE = get(ENV, "V8_SMOKE", "0") == "1"

# Headline national results — the V9 Table 2 columns, same names/order as test/runtests.jl.
const METRICS = [
    ("Real HousCon", "RealHou"),
    ("Real Invest",  "RealInv"),
    ("Export vol",   "ExpVol"),
    ("Import vol",   "ImpVolUsed"),
    ("Real GNE",     "RealGNE"),
    ("Real GDP",     "RealGDP"),
    ("Employment",   "AggEmploy"),
    ("CPI",          "CPI"),
]
const COAL = 5

agg6, params = cached_pipeline(6)

# ── guards: a sweep that silently varies nothing is worse than no sweep ─────────────────────
# `build_model!` reads these keys with a `haskey` fallback to a hardcoded default. If a key were
# absent, scaling it would change NOTHING and the sweep would report false robustness. Likewise
# an all-zero vector (possible for CET) is immune to scaling. Both must be fatal, not skipped.
groups = SMOKE ? GROUPS[1:1] : GROUPS
factors = SMOKE ? FACTORS[1:1] : FACTORS

println("="^78)
println("V8 — elasticity sensitivity", SMOKE ? "  [SMOKE MODE — plumbing check only]" : "")
println("="^78)

for (key, label) in groups
    haskey(params, key) || error(
        "V8: params has no key \"$key\" — build_model! would silently fall back to a " *
        "hardcoded default and this sweep would vary nothing. Fix prepare_parameters.jl.")
    v = params[key]
    v isa AbstractVector || error("V8: params[\"$key\"] is a $(typeof(v)), expected a vector.")
    any(!iszero, v) || error(
        "V8: params[\"$key\"] is all zeros — scaling it is a no-op, so this group cannot be " *
        "swept by a multiplicative factor. Use an additive perturbation for it instead.")
    @printf("  %-6s %-42s n=%3d  min=%+.4f  max=%+.4f\n", key, label, length(v),
            minimum(v), maximum(v))
end
println()

# ── metric extraction (same construction as test/runtests.jl) ───────────────────────────────
# Weights and benchmark come from the UNPERTURBED params on purpose: every sweep point must be
# measured against one fixed yardstick, or the comparison is not like-for-like.
const VTOT0 = parent(params["VTOT"])
const BMK0 = benchmark_levels(params)

function headline(res)
    out = Dict{String,Float64}()
    natmacro = res.values["NatMacro"]
    for (lbl, mac) in METRICS
        out[lbl] = 100 * (natmacro[findfirst(==(mac), MAINMACROS)] - 1.0)
    end
    xtot = res.values["xtot"]
    tw = sum(VTOT0[COAL, d] for d in 1:size(xtot, 2))
    out["Coal output"] =
        100 * sum(VTOT0[COAL, d] * (xtot[COAL, d] - 1) for d in 1:size(xtot, 2)) / tw
    xlab = res.values["xlab_o"]
    bl = BMK0["xlab_o"]
    out["Coal employ"] = 100 * (sum(xlab[COAL, d] for d in 1:size(xlab, 2)) /
                                sum(bl[COAL, d] for d in 1:size(bl, 2)) - 1)
    return out
end

const LABELS = vcat([l for (l, _) in METRICS], ["Coal output", "Coal employ"])

function solve_point(p, tag)
    println("── $tag")
    t0 = time()
    res = run_model!(agg6, p, COALPRICE_REFERENCE; h0 = 0.25, tol = 1e-8,
                     gdp = false, verbose = false)
    dt = time() - t0
    res.solved || error("V8: sweep point '$tag' did not solve (residual $(res.residual)) — " *
                        "the sweep is incomplete and its range is not interpretable.")
    @printf("   solved in %.0fs, residual %.2e\n", dt, res.residual)
    return headline(res)
end

# ── run the sweep ───────────────────────────────────────────────────────────────────────────
results = Dict{String,Dict{String,Float64}}()
results["centre"] = solve_point(params, "centre point (unperturbed elasticities)")

for (key, _) in groups, f in factors
    tag = "$key x$f"
    p2 = copy(params)                       # shallow copy is sufficient: build_model! only reads
    p2[key] = params[key] .* f
    @assert p2[key] != params[key] "V8: perturbation '$tag' did not change params[\"$key\"]"
    results[tag] = solve_point(p2, tag)
end

# ── report ──────────────────────────────────────────────────────────────────────────────────
points = collect(keys(results))
sweep_pts = filter(!=("centre"), points)

println()
println("="^78)
println("V8 RESULTS — % change vs benchmark, COALPRICE_REFERENCE")
println("="^78)
@printf("%-14s %10s %10s %10s %10s   %s\n", "metric", "centre", "min", "max", "spread", "sign")

unstable = String[]
rows = String[]
for lbl in LABELS
    c = results["centre"][lbl]
    vals = [results[p][lbl] for p in sweep_pts]
    lo, hi = minimum(vals), maximum(vals)
    lo, hi = min(lo, c), max(hi, c)
    stable = (lo > 0 && hi > 0) || (lo < 0 && hi < 0) || (abs(lo) < 5e-2 && abs(hi) < 5e-2)
    stable || push!(unstable, lbl)
    @printf("%-14s %10.3f %10.3f %10.3f %10.3f   %s\n", lbl, c, lo, hi, hi - lo,
            stable ? "stable" : "FLIPS")
    push!(rows, @sprintf("%s,%.6f,%.6f,%.6f,%.6f,%s", lbl, c, lo, hi, hi - lo,
                         stable ? "stable" : "unstable"))
end

outdir = joinpath(@__DIR__, "..", "analysis", "output")
mkpath(outdir)
open(joinpath(outdir, "v8_sensitivity.csv"), "w") do io
    println(io, "metric,centre,min,max,spread,sign_stability")
    foreach(r -> println(io, r), rows)
end

println()
if !isempty(unstable)
    println("UNRESOLVED — sign is not stable across the sweep for: ", join(unstable, ", "))
    println("Per VV_PLAN.md §V8 this is a FINDING, not a gate failure: report these results as")
    println("unresolved under elasticity uncertainty rather than as signed conclusions.")
else
    println("Sign is stable across every sweep point for all $(length(LABELS)) headline metrics.")
end

println()
@printf("✅ V8 sweep complete: %d points, all solved, all perturbations material.\n",
        length(points))
SMOKE && println("   (SMOKE MODE — 2 points only; rerun without V8_SMOKE for the full sweep.)")
println("   Table written to analysis/output/v8_sensitivity.csv")
