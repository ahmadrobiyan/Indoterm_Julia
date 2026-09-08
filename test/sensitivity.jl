"""
V8 — elasticity sensitivity analysis (VV_PLAN.md §V8; design in V8_ELASTICITY_SENSITIVITY.md).

Every result this project reports is a point estimate from ONE set of elasticities baked into
the 2016 HAR data. This gate asks which conclusions survive that uncertainty.

METHOD (Phase A — GTAP-style, group-level, one-at-a-time).
Six elasticity groups x {x0.5, x1.5}, plus the unperturbed centre point = 13 solves against
`COALPRICE_REFERENCE`, the project's external-validation scenario, so every number here is
directly comparable to the V9 baseline.

WHERE THE PERTURBATION MUST BE APPLIED — this is the whole correctness argument.
An elasticity is consumed TWICE: `prepare_parameters!` uses it to CALIBRATE the CES share and
scale parameters from the benchmark data (`ces_calibrate`), and `build_model!` uses it again
in the equations. Scaling `params[key]` after the fact changes the second use but not the
first, leaving alpha/gamma calibrated for the OLD sigma. The CES can then no longer reproduce
the base year, and `run_model!`'s benchmark gate rejects the point:

    benchmark solve failed under this closure (||F||inf = 20.108 > 1.0e-8)

Observed for P015 x0.5 on 2026-09-08 — the first perturbation ever attempted. Four of the six
groups are calibration-coupled: P015 (ALPHA/GAMMA_ARMINT, prepare_parameters.jl:295), P028
(primary-factor CES, :98), SLAB (labour CES, :87), SCET (CET make nest, :113). Only SMAR and
P018 are stored without feeding any calibration, so those two alone would have survived an
in-place perturbation — which is why a smoke run over one group cannot clear this design.

The fix is to perturb `agg[key]` — the PRE-calibration input — and re-run
`prepare_parameters!`, which is pure (it never writes to `agg`) and costs ~1.8 s against a
~126 s solve. sigma and the share parameters then move together, the benchmark still
replicates, and the sweep measures the elasticity rather than a calibration inconsistency.

THERE IS NO PASS/FAIL ON THE ECONOMICS. Per VV_PLAN.md the deliverable is a table: centre
value, range across the sweep, and whether the SIGN is stable. A conclusion whose sign flips
inside the sweep is reported UNRESOLVED — that is the finding, not a defect. The gate's own
verdict reports only whether the sweep itself is sound: every key present, every perturbation
material and consistently recalibrated, every point solved.

COST. ~2 min per point at 6 regions => ~30 min for the full 13. Set V8_SMOKE=1 for a 2-point
smoke run (one group, one factor) to validate the plumbing first.

Run:  julia --project=. test/sensitivity.jl
      V8_SMOKE=1 julia --project=. test/sensitivity.jl
"""

using Printf

include(joinpath(@__DIR__, "pipeline_cache.jl"))

# ── sweep definition ────────────────────────────────────────────────────────────────────────
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
# The perturbation targets `agg`, so the key must exist THERE. If it were absent,
# prepare_parameters! would fall back to a hardcoded default (zeros, or fill(5.0) for P015)
# and the sweep would vary nothing while reporting false robustness. Fatal, not skipped.
groups  = SMOKE ? GROUPS[1:1]  : GROUPS
factors = SMOKE ? FACTORS[1:1] : FACTORS

println("="^78)
println("V8 — elasticity sensitivity", SMOKE ? "  [SMOKE MODE — plumbing check only]" : "")
println("="^78)

for (key, label) in groups
    haskey(agg6, key) || error(
        "V8: agg has no key \"$key\" — prepare_parameters! would fall back to a hardcoded " *
        "default and this sweep would vary nothing. Fix the aggregation step.")
    v = vec(parent(agg6[key]))
    any(!iszero, v) || error(
        "V8: agg[\"$key\"] is all zeros — scaling it is a no-op, so this group cannot be " *
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

# ── perturb-then-recalibrate ────────────────────────────────────────────────────────────────
function same_value(a, b)
    try
        a === b && return true
        if a isa AbstractArray && b isa AbstractArray
            return size(a) == size(b) && all(a .== b)
        end
        return isequal(a, b)
    catch
        return false
    end
end

"""
    perturbed_params(key, factor) -> (params2, ncoupled)

Scale `agg[key]` by `factor` and re-derive the whole parameter set from it. `agg6` itself is
never mutated: the copy is shallow because `prepare_parameters!` only READS `agg`, and the one
entry that changes is REPLACED by a fresh array rather than scaled in place.

`ncoupled` counts the other `params` entries that moved — direct evidence that recalibration
actually happened. For P015/P028/SLAB/SCET it must be > 0 (the CES alpha/gamma were re-fitted);
for SMAR/P018 it is 0, and that is correct — they feed no calibration.
"""
function perturbed_params(key::AbstractString, factor::Float64)
    agg2 = copy(agg6)
    agg2[key] = vec(parent(agg6[key])) .* factor
    params2 = prepare_parameters!(agg2)
    ncoupled = count(k -> k != key &&
                          !same_value(get(params, k, nothing), get(params2, k, nothing)),
                     collect(keys(params2)))
    return params2, ncoupled
end

# `solved` is true ONLY when the continuation reached t = 1 (run_model!.jl:372). `residual` is
# measured at the last ACCEPTED point, which is a genuine solution of the system at
# `t_reached` — so a branch that stalls halfway reports a small residual anyway. Reporting the
# residual alone therefore reads as success: `P028 x0.5` stalled with residual 3.77e-09, well
# under the 1e-8 tolerance. Report `t_reached`, which the ScenarioResult docstring calls "the
# useful number for diagnosing a fold".
#
# A stalled point returns `nothing` rather than aborting: one 30-minute run should report every
# problem point, not just the first. Completeness is enforced once, at the end.
function solve_point(p, tag; header::Bool = true)
    header && println("── $tag")
    t0 = time()
    res = run_model!(agg6, p, COALPRICE_REFERENCE; h0 = 0.25, tol = 1e-8,
                     gdp = false, verbose = false)
    dt = time() - t0
    if !res.solved
        @printf("   INCOMPLETE after %.0fs — reached t = %.4f of 1.0 in %d step(s), %d rejection(s)\n",
                dt, res.t_reached, res.nsteps, res.nrejects)
        @printf("   residual %.2e there — a genuine solution at t < 1, not a stalled iterate\n",
                res.residual)
        return nothing
    end
    @printf("   solved in %.0fs, residual %.2e\n", dt, res.residual)
    return headline(res)
end

# ── run the sweep ───────────────────────────────────────────────────────────────────────────
results = Dict{String,Dict{String,Float64}}()
incomplete = String[]

centre = solve_point(params, "centre point (unperturbed elasticities)")
centre === nothing && error(
    "V8: the CENTRE point did not reach t = 1, so there is no baseline to compare against " *
    "and nothing in this sweep is interpretable. Investigate before rerunning.")
results["centre"] = centre

# Recalibration counts per point, kept for the provenance block in the report. A count of 0
# means the elasticity is used only inside the equations; a nonzero count names the share/
# scale blocks that had to be re-fitted so that calibration and equations agree (see the
# perturbed_params docstring — perturbing params in place instead would desynchronise them).
coupled = Dict{String,Int}()

for (key, _) in groups, f in factors
    tag = "$key x$f"
    println("── $tag")
    t0 = time()
    p2, ncoupled = perturbed_params(key, f)
    @assert vec(parent(p2[key])) != vec(parent(params[key])) "V8: perturbation '$tag' did not change params[\"$key\"]"
    coupled[tag] = ncoupled
    @printf("   recalibrated in %.1fs — %d coupled parameter block(s) re-fitted\n",
            time() - t0, ncoupled)
    r = solve_point(p2, tag; header = false)
    r === nothing ? push!(incomplete, tag) : (results[tag] = r)
end

# ── report ──────────────────────────────────────────────────────────────────────────────────
points = collect(keys(results))
sweep_pts = filter(!=("centre"), points)

println()
println("="^78)
println("V8 RESULTS — % change vs benchmark, COALPRICE_REFERENCE")
println("="^78)
if !isempty(incomplete)
    @printf("PARTIAL — %d of %d sweep points did not reach t = 1 and are absent from the\n",
            length(incomplete), length(groups) * length(factors))
    println("range below. Treat min/max as a LOWER BOUND on the true spread: the omitted")
    println("points are the ones the model found hardest, not a random subset.")
    println()
end
@printf("%-14s %10s %10s %10s %10s  %-12s %-12s %s\n",
        "metric", "centre", "min", "max", "spread", "min at", "max at", "sign")

unstable = String[]
rows = String[]
for lbl in LABELS
    c = results["centre"][lbl]
    # The extremum must be attributed to the point that attains it: a spread is only
    # interpretable if you know which elasticity moved it. Ties resolve to the first point
    # in sorted order, which is deterministic across runs.
    ilo = argmin([results[p][lbl] for p in sweep_pts])
    ihi = argmax([results[p][lbl] for p in sweep_pts])
    lo, hi = results[sweep_pts[ilo]][lbl], results[sweep_pts[ihi]][lbl]
    lo_at, hi_at = sweep_pts[ilo], sweep_pts[ihi]
    lo < c || (lo, lo_at = c, "centre")
    hi > c || (hi, hi_at = c, "centre")
    stable = (lo > 0 && hi > 0) || (lo < 0 && hi < 0) || (abs(lo) < 5e-2 && abs(hi) < 5e-2)
    stable || push!(unstable, lbl)
    @printf("%-14s %10.3f %10.3f %10.3f %10.3f  %-12s %-12s %s\n", lbl, c, lo, hi, hi - lo,
            lo_at, hi_at, stable ? "stable" : "FLIPS")
    push!(rows, @sprintf("%s,%.6f,%.6f,%.6f,%.6f,%s,%s,%s", lbl, c, lo, hi, hi - lo,
                         lo_at, hi_at, stable ? "stable" : "unstable"))
end

outdir = joinpath(@__DIR__, "..", "analysis", "output")
mkpath(outdir)
open(joinpath(outdir, "v8_sensitivity.csv"), "w") do io
    println(io, "metric,centre,min,max,spread,min_at,max_at,sign_stability")
    foreach(r -> println(io, r), rows)
end

# The point-level matrix is the actual product of an ~85-minute run; the summary above is a
# lossy projection of it. Persist it in long form so that later questions ("which elasticity
# drives this spread?", "is this group inert?") are answerable without re-solving. Incomplete
# points are recorded with their t_reached rather than omitted silently.
open(joinpath(outdir, "v8_sensitivity_points.csv"), "w") do io
    println(io, "point,metric,value,coupled_blocks,status")
    for lbl in LABELS
        @printf(io, "centre,%s,%.6f,0,solved\n", lbl, results["centre"][lbl])
    end
    for p in sort(sweep_pts), lbl in LABELS
        @printf(io, "%s,%s,%.6f,%d,solved\n", p, lbl, results[p][lbl], get(coupled, p, -1))
    end
    for p in incomplete, lbl in LABELS
        @printf(io, "%s,%s,NaN,%d,incomplete\n", p, lbl, get(coupled, p, -1))
    end
end

# An elasticity group whose every point reproduces the centre exactly is INERT in this
# aggregation: it contributes nothing to the reported range, so the range must not be read as
# though it covered that group. The 6-region make matrix is fully diagonal, which makes SCET
# inert by construction (no output mix to transform) — flag it rather than let it pad the
# apparent breadth of the sweep.
inert = String[]
for (key, _) in groups
    pts = filter(p -> startswith(p, key * " "), sweep_pts)
    isempty(pts) && continue
    all(p -> all(lbl -> results[p][lbl] ≈ results["centre"][lbl], LABELS), pts) &&
        push!(inert, key)
end
if !isempty(inert)
    println()
    println("INERT in this aggregation (every point reproduces the centre): ",
            join(inert, ", "))
    println("These groups widen no interval. The sweep effectively covers ",
            length(groups) - length(inert), " of ", length(groups), " elasticity groups.")
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
println("   Table written to analysis/output/v8_sensitivity.csv")
println()

# The gate's own verdict is about the SOUNDNESS OF THE SWEEP, never about the economics: a
# sign flip above is a finding, reported UNRESOLVED, and must not trip this. An incomplete
# sweep is different — its range is a lower bound, so it may not be certified.
if !isempty(incomplete)
    for tag in incomplete
        println("   incomplete point: $tag")
    end
    error("V8: $(length(incomplete)) of $(length(groups) * length(factors)) sweep points did " *
          "not reach t = 1. The reported range understates the true spread, so this sweep " *
          "cannot be certified. Each stalled point printed its t_reached above — that is the " *
          "number to diagnose. This is a continuation/robustness limit at those elasticities, " *
          "not an economic result.")
end

@printf("✅ V8 sweep complete: %d points, all reached t = 1, all perturbations material.\n",
        length(points))
SMOKE && println("   (SMOKE MODE — 2 points only; rerun without V8_SMOKE for the full sweep.)")
