"""
V8 — elasticity sensitivity sweep (Phase A: one-at-a-time, group-level).

**Status: not executed yet.** This file implements the design already written up in
`V8_ELASTICITY_SENSITIVITY.md` (2026-08-01) and is being handed off for someone (human or agent)
to run and interpret — it has NOT been run in this session. Read this whole header before
running it; it explains what to expect and what NOT to conclude from a first run.

## The question

Every result this project has reported (V1's Walras check, V4's aggregation consistency, V9's
coalprice comparison, etc.) is a point estimate from ONE fixed set of point elasticities baked
into the 2016 HAR data. Standard CGE/GTAP practice is to check whether headline conclusions
survive plausible uncertainty in those elasticities. This script asks: does any of this
project's headline finding (from `COALPRICE_REFERENCE`, the model's own external-validation
scenario — see V9 in VV_PLAN.md) flip SIGN when a substitution/CET elasticity group is scaled
±50% from its calibrated value?

## What it does

13 solves of `COALPRICE_REFERENCE` against the same cached 6-region database
(`cached_pipeline(6)`, already built and reused by every other gate — no rebuild needed per
sweep point):

1. **Center point** — unperturbed `params`, i.e. the same run V9 already reports.
2. **12 sweep points** — one at a time, each of the 6 elasticity groups (`SLAB`, `P028`, `P015`,
   `SMAR`, `SCET`, `P018`; see the table in `V8_ELASTICITY_SENSITIVITY.md` for what each one is
   and where it's consumed) scaled by ×0.5 and ×1.5, with the other five held at their
   calibrated value. This is the GTAP-style scalar-multiplier design, not a per-sector sweep —
   6 groups × 25 sectors × 2 directions (300 runs, ~20-30h) was explicitly rejected in the
   design doc as intractable and not the right question.

For each of the 13 runs it reads the same 8 national headline metrics V9 already validates
against `draftreport.pdf` Table 2 (`RealHou`, `RealInv`, `ExpVol`, `ImpVolUsed`, `RealGNE`,
`RealGDP`, `AggEmploy`, `CPI`, via `res.values["NatMacro"]` / `MAINMACROS` — same accessor
`verify_coalprice_reference.jl` uses) and prints a range table: center value, min, max across
the 12 sweep points, and whether the sign is stable across that range.

## What it does NOT do

- It does not touch `src/` — the six elasticity keys are already read out of `params` at
  `build_model!`/`build_model_full!` call time with no other plumbing needed (confirmed in the
  design doc). A sweep point is `copy(params)` plus one key overwrite.
- It is not a pass/fail gate. Per `VV_PLAN.md` §V8's own framing, the deliverable is a table,
  not a verdict — a conclusion whose sign flips within the sweep is reported **unresolved**,
  which is itself the finding this exists to surface, not a bug in the script.
- Phase B (finer grid / Gaussian quadrature around a flagged group) is explicitly NOT
  implemented here — only run it if Phase A below flags something, per the design doc's
  escalation path.

## Cost — read before running

~4-6 minutes wall-clock per solve at 6 regions (`build_model_full!` ~86-178s + homotopy ~120-
150s), so **13 solves ≈ 55-85 minutes** for this script end to end. This is the "high cost, N
solves" reason V8 is sequenced last (Phase 3) in `VV_PLAN.md` — every cheaper gate (V1-V7, V9)
is also a *precondition* for this one meaning anything: a sensitivity sweep on a model that
isn't already Walras-consistent, homogeneous, path-independent, and aggregation-faithful would
just be characterizing the uncertainty of a possibly-broken model.

Run with:
```bash
julia --project=IndotermJulia IndotermJulia/test/sensitivity.jl
```
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using Printf

# ── the 6 elasticity groups (see V8_ELASTICITY_SENSITIVITY.md for the full trace
#    from HAR code → params key → consumption site) ─────────────────────────────
const ELASTICITY_GROUPS = [
    ("SLAB", "Labour-type CES (occupation mix within a sector)"),
    ("P028", "Primary-factor CES (labour/capital/land)"),
    ("P015", "Armington domestic/import CES"),
    ("SMAR", "Margin-sourcing CES across provenance"),
    ("SCET", "CET output-mix (multi-product firms)"),
    ("P018", "Export demand elasticity"),
]

const SWEEP_FACTORS = (0.5, 1.5)

# ── same 8 national headline metrics V9 validates against draftreport.pdf Table 2 ──
const HEADLINE_METRICS = [
    ("RealHou",    "Real HousCon"),
    ("RealInv",    "Real Invest"),
    ("ExpVol",     "Export vol"),
    ("ImpVolUsed", "Import vol"),
    ("RealGNE",    "Real GNE"),
    ("RealGDP",    "Real GDP"),
    ("AggEmploy",  "Employment"),
    ("CPI",        "CPI"),
]

"""
    sweep_params(params, key, factor) -> Dict{String,Any}

One perturbed copy of `params` with `params[key]` scaled by `factor`. Every other entry is
read, never mutated, by `build_model!`/`build_model_full!`, so a shallow `copy` is sufficient
(same reasoning as the design doc's sketch — this is the implementation of it).
"""
function sweep_params(params::Dict{String,Any}, key::String, factor::Float64)
    p2 = copy(params)
    p2[key] = params[key] .* factor
    p2
end

"""
    run_point(agg, params, label; verbose=true) -> Dict{String,Float64}

Run `COALPRICE_REFERENCE` against `params` and extract the 8 headline % changes from
`NatMacro`. Returns `nothing` (with a printed warning, not an exception) if the homotopy did
not reach `t=1` — a partial path is a different experiment, not a scaled-down version of the
sweep point, exactly as `verify_coalprice_reference.jl` treats it.
"""
function run_point(agg, params::Dict{String,Any}, label::AbstractString; verbose::Bool = true)
    res = run_model!(agg, params, COALPRICE_REFERENCE; h0 = 0.25, tol = 1e-8, gdp = false,
                      verbose = verbose)
    if !res.solved
        println("  ⚠ $label did NOT reach t=1 (got to t=$(round(res.t_reached; sigdigits=6))) " *
                "— excluded from the range table, reported separately.")
        flush(stdout)
        return nothing
    end
    natmacro = res.values["NatMacro"]
    pctchg(name) = 100 * (natmacro[findfirst(==(name), MAINMACROS)] - 1.0)
    result = Dict(key => pctchg(key) for (key, _) in HEADLINE_METRICS)
    println("  $label results: ", join(["$k=$(round(result[k]; digits=4))" for (k,_) in HEADLINE_METRICS], ", "))
    flush(stdout)
    result
end

function run_sensitivity_sweep()
    agg6, params = cached_pipeline(6)

    println("="^78)
    println("V8 — elasticity sensitivity sweep, Phase A (one-at-a-time, group-level)")
    println("="^78)
    println("13 solves of COALPRICE_REFERENCE: 1 center point + 6 groups × {×0.5, ×1.5}.")
    println("Expected wall-clock: ~55-85 minutes total.\n")

    results = Dict{String,Dict{String,Float64}}()

    println("--- center point (unperturbed elasticities) ---")
    center = run_point(agg6, params, "center")
    center === nothing && error("center point did not solve — fix the base scenario before " *
                                 "attempting the sweep, a broken center point makes every " *
                                 "sweep point meaningless")
    results["center"] = center

    for (key, desc) in ELASTICITY_GROUPS, factor in SWEEP_FACTORS
        label = "$key × $factor"
        println("\n--- $label  ($desc) ---")
        p2 = sweep_params(params, key, factor)
        r = run_point(agg6, p2, label)
        r !== nothing && (results["$key×$factor"] = r)
    end

    println("\n" * "="^78)
    println("RANGE TABLE — center vs. min/max across all successful sweep points")
    println("="^78)
    sweep_labels = [k for k in keys(results) if k != "center"]
    @printf("%-14s %10s %10s %10s   %s\n", "metric", "center", "min", "max", "sign stable?")
    println("-"^68)
    unresolved = String[]
    for (key, label) in HEADLINE_METRICS
        c = center[key]
        vals = [results[l][key] for l in sweep_labels]
        lo, hi = minimum(vals), maximum(vals)
        # near-zero center values have no meaningful sign to test, same convention V9 uses
        signless = abs(c) < 0.05
        stable = signless || (sign(lo) == sign(hi) == sign(c))
        stable || push!(unresolved, label)
        @printf("%-14s %10.3f %10.3f %10.3f   %s\n", label, c, lo, hi, stable ? "✅" : "❌ UNRESOLVED")
    end

    println("\n" * "="^78)
    println("VERDICT")
    println("="^78)
    if isempty(unresolved)
        println("""
        No headline metric flips sign across the ±50% one-at-a-time sweep. Per VV_PLAN.md §V8
        there is no pass/fail here — this is the "clean" outcome, meaning the reported V9
        conclusions do not appear to be an artifact of one 2016 elasticity choice.""")
    else
        println("""
        UNRESOLVED: $(join(unresolved, ", ")). At least one headline metric's sign is sensitive
        to a ±50% elasticity perturbation. This is not a bug — per the design doc, this is
        exactly the failure mode a point-elasticity result silently hides. Consider Phase B
        (finer grid on the flagged group(s) only) before citing the affected metric's sign.""")
    end
    (; center, results, unresolved)
end

# Runs when this file is executed directly (`julia --project=IndotermJulia
# IndotermJulia/test/sensitivity.jl`), same convention as every other test/verify_*.jl gate.
# This has NOT been run in this session — see the header's cost estimate (~55-85 min) before
# invoking it.
run_sensitivity_sweep()
