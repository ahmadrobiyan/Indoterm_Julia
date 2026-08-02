# V8 — Elasticity sensitivity analysis: technical design

**Status: 🔴 not started.** This is a design document only — no code has been written. VV_PLAN.md
§V8 (lines 563-580) states the question and the one-line method; this file works out how to
actually implement it against this codebase, which parameters are the targets, where they live,
what a sweep costs, and what the deliverable looks like.

## The question, restated precisely

Every number this project has produced or validated (V1's Walras check, V2's homogeneity, V3's
path independence, V4's aggregation consistency, V9's coalprice comparison) is a point estimate
from ONE set of point elasticities — values baked into the 2016 HAR data and never revisited.
Standard CGE/GTAP practice is systematic sensitivity analysis precisely because a substitution
elasticity is usually the least-certain number in the model. V8 asks: which of this project's
conclusions are robust to that uncertainty, and which are artifacts of one elasticity choice?

## Where the elasticities actually live

Traced from `src/build_reg0!.jl` (source) through `src/build_premod!.jl` (packaging) to
`src/prepare_parameters.jl` (final `params` dict) to `src/build_model!.jl` (consumption):

| HAR code | `params` key | Meaning | Indexed by | Consumed at |
|---|---|---|---|---|
| `SIGMALAB` | `"SLAB"` | Labour-type CES (occupation mix within a sector) | commodity (`na=25`) | `build_model!.jl:301,406` → `E_xlab!` |
| `SIGMAPRIM` | `"P028"` | Primary-factor CES (labour/capital/land) | commodity (`na=25`) | `build_model!.jl:302,407` → `E_xlab_o!`, `E_pcap!`, `E_plnd!` |
| `1ARM`/`ARMSIGMA` | `"P015"` | Armington domestic/import CES | commodity (`na=25`) | `build_model!.jl:300,405` → `E_xint!`, `E_xhou!`, `E_xinv!` |
| `SIGMAMAR` | `"SMAR"` | Margin-sourcing CES across provenance | margin commodity (`nm`) | `build_equations.jl:919` directly |
| `SIGMAOUT` | `"SCET"` | CET output-mix (multi-product firms) | commodity (`na=25`) | `build_equations.jl:961` directly, and `build_model!.jl` |
| `EXP_ELAST` | `"P018"` | Export demand elasticity | commodity (`na=25`) | `build_model!.jl:304,500` → `E_xexpd!` |

Two things this confirms and that the sweep design depends on:

1. **All six are national-level, commodity-indexed vectors** (`SMAR` by margin commodity) — they
   are set once in `build_reg0!.jl`/`build_premod!.jl`, survive region aggregation untouched
   (`aggregate_model!.jl:140,228` just carries them through), and land in `params` as plain
   `Vector{Float64}` of length `na` (or `nm` for `SMAR`). There is no region axis to worry about.
2. **They are read out of `params` at `build_model!.jl` call time**, with a `haskey` fallback to
   a hardcoded default (`fill(0.5, na)` etc. — those defaults exist for missing-data robustness,
   not because they're meant to be used). This means a sweep does **not** need to touch the
   Steps 0-4 data pipeline (`build_reg0!`→`aggregate_model!`, ~250s, already cached by
   `cached_pipeline`) — it only needs to hand `build_model!`/`build_model_full!` a `params` dict
   with one of these six vectors scaled, everything else identical.

## Cost model, from measurements already on record

From `test/scale_probe.jl` and the V1/V3 run logs (`build_model_full!` timings) at 6 regions:
- `build_model_full!`: ~86-178s per call (varies with cache warmth).
- Homotopy solve to `t=1`: ~120-150s for a well-behaved scenario (3-4 steps × 40-65s), longer
  if the path needs more steps.
- **≈4-6 minutes wall-clock per sweep point** at 6 regions. This is why VV_PLAN tags V8 "high
  cost, N solves" and sequences it last (Phase 3).

A naive design — one-at-a-time on all 6 elasticity groups, per-sector (25 sectors), ±50%, both
directions — is 6 × 25 × 2 = 300 runs ≈ 20-30 hours. **Not the right design.** GTAP-style
systematic sensitivity varies elasticities as a **uniform scalar multiplier per group**, not
per-sector, precisely to keep the sweep tractable while still asking the right question ("is
the CONCLUSION sensitive to how confident we are in Armington substitution overall").

## Recommended design

**Phase A — one-at-a-time, group-level (the method VV_PLAN already specifies).**
6 elasticity groups × {×0.5, ×1.5} = **12 runs**, each against the same shocked scenario
(`COALPRICE_REFERENCE` — already the project's external-validation scenario, so results are
directly comparable to the V9 baseline) plus the unperturbed run as the center point = **13
solves ≈ 55-85 minutes**. Feasible in one session, same order of cost as V4's Stage 3.

```julia
# sketch — not implemented
function sweep_group(agg, params, key::String, factor::Float64, scenario::Scenario)
    p2 = copy(params); p2[key] = params[key] .* factor
    run_model!(agg, p2, scenario; tol=1e-8, gdp=true)
end
```
`params` is a flat `Dict{String,Any}`; a shallow `copy` plus one key overwrite is sufficient
since every other entry is read, never mutated, by `build_model!`/`build_model_full!`.

**Phase B (only if Phase A shows a fragile conclusion) — refine around the flagged group.**
Finer grid (e.g. ±10/25/50%) or Gaussian quadrature on just that elasticity, not all six.
This is the "upgrade if results warrant" escalation VV_PLAN already allows for.

## What "pass" means here — there is no pass/fail

Per VV_PLAN.md's own framing (line 577-579): **the deliverable is a table, not a verdict.**
For each headline result already reported elsewhere in this project (e.g. V9's national Real
GDP −0.09%, Employment −0.27%, or V1's national income/expenditure reldiff), record:
- the center-point value (unperturbed elasticities — already computed, e.g. in V9/V1 logs),
- the range across all 12 sweep points,
- **whether the sign is stable across the whole range.**

A conclusion whose sign flips within the sweep is reported as **unresolved**, not wrong — that
is itself the finding, and exactly the failure mode point-elasticity results can silently hide.

## Files (once implemented)

- `test/sensitivity.jl` — the gate: builds the 12+1 param variants, runs each against
  `COALPRICE_REFERENCE` (reusing `cached_pipeline(6)` — no rebuild needed per point), collects
  the headline metrics above, prints the range table.
- No changes needed to `src/` — the six keys are already read from `params` with no other
  plumbing required, which is exactly what makes this cheap to implement once reached.

## Why this is sequenced last

Every other gate (V1-V7, V9) is cheaper AND is a precondition for V8 meaning anything: if the
model isn't Walras-consistent (V1), homogeneous (V2), path-independent (V3), or aggregation-
faithful (V4), a sensitivity sweep would be characterizing the uncertainty of a possibly-broken
model. V8 answers "how much should we trust this number", which only matters once the prior
gates have established "this number is what the model actually says."
