"""
Experiment — labour productivity +2%/yr vs +3%/yr, national KPI comparison.

**Status: not executed yet.** Handed off, same convention as `test/sensitivity.jl` — read this
header, then run it. It has NOT been run in this session.

## The question

`TERM_CMF_REFERENCE` (`src/scenarios.jl`) IS the model's own reference labour-productivity
experiment: `blabnat = -3` (a 3% productivity GAIN — `blabnat` is labour-AUGMENTING TECHNICAL
CHANGE, not labour supply, see that scenario's docstring) with `delUnity = 1` (releases the
capital block for one year of accumulation; **load-bearing** — at `delUnity = 0` the branch
folds at λ≈0.9908083, an artifact of the omission, not a model property). This experiment asks
a comparative question the project has not yet run: how do national KPIs move under a smaller,
2%/yr productivity gain, and does the response scale roughly proportionally (2%→3% should move
KPIs by ~1.5×) or does something bend?

## Scenario definitions

`BLABNAT_3PCT` is just `TERM_CMF_REFERENCE` reused as-is (no need to redefine it — redefining a
scenario that already exists is exactly the kind of duplication `scenarios.jl`'s header warns
against). `BLABNAT_2PCT` is the same closure and the same `delUnity=1` companion shock, with only
the `blabnat` magnitude changed 3→2, following `TERM_CMF_REFERENCE`'s exact structure so the two
runs differ in only the one input being tested.

## KPIs

`SELMACROS` (`src/build_macros!.jl`) — this project's own "important national KPIs" shortlist:
`RealGDP, RealGNE, RealHou, RealInv, RealGov, ExpVol, ImpVolUsed, AggEmploy, realwage_io, CPI,
AggCapStock`. Read via `res.values["NatMacro"]` / `MAINMACROS`, the same accessor
`verify_coalprice_reference.jl` and `sensitivity.jl` use.

## What it does NOT do

- Does not touch `src/` or redefine `TERM_CMF_REFERENCE` — only adds one sibling `Scenario`
  locally in this file.
- Does not test regional detail — this is a national-KPI comparison only. If a regional
  breakdown of this experiment is wanted later, run
  `print_regional_confidence_report(res.values, bmk, 6)` on each result first, per the V4
  disclosed-limitation guidance in `VV_PLAN.md` and `AGENTS.md`, before citing any region-level
  sign.
- Ratio column divides the 3% run's % change by the 2% run's — read it as a proportionality
  check, not a formal elasticity estimate; a KPI near zero in the 2% run can produce a large or
  undefined ratio without that meaning much (flagged in the table, not hidden).

## Cost

2 solves × ~4-6 min (6-region `build_model_full!` ~86-178s + homotopy ~120-150s) ≈
**10-15 minutes** total — far cheaper than V8's 13-solve sweep, no rebuild of the data pipeline
needed (`cached_pipeline(6)` is reused).

Run with:
```bash
julia --project=IndotermJulia IndotermJulia/test/experiment_labor_productivity.jl
```
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using Printf

const BLABNAT_3PCT = TERM_CMF_REFERENCE   # already exactly this experiment — reused, not redefined

const BLABNAT_2PCT = Scenario(
    name   = "labour productivity +2%/yr — blabnat = -2%, delUnity = 1",
    source = "derived from origin/TERM.CMF:112-113 (blabnat magnitude changed 3% → 2%)",
    swaps  = TERM_CMF_SWAPS,
    shocks = [
        "blabnat"  => pct(-2),   # 2% productivity GAIN (blabnat is labour-augmenting, not supply)
        "delUnity" => 1.0,       # load-bearing — see TERM_CMF_REFERENCE's docstring
    ],
    notes = """
  Sibling of TERM_CMF_REFERENCE with only the blabnat magnitude changed (3% → 2%). Same closure,
  same delUnity=1 companion shock, same numeraire (:gdppi, the Scenario default — TERM_CMF_REFERENCE
  does not override it). Built for test/experiment_labor_productivity.jl's 2%-vs-3% KPI comparison.
  Quote no number from this scenario without naming the closure, same caveat as the 3% reference.""",
)

function run_kpis(agg, params::Dict{String,Any}, sc::Scenario)
    res = run_model!(agg, params, sc; h0 = 0.25, tol = 1e-8, gdp = false)
    if !res.solved
        println("  ⚠ '$(sc.name)' did NOT reach t=1 (got to t=$(round(res.t_reached; sigdigits=6))) " *
                "— no KPI table below is meaningful, fix the solve first.")
        return nothing
    end
    natmacro = res.values["NatMacro"]
    pctchg(name) = 100 * (natmacro[findfirst(==(name), MAINMACROS)] - 1.0)
    Dict(k => pctchg(k) for k in SELMACROS), res
end

function run_experiment()
    agg6, params = cached_pipeline(6)

    println("="^78)
    println("Experiment — labour productivity +2%/yr vs +3%/yr, national KPI comparison")
    println("="^78)

    println("\n--- run 1: blabnat = -3% (TERM_CMF_REFERENCE) ---")
    kpi3, res3 = run_kpis(agg6, params, BLABNAT_3PCT)
    kpi3 === nothing && error("3% run did not solve — fix before comparing")

    println("\n--- run 2: blabnat = -2% ---")
    kpi2, res2 = run_kpis(agg6, params, BLABNAT_2PCT)
    kpi2 === nothing && error("2% run did not solve — fix before comparing")

    println("\n" * "="^78)
    println("NATIONAL KPI COMPARISON")
    println("="^78)
    @printf("%-14s %10s %10s %10s   %s\n", "KPI", "+2%/yr", "+3%/yr", "ratio", "")
    println("-"^62)
    for k in SELMACROS
        v2, v3 = kpi2[k], kpi3[k]
        note = abs(v2) < 0.02 ? "(2% run ≈0 — ratio not meaningful)" : ""
        ratio = abs(v2) < 0.02 ? NaN : v3 / v2
        @printf("%-14s %10.3f %10.3f %10.3f   %s\n", k, v2, v3, ratio, note)
    end

    println("""

    Read: a ratio near 1.5 means that KPI's response scales roughly proportionally with the
    size of the productivity shock (2%→3% is a 1.5× larger shock). A ratio well away from 1.5
    (and not flagged near-zero above) means the response is disproportionate to shock size for
    that KPI — worth a second look before citing either run's number as "roughly what a policy
    of size X would produce" by linear scaling from the other.
    """)
    (; kpi2, kpi3, res2, res3)
end

# Runs when this file is executed directly, same convention as sensitivity.jl and the other
# test/verify_*.jl gates. This has NOT been run in this session — see the header's ~10-15 min
# cost estimate before invoking it.
run_experiment()
