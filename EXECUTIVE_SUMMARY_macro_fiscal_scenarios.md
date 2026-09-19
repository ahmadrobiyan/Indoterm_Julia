# Executive summary — Macro & fiscal policy scenarios for INDOTERM

**Model:** INDOTERM (TERM-family, ORANI-G-derived multi-region CGE of Indonesia),
Julia port at 25 sectors × 6 island groups. All figures are % deviations of the
national aggregates from the model's own benchmark, computed after the full
scenario is applied (`t = 1`).

---

## Bottom line

Three economy-wide policy scenarios were built and all three solve to a converged
solution in 3 homotopy steps with zero rejections:

| scenario | shock | solved | ‖F‖∞ |
|---|---|---|---|
| `GOVSPEND_EXPANSION` | government demand **+10%** | ✅ | 3.3e-9 |
| `IMPORTPRICE_DOWN` | world import price **−10%** | ✅ | 4.0e-9 |
| `LABPROD_LONGRUN` | labour productivity **+3%** | ✅ | 5.6e-9 |

---

## Results — read against the baseline, not as raw numbers

Every scenario carries `delUnity = 1` (advance one year), which on its own is a
substantial event: capital accumulation lifts the capital stock ~5.6% and real GDP
~2.4% with **no policy at all**. So the raw scenario numbers are mostly baseline.
The right column is what the policy itself contributes.

**Baseline** (`delUnity = 1`, no policy shock): RealGDP **+2.4277%**, CPI +0.1891%,
employment 0.0000%, capital stock **+5.5646%**.

| scenario | RealGDP | CPI | Employ | CapStock | **policy effect on RealGDP** | **on CPI** |
|---|---|---|---|---|---|---|
| baseline | +2.4277% | +0.1891% | 0.0000% | +5.5646% | — | — |
| gov +10% | +2.3245% | −0.4833% | 0.0000% | +5.5646% | **−0.10 pp** | −0.67 pp |
| imp −10% | +2.3286% | −1.6647% | 0.0000% | +5.5646% | **−0.10 pp** | −1.85 pp |
| lab +3% | +4.0597% | +0.1930% | −0.0000% | +5.5646% | **+1.63 pp** | +0.00 pp |

Read as: **a 10% government spending increase and a 10% fall in world import
prices each leave real GDP essentially unchanged (−0.1%), while a 3% labour
productivity gain raises it by 1.6%.**

That is exactly what this closure should produce, and it is the single most
important thing to understand about these results. Under the long-run closure,
national employment is fixed and the capital stock is set by the accumulation
equation, so **real GDP is supply-determined** — demand shocks redistribute and
move prices, not output. The +1.63 pp from the productivity shock is a supply
shock and is consistent with labour's factor share (0.5145 × 3 ≈ 1.54 pp, the
project's own V5 decomposition).

Consequently the interesting action is in **prices and composition**, not GDP:
import-price falls show up almost entirely in the CPI (−1.85 pp), and government
spending shows up in the composition of demand and a milder CPI effect (−0.67 pp).
Capital stock is *identical* in all three scenarios (+5.5646%) because it is
driven entirely by the accumulation step, not by which shock is applied.

### Internal consistency

| scenario | national income-vs-expenditure reldiff | worst regional change-consistency |
|---|---|---|
| gov +10% | 0.86% | 2.80% (MalukuPapua) |
| imp −10% | 2.28% | 1.81% (Java) |
| lab +3% | 0.85% | 0.52% (Kalimantan) |

All three national gaps sit well inside the project's 5% sanity bound, and the
regional dispersion is the same order as the ~3.8% the project documents for its
own reference scenario (V1 — regional income ≠ expenditure is normal
multi-region behaviour; the identity that matters is the national one).

---

## The closure (what is fixed, what is free)

All three scenarios use one closure, `TERM_LR_SWAPS_MATCHED` from
`src/closures.jl`, with the GDP price index as numeraire and `delUnity = 1`:

| swap | effect |
|---|---|
| `xhouhtot = fhou` | regional household consumption follows wage income |
| `houslack = shrBoTnom2` | fixes the nominal balance of trade / GDP |
| `fgovtot = fgovtot3` | regional real government follows regional real GDP |
| `xcap = faccum` | capital accumulation ON — `xcap` becomes endogenous |
| `finv1 = finv4` | dynamic investment rule |
| `flabsup_id = xlab_id` | **fixes national employment**; the real wage adjusts |

**Fixed (exogenous):** capital is *not* fixed — it is mobile; land `xlnd` is fixed;
national employment `xlab_id` is fixed; the trade balance and government-revenue
side are closed by the swaps above; the exchange rate is freed and the GDP price
index is pinned as numeraire.

**Free (endogenous):** `xcap`, `finv1`, `xhouhtot`, `houslack`, `fgovtot`,
`flabsup_id`, `phi`, the real wage, and all remaining prices and quantities.

A full instrument-by-instrument statement, including which `.CMF` each scenario
borrows from and where it deviates, is in
[`MACRO_FISCAL_SCENARIOS.md`](MACRO_FISCAL_SCENARIOS.md).

---

## Key methodological finding: the conventional short-run closure had to be rejected

The government-spending scenario was first built on `COALPRICE_SWAPS` — the
project's short-run comparative static, and the natural home for an economy-wide
demand shock. Under it:

- the **benchmark still solves exactly** (‖F‖∞ = 1.4e-9, "start point already
  satisfies F(x)=0"), so the closure is *valid* and the shock is correctly applied;
- but **every shocked point grinds** — the corrector stalls around ‖F‖∞ ≈ 1e-5
  with the damped-CGNR fallback maxed at 200 iterations per step — at **+10% and
  +2% alike**, and never reaches tolerance.

Because the failure is identical at a fifth of the magnitude, it is a property of
the closure rather than of the shock size, matching the project's own recorded
finding (V6) for economy-wide shocks under rigid-factor closures. All three
scenarios therefore use the long-run closure, which reaches `t = 1` in 3 steps
with 0 rejections. **The short-run variant is unreachable in this port and must
not be reported as an alternative result.**

---

## Caveats

1. **No external validation.** Unlike the coal-price scenario (V9), these three
   have no published GEMPACK result to check against. They are internally
   consistent (identity table above) — the benchmark re-solve gate inside
   `run_model!` enforces that at every run — but they are **not** externally
   validated. Do not present them as calibrated forecasts.
2. **These are long-run, supply-side results.** They say little about short-run
   multipliers, which is where fiscal-expansion questions are usually asked. The
   short-run closure does not currently solve here.
3. **The GDP effects are small by construction, not by finding.** A near-zero real
   GDP response to demand shocks is what a fixed-employment, supply-determined
   closure delivers. Reporting only the GDP column would badly misrepresent these
   scenarios; the CPI and composition results carry the policy signal.
4. **No financing leg.** Government spending is not financed by any tax or
   transfer instrument — the deficit absorbs it. This is not a balanced-budget
   experiment.
5. **Numeraire differs from LONGRUN.CMF** for the import-price scenario
   (`:gdppi` here vs. its exchange-rate numeraire). Real quantities are unaffected;
   the CPI column is only comparable against a run using the same numeraire.
6. **Regional detail is resolution-sensitive.** Use
   `print_regional_confidence_report` before quoting regional *signs* — V4's
   disclosed limitation applies to any scenario.

---

## Deliverables

| file | what |
|---|---|
| `src/scenarios.jl` | the three `Scenario` constants, with full provenance and closure notes |
| `src/IndotermJulia.jl` | exports `GOVSPEND_EXPANSION`, `LABPROD_LONGRUN`, `IMPORTPRICE_DOWN` |
| `MACRO_FISCAL_SCENARIOS.md` | instrument-by-instrument statement of every closure + all deviations from the source `.CMF` files |

## Reproducing

```julia
using IndotermJulia
include("test/pipeline_cache.jl")

agg, params = cached_pipeline(6)
r = run_model!(agg, params, GOVSPEND_EXPANSION; verbose = true)  # describe() prints first
```

`describe(sc)` prints the closure, every shock and the notes before each run, so a
log always records which closure produced a number.
