# Macro & fiscal policy scenarios — statements of closure

Three constructed policy scenarios for the INDOTERM Julia port, living in
`src/scenarios.jl` as `GOVSPEND_EXPANSION`, `LABPROD_LONGRUN` and
`IMPORTPRICE_DOWN`. This file states, for each one, **what is fixed and what is
free** — because in a CGE model the closure *is* the experiment: the same shock
under two closures is two different experiments, and this project has already
published one wrong number by quoting a result without naming its closure
(`PLAN.md` §N; `closures.jl` header).

None of the three is a transcription of a `.CMF` file. They are constructed, and
each `.CMF` whose lever or magnitude they borrow is named below, together with
every point where they depart from it. There is no external reference result for
any of them — the only externally validated scenario in the project remains
`COALPRICE_REFERENCE` (V9).

---

## The closure: long run (`TERM_LR_SWAPS_MATCHED`)

All three scenarios use one closure, a pre-existing tested constant from
`src/closures.jl`:

| swap | effect |
|---|---|
| `xhouhtot = fhou` | regional household consumption follows wage income |
| `houslack = shrBoTnom2` | fixes the nominal balance of trade / GDP |
| `fgovtot  = fgovtot3` | regional real government follows regional real GDP |
| `xcap     = faccum` | switches ON capital accumulation — `xcap` becomes endogenous |
| `finv1    = finv4` | dynamic investment rule |
| `flabsup_id = xlab_id` | **fixes national employment**; the real wage carries the adjustment |

**Fixed (exogenous).** Land `xlnd` — no closure in the project ever swaps it.
National employment by occupation `xlab_id`. The new exogenous sides of the swaps:
`fhou`, `shrBoTnom2`, `fgovtot3`, `faccum`, `finv4`. Plus the automatic exogenous
set: technical change (`acap`, `alnd`, `atot`, `atrad`, `atradmar`, `blab*`,
`bprim*`, `bint_scd`, `srctwist`), tax shifters (`delPTXRATE`,
`tuser_su`/`_sud`/`_ud`), trade shifters (`fpexp`, `fqexp`, `fpexp_d`, `fqexp_d`,
`pfimp`, `natfpexp`, `natfqexp`), government demand (`fgov`, `fgov_s`, `fgovgen`),
`nhou`, `invslack`, `capslack`, `emptrend`, `delUnity`, `delfwage_o`, `frnorm`,
`frnorm_id`, `gtrend` — minus the ones the swaps moved: `xcap`, `finv1`,
`xhouhtot`, `houslack`, `fgovtot` and `flabsup_id` are no longer fixed.

**Free (endogenous).** `xcap` (capital reallocates across sectors), `finv1`,
`xhouhtot`, `houslack`, `fgovtot`, `flabsup_id`, `phi`, and every remaining price
and quantity in the model.

**Numeraire.** The GDP price index — `swap phi = NatMacro("GDPPI")` applied, so
`phi` is free and `NatMacro("GDPPI")` is pinned at 1.0 (`numeraire = :gdppi`).

**`delUnity = 1` is mandatory with this closure.** Without it the accumulation
equations weld all 150 `xcap` cells back to their benchmark values and the
capital-mobility swaps are inert — the pinned-closure failure that produced the
spurious fold at λ\* = 0.9908083 (`TERM_CMF_REFERENCE`'s notes). All three
scenarios carry it.

---

## The three scenarios

| | `GOVSPEND_EXPANSION` | `IMPORTPRICE_DOWN` | `LABPROD_LONGRUN` |
|---|---|---|---|
| **Policy** | economy-wide government spending **+10%** | world import price **−10%** | labour productivity **+3%** |
| **Instrument** | `fgovgen` | `pfimp`, all 25 commodities | `blabnat` |
| **Target form** | `logpct(10)` = log(1.10) | `pct(-10)` = 0.90 each | `pct(-3)` = 0.97 |
| **Closure** | `TERM_LR_SWAPS_MATCHED` | `TERM_LR_SWAPS_MATCHED` | `TERM_LR_SWAPS_MATCHED` |
| **Numeraire** | `:gdppi` | `:gdppi` | `:gdppi` |
| **Shocks** | 2 (`fgovgen`, `delUnity`) | 26 (25 `pfimp`, `delUnity`) | 2 (`blabnat`, `delUnity`) |
| **Borrows from** | `fund1.cmf` (lever, widened) | `LONGRUN.CMF` (shock) | `TERM.CMF:112` (shock) |

### `GOVSPEND_EXPANSION` — government spending +10%

`fgovgen` is the scalar economy-wide term of `E_xgov`
(`build_equations.jl:684`):

```
xgov(c,s,d) = XGOV0(c,s,d) · exp( fgovtot(d) + fgov(c,s,d) + fgov_s(c,d) + fgovgen )
```

so `logpct(10)` raises real government demand 10% in every commodity, source and
region, leaving its **composition** unchanged — the `fgovtot`/`fgov_s` terms stay
at 0. `fund1.cmf` asks the composition/targeting question with
`fgov_s("construction","SCC")`; this is its economy-wide counterpart.

**`logpct`, not `pct`.** `fgovgen` is declared unbounded (`build_model!.jl:206`) —
an additive, benchmark-0 shifter entering the exponent as a bare summand.
`pct(10)` = 1.10 would apply `exp(1.10) − 1` ≈ **+200%** and would still solve
without complaint.

**No financing leg.** There is no financing swap: government revenue and the
budget deficit absorb the expansion, and no tax or transfer instrument is moved.
Read the result as the demand-side incidence of the spending, not a
balanced-budget exercise. (`fund1.cmf` pairs its spending rise with a consumption
cut — `fgovtot = -11.5`, commented out. Nothing here does.)

### `IMPORTPRICE_DOWN` — world import price −10%

`pfimp` is the commodity-level world price of imports in foreign currency, indexed
by commodity and benchmarked at 1.0 (`build_model!.jl:85`), so the uniform −10% is
`pct(-10)` on each of the 25 aggregated commodities. `pfimp` is indexed by
commodity only — there is no national import-price scalar (unlike exports, which
have `natfpexp`/`natfqexp`) — so a uniform shock must be one target per commodity.
That makes 26 shocks, which also means `run_model!`'s arclength fallback is
unavailable (it traces a single variable).

**Numeraire deviates from `LONGRUN.CMF`.** LONGRUN.CMF leaves `phi` exogenous (an
exchange-rate numeraire). This scenario uses `:gdppi`, for consistency with the
other two macro scenarios. Real quantities are unaffected by that choice — only
the nominal yardstick — but the CPI column is comparable only against a run using
the same numeraire, exactly as `test/reference/draftreport_coalprice.md` records
for `COALPRICE_REFERENCE`.

### `LABPROD_LONGRUN` — labour productivity +3%

`blabnat` drives `alab_o`, labour-**augmenting** technical change
(`TERM.TAB:460`, `:468`), so `-3` is a productivity **gain** and real GDP must
rise. It is not a labour-supply cut.

**How this differs from `TERM_CMF_REFERENCE`.** Same two shocks, same
capital/investment/consumption/government swaps, same GDPPI numeraire — the
entire difference is the labour axis:

- `TERM_CMF_REFERENCE` uses `TERM.CMF:110 delfwage = flabsup_id`, the national
  real-wage adjustment mechanism (a third variant, neither textbook long-run nor
  short-run — see `closures.jl`).
- This scenario uses `flabsup_id = xlab_id`: national employment **fixed**.

The same cell reads differently under the two. Do not quote one against the other
without naming the closure.

---

## The short-run closure was tried and rejected — empirically

`COALPRICE_SWAPS` (the draft report's short-run comparative static: `xcap`/`xlnd`
fixed in place, real wages fixed, national labour-supply mechanism off) is the
natural home for an economy-wide demand shock, and is what the constructed
`HILIRISASI_*` scenarios use. It was the first thing tried for
`GOVSPEND_EXPANSION`. The result:

- **The benchmark still solves exactly** — ‖F‖∞ = 1.4e-9, "start point already
  satisfies F(x)=0". So the closure is valid and the shock variable is genuinely
  exogenous: the problem is not a modelling error.
- **Every shocked point grinds.** The corrector stalls around ‖F‖∞ ≈ 1e-5 with the
  damped-CGNR fallback maxed at 200 iterations per step, at both **+10%** and
  **+2%** government demand, and never reaches tol. Each Newton iteration costs
  seconds, so the branch is not reachable in practical time.

Because the failure is identical at +2% and +10%, it is a property of the closure,
not of the shock size — the same "grinding" mode the project records for
economy-wide shocks under a rigid-factor closure (V6's `TERM_SR_SWAPS` note). The
long-run closure reaches t = 1 in 3 steps with 0 rejections. **Report the
short-run variant as unreachable in this port, not as an alternative result.**

---

## Deviations from the source `.CMF` files — read before comparing

1. **`LONGRUN.CMF`'s closure is not used, and cannot be.** Its long-run block opens
   with `swap flabsup = flab_io;`. There is **no variable `flabsup`** in this
   model: `TERM.TAB:1983–1996` (Excerpt 38) declares `flabsupA`, `flabsupB` and
   `flabsup_id`, and the plain name is a pre-rename leftover. The port declares no
   such variable, so `apply_swaps!` would raise on an undeclared name.
   `IMPORTPRICE_DOWN` therefore keeps LONGRUN.CMF's **shock** but uses the port's
   own long-run closure.

   LONGRUN.CMF also states `NatMacro("RealInv") = invslack` twice (lines 39 and 52,
   the same pair both ways), which cancels — so its net long-run closure does *not*
   fix national investment. Worth knowing before treating that file as a
   specification.

2. **All three numeraire choices are `:gdppi`**, which for `IMPORTPRICE_DOWN`
   deviates from LONGRUN.CMF's exchange-rate numeraire (see above).

3. **`GOVSPEND_EXPANSION` has no financing leg** (see above).

4. **No external validation.** Unlike `COALPRICE_REFERENCE` (V9), these have no
   published GEMPACK result to check against. They are internally consistent — the
   benchmark re-solve gate inside `run_model!` enforces that — but not externally
   validated.

---

## Running them

```julia
using IndotermJulia
include("test/pipeline_cache.jl")

agg, params = cached_pipeline(6)
r = run_model!(agg, params, GOVSPEND_EXPANSION; verbose = true)  # describe() prints first
```

`describe(sc)` prints the closure, every shock and the notes before each run, so a
log always records which closure produced a number. For region-level results, run
`print_regional_confidence_report` first — V4's disclosed limitation applies to any
scenario (near-zero BaliNusa cells are resolution-sensitive).
