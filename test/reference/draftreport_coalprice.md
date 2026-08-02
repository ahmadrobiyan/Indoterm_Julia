# Reference target for V9 — the authors' coal-price simulation

Extracted 2026-07-30 from `origin/draftreport.pdf`, "Updating the IndoTERM database"
(Horridge & Roos, CoPS Victoria University; Yusuf & Suryana, Universitas Padjadjaran).

This is the **only external reference target the project has**. Everything else in `test/`
checks the model against itself. See `VV_PLAN.md` V9.

## The simulation

- **Shock:** `shock fpexp_d("coal") = 50` — a 50% rise in the export price of Indonesian
  coal. Source: `origin/coalprice.CMF:116`.
- **Closure:** short-run comparative static. `coalprice.CMF` active swaps:
  | line | swap | meaning |
  |---|---|---|
  | `:79` | `xhouhtot = fhou` | regional consumption follows wage income |
  | `:83` | `houslack = natfhou` | fix the national propensity to consume |
  | `:108` | `flab_i = flabsupA` | switch off regional wage differentials |
  | `:109` | `labslack = flabsup_id` | switch off the national labour supply mechanism |
- `xcap`, `xlnd`, government demand and real wages stay exogenous. `delUnity` is **not**
  shocked (`:114` is commented out) — no capital accumulation, no year advanced.
- **Numeraire:** `phi` (the exchange rate) stays exogenous; the `swap phi = NatMacro("GDPPI")`
  at `:72` is commented out. ⚠️ This differs from `initialize_model!.jl:207-218`, which applies
  that swap unconditionally. Real quantities are unaffected; **the CPI column below is not
  comparable until the numeraire matches.**
- **Solution method:** `method = euler; steps = 8` — an 8-step Euler linearisation. This port
  solves the levels system exactly. The two are not expected to agree to more than a few
  significant figures, and demanding otherwise would be a category error.

## Aggregation, and why the regional rows *are* usable

The report is **18 regions × 30 sectors**; this port is locked at **6 × 25**. But the report's
18 regions **nest exactly inside our 6**, so its regional rows aggregate up to ours:

| our REG6 | report rows |
|---|---|
| 1 Sumatra | 1 NorthSumatra, 2 CentrSumatra, 3 Jambi, 4 SumSel, 5 SouthSumatra |
| 2 Java | 6 DKI, 7 NthWestJava, 8 EastJava |
| 3 Kalimantan | 9 KalBar, 10 KalTeng, 11 KalSel, 12 KalTim, 13 KalUt |
| 4 Sulawesi | 14 Sulawesi |
| 5 BaliNusa | 15 Bali, 16 NusaTeng |
| 6 MalukuPapua | 17 Maluku, 18 Papua |

Aggregating the report's percentage changes needs weights (regional shares of the relevant
aggregate), which we have in our own database only approximately — so treat aggregated
regional comparisons as indicative and the **national row as the primary target**.

The sectoral aggregation does **not** nest as cleanly; Table 5 of the report gives the full
185 → 30 mapping if a sectoral comparison is ever attempted.

## Table 2 — regional and national macro variables, % change

| # | Region | Real HousCon | Real Invest | Export vol | Import vol | Real GNE | Real GDP | Employment | CPI |
|---|---|---|---|---|---|---|---|---|---|
| 1 | NorthSumatra | 0.84 | -0.39 | -4.74 | 1.89 | 0.37 | -0.23 | -0.51 | 1.37 |
| 2 | CentrSumatra | 0.46 | -0.63 | -4.58 | 1.21 | 0.08 | -0.40 | -0.89 | 1.19 |
| 3 | Jambi | 1.32 | 1.13 | -2.26 | 2.88 | 1.12 | 0.03 | -0.05 | 1.63 |
| 4 | SumSel | 1.81 | 2.24 | -0.34 | 3.10 | 1.78 | 0.26 | 0.44 | 1.86 |
| 5 | SouthSumatra | 1.03 | -0.14 | -1.27 | 2.11 | 0.57 | -0.14 | -0.33 | 1.48 |
| 6 | DKI | 0.70 | -0.17 | -5.68 | 1.48 | 0.34 | -0.30 | -0.65 | 1.33 |
| 7 | NthWestJava | 0.33 | -0.53 | -4.87 | 0.82 | 0.05 | -0.50 | -1.02 | 1.19 |
| 8 | EastJava | 0.61 | -0.44 | -5.72 | 1.41 | 0.24 | -0.35 | -0.75 | 1.30 |
| 9 | KalBar | 1.00 | -0.23 | -7.05 | 2.61 | 0.54 | -0.16 | -0.35 | 1.71 |
| 10 | KalTeng | 2.21 | 3.30 | -1.08 | 4.65 | 2.22 | 0.51 | 0.83 | 2.31 |
| 11 | KalSel | 6.07 | 15.32 | 5.91 | 10.90 | 7.90 | 2.29 | 4.65 | 3.84 |
| 12 | KalTim | 9.17 | 19.29 | 4.56 | 13.36 | 12.13 | 3.02 | 7.71 | 4.98 |
| 13 | KalUt | 4.80 | 10.69 | 2.65 | 8.53 | 5.97 | 1.69 | 3.40 | 3.21 |
| 14 | Sulawesi | 1.21 | 0.01 | -6.59 | 3.04 | 0.71 | -0.03 | -0.15 | 1.73 |
| 15 | Bali | 0.96 | -0.11 | -5.59 | 2.32 | 0.51 | -0.16 | -0.40 | 1.54 |
| 16 | NusaTeng | 1.34 | 0.05 | -4.06 | 2.85 | 0.75 | 0.04 | -0.02 | 1.75 |
| 17 | Maluku | 1.42 | 0.22 | -5.98 | 2.87 | 0.82 | 0.08 | 0.06 | 1.76 |
| 18 | Papua | 1.23 | -0.21 | -2.83 | 2.50 | 0.62 | -0.02 | -0.13 | 1.63 |
| **19** | **National** | **1.10** | **0.90** | **-3.41** | **2.00** | **0.90** | **-0.09** | **-0.27** | **1.54** |

Table 2 extracted cleanly (`pdftotext -layout`) and every row has all 8 entries. **This is the
primary comparison target.**

## Table 3 — national sectoral variables, % change

⚠️ **Provenance warning, read before using these numbers.** `pdftotext -layout` misaligned this
table: the sector-name column is a separate text frame and came out offset by three rows, which
put Coal's +51.90 price change against `PetrolLNG`. The columns below were recovered with
`pdftotext -table`, which assigns cells by x-position, and the alignment was then confirmed
against three independent economic anchors:

- **Coal price +51.90** — must belong to Coal; the shock is +50% on coal.
- **Electricity price +9.95** and **Cement +10.29** — the two most coal-intensive users.
- **OthNMmetlPrd +5.06** — kiln-fired, third most coal-intensive.

Under the misaligned `-layout` reading those same numbers landed on VehicEquip and Transport,
which have no coal-cost story. Three anchors agreeing is why the table below is trusted.

**The Imports and Exports columns are NOT recovered and are deliberately omitted.** The report
blanks a cell "where imports or exports are very small", and neither extraction mode assigns
the surviving values to rows consistently — `-layout` and `-table` disagree about them. Do not
reconstruct them from the text; re-read them from the PDF if they are ever needed.

| # | Sector | Output | Employment | Price |
|---|---|---|---|---|
| 1 | Agriculture | -0.62 | -0.84 | 0.86 |
| 2 | OthServices | 0.16 | 0.31 | 1.69 |
| 3 | Forestry | -0.41 | -0.72 | 0.33 |
| 4 | Fishing | 0.19 | 0.35 | 2.37 |
| 5 | **Coal** | **2.21** | **14.99** | **51.90** |
| 6 | OilGas | -0.13 | -0.64 | -0.02 |
| 7 | OtherMining | -0.11 | -0.27 | 0.97 |
| 8 | PetrolLNG | -0.18 | -1.02 | 0.72 |
| 9 | FoodProds | -0.55 | -1.30 | 0.88 |
| 10 | TCF | -2.28 | -3.87 | 0.66 |
| 11 | WoodPrd | -1.08 | -1.90 | 1.00 |
| 12 | PaperPulp | -4.68 | -8.33 | 1.48 |
| 13 | OthNMmetlPrd | -1.10 | -1.23 | 5.06 |
| 14 | BasChemicals | -0.81 | -7.89 | -0.03 |
| 15 | Fertilizer | -1.49 | -5.30 | 0.65 |
| 16 | ChemRubrPlas | -1.08 | -2.76 | 0.64 |
| 17 | OthManufact | -2.10 | -2.66 | 1.29 |
| 18 | Cement | -0.72 | -1.32 | 10.29 |
| 19 | BasIronSteel | -5.92 | -11.36 | 1.92 |
| 20 | MetalProds | -1.91 | -4.75 | 0.73 |
| 21 | VehicEquip | -1.02 | -1.84 | 0.70 |
| 22 | Electricity | -1.85 | -7.58 | 9.95 |
| 23 | Construction | 0.79 | 1.53 | 1.89 |
| 24 | Trade | -0.17 | -0.23 | 1.47 |
| 25 | Transport | 0.56 | 1.69 | 1.72 |
| 26 | RealEstate | 0.18 | 2.71 | 3.27 |
| 27 | GeneralGov | -0.07 | -0.10 | 1.62 |
| 28 | GovEducSvc | 0.05 | 0.05 | 1.56 |
| 29 | GovHealthSvc | 0.13 | 0.16 | 1.48 |
| 30 | OthGovSvc | -0.10 | -0.11 | 1.71 |

Note the sector ordering is odd as printed — `OthServices` sits at position 2, between
Agriculture and Forestry. Both extraction modes agree on that ordering, so it is reproduced
as printed rather than "corrected".

## Qualitative claims the report makes — cheap sign tests

These are worth asserting before any magnitude is compared, because a translation error large
enough to matter will usually break one of them:

1. **"exports fell for all sectors but coal"** — exactly one sector has a positive export
   response.
2. **Prices rose for nearly all Indonesian-produced goods**; imports rose for all goods.
3. **Dutch disease.** Terms of trade improve; Indonesia consumes more (Real GNE +0.90) while
   producing slightly less (Real GDP **−0.09**). The opposite signs on GNE and GDP are the
   single most diagnostic pair in Table 2.
4. **Trade balance improves by ≈0.2% of GDP.**
5. Coal provinces (KalTim, KalSel, KalUt, KalTeng, SumSel, Jambi) gain; the most
   trade-exposed non-coal provinces lose most employment and real GDP.
6. Non-coal provinces still raise household consumption — the report attributes this to
   capital income, of which each province spends a constant share of the national total,
   rather than to wage income, which is spent where earned.

## Implementation status (updated 2026-07-31)

Both gaps named in the original version of this file are closed:

| gap | closed by |
|---|---|
| (a) a `COALPRICE` closure and scenario | `COALPRICE_SWAPS` in `src/closures.jl`; `COALPRICE_REFERENCE` in `src/scenarios.jl` |
| (b) an opt-out for the unconditional GDPPI numeraire swap | `initialize_model!(…; numeraire = :gdppi \| :exrate)`; `Scenario` carries the choice, so it cannot be set at the call site |

The comparison itself is `test/verify_coalprice_reference.jl`.

⚠️ The numeraire warning at the top of this file is therefore **resolved** — the CPI column
*is* now comparable, because `COALPRICE_REFERENCE` sets `numeraire = :exrate` to match
`coalprice.CMF:72`. The warning is left in place above because it documents why the column
would otherwise have been silently wrong.

A third gap this file did not anticipate, found only by building the closure: **`E_natfhou!`
was an empty stub.** `natfhou` therefore appeared in no equation, and `solve_newton!`'s
orphan-pinning step squared the system anyway — so every closure tried before this one ran
correctly and the omission was invisible. `:83 swap houslack = natfhou` makes natfhou
exogenous, the orphan vanished, and the model came out one variable short of square.
TERM.TAB:2053 is now implemented. The general lesson is recorded in `VV_PLAN.md` route 1b: an
orphan variable is a missing equation until proven otherwise.

One translation point this file did not anticipate: the shock is **`logpct(50)` = log(1.5),
not `pct(50)` = 1.5**. `fpexp_d` is declared unbounded (`build_model!.jl:214`), so it is a
benchmark-0 additive shifter, and `E_xexpd` uses it in log space beside
`pfexp = log(ppur) − log(phi)`. Passing `pct(50)` would have applied ≈ **+348%** and solved
without complaint — Table 3's Coal price row (+51.90) is the cheapest check that the units
are right.
