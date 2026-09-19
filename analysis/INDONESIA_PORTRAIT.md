# Preliminary portrait of the Indonesian economy from INDOTERM (2016 base)

**Status:** exploratory research note, not a validated result. Every number below is either
(a) a *structural statistic* read directly from the model's 2016 base data, or (b) a
*scenario result* already validated elsewhere in this project and cited to its gate. Nothing
here is a forecast, and nothing here is a new solve.

**Data basis:** `data/national_data.csv` (185-commodity ORANI-G national core) and
`data/regsupp_data.csv` (regional shares `R001`–`R006`, `TSPD`, distances). Aggregations to
25 sectors use `src/aggregation_data.jl` `SEC_MAP_185_to_25`; 34 provinces → 6 island groups
use `REG_MAP_34_to_6`. Linkage indices from `analysis/output/linkage_25sector.csv`; GVC
measures from `analysis/gvc_accounting.jl`; commodity import/DVA screens from
`analysis/tkdn_screening.csv`.

**Unit note.** The base data are in IDR billion. Total expenditure-side GDP in the base is
**12,915,060** (≈ IDR 12,915 trillion), consistent in order of magnitude with official 2016
nominal GDP. This is a 2016 vintage: **pre-nickel-downstreaming, pre-coal-cycle**, so it is a
structural map, not a current-market snapshot.

---

## 1. Macro shape (2016 base)

| Aggregate | Value | Share of GDP |
|---|---|---|
| Household consumption C | 7,554,347 | **58.5 %** |
| Investment I | 3,635,068 | **28.1 %** |
| Government G | 1,611,680 | **12.5 %** |
| Exports X | 2,356,298 | **18.2 %** |
| Imports M | 2,276,289 | **17.6 %** |
| **GDP (expenditure)** | **12,915,060** | 100 % |

Factor income (income side, ≈12,359,554 — see caveat): labour 6,331,821 (**51.2 %**),
capital 5,565,891 (**45.0 %**), land 461,843 (3.7 %).

- **Labour-income dominated, consumption-driven.** Growth shocks transmit mainly through
  household income and consumption (58.5 % of demand), not through investment.
- **Open but not hyper-open.** Trade is ~18 % of GDP on each side; the economy is roughly
  balanced on trade in the base.
- **Caveat (documented):** the income-side and expenditure-side national GDP do not fully
  agree (the model carries a ~4 % national gap, and per-region GDP income≠expenditure up to
  ~36 %). This is a **defect in the shipped 2016 data**, not the translation — see `VV_PLAN.md`
  V1 and `PLAN.md` §N. Use *shares and signs*, not level identities, and never report
  regional income=expenditure.

---

## 2. Key sectors — linkage classification (Rasmussen–Hirschman, 25-sector mean = 1.0)

BL = backward linkage (power of dispersion), FL-Ghosh = forward linkage (supply-side).
Source: `analysis/output/linkage_25sector.csv`.

**Key sectors (above average both ways):**
| Sector | BL | FL-Ghosh | Output | Read |
|---|---|---|---|---|
| **Utilities** | **1.58** | 1.48 | 721,613 | Highest BL in the economy — everything buys power |
| **FoodProc** | 1.17 | 0.88 | **1,917,598** | Largest output; strongly backward (agri pull) |
| **Chemicals** | 1.05 | 1.16 | 1,562,805 | Classic key sector; deep forward into manufacturing |
| **BasMetals** | 1.14 | 1.03 | 272,640 | Key sector, but its pull is **construction**, not machinery |
| **NonMetalPrd** | 1.14 | 1.25 | 190,633 | Cement/glass — coal-intensive, construction-facing |
| **WoodPaper** | 1.08 | 1.18 | 431,984 | Forward into construction and packaging |

**Backward-oriented (pull the economy, don't supply it):** Construction (1.07 BL), Transport,
HotelsRest, Apparel, OthServices.

**Forward-oriented (supply the economy, don't pull it):** **OthMining** (FL 1.30, upstreamness
2.27), **OilGas** (FL 1.51 — the single highest forward linkage), Crops, Forestry, Coal, Livestock.

**Weakly linked both ways (policy caution):** **MetalMach** (BL 0.96, FL 0.86),
**TranspEquip** (BL 0.91, FL 0.83), Trade, FinBusSvc, OthManuf, BevTobText.

> **Contrarian-but-important finding (already documented in `linkage_analysis.md`):** the
> "metals + machinery + transport equipment" cluster is **not** one high-multiplier block. Only
> **BasMetals** is an above-average key sector, and its key-sector status is anchored in
> **construction demand and mining inputs**, not in the machinery chain. The real internal chain
> is directional and narrow: **BasMetals → MetalMach → TranspEquip** (MetalMach supplies 5.5 %
> of TranspEquip's inputs and takes 5.1 % of its own from BasMetals), forming a corridor, not a
> symmetric cluster.

---

## 3. Regions — 34 provinces collapsed to the 6 validated island groups

Regional detail exists for all **34 provinces** (`R001`–`R006` in the base data). The locked
*validation* model is the 6 island groups below.

### 3.1 Island weights and character

| Island | Output share | Resource* share | Manufacturing share | Agriculture share |
|---|---|---|---|---|
| **Java** | **61.5 %** | 3.0 % | **60.3 %** | 14.8 % |
| **Sumatra** | 19.4 % | 9.0 % | 43.5 % | 26.6 % |
| **Sulawesi** | 6.3 % | 11.9 % | 29.8 % | 36.5 % |
| **Kalimantan** | 5.7 % | **27.4 %** | 26.4 % | 17.2 % |
| **MalukuPapua** | 3.8 % | **56.8 %** | **5.0 %** | 21.3 % |
| **BaliNusa** | 3.4 % | 18.4 % | 19.4 % | 32.7 % |

*resource = Coal + OilGas + OthMining share of the island's gross output.

### 3.2 Province-level structure (top provinces)

| Province | Island | Output % | Resource % | Manuf % |
|---|---|---|---|---|
| JaBar | Java | 19.3 | 2.5 | 71.8 |
| JaTim | Java | 16.8 | 3.6 | 56.7 |
| JaTeng | Java | 9.8 | 3.7 | 52.3 |
| DKI | Java | 8.2 | 1.3 | 42.1 |
| Banten | Java | 6.5 | 3.7 | 74.0 |
| SumUt | Sumatra | 3.9 | 3.4 | 49.9 |
| RiauProv | Sumatra | 3.8 | 9.1 | 41.7 |
| **PapuaProv** | MalukuPapua | 2.5 | **70.4** | 1.8 |
| **KalTim** | Kalimantan | 2.2 | **31.1** | 25.3 |
| **BaBel** | Sumatra | 0.95 | **58.6** | 13.9 |
| **MalUt** | MalukuPapua | 0.52 | **50.8** | 7.0 |
| **KalBar** | Kalimantan | 1.27 | **43.5** | 17.6 |

**The core regional fact:** the economy is a **Java-centred manufacturing/consumption core
supplied by resource- and agriculture-exporting peripheries.**
- **Java = 61.5 % of output, 60 % manufacturing** — JaBar/Banten are the factory belt (72–74 %
  manufacturing); DKI is services/consumption-heavy.
- **MalukuPapua is 57 % resource** (PapuaProv alone 70 % — copper/gold/gas) yet only 5 %
  manufacturing: a pure extraction periphery.
- **Kalimantan is the coal island** but is *more diversified* than its reputation (27 % resource,
  26 % manufacturing — mainly wood/paper and metals processing).
- **Sumatra is the most balanced** large island (43 % manufacturing, 27 % agriculture, palm
  oil + coal + gas).

---

## 4. Global supply-chain risk (import dependence and GVC position)

Source: `analysis/gvc_accounting.jl` (VS = imported content of exports, Hummels–Ishii–Yi) and
`analysis/tkdn_screening.csv` (commodity import shares).

**National aggregate VS = 13.29 %** — Indonesia's exports are ~87 % domestically value-added,
low GVC participation by regional standards. But that average hides the exposed sectors.

### 4.1 Sector import dependence (imports as share of total use)

| Sector | Import share of use | Export/output | VS (import content of exports) |
|---|---|---|---|
| **BasMetals** | **50.0 %** | 36.6 % | 13.4 % |
| **MetalMach** | **46.4 %** | 20.0 % | **28.6 %** |
| **Chemicals** | 28.8 % | 18.2 % | 18.3 % |
| **TranspEquip** | 24.8 % | 13.3 % | 19.1 % |
| **OilGas** | 22.2 % | 21.9 % | 4.5 % |
| Apparel | 16.9 % | **38.1 %** | **25.9 %** |
| WoodPaper | 11.8 % | 22.7 % | 15.3 % |
| OthManuf | 11.6 % | 29.6 % | 25.3 % |
| OthMining | 4.5 % | 13.4 % | 8.7 % |
| **Coal** | 3.4 % | **65.9 %** | 8.7 % |

### 4.2 The two opposite exposures

Two very different risk profiles dominate:

1. **Import-dependent manufacturing (supply-chain fragility).** BasMetals, MetalMach,
   Chemicals, TranspEquip, Apparel. MetalMach is the sharpest: ~46 % import share and **28.6 %
   of its exports are imported content** — the deepest GVC participant in the economy.
   Commodity level (TKDN), the extreme cases are near-total import dependence:
   **Aircraft 96.5 %, Varnish 91.7 %, CoarseSalt 89.5 %, PrimaryMover 84.0 %,
   BasIronSteel 68.8 %, Soy 65.4 %, Electronics 60.0 %, Textile 45.0 %, Plastics 46.3 %,
   Fertilizer 27.9 %.** These are the items where a trade disruption bites first.

2. **Export-dependent resources (demand/price fragility).** **Coal (65.9 % of output
   exported)**, Apparel (38.1 %), BasMetals (36.6 %), OthManuf (29.6 %), WoodPaper (22.7 %),
   OilGas (21.9 %). These sectors are *not* import-fragile (Coal VS is only 8.7 %) — their risk
   is a **world-price / external-demand** shock, not a sourcing shock.

3. **The energy-balance anomaly:** Indonesia exports coal/gas but is a **net crude importer**
   (CrudeOil import share **43.1 %**). This is the classic dual exposure — the same price move
   (e.g. a coal boom) helps the coal periphery and hurts the oil-importing core.

**Siting fact for the metals chain:** the upstream ore base is geographically specific —
the model's `OthMining` export base is **72 % MalukuPapua / 18 % BaliNusa** (copper concentrate),
while BasMetals/MetalMach capacity is **~83 % in Java** (JaBar/Banten/JaTim). The ore is in the
east; the smelter and the machinery are in Java. That geographic mismatch is the structural
tension any *hilirisasi* (downstreaming) policy has to cross.

---

## 5. Domestic supply-chain risk (inter-provincial)

Sources: `R001` regional output shares, `regsupp_data.csv` distances, `SCALE34_PLAN.md` §1.3.

- **The inter-provincial trade graph is complete.** In the model *every* province trades every
  commodity with *every other* province (`xtrad[a,s,r,d]` exists for all r,d), and inter-regional
  trade is a large share of each region's use. This is why the 34-region Jacobian has no small
  separator and why factorizing it is hard — it is also the economic statement that **provinces
  are tightly coupled, not a set of independent islands.**
- **Concentration risk is severe and sector-specific:**
  - **Coal:** HHI 0.44 — **90 % of output in three provinces** (KalTim, KalSel, SumSel). A coal
    shock is not a national shock; it is a Kalimantan/South-Sumatra shock with fiscal spillovers.
  - **MetalMach:** HHI 0.45 — **83 % in JaBar/Banten/JaTim.** Any disruption to the Java
    industrial belt takes out most of the machinery supply chain at once.
  - Apparel HHI 0.27 (75 % in Java), BasMetals HHI 0.15 (59 % Java), Chemicals HHI 0.14 (56 % Java).
- **The domestic multiplier that stays home.** GVC "DVA-from-other-domestic-sectors" (the value
  added a sector pulls from *other* domestic industries — the thing an onshore-deepening policy
  is trying to raise) is highest for **FoodProc 52.1 %, Utilities 48.6 %, BasMetals 47.1 %,
  NonMetalPrd 44.7 %, Construction 43.3 %** — and notably **low for MetalMach (29.7 %)** and
  **TranspEquip (25.5 %)**, consistent with assembly-type, import-fed activity.
- **Distance matters but is not modelled as a cost shock instrument here:** `DISTGONE` records
  average origin→destination distances per commodity (rice ≈ 0.2–3.1 thousand-km units across
  province pairs; a full 34×34 matrix exists). Eastern provinces (PapuaProv, MalUt) sit at the
  long end of nearly every commodity's distance distribution.

**Net read:** domestic risk is a **concentration risk**, not a connectivity risk. The system is
well connected; what is fragile is that critical sectors sit in one or two places.

---

## 6. Vulnerabilities, with the project's own shock evidence

| Vulnerability | Evidence | Status |
|---|---|---|
| **Coal-price exposure is a regional, not national, event** | +50 % coal export price → national real GDP **−0.205 %**, CPI **+1.42 %**, exports **−3.88 %**, coal output +2.35 %; but **Kalimantan real GDP +5.2 %** vs every other island **<0.4 %** | Validated (V9 vs published GEMPACK; V4 sign-robust at 6 and 12 regions) |
| **Dutch-disease pattern** | Terms of trade improve; real GNE +0.90 % while real GDP −0.09 % in the reference run | Validated (V9 qualitative claims) |
| **Short-run closure fragility** | The short-run closure **folds** under the coal shock at λ*≈0.9998 across three independent attempts (200 clean steps) | Disclosed limitation (V6) — closure ordering "not established", not a bug |
| **Regional sign fragility near zero** | Region-level *signs* for a small set of near-zero variables (concentrated in BaliNusa `xinvi`/`xinv`) flip with resolution | Disclosed limitation (V4) — screen with `print_regional_confidence_report` |
| **Magnitudes are elasticity-conditional** | All 10 headline metrics sign-stable across a ±50 % elasticity sweep, but the ranges are **lower bounds** (1 of 13 points folded and is missing) | Disclosed limitation (V8) |
| **Import-dependence of the industrial belt** | BasMetals/MetalMach/Chemicals/TranspEquip 25–50 % import share | Structural (base data) |

---

## 7. Can we use all 34 regions now? — **Yes for structure, 34-province shocks `non-gating` (2026-09-17).**

Short answer: **I can read, aggregate, and report structure at all 34 provinces today. 34-province *shocked* impacts are closed `non-gating` by user decision — 6 island groups stays the publication model.**

What works at 34 regions (verified in `logs/`):
- The data pipeline: `data/cache_34reg.jls` builds cleanly (34 regions × 25 sectors, no flattening).
- The full model builds: **2.4 M variables**, ~98 s.
- The **benchmark solve converges**: ‖F‖∞ = **8.7e-9**, via the Schur condensation solver.
- Factorization fits: **234 M LU nonzz, ~15.7 GB peak** on a 23.8 GB machine.
- All regional base data (`R001`–`R006`, `TSPD`, distances) is available at 34 provinces — which
  is what this note's regional tables are built from.

What does **not** yet work at 34 regions:
- **Shock convergence.** A 0.1 % shock (`blabnat 1.0→0.999`) converges to t=0.35 then **stalls
  at t=0.6** (`maxit`, ‖F‖∞ = 7.2e-5). A 1 % shock stalls outright (line search degenerates over
  20 backtracks). So **34-region scenario numbers would be unreliable and must not be quoted.**
- Consequently the model's own rule stands: **6 island groups is the locked, validated
  publication model; 12 and 34 regions are test instruments** (`AGENTS.md`, `VV_PLAN.md` V4).

Practical rule: **report 34-province *structure* (shares, linkages, exposure); report
6-island *impacts* (scenario results),** and say which is which.

---

## 8. Where to go deeper (ranked by value per effort)

1. **Provincial exposure mapping** — combine §3 province weights with §4 sector import/export
   exposure to produce a province × risk matrix (e.g. Banten's 74 % manufacturing × MetalMach's
   46 % import share). Cheap, purely from base data, and directly policy-relevant.
2. **A second scenario** under the validated 6-region closure (e.g. oil/gas-price shock, or an
   ore-export restriction) to complement the coal story — `analysis/gvc_policy.jl` already
   defines a 4-leg *hilirisasi* decomposition (ban / smelter-in-Java / smelter-at-the-ore / full);
   the interaction term is the actual policy question and it has not been run to completion.
3. **34-region shock convergence** (Stage 2b.4) — **closed `non-gating` 2026-09-17** (line-search / conditioning wall remains, but not required for validation).
4. **IO-2020 rebase** (`IO2020_REBASE_PLAN.md`) — moves the portrait from 2016 (pre-nickel-boom,
   pre-coal-cycle) to a current structure. Would materially change §3/§4 for nickel and coal.

---

## Caveats to carry into any write-up

- 2016 base: **pre-nickel-downstreaming, pre-coal-cycle.** Not a current-market description.
- Regional *impacts* are at 6 island groups only; provincial claims are out of scope until
  Stage 2b.4 passes.
- V4 (near-zero regional signs, notably BaliNusa), V6 (short-run folds), V8 (elasticity ranges
  are lower bounds) — all disclosed, all load-bearing on how numbers may be stated.
- Income-side ≠ expenditure-side GDP in the shipped data (national ~4 % gap, regional up to
  ~36 %): report shares and signs, not level identities.
