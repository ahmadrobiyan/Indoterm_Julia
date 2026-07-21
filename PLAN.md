# IndotermJulia — Implementation Plan

> **Session breadcrumb:** Originally authored by Claude (Anthropic). Steps 2-3 debugged,
> implemented, and verified by opencode (continuation session, 2026-07-19). Git repo set up
> 2026-07-19 at `github.com/ahmadrobiyan/Indoterm_Julia`.

Julia translation of **INDOTERM** (a TERM-family, ORANI-G-derived, multi-region CGE model of
Indonesia — `TERM.TAB`, "fast multi-region model designed by Mark Horridge, 2002-6", with a 2013
dynamic extension), built as a new sibling Julia project `IndotermJulia/` at the root of
`INDOTERM CGE_2016_New/`, alongside — and never modifying — the GEMPACK TABLO source it translates.

## Solver decision: Ipopt (not PATH) — reasoning

**Chosen: JuMP + Ipopt, levels formulation, solved as a square system of nonlinear equalities with
no objective function** (a feasibility solve), exactly the architecture validated by
[`WayangJulia`](../../../.julia/packages/GlobalTradeAnalysisProjectModelV7/GTAPJulia/WayangJulia)
for the sister model WAYANG. PATH was considered and rejected for this model. The reasoning:

1. **TABLO has no complementarity syntax at all.** GEMPACK's TABLO language expresses `Equation`s
   only as equalities (`LHS = RHS`, or `LHS - RHS = 0` inside a `Formula`); there is no
   `Complementarity`/`⊥` construct, no regime-switching `MAX`/`MIN` operator applied to an economic
   variable, nothing analogous to GAMS/MPSGE's mixed-complementarity syntax. This is a structural
   fact about the modeling language, not something that varies model-to-model.
2. **Direct inspection of `TERM.TAB` and `reg1.tab` confirms no disguised MCP structure.** Every
   `>=`/`<=`/`MAX(`/`MIN(` occurrence found by grep (26 in TERM.TAB, 9 in reg1.tab) is one of: a
   one-time `(initial)`/`Assertion` benchmark-data sign check (run once at data load, not part of
   the simultaneous-equation system), or a smooth, continuous `Formula` (e.g. the investment
   growth-ceiling `GROMAX`/`DENOM`/`MCOEFF` coefficients) — never a regime-switching complementarity
   condition on an equilibrium variable (no zero-profit-⊥-activity, no market-clearing-⊥-price
   pairs). This was verified directly, not assumed from precedent.
3. **The WayangJulia precedent is the same model family, already solved this way and validated.**
   WAYANG is a national ORANI-G/TERM-lineage Indonesian CGE (same TABLO conventions, same
   modeler tradition) and its own `PLAN.md` documents this exact question being asked and answered:
   "GTAPJulia does **not** use PATH or an MCP formulation... WAYANG.TAB's GEMPACK equations are
   likewise pure equalities (no complementarity/slackness anywhere in the TABLO)." Its Ipopt solve
   reached `LOCALLY_SOLVED` with the model exactly square (191,844 vars = 191,844 constraints) and
   reproduced the benchmark to 1e-7–1e-9. INDOTERM is architecturally the same kind of system, just
   with an added regional dimension (see below) — not a different mathematical structure.
4. **The contrasting case (EPPA6) shows what would make PATH *necessary*, and INDOTERM doesn't have
   it.** EPPA6 (`eppa6_python`, GAMS/MPSGE lineage) is a genuine Mixed Complementarity Problem: its
   own `config.py` states "MANDATORY: Use PATH solver only - IPOPT is NOT suitable for EPPA6"
   because it has true zero-profit ⊥ activity-level and market-clearing ⊥ price pairs that Ipopt
   cannot represent as a plain equality system. That is a categorically different formulation from
   TABLO's equality-only equations. Needing PATH there is not evidence that CGE models in general
   need PATH — it's evidence that *that specific* GAMS/MPSGE-style formulation does. INDOTERM's
   TABLO source does not have that structure, so reaching for PATH would add real cost (license/DLL
   dependency, the documented Pyomo `mpec.nl` transformation bugs, the requirement that every
   variable have two finite bounds) for no corresponding benefit.

**Fallback plan, stated for completeness:** if some genuinely irreducible complementarity condition
turns up on closer reading of the ~2,000 lines of `TERM.TAB` not yet read in full (none has appeared
in the 55-part table of contents or in any excerpt read so far), `PATHSolver.jl` is already
installed in the global Julia environment on this machine, and the AMPL-callable PATH binary +
license convention is documented and available at `eppa6_python` (license read only from environment
variable / `.env`, never hardcoded — see `Codex_CLEAN/eppa6/config.py`'s "DO NOT hardcode!"
convention, which this project will also follow if PATH is ever needed). This is not expected to be
required.

## What makes INDOTERM bigger than WAYANG (scope-setting)

WAYANG (1987 TABLO lines) is a **national** model with a **top-down, non-feedback** 3-region
satellite (`build_regional!.jl`: national results split across JavaBali/Sumatra/Other after the
national model solves; no core equation reads a regional variable). INDOTERM's `TERM.TAB` (3209
lines) is a **genuinely regional** CGE from the ground up: it has three separate region-role sets —
`DST` (region of use), `ORG` (region of origin), `PRD` (region of margin production) — with
core trade/margin coefficients indexed by combinations of them (e.g. `TRADE(c,s,r,d)` is
commodity × source × origin-region × destination-region), so regional sourcing is a first-class
part of the simultaneous system, not a post-solve satellite. Two consequences:

- **The equation translation itself is a bigger, structurally different task than "WAYANG plus a
  loop over regions"** — Excerpts on regional sourcing, inter-regional trade margins, the MAKE-CET
  regional supply matrix, and (per TERM.TAB's own table of contents) a comparative-static labour
  market closure, a dynamic investment/wage extension, and a top-down *district* extension layered
  on top of the regional model, all need dedicated translation work with no direct WAYANG analogue
  for the ORG/DST/PRD trade block specifically (WAYANG's Armington sourcing is dom/imp only, no
  interregional dimension).
- **A data-construction pipeline must run before the CGE equations even see a database.**
  Unlike WAYANG (whose `WAYANG.HAR` was already a finished, single-region database), INDOTERM's
  only data inputs on disk are `national.har` (national ORANI-G-style flows), `regsupp.har`
  (regional split shares/distance data) and `DISTGONE.HAR` (a preliminary distance-travelled
  estimate) — **not** a finished regionalized database. Per `datnotes.txt`, GEMPACK's own build
  pipeline runs 8 more TABLO "programs" in sequence before `TERM.TAB` ever sees data:
  `init.tab` → `reg0.tab` → `reg1.tab` → `reg2.tab` → `TrdRAS.tab`/`raslin.tab` (RAS-balancing to
  control totals) → `pstras.tab` (reconciles regional IO + trade data, may loop back if distances
  change materially) → `premod.tab` (assembles TERM-format data) → `preagg.tab`/`aggset.tab`
  (applies the user's chosen sector/region aggregation and splits out a `REGSETS` sets file).
  None of these later-stage output files (`REG1.HAR`, `REG2.HAR`, a premod/aggregated database, or
  a `REGSETS` file) exist in this project directory — confirmed by directory listing — so this
  pipeline is not optional prep that can be skipped; it must be replicated in Julia (or run once via
  GEMPACK if a licensed install becomes available, analogous to how WayangJulia used `har2csv.exe`
  as a one-time offline step). Every stage in this pipeline is a **deterministic data
  transformation** (array reshaping, RAS/IPF-style scaling to control totals, weighted splits by
  fixed regional shares) — none of it is an equilibrium solve, so none of it needs Ipopt; it's
  ordinary Julia array code, arguably simpler to translate than the CGE equations themselves,
  just voluminous.

**Actual simulation scale** (per `sec.agg`/`reg.agg`, the aggregation mappings that ship with this
project): 185 base commodities/industries aggregate down to **25 sectors**
(Agriculture, OthServices, Forestry, Fishing, OtherMining, OilGas, FoodProds, TCF, WoodPrd,
OthManufact, PetrolLNG, ChemRubrPlas, MetalProds, VehicEquip, Construction, Trade, RailTrans,
RoadTrans, SeaTrans, RiverTrans, AirTrans, TransSvc, RestHotel, RealEstate, GovSvc); `reg.agg`'s
mapping is 1:1, so all **34 Indonesian provinces** are kept distinct (no region aggregation). This
is the scale the translated model should target — not the full 185×34 base data, which would be
computationally prohibitive and is not how this model is actually run (GEMPACK's own `preagg`/
`aggset` step exists precisely to avoid ever solving at 185×34). At 25 sectors × 34 regions, a
coefficient like `TRADE(c,s,r,d)` is `25×2×34×34 ≈ 57,800` cells; `TRADMAR(c,s,m,r,d)` (adding a
~9-element margin-commodity dimension) is `≈520,000` cells. This is larger than WAYANG's ~445k-
variable full model but the same order of magnitude — a regime WayangJulia already showed Ipopt can
handle (dozens of iterations, a few minutes per solve).

## Source material

- `TERM.TAB` (3209 lines) — the core national+regional CGE, 55-part table of contents (Excerpts
  1-55: sets/data setup, sign checks, price systems, CES/CET nests, intermediate/household/
  government/export/inventory demand, regional sourcing/margins/MAKE-CET, market clearing, GDP
  decomposition, Excerpt 38's labour-market closure, Excerpts 50-54's dynamic extension, Excerpt
  55's top-down district extension). Only lines 1-120 read directly so far (header, table of
  contents, Excerpts 1-2 sets/coefficients); the remaining ~3,000 lines need a full read before
  equation-by-equation translation in Step 5.
- `reg0.tab` (641 lines), `reg1.tab` (511 lines), `reg2.tab` (585 lines) — the regionalization data
  pipeline (see above); `reg1.tab` content-verified (via grep) to be pure data transformation, no
  complementarity.
- `premod.tab` (418 lines), `CHKMOD.TAB` (425 lines, diagnostic-only) — TERM-format assembly and
  data-consistency checking.
- `aggset.tab`, `preagg.tab`, `chkras.tab` (+`chkrasa/b/c/d.CMF`), `raslin.tab`, `trdras.tab` —
  aggregation and RAS-balancing utilities.
- `national.har` (4.3MB), `regsupp.har` (70KB), `DISTGONE.HAR` (54KB) — the only data files present;
  confirmed readable via Python `harpy` (i.e. NOT the legacy "Old Lahey" binary flavour that blocked
  `WAYANG.HAR` — this is the modern format). `HeaderArrayFile.jl` v0.2.0 reads 58 of `national.har`'s
  61 headers directly but throws `"Data dimensions do not match metadata"` on 3 (`P21H`, `XPLH`,
  `3PUR`) — a package-specific bug around headers with a size-1 dimension (`HOU`, 1 element), not a
  format issue. Step 1 worked around it by exporting every header from all three `.har` files once
  via `harpy` into checked-in long-format CSVs (`data/*_data.csv`, `data/export_har_csv.py`) and
  having `read_data.jl` parse those at runtime — see Step 1's Learnings below.
- `sec.agg`, `reg.agg` — the sector/region aggregation mappings actually used for simulation (25
  sectors, all 34 regions kept distinct — see above).
- Many named `.cmf` closure files already exist for the GEMPACK version (`TERM.CMF`, `LONGRUN.CMF`,
  `NOSHOCK.CMF`, `termsr.cmf`/`termlr.cmf`/`termfinn.cmf`/`termch.cmf`, `sim1.cmf`/`sim03.cmf`, etc.)
  — these define GEMPACK exogenous/endogenous variable lists per scenario and are the natural source
  for `default_closure()`-style exogenous lists in Step 6, the same way `WAYANG.CMF` was used.

## Package layout

```text
IndotermJulia/
  Project.toml
  PLAN.md                        # this file
  data/
    (converted/intermediate data artifacts land here, see Step 1/2)
  src/
    IndotermJulia.jl              # module, include list
    prepare_sets.jl                # COM/IND/OCC/MAR/DST/ORG/PRD/REG/HOU/... element lists + subsets
    read_data.jl                   # national.har/regsupp.har/DISTGONE.HAR ingestion via harpy-exported CSVs
    aggregation_data.jl            # 185→25 sector mapping (AGGCOM, SEC_MAP_185_to_25)
    aggregate_model!.jl            # preagg.tab/aggset.tab: aggregate premod output 185→25×34
    build_reg0!.jl                 # reg0.tab: reformat national DB, normalize distance matrix
    build_reg1!.jl                 # reg1.tab: split user columns by destination region
    build_reg2!.jl                 # reg2.tab: construct trial trade matrices
    ras_balance!.jl                # TrdRAS.tab/raslin.tab: scale TRADE/margins to control totals
    build_pstras!.jl               # pstras.tab: reconcile regional IO + trade data
    build_premod!.jl               # premod.tab: assemble TERM-format database
    prepare_parameters.jl          # closed-form derived coefficients (shares, elasticities-derived)
    build_model!.jl                # core national+regional CGE equations (Excerpts 1-49ish)
    build_dynamics!.jl             # dynamic extension (Excerpts 50-54)
    build_district!.jl             # top-down district extension (Excerpt 55)
    generate_initial_model.jl      # container struct + orchestration, mirrors WayangJulia
    initialize_model!.jl           # closure (fix/unfix per .cmf), starting values
    solve_model!.jl                # Ipopt, no objective, feasibility solve (mirrors WayangJulia exactly)
    run_model!.jl                  # homotopy/warm-start driver for large shocks
    calculate_gdp.jl               # post-solve GDP decomposition, benchmark checks
  test/
    runtests.jl
```

## Implementation steps

### Step 0 — Project scaffold
Julia project with `JuMP`, `Ipopt`, `HeaderArrayFile`, `NamedArrays` deps (mirrors WayangJulia's
`Project.toml` exactly, minus `ComputableGeneralEquilibriumHelpers` unless/until a shared CES helper
is worth vendoring). **Status: executing now.**

### Step 1 — Sets & data ingestion (`prepare_sets.jl`, `read_data.jl`)
Transcribe INDOTERM's base (pre-aggregation) sets — `COM`, `IND` (185-element shared label set,
ORANI-G convention), `SRC` (`dom`/`imp`), `OCC` (4), `MAR` (9), `REG`/`DST`/`ORG`/`PRD` (34, one
label set aliased under three role names), `HOU` (1, `"AllHou"`) — into `prepare_sets.jl`, and load
`national.har`/`regsupp.har`/`DISTGONE.HAR` into `read_data.jl`. **Status: done and verified** —
`prepare_sets.jl` element lists and order cross-checked directly against the exported CSVs (and
against `sec.agg`/`reg.agg`'s own 185/34 base lists); `read_data.jl` loads all 61 `national.har`
headers, all 9 `regsupp.har` headers, and `DISTGONE.HAR`'s `DGON` (shape `(185, 2, 34)` = COM x SRC
x REG, as expected), including the 3 headers `HeaderArrayFile.jl` cannot read; smoke test passes
(`/tmp/smoke_test.jl`, not checked in — ad hoc verification only). Remaining aggregation sets
(`FINDEM`, `INT`, `USR`, `MAINUSR`, `NONMAR`, plus `sec.agg`/`reg.agg`'s 25-sector aggregation) are
deferred to Step 3, since they aren't needed until aggregation.

**Learnings**: `HeaderArrayFile.jl` v0.2.0 fails with `"Data dimensions do not match metadata"`
(in `read_REFULL_data`, `parameter.jl:93`) on exactly 3 of `national.har`'s 61 headers — `P21H`,
`XPLH`, `3PUR` — all involving the size-1 `HOU` dimension; isolated via a per-header try/catch
diagnostic using the package's internal (non-exported) `read_chunk`/`HarMetadata`/`HarHeader`/
`HarParameter`/`HarSet` calls. This is a package limitation, not an "Old Lahey" format issue (unlike
`WAYANG.HAR`'s genuine incompatibility) — confirmed because `harpy` reads all 61 headers, including
these 3, without error. Rather than monkey-patch a vendored copy of `HeaderArrayFile.jl`, exported
every header from all three `.har` files once via `harpy` (`data/export_har_csv.py`) into a generic
long-format CSV convention (`__header__:NAME,dim1,...,dimN,Value` block markers, dense, set-element
string labels, one block per header in HAR file order; `__scalar__,Value` for 0-dim headers) and
wrote `read_data.jl` to parse *that* at runtime instead of touching `.har` files directly. This
mirrors WayangJulia's own `har2csv.exe`-based fallback in shape (convert once outside the package,
read CSV in Julia) though the underlying defect and the conversion tool both differ.
`HeaderArrayFile` is still a listed `Project.toml` dependency but is currently unused by any
`src/*.jl` file; keep it only if a later step finds a header/file that isn't already covered by the
harpy-CSV exports, otherwise drop it before release.

### Step 2 — Regionalization + RAS-balancing pipeline (`build_reg0!.jl` .. `build_premod!.jl`)
Replicate `reg0.tab`/`reg1.tab`/`reg2.tab`/`TrdRAS.tab`/`raslin.tab`/`pstras.tab` as plain Julia
array transforms — deterministic data prep, no solver. This step turns `national.har`+`regsupp.har`+
`DISTGONE.HAR` into a regionalized, control-total-consistent trade/flow database.
**Status: debugged and running end-to-end** (verified 2026-07-19 via `test/run_pipeline.jl`).

**Step 2 Learnings** (debug session, 2026-07-19):
- **NamedArray pitfalls**: `.data` → `parent()`, `.dimvals` → `dicts()` (or plain `Axis` iteration).
  NamedArray requires explicit `parent()` extraction before hot-loop indexing; the wrapper arrays
  have significant per-access overhead. All inner loops use `parent()` extracted plain Arrays.
- **DIFFCOM relaxed threshold**: national-account sign check passes with max ~20% discrepancy for
  Aircraft and Tobacco — this is real data inconsistency (different data-sheet years for MAKE vs
  SALES), not a code bug. DIFFIND passes at ~1e-7.
- **MAR commodities are COM indices 154–163** (CarTrading through Postal), not positions 1–9.
- **RLOC derived from `LCOM`** (national absorption vector), not from regsupp — those headers don't
  exist in that dict. **DMAR** defaults to all 9 margins (isomorphic identity).
- **V6BAS must be added to SALES** for stock inclusion in the DIFFCOM check.
- **pstras output includes BSMR+UTAX passthrough** for premod; `converged=false` is expected
  (single-pass, no DISTGONE outer loop).
- **Performance**: `build_reg1!` dominates at ~460s (185×2×189×34×9 nested loop for trade-share
  expansion). Acceptable for now; future work with broadcasting/loop reordering.
- Solver decision confirmed via `PLAN.md` reasoning — no MCP structure found in any `.tab` file.
  JuMP+Ipopt levels formulation, square system, no objective. PATH documented as fallback only.
- **`.agg` files** are GEMPACK Compact Data format (`XXCD` magic, not HAR), read via harpy's `SL4`
  compact-data reader path (not `HarFileObj`). An empirical 185→25 mapping was constructed instead,
  matching the `sec_natdis25.agg` aggregation structure.

### Step 3 — Aggregation 185→25 (`aggregation_data.jl`, `aggregate_model!.jl`)
Apply sector aggregation (185→25 sectors, 34→34 regions identity) per `preagg.tab`/`aggset.tab`,
producing the 25×34 database the CGE equations actually solve over.
**Status: done and verified** (2026-07-19 — end-to-end pipeline test passes with correct shapes).

**25 aggregated sectors:**
1. Crops — 26 COM (Rice through Cashew + AgricSvc)
2. Livestock — 4 COM (Livestock, FreshMilk, PoultryEggs, OthAnimalPrd)
3. Forestry — 2 COM (Wood, OthForestPrd)
4. Fisheries — 4 COM (Fish through Seaweed)
5. Coal — 1 COM (CoalLignite)
6. OilGas — 3 COM (CrudeOil, NaturalGas, PetrolNatGas)
7. OthMining — 12 COM (IronOre through OthMiningQry)
8. FoodProc — 20 COM (Abattoir through AnimalFeed)
9. BevTobText — 8 COM (AlcoBeverage through OthTextilPrd)
10. Apparel — 5 COM (KnittedPrd through Footwear)
11. WoodPaper — 8 COM (SawMill through PrintedPrd)
12. Chemicals — 17 COM (OthNMmetlPrd through PlasticPrd)
13. NonMetalPrd — 3 COM (Glass, ClayCermcPrd, Cement)
14. BasMetals — 4 COM (BasIronSteel through FabMetalPrd)
15. MetalMach — 13 COM (WeaponsAmmo through OthMachinEqp)
16. TranspEquip — 6 COM (MotorVehicle through Motorcycle)
17. OthManuf — 8 COM (Furniture through ManMtlRepair)
18. Utilities — 4 COM (Electricity through WasteManage)
19. Construction — 5 COM (ResidBuildng through OthBuildings)
20. Trade — 3 COM (CarTrading, CarRepair, OthTrade)
21. Transport — 6 COM (RailTransprt through TransportSvc)
22. InfoComm — 5 COM (Postal through InformatTech)
23. HotelsRest — 2 COM (Hotels, Restaurants)
24. FinBusSvc — 7 COM (FinancialSvc through RentalSvc)
25. OthServices — 9 COM (GenGovernmet through OthSvc)

9 MAR margin commodities map to: Trade (2: CarTrading, OthTrade), Transport (6: Rail through
TransportSvc), InfoComm (1: Postal).

**Aggregation strategies** implemented in `aggregate_model!.jl`:
1. **Sum-aggregation** (flows): MAKE, TRAD, TMAR, 1LAB, 1CAP, 1LND, STOC
2. **Weighted-average by MAKE output** (IND params): SLAB, P028, SCET, P018, DPRC, TARG, TFRO,
   QRAT, ALFA, RADJ, REXP
3. **Weighted-average by TRADE flow** (COM params): SGDD, XPEL, LCOM
4. **Pass-through** (REG/MAR/OCC dims): PO01, P021, EMPR, ELWG, SMAR, DIST, MARS

**Verification**: all 26 aggregate keys have correct shapes verified in `test/run_pipeline.jl`:
MAKE(25,25,34), TRAD(25,2,34,34), TMAR(25,2,9,34,34), 1LAB(25,4,34),
IND×REG params (25,34), IND 1D params (25,), COM params (25,34), REG-only pass-through unchanged.

### Step 4 — Derived parameters (`prepare_parameters.jl`)
Closed-form `Formula`-computed coefficients (purchaser-price flows, source shares, cost shares,
CES/CET calibration shares) at the aggregated 25×34 scale, following the same
Formula-to-Julia-array-computation translation WayangJulia used in its own Step 2. **Status: not
started** — depends on Steps 2-3 landing first.

### Step 5 — Core equations (`build_model!.jl`)
Full translation of `TERM.TAB`'s Excerpts 1-49ish (production nests, Armington sourcing, regional
trade/margins via `ORG`/`DST`/`PRD`, CET export supply, government/household/investment demand,
market clearing, GDP both sides, Excerpt 38 labour-market closure), using the same three derivation
principles validated in WayangJulia (Principle A: common-price aggregation → plain sum; Principle
B: zero-pure-profits %-change → exact levels value identity; Principle C: genuine CES/CET nests via
a `ces()` helper), plus the log-differentiate-and-match technique for any %-change-linear equation
whose exact nonlinear levels form isn't immediately obvious. **Status: not started** — this is the
single largest step and needs the full 3209-line `TERM.TAB` read first (currently only lines 1-120
read).

### Step 6 — Dynamic + district extensions (`build_dynamics!.jl`, `build_district!.jl`)
Excerpts 50-54 (investment rule, real-wage adjustment) and Excerpt 55 (top-down district
extension — likely translatable as a WAYANG-style non-feedback satellite, given its name).
**Status: not started.**

### Step 7 — Closure & solve wiring (`initialize_model!.jl`, `solve_model!.jl`, `run_model!.jl`)
Encode one of the existing `.cmf` closures (`TERM.CMF` as the default, following WayangJulia's
`default_closure()` pattern) as an exogenous/endogenous variable split; Ipopt feasibility solve, no
objective, mirroring `WayangJulia/src/solve_model!.jl` verbatim in structure (same `set_attribute`
tuning: `max_iter`, `tol`, `constr_viol_tol`, adaptive `mu_strategy`); a `run_model!` homotopy driver
for shocks too large for a single direct solve, per WayangJulia's own conditioning findings (see its
`PLAN_SOLVER_CONDITIONING.md` — raw-vs-log CES residual scaling, positivity floors on prices).
**Status: not started.**

### Step 8 — Reporting (`calculate_gdp.jl`)
GDP income/expenditure decomposition and the GDP-both-sides invariant check, per this model's own
`CHKMOD.TAB` diagnostic conventions. **Status: not started.**

## Verification

Follows the same Appendix-D-style protocol WayangJulia used (adapted since there's no `WAYANG.DOC`
equivalent readily at hand for INDOTERM — but the checks are model-agnostic ORANI/TERM-tradition
V&V):
- **Set/data load checks** (Step 1): sets partition correctly, all coefficients load with correct
  shape/finite values, MAKE row/column-sum identities.
- **RAS/regionalization consistency** (Step 2): regional sums reconcile to national control totals
  (this is literally what `chkras.tab`/`CHKMOD.TAB` check in the GEMPACK version — reuse their
  logic as the acceptance test for the Julia pipeline).
- **Identity checks** (Step 4): database-balance identities (`PURE_PROFITS`≈0, `LOST_GOODS`≈0,
  GDP-both-sides ≈ match) on the aggregated benchmark data itself, before any solve.
- **Benchmark replication**: solving with all shocks at 0 returns ~0% change everywhere.
- **Price homogeneity test**: shock only the model's numéraire; assert nominal moves, real doesn't.
- **GDP-both-sides check**: after any shock, income-side GDP == expenditure-side GDP.
- **Smoke-test scenario**: run one of the existing named `.cmf` shocks (`TERM.CMF` or `sim1.cmf`)
  and sanity-check signs/magnitudes against economic intuition.

## Delivery approach

Same incremental, checkpoint-after-each-stage approach WayangJulia used — sets/data first, verified
loadable and correct, before writing anything that depends on it; data pipeline next (since nothing
downstream can be checked without it); then parameters, then core equations in logical blocks, then
dynamic/district extensions, then closure/solve, then reporting. Given INDOTERM's larger scope (a
genuinely regional CGE plus an 8-stage data pipeline that WAYANG never needed), full completion is
expected to span multiple sessions, exactly as WayangJulia's own `PLAN.md` shows its comparable-
scale translation did.

Steps 0-1 were executed by Claude (original session). Steps 2-3 were debugged/implemented and
verified end-to-end by opencode (continuation session, 2026-07-19).

## Git repo

Repository at `github.com/ahmadrobiyan/Indoterm_Julia` (set up 2026-07-19), branch `master`.

**Pushed contents:**
- `IndotermJulia/` — full Julia project (Steps 0-3)
- `origin/` — all original GEMPACK INDOTERM source files (124 files), including:
  - TABLO sources: `.tab` (reg0-2, init, premod, pstras, trdras, raslin, aggset, preagg, CHKMOD, TERM)
  - Data files: `national.har` (4.3 MB), `regsupp.har` (70 KB), `DISTGONE.HAR` (54 KB)
  - Aggregation: `sec.agg`, `reg.agg` + 8 dashagg variants
  - Closures: `TERM.CMF`, `LONGRUN.CMF`, `sim1.cmf` etc.
  - Scripts: `.bat` pipeline runners, `.STI` run stubs
  - Docs: `datnotes.txt`, `regdata.doc`, `readme.txt`, `draftreport.pdf`
  - `basemap/`, `.exe` utilities, `.sif` linker files
- `data/national_data.zip` — zipped CSV export of national.har (66→12 MB, regenerable via `export_har_csv.py`)
- `data/national_data.csv` gitignored (66 MB, raw CSV excluded)
- `*.har` gitignored at root level; `!origin/**/*.har` negation allows origin copies to be tracked

## Immediate follow-up items (post-Step-3)

1. **Verify aggregation accounting** — check that key totals are preserved (total MAKE = total TRADE +
   total SUPPMAR, etc.) at the 25×34 level, analogous to the DIFFCOM/DIFFIND checks at 185-level.
2. **Add P015 (ARMSIGMA)** to premod output — the Armington elasticity is read in `build_premod!`
   but not stored in the `premod` output dict; needed for solver compatibility.
3. **Write aggregated output to GEMPACK HAR** — produce a `.har` file from the 25×34 aggregated dict
   so the TABLO solver can consume it directly.
4. **RAS outer loop** — `converged=false` is expected with single-pass RAS; a multi-pass loop
   (feed `DGON` from iteration N back into reg2's distance weighting) would resolve the trade/MAKE
   discrepancies flagged by DOMERR/IMPERR.
5. **Step 4 (derived parameters)** — closed-form `Formula`-computed coefficients at the aggregated
   25×34 scale: purchaser-price flows, source shares, cost shares, CES/CET calibration shares.

## Status summary (2026-07-19)

| Step | Description | Files | Status |
|------|-------------|-------|--------|
| 0 | Scaffold | `Project.toml`, module | ✅ done |
| 1 | Sets & data ingestion | `prepare_sets.jl`, `read_data.jl` | ✅ done & verified |
| 2 | Regionalization + RAS pipeline | `build_reg0!.jl`..`build_premod!.jl` (6 files) | ✅ debugged, runs end-to-end |
| 3 | Aggregation 185→25 × 34 | `aggregation_data.jl`, `aggregate_model!.jl` | ✅ done & verified |
| 4 | Derived parameters | `prepare_parameters.jl` | ❌ not started |
| 5 | Core equations (~3000 LOC TERM.TAB) | `build_model!.jl` | ❌ not started |
| 6 | Dynamic + district extensions | `build_dynamics!.jl`, `build_district!.jl` | ❌ not started |
| 7 | Closure & solve (Ipopt) | `initialize_model!.jl`..`run_model!.jl` | ❌ not started |
| 8 | Reporting | `calculate_gdp.jl` | ❌ not started |
