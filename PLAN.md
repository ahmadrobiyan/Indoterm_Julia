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

## Pipeline data-flow modifications (Step 4 enablers)

Step 4 requires the full set of TERM.TAB flow coefficients at the aggregated 25×34 scale. The following
modifications were made to existing pipeline code to carry all necessary arrays forward:

- **`build_pstras!.jl`**: Added `1PTX` (production tax, from reg1 FACT[4]) to output.
- **`build_premod!.jl`**: Added pass-through of `BSMR` (USE delivered), `UTAX` (commodity taxes),
  `2PUR` (investment by COM×IND×DST), `STOK` (inventories), `1PTX` (production tax) from pstras.
- **`aggregate_model!.jl`**: Added `_agg_use()` for BSMR/UTAX (aggregates both COM and USR-IND dims),
  `_agg_make()` for 2PUR, and standard `_agg_first()` for STOK/1PTX.

## Step 4 implementation details

`prepare_parameters!` computes all TERM.TAB Coefficient-Formula blocks at the aggregated 25×34 scale.
Implemented formula blocks by excerpt:

| Excerpt | Coefficients computed |
|---------|----------------------|
| 7 | `PUR`, `PUR_S`, `PUR_CS`, `SRCSHR`, `PUR_D`, `TAXRATE` |
| 10 | `LAB_O` |
| 11 | `PRIM`, `PRIMCOST` |
| 12 | `VARCST`, `VCST`, `VTOT`, `PTXRATE`, `COSTMAT` |
| 13 | `HOUPUR`, `HOUPUR_C`, `BUDGSHR`, `EPSH`, `EPSAVE`, `BLUX`, `SLUX`, `HOUSHR` |
| 14 | `INVEST_I`, `INVEST_C` |
| 17 | `USE_U`, `USE_I` |
| 18 | `LOCUSE`, `LOCUSE_S`, `LOCUSE_SD`, `IMPSHR` |
| 19 | `DELIVRD`, `BASSHR`, `MARSHR` |
| 20 | `DELIVRD_R` |
| 21 | `TRADMAR_CS`, `SUPPMAR_P/D/RD/R` |
| 22 | `MAKE_C`, `MAKE_I`, `MAKESHR1`, `MAKESHR2`, `MAKE_D` |
| 23 | `TRDIAG`, `TRADE_D`, `TRADE_R`, `TRADE_RD` |
| 25 | `PRIM_I`, `PRIMSHR` |
| 27 | `LAB_I`, `LAB_IO`, `LND_I`, `CAP_I`, `SLAB_I` |
| 28 | `GDPINCSUM`, `GDPINC` |
| 29 | `GDPEXPSUM`, `GDPEXP` |
| 31 | `TRADE_CR`, `IMPUSED_C`, `IMPLANDED_C` |
| 34 | `PRIM_D` |
| 35 | `NATVTOT`, `LAB_OD`, `CAP_D` |
| 40 | `ROWDEM`, `EXPSHR` |
| 41 | `CHECKA`, `CKRATA`, `CHECKB`, `CKRATB` |
| 44 | `VMAINUSE`, `LOCSHR`, `LOCSHR_R` |

**Total: 75 derived parameters → 85** (MAKE, BSMR, UTAX, 2PUR, STOK, 1PTX pass-through added for equation
functions; verified with mock data on 2026-07-21).

## Step 5 — Core equations (executing)

### Approach
- **Linearized %-change form** (all equations are linear, following TERM.TAB's own %-change convention)
- **JuMP array-constraint syntax** (`@constraint(m, [c=1:na, s=1:ns, ...], ...)`) for per-block compilation
  efficiency — one compilation per equation block instead of per individual constraint
- **Index dictionaries** (`TRADE_idx`, `USE_usc_idx`, etc.) for sparse sums over conditionally-zero flows
- **JuMP + Ipopt** (feasibility solve, no objective), matching the WayangJulia architecture already
  validated for national-level ORANI-G/TERM models

### Completed
| Module | Equations | Line count |
|--------|-----------|------------|
| `build_model!.jl` | ~168 variable arrays (4564 JuMP vars at 3×5 scale), `build_model_full!` orchestrator | ~560 |
| `build_equations.jl` | ~130 equation functions covering Excerpts 6–29 (basic prices → GDP expenditure-side), 10 index-dictionary setup functions | ~1180 |

### Equation blocks implemented
1. **Excerpt 6** — Import prices (`pimp`, `pfimp`, `phi`)
2. **Excerpt 7** — Basic + purchaser prices (`pbasic`, `ppur`, `ppur_s`, `phou`, `tuser`)
3. **Excerpt 8** — Intermediate demands (`xint`, `xint_s`, `aint_s`, `pint`)
4. **Excerpt 9** — Factor demands (`xlab`, `plab_o`, `wlab_o`, `xlab_o`, `pcap`, `plnd`, `pprim`, `xprim`, `aprim`, `alab_o`, `wprim`)
5. **Excerpt 11** — Production costs (`pvar`, `pcst`, `delPTX`, `ptot`)
6. **Excerpt 13** — Household demands (`xsub`, `xlux`, `xhouh_s_agg`, `alux`, `asub`, `wlux`, `phouhtot`, `whouhtot`, `xhoutot`, `phoutot`)
7. **Excerpt 14** — Investment demands (`xinvi`, `pinvest`, `pinvitot`)
8. **Excerpt 15** — Capital accumulation (`gret`, `xinvitot`, `finv2`)
9. **Excerpt 16** — Government demands (`xgov`, `xgov_s`, `fgovtot2`, `fgovtot3`)
10. **Excerpt 17** — Export/total demand (`pfexp`, `xexpd`, `xexp`, `xexp_s`, `xstocks`, `xint_i`, `xuse`)
11. **Excerpt 19** — Margins (`pdelivrd`, `xtradmar`, `psuppmar_p`)
12. **Excerpt 20** — Regional sourcing (`puse`, `xtrad`)
13. **Excerpt 21** — Margin supply (`xsuppmar_p`, `psuppmar_p`, `xsuppmar`, `xsuppmar_d`, `xsuppmar_rd`)
14. **Excerpt 22** — MAKE/CET (`xmake`, `xtotA_B`, `xcomA_B`, `pmake`)
15. **Excerpt 23** — Market clearing (`xtrad_d`, `xtrad_r`, `pdomA_sum`)
16. **Excerpt 25** — Factor-market closure (`pfin`, `xfina`, `xfinb`, `xfinc`, `xfind`, `wfin`)
17. **Excerpt 26** — Tax revenue (`delTAXint`, `delTAXhou`, `delTAXinv`, `delTAXgov`, `delTAXexp`)
18. **Excerpt 27** — Labour/occupation closure (`wlnd_i`, `wcap_i`, `wprim_i`, `xlnd_i`, `xcap_i`, `plab`, `plab_i`, `realwage_i`, `xlab_i`, `wlab_i`, `rlab_i`, `plab_io`, `xlab_io`, `realwage_io`, `wlab_io`, `rlab_io`)
19. **Excerpt 28** — GDP income-side (`delGDPINCa`, `delGDPINCb`, `delGDPINCc`, `delGDPINCd`, `delGDPINCe`, `wgdpinc`)
20. **Excerpt 29** — GDP expenditure-side (`delXGDPEXPa_setup`, `delXGDPEXPb`, `delXGDPEXPc`, `xgdpexp`, `delPGDPEXPa`, `delPGDPEXPb`, `delPGDPEXPc`, `pgdpexp`, `wgdpexp`, `wgdpdiff`, `xgne`, `pgne`, `wgne`, `delINDTAX`, `delBUDG1`, `delBUDG2`, `delVGDPEXP`)
21. **Excerpt 38** — Labour-market closure (`labslack`, `flab_i`, `realwage`)
22. **Excerpt 39** — Household closure (`fhou`, `natfhou`)

### Bugs found and fixed during implementation
1. **Hard-coded `nm = 9`** in `build_model_full!` — `build_model!.jl:370` had `nm = 9` regardless of actual TMAR dimensions; changed to `haskey(agg, "TMAR") ? size(agg["TMAR"], 3) : 9`.
2. **Generator variable `m` conflicted with JuMP model `m`** — `E_pdelivrd!` used `sum(MARSHR[c,s,m,r,d] * ... for m in 1:nm)` where `m` is both the JuMP model and the loop variable; changed to explicit `for` loop with `mh`.
3. **Function name mismatch** — `build_model!.jl` called `E_xtradmar!` but function was `E_xtradmar_na!`; called `E_delXGDPEXPa!` but function was `E_delXGDPEXPa_setup!`.
4. **Missing variable `rlab_io`** — referenced by `E_rlab_io!` but not declared in scaffold nor added to `vars` dict.
5. **Missing `params["MAKE"]`** — `E_xtotA_B!` and `E_xcomA_B!` expected `params["MAKE"]` but it wasn't stored by `prepare_parameters!`.
6. **Removed redundant `E_nhouh!`/`E_xhou_s!` calls** — single-household-type simplification (uses `E_xhouh_s_agg!` instead).

### Current performance (3×5 mock, 2026-07-21)
| Metric | Value |
|--------|-------|
| Variables | 4564 |
| Constraints | 3372 |
| Build time (incl. compilation) | ~46s |
| Precompilation time | ~35s |
| Runtime (after precompile) | ~50ms |

The per-equation per-iteration `@constraint` calls (for nested loops with conditionals) dominate build
time; JuMP's array-constraint syntax (`@constraint(m, [i=1:n], ...)`) is used wherever the equation
is unconditional. Remaining scalar loops are unavoidable due to conditional structure (e.g. "sum only
where coefficient > 0").

### Next steps
1. Test at full 25×34 scale with mock data (sub-second data generation, verify build completes)
2. Test at full 25×34 scale with real aggregated pipeline output
3. Add stub/missing closure equations (Excerpts 38–39 already stubbed)
4. Bind into a solvable system (fix shock variables, set up Ipopt solve)
5. Run a small shock (e.g. −5% export shift) and verify against expected signs
6. Benchmark replication (all shocks at 0 → ~0% change)
7. Price homogeneity test

## Immediate follow-up items

1. **Step 5 closure** — bind equations into solvable Ipopt system; handle remaining Excerpts 30–37, 40–55.
2. **Test at full 25×34 scale** — run `prepare_parameters!` and `build_model_full!` with real pipeline output.
3. **Write aggregated output to GEMPACK HAR** — produce a `.har` file from the 25×34 aggregated dict
   for cross-checking against the original TABLO.
4. **RAS outer loop** — `converged=false` is expected with single-pass RAS; a multi-pass loop
   (feed `DGON` from iteration N back into reg2's distance weighting) would resolve trade/MAKE
   discrepancies flagged by DOMERR/IMPERR.
5. **Add P015 (ARMSIGMA)** to premod output — the Armington elasticity is read in `build_premod!`
   but not stored; needed for solver compatibility.

## Status summary (2026-07-21)

| Step | Description | Files | Status |
|------|-------------|-------|--------|
| 0 | Scaffold | `Project.toml`, module | ✅ done |
| 1 | Sets & data ingestion | `prepare_sets.jl`, `read_data.jl` | ✅ done & verified |
| 2 | Regionalization + RAS pipeline | `build_reg0!.jl`..`build_premod!.jl` (6 files) | ✅ debugged, runs end-to-end |
| 3 | Aggregation 185→25 × 34 | `aggregation_data.jl`, `aggregate_model!.jl` | ✅ done & verified |
| 4 | Derived parameters | `prepare_parameters.jl` | ✅ done & verified (85 params) |
| 5 | Core equations (~3000 LOC TERM.TAB) | `build_model!.jl`, `build_equations.jl` | ✅ ~130 equation functions (Excerpts 6–29), builds at 3×5 (4564 vars, ~46s) |
| 5a | Scale to 25×34 | — | ✅ done, 2,436,062 vars / 1,370,212 cons |
| 5b | Ipopt solve | `initialize_model!.jl` | 🔄 base static closure wired + verified; `solve_model!.jl` not started |
| 6 | Dynamic + district extensions | `build_dynamics!.jl`, `build_district!.jl` | ❌ not started |
| 7 | Closure & solve (Ipopt) | `initialize_model!.jl`..`run_model!.jl` | ❌ not started |
| 8 | Reporting | `calculate_gdp.jl` | ❌ not started |

## Session continuation (2026-07-22)

Picked back up via Claude (VS Code extension session). Verified the state this file claims still
matches reality: Julia 1.12.6 + all four deps instantiated (`Ipopt` v1.15.0, `JuMP` v1.30.1,
`NamedArrays` v0.10.5, `HeaderArrayFile` v0.2.0); `build_model!.jl` (562 lines) and
`build_equations.jl` (1178 lines) line counts match this file's own claims exactly. Re-ran
`test/run_pipeline.jl` (Steps 0-4) in the background to reconfirm end-to-end.

**Found: Step 4/5 work is uncommitted.** `git log` shows the last commit is `de84935` ("update
PLAN.md: git repo setup, origin/, national_data.zip"); `git status` shows `build_equations.jl`,
`build_model!.jl`, `prepare_parameters.jl` as untracked and `PLAN.md`, `src/IndotermJulia.jl`,
`src/aggregate_model!.jl`, `src/build_premod!.jl`, `src/build_pstras!.jl`, `test/run_pipeline.jl` as
modified-but-unstaged — i.e. everything this file documents as "Step 4/5 done" only exists in the
working tree, not in any commit. Not committed yet pending user confirmation (nothing here has been
asked for this session).

Added `TODO.md` alongside this file as the short, checkable continuation list (this file stays the
narrative/rationale/learnings record, per its own existing style). See `TODO.md` for the concrete
next actions — the immediate ones are: build at full 25×34 scale (Step 5a, not yet attempted — only
the 3×5 mock scale has been tested), audit variable/equation coverage for gaps before wiring the
Ipopt solve, then proceed to closure wiring (Step 5b) and the verification protocol already
documented above.

**Pipeline re-run found and fixed a real bug (not just a stale-claim issue).** The background
`test/run_pipeline.jl` run **failed**: `MethodError: no method matching Float64(::Vector{Float64})` in
`prepare_parameters.jl:29`. Root cause: `P021`, the Frisch (marginal-budget-share) parameter, is
region-specific in the actual data — `build_premod!.jl:248` stores it as
`NamedArray(FRISCH_d, REG, (:DST,))`, one value per of the 34 regions — but `aggregate_model!.jl`'s
generic REG-only pass-through block (`for k in ["PO01", "P021", "EMPR", "ELWG", "SMAR"]`, line 216)
unwraps every `NamedArray` in that list to a plain `Vector`, and `prepare_parameters.jl` had assumed
`P021` was a single national scalar. Note that even the pre-existing "success" branch of that ternary
(`vec(parent(agg["P021"]))[1]`) was silently wrong in the same way — it would have taken only region
1's Frisch value and applied it to all 34 regions, rather than crashing. Fixed by keeping `P021`/
`FRISCH` as an `nr`-length vector end-to-end and indexing `FRISCH[d]` in the `BLUX`/`SLUX` loop.
Re-ran the full pipeline after the fix: Steps 0-4 now pass end-to-end, `prepare_parameters!` derives
83 parameters. **This means Step 2's "verified 2026-07-19" claim above was true for the pipeline
shape/plumbing at the time, but the newer `prepare_parameters.jl` code (Step 4, added after that
verification) had never actually been run end-to-end until now — worth remembering that "verified"
timestamps on this file only cover what existed as of that date.**

**Excerpt 49 cross-check completed — no missing core equations.** Read `TERM.TAB` Excerpts 30-49 in
full to settle the open question of whether anything beyond Excerpts 1-29/38-39 is needed for a square
core solve: **no.** Excerpts 30-37 and 40-48 are entirely national/regional aggregation, contribution
decomposition (Keller decomposition, terms-of-trade contributions), and `Assertion`/`Write ... to file
SUMMARY` diagnostics — one-way reads *from* the core solution, nothing the core equations depend on.
Cross-checked Excerpt 49's own `Substitute`/`Backsolve` variable list (GEMPACK's condensation
directives) against `build_equations.jl`'s ~127 `E_*!` functions: every name resolves once naming
mismatches and the single-household simplification are accounted for (`xtradmar` → implemented as
`E_xtradmar_na!`, matching the historical bug-fix note above; `xhou_s`/`xhouh_s` → both defined by
`E_xhouh_s_agg!`). The remaining names not found anywhere in `build_model!.jl`/`build_equations.jl`
(`contCPI, continccom, contincind_d, contMainMacro, contnatxtot, contxprim_i, xrowdem, xrowdem_d`) are
all Excerpt 30-40 reporting variables, confirmed out of scope. **Step 5's equation set is complete
against this cross-check** — remaining Step 5 work is Step 5a (full-scale build test) and Step 5b
(solve wiring), not new equations.

**Step 5a — full 25×34 scale build: done, one real bug found and fixed.** Added
`test/run_full_model.jl` to run the whole pipeline (Steps 0-4) and then `build_model_full!` at the
real aggregation scale, not the 3×5 mock used during Step 5's original implementation. First attempt
failed: `KeyError: key "INVEST_C" not found` in `build_model!.jl:49`. Traced it back through every
pipeline stage (`reg1`→`reg2`→`ras_balance!`→`build_pstras!`→`build_premod!`→`aggregate_model!`) by
comparing each stage's printed key list: `2PUR` (the investment-by-commodity-by-industry matrix)
is present in `reg1`, `reg2`, and `ras`'s output dicts (it's carried through unchanged, never
RAS-balanced), but **`build_pstras!.jl`'s output dict never included a `"2PUR"` key at all** — the
one stage that silently dropped it. Every downstream stage handled the missing key "correctly" in
isolation (`build_premod!.jl`'s `haskey(pstras,"2PUR") ? ... : nothing`, `aggregate_model!.jl`'s
`premod["2PUR"] !== nothing` guard, `prepare_parameters.jl`'s `V2PUR !== nothing` guard around the
`INVEST_I`/`INVEST_C` block) — which is exactly why this went unnoticed until a consumer
(`build_model!.jl`) read `params["INVEST_C"]` unconditionally. Fixed with a one-line pass-through
addition to `build_pstras!.jl`, mirroring the existing `BSMR`/`UTAX` pattern. After the fix:

- Pipeline: `agg keys: 30, params: 85` (up from 29/83 — the two new `INVEST_I`/`INVEST_C` keys).
- Full model build: **98.0s**, **2,376,562 variables**, **1,368,512 constraints**, 169 entries in
  the `vars` dict.
- The ~1.008M vars-minus-constraints gap is **expected, not a bug** — this is the normal pre-closure
  state of any GEMPACK-style CGE model: the full variable set always exceeds the core equation count,
  and a *closure* (Step 5b) designates exactly that many variables exogenous (fixed at their
  benchmark value) to square the system. `TERM.CMF`'s exogenous-variable list should specify close to
  1,008,050 variables; this is worth a sanity-check once `initialize_model!.jl` is written, as a large
  mismatch there would indicate a missed equation or an extra undeclared variable.

Also fixed a second, unrelated bug uncovered by the same test run: `test/run_full_model.jl` itself
called `num_constraints(m, F, S; count_variable_in_set_constraints=false)` per constraint-type pair,
but that keyword isn't accepted by that method signature in JuMP 1.30.1 — only the no-argument
`num_constraints(m; count_variable_in_set_constraints=false)` form supports it. Not a translation
bug, just a test-script API mismatch; fixed by switching to the simpler call.

**Step 5b prep — Exogenous-list cross-check found and fixed a real equation gap
(`srctwist`/`avesrctwist`).** Before writing `initialize_model!.jl`, cross-checked `TERM.CMF`'s
default-closure `Exogenous` block (48 names, lines 21-68) against all 169 JuMP variable names
declared in `build_model!.jl`. This is a deeper check than the Excerpt 49 cross-check above — that
one only verified every `Substitute`/`Backsolve` *equation name* has a corresponding `E_*!`
function; this one checks that every variable a closure or equation actually *needs* is modeled at
all. 7 of the 48 names came up unmatched. Traced each by grepping `TERM.TAB` for the name and
checking which `! Excerpt N of TABLO input file: !` marker it falls under: `delfwage_o` (line
2906), `delUnity` (2696), `emptrend` (2878), `frnorm`/`frnorm_id` (2724-2725), and `gtrend` (2732)
all land inside Excerpts 50-54 (2664-2920) — confirmed Step 6 dynamic scope (investment rule,
real-wage adjustment), correctly out of Step 5. `srctwist` did not — it traced to Excerpt 20 (line
986), core scope.

Read `TERM.TAB` lines 960-1005 to understand why: it contains *two* candidate `E_xtrad`
formulations back to back, which looks like a duplicate-declaration problem until you track TABLO's
`!...!` comment delimiters character-by-character rather than line-by-line. The first formulation
(lines 992-996) uses `srctwist(c,s,r,d)`/`avesrctwist(c,s,d)` (declared as `Variable`s at 985-987,
with their own defining equation `E_avesrctwist` at 988-990) and the `SIGMADOMDOM(c)` CES
elasticity, ending in the real `Substitute xtrad using E_xtrad;` directive — this is live, compiled
code. The second formulation (lines 998-1005, introducing a `twistsrc(i,s,k)` variable instead)
looked equally live at a glance, but line 998 opens a TABLO comment (`! alternative form with
"twists":`) with no closing `!` on that line — the comment doesn't actually close until the lone
trailing `!` at the very end of line 1005. Everything in between, including the `twistsrc` variable
declaration and the second `Equation E_xtrad`, is dead documentation that was never compiled. So the
first formulation is unambiguously canonical.

Compared against it, `build_equations.jl`'s `E_xtrad!` (`xtrad-atrad == xuse-(pdelivrd+atrad-puse)`)
was missing three things: the `srctwist`/`avesrctwist` shift terms (and `E_avesrctwist` didn't exist
at all — neither variable was referenced anywhere in the Julia code); and the `SIGMADOMDOM(c)`
elasticity coefficient, implicitly using 1.0 in its place. The elasticity's aggregated value
(`SGDD`) turned out to already be computed correctly by `build_premod!.jl`/`aggregate_model!.jl` and
read into a local `SGDD` binding at `prepare_parameters.jl:26` — but, like the `P021` and `2PUR`
bugs earlier this session, never added to that function's output dict, so `params["SGDD"]` didn't
exist. This class of bug (data computed correctly upstream, silently dropped by a missing dict-key
write, masked by a downstream `haskey`/`!== nothing` guard) has now recurred three times
(`P021`/FRISCH, `2PUR`/investment, `SGDD`/SIGMADOMDOM) — worth treating as a standing suspicion
whenever a parameter that should obviously matter seems to have no effect. Grepping for the same
pattern turned up **five more instances**: `SLAB`, `P028`, `SMAR`, `PO01`, `SCET`, `P018` are all
read into local bindings in `prepare_parameters.jl` (lines 24-38) and never written to its output
dict either — confirmed dead by grep (each name appears exactly once in the file, at its own
declaration). `build_model!.jl` silently compensates with hardcoded placeholder elasticities
(`sigmalab`/`sigmaprim`/`sigmaout` all `fill(0.5, na)`, `sigmadomimp = fill(5.0, na)`) instead of the
real, data-derived values. Logged as a follow-up in TODO.md's "Smaller loose ends" section (folded
into the pre-existing `P015`/`ARMSIGMA` item, since it's the same fix shape) — lower urgency than
`srctwist` because it changes solved *magnitudes* on a real shock, not the variable/constraint
*count*, so it won't be caught by the benchmark-replication (all-zero-shock) test.

Fixed the confirmed gap: added `p["SGDD"] = SGDD` to `prepare_parameters.jl`; added
`srctwist[1:na,1:ns,1:nr,1:nr]` and `avesrctwist[1:na,1:ns,1:nr]` `@variable`s to `build_model!.jl`;
added `E_avesrctwist!` to `build_equations.jl` (mirrors `E_puse!`'s
`DELIVRD_R[c,s,d] * lhs == sum_r(DELIVRD[c,s,r,d] * rhs)` pattern, i.e. TABLO's `ID01(x)*lhs = ...`
convention is just `x` used directly as the coefficient); rewrote `E_xtrad!` to include the
`srctwist`/`avesrctwist` terms and the `SGDD[c]` coefficient. `srctwist` itself is left with no
defining equation on purpose — `TERM.CMF` lists it `Exogenous` in the base closure, so it's meant to
be fixed (at 0, for a no-shock benchmark) by Step 5b's `initialize_model!.jl`, exactly like any other
closure-exogenous variable; `avesrctwist` is fully determined by `E_avesrctwist` given `srctwist`.

Verified with a full pipeline + full model build re-run (after also killing a stale julia.exe from
an earlier background run in this session that was still resident at ~1.3GB on this 7.8GB-RAM
machine and had caused a transient `OutOfMemoryError` on the first retry attempt — unrelated to the
code change, confirmed by the crash site being an unrelated pre-existing variable declaration
several lines away from anything just edited). Second attempt succeeded: **61.6s build time,
2,436,062 variables (+59,500 = 25×2×34×34 + 25×2×34, exactly the new `srctwist`+`avesrctwist`
sizes), 1,370,212 constraints (+1,700 = 25×2×34, i.e. every `E_avesrctwist` instance active, none
skipped by the `DELIVRD_R > 0` guard), 171 `vars` dict entries (+2)** — matches hand-calculated
expectations exactly, confirming the fix is wired correctly end to end.

## Step 5b — base static closure + touched/untouched diagnostic (session continuation)

**`initialize_model!.jl` written and wired in.** Encodes `TERM.CMF`'s default `Exogenous` list (48
names) as a fix/free split over the JuMP variables from `build_model!`/`build_model_full!`:
`BASE_CLOSURE_SCALARS` (10 names) and `BASE_CLOSURE_ARRAYS` (32 names) get `JuMP.fix(...; force=true)`
at a shock value (default 0.0 — benchmark replication); everything else stays free for Ipopt.
`DYNAMIC_ONLY_CLOSURE_NAMES` (6 names: `delfwage_o`, `delUnity`, `emptrend`, `frnorm`, `frnorm_id`,
`gtrend`) documents the Excerpt 50-54 names that aren't modeled yet, so the omission is greppable
rather than silent. Deliberately does not implement `TERM.CMF`'s numéraire swap
(`phi = Natmacro("GDPPI")`, needs Excerpt 30-40 reporting that's out of Step 5's scope) — `phi` stays
fixed as the default numéraire. Every declared variable also gets a `0.0` start value, the correct
guess for a %-change-at-benchmark formulation. Wired into `IndotermJulia.jl`'s include/export list;
`test/run_full_model.jl` extended to call it and report fixed/free counts.

**Built a touched/untouched constraint-term introspection diagnostic to verify the closure is sized
correctly** (`test/diagnose_gap.jl`): walks `list_of_constraint_types(m)` filtered to affine
constraints, collects every `VariableRef` referenced by any constraint's `.func.terms`, and reports
any *free* (non-fixed) variable that never appears — i.e. has no defining or using equation at all.
This is a stronger check than a raw vars-minus-constraints count, since it locates *which* variables
are the problem, not just how many. First run: **13,877 untouched free variables** — `wlab_o=>850,
xsuppmar_d=>10404, psuppmar_p=>1258, xsuppmar_rd=>306, xtrad_d=>162, fgret=>850, fhou2=>34,
plab_id=>4, wlab_id=>4, rlab_id=>4, natfhou=>1`. None of this was expected from the closure alone, so
each entry needed individual root-causing — this uncovered two more real bugs, described below.

**Bug: `LAB_O` identically zero — traced to a leftover placeholder stub, not a data gap.**
`wlab_o=>850` traced to `E_wlab_o!`/`E_plab_o!`/`E_wprim!`'s shared `LAB_O[i,d] > 1e-10` guard never
passing. `params["LAB_O"]` derives from `agg["1LAB"]`, which derives from `build_pstras!.jl`'s
`V1LAB_iod` — which had literally been left as `for o in 1:NO; V1LAB_iod[i,o,d] = 0.0; end`, a
"We'll recompute below" comment that was never acted on. The correct computation needs `reg0`'s
national labour-occupation shares (`OSHR`, IND×OCC) applied to `reg1`'s regional labour factor
payment (`FAC_r[i,g_lab,d]`), but `OSHR` had never been added to `reg1`'s own output dict either.
**Fixed**: added `"OSHR" => reg0["OSHR"]` to `build_reg1!.jl`'s output dict; replaced the zero-stub
in `build_pstras!.jl` with `V1LAB_iod[i,o,d] = FAC_r[i,g_lab,d] * OSHR[i,o]`. Re-verified: untouched
count 13,877 → 13,027 (−850, `wlab_o` fully gone), every other count unchanged — confirms no side
effects on the rest of the model.

**Bug (by far the most severe found in this project so far): 10 of 11 lookup-dict `*_setup!`
functions in `build_equations.jl` were defined and exported but never called anywhere.** Many
equations depend on a module-level `const XXX_idx = Dict{Tuple{...},Float64}()` that a dedicated
`*_setup!()` function must populate *before* the equations that read it via `get(dict, key, 0.0)` are
built — the equations never call their own setup function. Grepping every call site in
`build_model!.jl` (the only place any `E_*!`/`*_setup!` function is invoked) found only one of these
11 functions ever called (`E_delXGDPEXPa_setup!`, which builds constraints directly and doesn't share
this shape). The other 10 — `E_plab_o_setup!`, `INVEST_setup!`, `USE_IS_setup!`, `USE_usc_setup!`,
`TRADMAR_setup!`, `SUPPMAR_setup!`, `SUPPMAR_D_setup!`, `TRADE_setup!`, `PUR_src_setup!`,
`TAX_PUR_setup!`, `STOCKS_setup!` — were dead code, meaning every equation reading their dicts
(`E_plab_o!`, `E_wlab_o!`, `E_wprim!`, `E_pinvitot!`, `E_xint_i!`, `E_xuse!`, `E_xsuppmar_p!`,
`E_psuppmar_p!`, `E_xsuppmar_d!`, `E_xsuppmar_rd!`, `E_xtrad_d!`, `E_xtrad_r!`, `E_xfind!`,
`E_delTAXint!`, `E_delTAXhou!`, `E_delTAXinv!`, `E_delXGDPEXPb!`, `E_delPGDPEXPb!`) had been silently
computing over an empty dict the whole time — a systemic silent-zero bug spanning labour prices,
margin supply, market clearing, investment allocation, tax revenue, and GDP-expenditure
decomposition, much wider-reaching than any single previous bug (`P021`, `2PUR`, `SGDD`, `LAB_O`).

This directly interacted with the `LAB_O` fix above: making `LAB_O` nonzero let
`E_plab_o!`/`E_wlab_o!`/`E_wprim!`'s guards start passing, but since `V1LAB_idx` (one of the 10 dead
dicts) was still empty, those equations began wrongly forcing `plab_o`/`wlab_o`/`wprim` to exactly 0
instead of being simply absent — "touched but wrongly valued" is invisible to the touched/untouched
diagnostic, so this had to be caught by reasoning about the interaction, not by re-running the
diagnostic alone. Also found in the course of the same investigation: `E_xinv_s!` (defines `xinv_s`,
consumed by `E_xfinb!`) was never called anywhere in `build_model_full!` at all — a different bug
shape (a wholesale missing equation-block call, not a missing setup call).

**Fixed**: read `prepare_parameters.jl` in full to find the correct already-computed source array for
each dict — `params["USE"]`, `params["TAX"]`, `params["PUR"]`, `params["INVEST"]`,
`params["STOCKS"]`, plus `TMAR`/`MARS` newly extracted from `agg` — and added all 8 remaining
`*_setup!()` calls (`E_plab_o_setup!` included) to `build_model_full!` in `build_model!.jl`, right
before the Excerpt 6 block, guarded with `haskey(params, ...)` where the source itself is
conditionally populated (`USE`/`TAX`/`PUR` only exist if `agg["BSMR"]`/`agg["UTAX"]` were present).
Added the missing `E_xinv_s!(m, vars, na, nr, params)` call to the Excerpt 14 block. Deliberately did
**not** resurrect `build_model!.jl`'s own pre-existing local `PUR_d = USE_d` computation (missing a
`+ TAX_d` term and never actually used by anything) — used the correct, already-computed
`params["PUR"]` instead.

**Bug: `STOK` (stock-change) pass-through gap, found during the same investigation.**
`params["STOCKS"]` (feeds `STOCKS_idx`, `GDPEXPSUM`, `CHECKA`/`CKRATA`,
`E_delXGDPEXPb!`/`E_delPGDPEXPb!`) was silently defaulting to `zeros(T,na,nr)` via
`prepare_parameters.jl`'s `agg["STOK"] !== nothing` guard, even though `reg1["STOK"]`/`reg2["STOK"]`
(regional industry stock-change data) existed correctly upstream — `build_pstras!.jl`'s output dict
simply never included a `"STOK"` key, the same silent-drop shape as the earlier `2PUR`/`UTAX` bugs.
**Fixed** by adding `"STOK" => haskey(reg1, "STOK") ? reg1["STOK"] : nothing,` to `build_pstras!.jl`'s
output dict, matching the existing `UTAX`/`2PUR` pass-through lines.

**Final verification (2026-07-22).** Re-ran `test/diagnose_gap.jl` with all three fixes (8 setup!
calls + `TMAR`/`MARS` extraction, `E_xinv_s!` call, `STOK` pass-through) applied together: untouched
free vars dropped from **13,027 → 3,575** (−9,452). `xsuppmar_rd` (previously 306) is completely
gone; `xsuppmar_d` fell from 10,404 → 1,258, now matching `psuppmar_p`'s count exactly — consistent
with both sharing the same sparse margin-flow index set (real zero-valued combinations in the 25×34
data, not a bug). Remaining breakdown (`xsuppmar_d=>1258, psuppmar_p=>1258, fgret=>850, xtrad_d=>162,
fhou2=>34, plab_id=>4, wlab_id=>4, rlab_id=>4, natfhou=>1`) is not yet individually root-caused but
looks like genuine data sparsity rather than new missing wiring — none of these changed shape when
the setup calls were added.

Also wrote `test/check_setup_dicts.jl` to numerically confirm every one of the 12 lookup dicts is
non-empty with a real nonzero sum after `build_model_full!` — e.g. `V1LAB_idx` (3,400 entries, 2,686
nonzero, sum 6.43e6, previously all zero), `TRADMAR_idx`/`SUPPMAR_idx`/`TRADE_idx`/`PUR_src_idx`/
`TAX_idx`/`INVEST_C_idx`/`USE_IS_idx`/`USE_usc_idx`/`STOCKS_idx` all populated similarly — and that
`agg["STOK"]` and `params["STOCKS"]` carry matching nonzero sums (33,956.96) end-to-end. This
confirms `plab_o`/`wlab_o`/`wprim` are no longer wrongly forced to zero, closing the gap the
touched/untouched diagnostic can't see on its own ("touched but wrongly valued").

**Recurring pattern, now the ninth documented instance across this project**: data computed or
available correctly at one pipeline stage, silently dropped by a missing dict-key write or missing
function call at the very next stage, masked by a downstream `haskey`/`!== nothing` guard that
degrades gracefully instead of erroring (`P021`/FRISCH, `2PUR`/investment, `SGDD`/SIGMADOMDOM,
`LAB_O`/OSHR, `STOK`/STOCKS, and now the ten `*_setup!` dicts). Worth treating as a standing
suspicion any time a variable seems to have no effect on a shock, or shows up fully/partially
untouched in the diagnostic, even when its own equation function looks correct in isolation.

**Bug (tenth+ instance, a new sub-shape): `E_ggro!`/`E_fgret!` (Excerpt 15) never called at all.**
Continuing to work through the remaining 3,575 untouched vars after the section above, `fgret=>850`
stood out: 850 = `na*nr` exactly, i.e. *100%* of the array was untouched, not a sparse subset — a
strong signal for "missing equation," not "genuine data sparsity." Grep confirmed `fgret` had zero
defining equations in `build_equations.jl`. Reading `TERM.TAB` lines 738-777 (with its own comment
that "normally `capslack` is exogenous and `fgret` endogenous," confirming the intended closure
direction) showed `ggro` *also* lacked its own equation — it only appeared as an input term inside
`E_xinvitot!`'s constraint, making it "touched" by the diagnostic even though it was underdetermined
by exactly `na*nr` equations. This is the same masking shape noted above for "touched but wrongly
valued," just for "touched but underdetermined" instead — the touched/untouched diagnostic alone
cannot tell these apart from a correctly-defined variable; only checking the equation count against
TABLO's source caught it. **Fixed**: added `E_ggro!` (`ggro[i,d] == finv1[i,d] + 0.33*(2*gret[i,d] -
invslack)`) and `E_fgret!` (`gret[i,d] == fgret[i,d] + capslack`) to `build_equations.jl`, wired both
into the Excerpt 15 block in `build_model!.jl` right after the existing `E_xinvitot!` call. All
referenced variables (`gret`, `ggro`, `finv1`, `invslack`, `fgret`, `capslack`) already existed as
JuMP variables from the original variable-declaration pass, so no new `@variable` calls were needed.
Re-ran `test/diagnose_gap.jl`: untouched dropped 3,575 → 2,725 (−850, exactly `na*nr`), `fgret`
completely gone from the breakdown, every other count unchanged.

**Bug: missing "_id" labour-aggregate family (Excerpt 27) — `E_plab_id!`, `E_realwage_id!`,
`E_xlab_id!`, `E_wlab_id!`, `E_rlab_id!` never called at all, plus a missing derived parameter.**
Same size-matching heuristic: `plab_id=>4, wlab_id=>4, rlab_id=>4` in the post-`ggro`/`fgret`
breakdown each equaled `no` (4 labour occupations) exactly — 100% of each array, not a sparse subset.
Grep confirmed the entire 5-equation "_id" family was missing from `build_equations.jl`.
`xlab_id`/`realwage_id` didn't appear in the untouched list at all, but only because they're
*consumed* (not defined) by an unrelated, pre-existing constraint at `build_equations.jl:1173-1174`
(`realwage_id[o] == 2.0*xlab_id[o] + flabsup_id[o]`, part of some other Excerpt-38-area block) — the
same "touched but not well-defined" masking as `ggro` above, just one layer more hidden since it
never shows up as a *count* at all. Reading `TERM.TAB` lines 1300-1370 confirmed the TABLO source
(`SLAB_ID(o,d)`-weighted sums over region `d`, mirroring the already-implemented `_i` family's sums
over industry `i`) and, cross-referencing `prepare_parameters.jl`, that the required coefficient
(`LAB_ID`/`SLAB_ID`) had never been computed at all — only the sibling `LAB_I`/`SLAB_I`
(per-industry) and `LAB_IO` (per-region) existed. **Fixed** in two steps: (1) added
`p["LAB_ID"]`/`p["SLAB_ID"]` to `prepare_parameters.jl` right after the existing `SLAB_I_arr` block
(`LAB_ID[o] = sum_d LAB_I[o,d]`; `SLAB_ID[o,d] = LAB_I[o,d]/LAB_ID[o]`, matching `LAB_I`/`SLAB_I`'s
own share-normalization pattern); (2) added the 5 equation functions to `build_equations.jl`,
mirroring the existing `E_plab_i!`/`E_realwage_i!`/`E_xlab_i!`/`E_wlab_i!`/`E_rlab_i!` functions
exactly but summing over region `d` with the `SLAB_ID` weight instead of summing over industry `i`
with `SLAB_I`/`V1LAB_idx` — since `SLAB_ID` is already a normalized share (sums to 1 over `d` for
fixed `o`), no extra division was needed in the constraint itself, unlike the `_i` family's
`lab_sum * lhs == sum(...)` form. Wired all 5 calls into `build_model_full!` in `build_model!.jl`,
placed between the existing `_i`-family and `_io`-family call blocks. All referenced variables
(`plab_id`, `xlab_id`, `wlab_id`, `rlab_id`, `realwage_id`) already existed as JuMP variables and were
already present in the `vars` dict — no new declarations needed. Re-ran `test/diagnose_gap.jl`:
untouched dropped 2,725 → 2,713 (−12, exactly `3*no`), `plab_id`/`wlab_id`/`rlab_id` completely gone
from the breakdown, every other count unchanged.

**Recurring pattern, now confirmed for an eleventh and twelfth time.** Both bugs above share a
distinct signature worth naming explicitly: when a variable's own defining equation is missing but
the variable is still referenced as an *input* to some other, unrelated equation, the touched/
untouched diagnostic reports it as fine — "touched" only means "appears somewhere in some
constraint," not "has its own equation." The reliable tell in both cases was a count that exactly
matched an index-set size (`na*nr` for `fgret`, `no` for `plab_id`/`wlab_id`/`rlab_id`) — i.e. *100%*
of an array was untouched, versus the genuine-data-sparsity entries (`xsuppmar_d`, `psuppmar_p`, etc.)
which are all partial counts within a larger, only-partially-touched array. This 100%-vs-partial
distinction is now the standard first check whenever a new untouched-variable entry needs
root-causing: compute what `na*nr`, `no*nr`, `no`, etc. would be for the plausible index sets, and if
the count matches exactly, assume missing equation before assuming data sparsity.

**Closing out the remaining 2,713 untouched vars, applying that same 100%-vs-partial test to each
entry.** `fhou2=>34` matched `nr` exactly — a thirteenth instance of the same bug shape. Grep
confirmed `build_equations.jl` had `E_fhou!` and `E_natfhou!` but no `E_fhou2!` at all, despite
`build_model!.jl` declaring `fhou2[1:nr]` as a real variable. `TERM.TAB` lines 2038-2053 (Excerpt 39)
show `E_fhou`/`E_fhou2` as a genuine TABLO equation chain — both share `whouhtot(h,d)` on the LHS,
which looks like over-determination until you notice `whouhtot` is already pinned by the pre-existing
`E_whouhtot!` (`whouhtot == phouhtot + xhoutot`, the nominal = price × real identity). Given that,
`E_fhou!` actually solves for `fhou` (propensity to consume from labour income) and `E_fhou2!` solves
for `fhou2` (propensity to consume from regional GDP) — each equation in the chain defining exactly
one new variable from ones already pinned by earlier equations, the standard GEMPACK/TABLO idiom.
**Fixed**: added `E_fhou2!` (`whouhtot[d] == wgdpexp[d] + fhou2[d] + houslack`, an exact mirror of
`E_fhou!` with `wgdpexp` swapped in for `wlab_io`) to `build_equations.jl`, wired the call into
`build_model_full!` right after `E_fhou!`. Re-ran `test/diagnose_gap.jl`: untouched dropped 2,713 →
2,679 (−34, exactly `nr`), `fhou2` completely gone from the breakdown.

**`natfhou=>1` is different: a deliberate, documented stub, not a bug.** `E_natfhou!` exists and is
called, but its body is completely empty (`function E_natfhou!(...) end`). Its TABLO source
(`NatMacro("NomHou") = natfhou + NatMacro("NomGDPexp")`) needs the `MainMacro`/`NatMacro` national
aggregation layer from Excerpts 30-40 — the same reporting/decomposition machinery already read in
full and confirmed out of Step 5's core-square-system scope (nothing in Excerpts 1-29/38-39 depends
on it for its own solution; it's a one-way read *from* the solved model, same category as the
un-ported `CHECKC`/`D`/`E` diagnostics). Left as an intentional, greppable stub rather than building
out the national-aggregation layer just for one diagnostic ratio — revisit only alongside Step 8
reporting or a future Excerpt 30-40 pass.

**The other three (`xtrad_d=>162`, `xsuppmar_d=>1258`, `psuppmar_p=>1258`) are confirmed genuine data
sparsity, not bugs.** `E_xtrad_d!` guards on `TRADE_D[c,s,r] > 1e-10`; the full array is
`na*ns*nr`=1,700, so 162 zero-valued commodity/source/region combinations is a small fraction of a
mostly-populated array — the opposite signature from every real bug found this session, all of which
were 100%-of-array. `xsuppmar_d`/`psuppmar_p` share the same sparse margin-flow index set and guard
shape, and neither changed count across three consecutive rounds of unrelated fixes — the kind of
stability that would be surprising for a live wiring bug still lurking underneath.

**Final state after this round: 2,679 untouched free vars**, of which only `natfhou=>1` remains as a
knowingly-deferred gap; everything else is either fixed or confirmed as real, expected data sparsity
in the underlying 25×34 dataset. This closes out the touched/untouched diagnostic's todo item from
Step 5b for now — remaining Step 5 work moves to `solve_model!.jl` and the verification protocol.
