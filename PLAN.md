# IndotermJulia — Implementation Plan

> **Session breadcrumb:** Originally authored by Claude (Anthropic). Steps 2-3 debugged,
> implemented, and verified by opencode (continuation session, 2026-07-19). Git repo set up
> 2026-07-19 at `github.com/ahmadrobiyan/Indoterm_Julia`.

Julia translation of **INDOTERM** (a TERM-family, ORANI-G-derived, multi-region CGE model of
Indonesia — `TERM.TAB`, "fast multi-region model designed by Mark Horridge, 2002-6", with a 2013
dynamic extension), built as a new sibling Julia project `IndotermJulia/` at the root of
`INDOTERM CGE_2016_New/`, alongside — and never modifying — the GEMPACK TABLO source it translates.

## Milestones & gates (read this first)

Each row is a go/no-go gate, not just a task description — "done" means the gate's check actually
passed, not that code was written. Detailed narrative for every row lives further down under the
matching `## Step N` heading; this table exists so "where are we" never requires reading the whole
file. **Update this table whenever a gate's status changes** — it is the single source of truth for
project state, everything below it is supporting narrative/evidence.

| # | Milestone | Gate (what "done" requires) | Status |
|---|-----------|------------------------------|--------|
| 0 | Project scaffold | `Project.toml` + module skeleton loads | ✅ done |
| 1 | Sets & data ingestion | `prepare_sets.jl`/`read_data.jl` load all `.har` headers, sets cross-checked against source CSVs | ✅ done & verified |
| 2 | Regionalization + RAS pipeline | `build_reg0!`..`build_premod!` run end-to-end, DIFFIND/DIFFCOM sign checks pass | ✅ done & verified |
| 3 | Aggregation 185→25×34 | `aggregate_model!` produces correctly-shaped 25×34 database | ✅ done & verified |
| 4 | Derived parameters | `prepare_parameters!` computes all CES/CET calibration shares at 25×34 scale | ✅ done & verified |
| 5 | Core equations — **levels form** (not %-change, see Course correction below) | Every excerpt (6-29, 27/38-39) converted to true levels; full model builds at 25×34 scale with no `UndefVarError`/shape errors | ✅ done — builds end-to-end, 2,436,077 vars / 1,396,876 cons |
| 5b | Base static closure | `initialize_model!.jl` fix/free split matches `TERM.CMF`'s default Exogenous list; touched/untouched diagnostic shows no unintended orphan variables | ✅ done (`test/diagnose_gap.jl`; the `natfhou` stub it left behind was closed 2026-07-31 — zero orphans now) |
| 5c | **Benchmark solve** (the "one final verify") | Solver converges on the zero-shock benchmark at 25×6, then replicates it | ✅ **equations VERIFIED** (2026-07-25): benchmark relative residual **1.86e-9**. Solver switched **Ipopt → sparse Newton** (Ipopt cannot solve a square system). See detail below |
| 5c-i | Make gate 5c self-executing | `test/solve_benchmark_6reg.jl` passes unaided (squaring + Newton moved from scratchpad into `src/`) | ✅ **done (2026-07-25)** — `src/solve_newton!.jl`; test passes end-to-end: square 88,599=88,599, converged 1 iter, scaled `‖F‖∞`=1.11e-9, max rel Δ vs benchmark 5.4e-7 |
| 5c-ii | Remove residual rank deficiency | Equilibrated Jacobian is full rank | ✅ **done (2026-07-25)** — deficiency **637 → 0** via zero-flow guards on `ppur_s`/`psuppmar_p`/`plnd`. Newton converges in **1 iter / 2.2 s** |
| 6 | Dynamic extension | Excerpts 50-54 (`build_dynamics!.jl`) built and solving | ✅ **done (2026-07-25)** — square 89,911=89,911, converges 1 iter, benchmark replication unchanged (rel Δ 5.42e-7) |
| 6b | District extension | Excerpt 55 (`build_district!.jl`) | ⛔ **N/A for this dataset** — the five required headers (`RGN`/`MRGN`/`SREG`/`LOCI`/`MVTO`) exist in none of the three `.har` files; INDOTERM's finest spatial unit is the 34 provinces. Deliberately not stubbed — see TODO.md |
| 6c | Macro reporting (Excerpts 31-32) | `MainMacro`/`NatMacro` built; defines `NatMacro("GDPPI")` for TERM.CMF's numeraire swap | ✅ **done (2026-07-25)** — `src/build_macros!.jl`; square 90,460=90,460, converges 1 iter, benchmark replication unchanged (5.42e-7) |
| 7 | Closure & solve wiring beyond base closure | `run_model!.jl` homotopy driver for shocks too large for a direct solve | ✅ **done (2026-07-29, §N.19e)** — TERM.CMF's reference simulation runs end to end under the swapped closure: `delUnity` 0→1 in 3 steps, then `blabnat` 1→0.97 by arclength in **22 steps / 0 rejections**, ‖F‖∞ = 6.73e-10. The old "CGNR stall" note was superseded by §F/§I and the CGNR path is gone |
| 7a | Pseudo-arclength continuation | `arclength_solve!` traces a branch through a limit point; every accepted point is a genuine solution | ✅ **done (2026-07-29, §N.19/§N.19c)** — `src/arclength.jl` |
| 7c | Scenario specified in ONE declarative place | Every `Shock` line in the source `.CMF` has a counterpart carrying its source line number; a missing shock is a visible absence | ✅ **done (2026-07-31)** — `src/scenarios.jl` is exported and `run_model!(agg, params, sc)` consumes a `Scenario` directly. Three are defined: `TERM_CMF_REFERENCE`, `TERM_CMF_NO_DELUNITY`, `COALPRICE_REFERENCE`. The struct also carries the **numeraire**, so a call site cannot silently pick the wrong one |
| 9 | Price homogeneity of degree zero | Scale numeraire `NatMacro("GDPPI")` by λ; prices scale by λ, real quantities invariant | ✅ **done (2026-07-31, §N.20)** — root cause was ours: `E_ppur_s!` pinned `ppur_s = 1.0` on 26 live dust-flow cells, which breaks degree-zero homogeneity by construction. Anchoring to benchmark source shares `SRCSHR` instead took the two worst cells from 3.7e-4/4.1e-4 to **2.4e-15/2.2e-15**, benchmark unregressed at 8.1e-10 |
| 10 | **External validation** — the port against someone else's numbers | Reproduce `origin/coalprice.CMF` + `draftreport.pdf` Table 2's national row in sign and magnitude | ✅ **done (2026-07-31)** — `test/verify_coalprice_reference.jl`. All 8 columns match in sign and within 0.5 pp; Dutch-disease wedge reproduced; coal-sector employment **+15.07 vs +14.99**. See `VV_PLAN.md` V9. **This is the only gate in the project that can distinguish "faithful translation of TERM" from "self-consistent model that is not TERM"** |

| 8 | Reporting | `calculate_gdp.jl` GDP income/expenditure decomposition + GDP-both-sides check | ✅ **done (2026-07-25)** — `src/calculate_gdp.jl` + `test/gdp_report_6reg.jl` pass. Model-side check `wgdpdiff = 1.2e-12` ✓. Levels check reports a **benchmark-data** imbalance (see below) |

**All static-model gates completed (2026-07-25).** The rank deficiency (637→0) is resolved,
dynamics and macros are verified. The `ras_balance!` performance regression is resolved and the
pipeline runs end-to-end in ~60s.

> **Status as of 2026-07-31 — the shock blockage is long resolved and the model is now
> externally validated.** The paragraph this note replaces said "a real shock has not yet been
> run"; that was true on 2026-07-25 and stopped being true on 2026-07-29 (gate 7). The current
> state:
>
> - Shocks run. TERM.CMF's reference simulation completes (§N.19e), and the published
>   `coalprice.CMF` simulation completes in **3 homotopy steps with 0 rejections**.
> - **The port reproduces a published GEMPACK result** (gate 10 above). That is the one thing
>   no amount of internal testing could establish.
> - Two model-equation defects were found and fixed on 2026-07-31, both of which had been
>   invisible for the whole project because each was masked by a defensive measure that
>   produced the right answer under every closure tried until then:
>   `E_ppur_s!`'s pinned price (gate 9) and `E_natfhou!`'s missing equation (gate 5b).
>   **Both were found by building a closure that had never been built before, not by any test.**
>   That is the strongest argument for running the remaining V&V gates in `VV_PLAN.md`.
>
> Remaining work is verification breadth (V1–V8) and the Excerpt 49 condensation needed to
> scale from 6 regions to 34 — not correctness blockers.

### Gate #5c RESOLVED (2026-07-25) — equations verified, solver replaced

Supersedes the 2026-07-23 "paused / MUMPS wall" detail retained below for history. Read this first.

**1. The equations are correct — proven twice, independently.**
- `test/check_residual_6reg.jl` evaluates every constraint at the benchmark point with no solver
  call: only **28 of 89,301** exceed 1e-6, worst **1.8e-4** — on equations carrying ~1e6
  coefficients, i.e. ~1e-10 relative. Pure floating-point roundoff.
- After row equilibration, benchmark `||F||_inf = **1.86e-9**` (relative).

  Ipopt's thrashing was **never** an equation bug. When a solve misbehaves, run the residual check
  first — it is fast, needs no solver, and settles correctness immediately.

**2. The system was structurally UNDER-DETERMINED; that is what stalled Ipopt.** Free variables ≠
equality constraints ⟹ non-unique KKT multipliers ⟹ `inf_du` → ~1e12. Localised with a
Hopcroft-Karp **bipartite maximum matching** (constraints × free vars; unmatched free vars name the
under-determined families, unmatched constraints the redundant rows — a Dulmage-Mendelsohn-style
diagnostic). Closure sequence:

  | Step | Free vs equations |
  |---|---|
  | start | 96,892 vs 89,301 → deficit **7,591** |
  | + `OMIT_NAMES` fix | → **1,597** |
  | + delTAX zero-flow guards | → **49** |
  | + pin 157 unmatched, delete 108 dead rows | **89,193 = 89,193 (square)** |

  Two real bug classes, both now fixed in `src/`:
  - **GEMPACK `Omit` left endogenous.** `Omit` DELETES a variable — it has no equation and sits at
    its neutral benchmark. `TERM.TAB:830` (`fxstocks`) and `:2575` (`ahou_s asuppmar bint_s`) were
    all left free. Fixed via `OMIT_NAMES` in `initialize_model!.jl` (ratio-type → 1.0, additive →
    0.0; never shockable — an `Omit` is structural, not a closure choice).
  - **Guarded zero-flow skips.** Blocks emitting `@constraint` only where a benchmark flow is
    non-zero (`PUR > 1e-10 || continue`) while the variable is declared *densely* leave every
    skipped cell an unknown with no equation. Fixed for `E_delTAX{int,hou,inv,gov,exp}!` by
    converting the guard to `if/else` and pinning the skipped cell at 0 (no tax base ⟹ no revenue
    change, for any shock).

  Note the 108 unmatched constraints were **not** Walras redundancy: their free-variable fingerprint
  is *empty* — they are all-fixed `0 == 0` rows created by the fixes themselves. `dependency_detector
  = "mumps"` does **not** drop them; they need explicit `JuMP.delete`.

**3. Ipopt cannot solve this model — supersedes the "Solver decision" section below.** Given the
exactly-square system, Ipopt **hung at `iter 0` for 56 minutes** at 5.6 GB without emitting a single
iteration line. Zero degrees of freedom is degenerate for an interior-point *optimizer*. Recognition
signature: log stops after `iter 0` and Ipopt prints `Too few degrees of freedom (n_x, n_c)`. Do not
wait it out. A square `F(x)=0` needs a square solver — the squaring work is the **prerequisite** for
Newton, not the cause of the trouble.

**4. Replacement: sparse Newton on `F(x)=0`.** Full build + solve ≈ **30 s** (vs. Ipopt's hang). No
new dependencies or licences, levels formulation untouched. `F` and the sparse Jacobian come from
the JuMP model via MOI's nonlinear evaluator — build `MOI.Nonlinear.Model()`, add each constraint
with `MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(rhs))` (accepts
affine/quadratic/nonlinear alike), then `MOI.Nonlinear.Evaluator(..., SparseReverseMode(), allidx)`.
Build the evaluator over **all** variables and select free *columns* afterwards — passing only free
vars errors, since fixed vars still appear in the expressions.

**5. Equation scaling is mandatory — the hidden root cause of most symptoms.** Jacobian row scales
span **2.7e-10 to 1.7e12** (22 orders of magnitude): value identities carry ~1e6 coefficients while
price equations are O(1). This single fact explains Ipopt's unreachable absolute `tol`, UMFPACK `lu`
throwing `SingularException`, and a wildly false rank estimate. Fix: divide each row by its largest
Jacobian entry at the benchmark and run residual, Newton step and convergence test in those units,
so `tol` is effectively *relative*.

  | Jacobian | measured rank deficiency |
  |---|---|
  | unscaled | 71,101 ← **artifact** (SuiteSparseQR's tolerance keys off the ~1e6 max column norm) |
  | row-equilibrated | 643 |
  | row + column equilibrated | **637 ← real** |

**6. Rank deficiency 637 → 0 — RESOLVED (2026-07-25).** All 637 null directions were numerically
empty columns, i.e. variables multiplying a *zero benchmark flow* (`ppur_s` 534, `psuppmar_p` 34,
`plnd` 15, plus ~54 that resolved once those were pinned — there was **no** residual Walras-type
dependency, contrary to the initial estimate). Two distinct mechanisms, both now guarded in `src/`:

- **Dense-but-uninformative equation** (`E_ppur_s!`, `E_plnd!`): the constraint IS emitted and DOES
  hold, but its derivative w.r.t. the variable is exactly 0 at a zero flow —
  `d/d(ppur_s) log(ppur_s*X + _TINY) = X/(ppur_s*X + _TINY) = 0` when `X = 0`. Fix: emit the
  constraint only where the benchmark flow is non-zero, and pin the variable at 1.0 otherwise.
- **Skipped equation** (`E_psuppmar_p!`): already had a `|| continue` guard, leaving the variable
  free with no equation. Fix: `if/else`, pinning at 1.0.

Result: equilibrated Jacobian is **full rank (88,599)**, `lu` no longer singular (the SPQR
least-squares fallback is not needed), Newton **converges in 1 iteration / 2.2 s**, and the drift off
benchmark fell from 6.7e-4 to **5.4e-7**. Note the model is now **naturally square from `src/`
alone**: the scratchpad squaring shrank from "pin 157 + delete 108 dead rows" to **pin 1, delete 0**,
and that single pin is `natfhou` — the documented stub already noted under gate 5b, not a new hack.

Diagnostic worth keeping: a numerically empty Jacobian column is invisible to structural (matching)
analysis, which sees the variable in an equation and calls it determined. Only the equilibrated
numerical rank catches it. Run both.

**7. Productionised (gate #5c-i, done 2026-07-25; solver updated 2026-07-28).** New
**`src/solve_newton!.jl`** exports `solve_newton!(m, vars; maxit, tol, verbose) -> NewtonResult`.
It owns the whole procedure: index the free variables, defensively square the system (pin *orphan*
variables that appear in no constraint — one at the time of writing, `natfhou`, **none since
2026-07-31, when its missing equation was written; see the superseded-stub note below**; delete
*dead rows* whose variables are all fixed — today none), build the MOI nonlinear evaluator, row-equilibrate, then run a
damped Newton with fraction-to-boundary and backtracking. The linear solver was changed from sparse
`lu`/`qr` to **regularized CGNR** (CG on the augmented system `[J; sqrt(μ)I]·d ≈ [−F; 0]`, μ=0.1) to
handle the ill-conditioned Jacobian without step explosion. If the system is not square it raises
with a diagnostic pointing at the matching procedure rather than silently solving the wrong problem.
`solve_model!.jl` (Ipopt) is retained for reference. `Project.toml` gains `SparseArrays` +
`LinearAlgebra`. `test/solve_benchmark_6reg.jl` now passes unaided:

```
squaring: pinned 1 orphan var(s) [natfhou], deleted 0 dead row(s)
system: 88599 equations, 88599 unknowns (square)
Solve: 42.8s  status=converged  iters=1  scaled ||F||_inf = 1.106e-9
Max |Δ| vs benchmark: 1.0e-6 absolute / 5.4e-7 relative
PHASE 1 GATE PASSED ✓
```

Note the test asserts on the **scaled** residual and on deviation from the **benchmark seed**
(snapshotted before the solve) — not on "ratio vars ≈ 1", which is wrong for the value-flow variables
that `bmk_levels` deliberately seeds at their benchmark flow.

### Gate #5c status detail (2026-07-23, paused per user direction) — HISTORICAL

Progress since the row above was first written:

1. **`INVALID_MODEL` bug — found and fixed.** `initialize_model!.jl` was fixing every `TERM.CMF`
   closure variable (`BASE_CLOSURE_SCALARS`/`BASE_CLOSURE_ARRAYS`) at a uniform benchmark of `0.0`.
   Correct for a %-change formulation, but wrong for the levels formulation this model now uses:
   17 of the 41 closure names are **ratio-type** variables (declared `>= 1e-6`, e.g. `phi`, `acap`,
   `blab`, `bprim`, `tuser_su`, `nhou`, `pfimp`) whose correct benchmark is `1.0`, not `0.0` — fixing
   them at `0.0` sent every `log(...)` reading them to `-Inf`, cascading through the model.
   **Fixed** by adding a `_is_ratio_type` helper (same `>= 1e-6`-bound convention
   `test/solve_benchmark.jl` already used for free variables) so closure variables now fix at `1.0`
   (ratio-type) or `0.0` (genuine shift/shock-type), matching each one's real benchmark value. Confirmed
   via `test/check_start_point.jl` (a diagnostic that evaluates every constraint at the fixed/start
   initial point, no solver call): **124,100 → 0 non-finite constraints** after the fix. The real
   Ipopt solve no longer terminates `INVALID_MODEL` either — it now runs genuine Newton iterations.
2. **New blocker found: MUMPS (Ipopt's default linear solver) cannot factor this model's KKT system.**
   Running `test/solve_benchmark.jl` for real now fails with `MUMPS returned INFO(1) =-13 - out of
   memory when trying to allocate N MB` (N ranging 5,288–38,784 MB depending on `mumps_mem_percent`),
   eventually `EXIT: Restoration Failed!` / `Status: OTHER_ERROR`. This is **not a genuine memory
   shortage**: the host has ~16.68GB free (confirmed via `Get-CimInstance Win32_ComputerSystem`/
   `Win32_OperatingSystem`), and a raw Julia allocation-loop test succeeded cleanly up to 20GB in one
   contiguous block. Lowering `mumps_mem_percent` from Ipopt's default (~1000%) to 50 did shrink the
   *requested* allocation (38,784 MB → 5,288 MB) but did **not** avoid the failure — it recurred at the
   smaller size too. Most likely explanation: **MUMPS's bundled build (via `Ipopt_jll`) hits an internal
   scaling limit — likely 32-bit integer indexing — for a KKT system this large** (2.4M variables /
   1.4M constraints, ~3.8M×3.8M combined with fill-in), a known practical ceiling for open-source MUMPS
   at multi-million-variable NLP scale, distinct from any correctness bug in the model. Only
   `Ipopt v1.15.0` is installed (`Pkg.status()` confirms no `HSL_jll`/Pardiso) — HSL's MA57/MA86/MA97,
   the standard alternative for CGE-scale KKT systems, needs a separate free-for-academic-use
   registration through STFC and building `HSL_jll`, an external step outside what can be done
   unilaterally in-session.
3. **User decision (2026-07-23): pause here.** Given the choice between installing HSL, investigating
   fill-in reduction, or pausing to report status, the user chose to pause active solver work and wait
   for further direction. **No further solver tuning, HSL installation, or fill-in investigation should
   be undertaken until the user asks again.** The three options, for reference when that direction comes:
   - **Install HSL (MA57)**: register (free, academic) at STFC, build `HSL_jll`, point Ipopt at
     `linear_solver="ma57"`. Likely the most robust real fix, but needs an external registration step
     the assistant can't complete alone.
   - **Investigate/reduce KKT fill-in**: profile which equation families contribute the most nonzeros
     to the Jacobian/Hessian and look for reformulations (e.g. eliminating dense sum-type constraints)
     that shrink fill-in enough for MUMPS to cope. Uncertain payoff, more invasive to the model.
   - **Pause** (chosen): stop here, document state, wait.
4. **Known follow-up, not yet done**: `solve_model!.jl`'s `print_level=8` default was documented as
   avoiding Ipopt's full-vector per-iteration dumps, but empirically it still produces them (`new
   vars[i]`, `curr_c[i]`, `curr_x[1][i]`, `final y_c/y_d/z_L/z_U` for all 2.4M variables) — one run
   produced a 1.5GB/34M-line log file. That docstring/default needs correcting in a future session
   (only `print_level>=11` was assumed to do this; 8 does too at this model's scale). Not urgent enough
   to act on alone — noted here so it isn't lost.

### Step 6 / Step 8 notes (2026-07-25)

**Dynamic extension — two real data bugs found and fixed on the way.**
1. `build_premod!.jl` was using **hard-coded placeholders** for every dynamic
   coefficient (`NATCAPSTOK = ones(NI)`, `DPRC = 0.08`, `RNORMAL = 0.1`, …) even though
   `DPRC`/`TARG`/`TFRO`/`QRAT`/`ALFA`/`RADJ`/`REXP` are **genuinely present in
   `national.har`** — `build_reg0!.jl` simply never forwarded them into `elast_dict`.
   Now forwarded and used, with the old constants kept only as documented fallbacks.
2. `CAPSTOK` (header `STOC`) really is absent from this dataset, so it is **derived**,
   not invented: `CAPSTOK = CAP / RNORMAL`, i.e. capital stock = rental ÷ gross rate of
   return. That makes benchmark `GROSSRET = CAP/CAPSTOK` equal `RNORMAL` exactly — the
   steady-state condition ORANIGRD assumes initially, which is also what lets the
   Excerpt 51 investment rule start from balance.
   The old `NATCAPSTOK = 1.0` placeholder implied a capital stock of one currency unit per
   industry, making `GROSSGRO = INVEST_C/CAPSTOK` astronomical and `MCOEFF` explode; the
   symptom was `finv4` drifting 0.95 off benchmark through `mratio → gro → finv4`. With
   the real derivation the drift is gone.
   Also fixed: `STOC` was dropped at aggregation — it is a value FLOW, so it now sums
   (like `1CAP`) in `aggregate_model!.jl`/`aggregate_regions!.jl` rather than being
   weighted-averaged with the rate parameters.

**GDP slot-ordering bug (same class as `E_pdomB`; found while writing Step 8).**
TERM.TAB indexes `delGDPINC` by **category name**, and
`GDPINCCAT = (Land, Labour, Capital, ProdTax, ComTax)` — so `E_delGDPINCb`, despite being
the *second* equation, writes `"Capital"` = **slot 3**, and `E_delGDPINCc` writes
`"Labour"` = **slot 2** (`TERM.TAB:1418-1421`). The port had translated the equation
letters positionally (a→1, b→2, c→3), silently swapping Labour and Capital relative to
`GDPINCSUM`, which `prepare_parameters.jl` fills in the correct category order.
`wgdpinc` sums over all categories, so **totals were right and benchmark replication
could never detect it** — only per-component reporting would. Fixed, with a regression
assert in `test/gdp_report_6reg.jl` that checks the Labour/Capital slots against
`LAB_IO`/`CAP_I` directly. The expenditure side was already correct
(`b→4, c→6, d→9, e→7, f→8`, and `FINDEM → cols [1,2,3,5]` skipping STOCKS).

**⚠️ Benchmark-database finding (not a translation defect): the regional GDP identity does not
hold in the source data.** Step 8's reporting compares income-side and expenditure-side GDP two
different ways, and they say different things — deliberately kept separate:

| check | what it measures | result |
|---|---|---|
| `gdp_change_consistency` (`wgdpdiff`) | **model**: do the two sides *move* together? | **1.2e-12 ✓** |
| `gdp_both_sides_check` (levels) | **data**: do the two sides *equal* each other? | **off by up to 36.1%** |

TERM's own `wgdpdiff = wgdpinc − wgdpexp` compares two bmk=1 ratio indices, so it is **zero at the
benchmark by construction and structurally blind to a levels mismatch** — which is why the levels
check had to be added separately rather than trusting the model's built-in diagnostic. At a zero
shock the levels check reduces to `GDPINCSUM` vs `GDPEXPSUM` as calibrated, i.e. it is a property of
the database. Nationally the gap is 1.07%; by region: Sumatra +1.2%, Java −1.1%, Kalimantan +18.8%,
Sulawesi −23.9%, BaliNusa +8.2%, MalukuPapua +36.1%. Same family as the imbalance that forces the
`K` constant in `E_pdomA_sum!` and shows up as `DIFFCOM_sc ≈ 0.20` in the pipeline diagnostics.
**Consequence:** national and Java/Sumatra results are on solid ground; results for the small
eastern regions should carry a caveat. The test reports this rather than asserting it — asserting
would amount to asserting the input data is perfect.

### Region-scaling measurement (2026-07-25) — Phase 3 condensation is REQUIRED

Measured with `test/scale_probe.jl` (pipeline run once, then build+solve+factorize per size).
All sizes converged in 1 Newton iteration and stayed exactly square, so this is purely a
resource question:

| regions | variables | jac nnz | **LU nnz** | **fill ratio** | LU time | peak RSS |
|---|---|---|---|---|---|---|
| 6  | 121,452 | 427,015   | 6,592,224   | 15.4 | 3.0 s  | 2.5 GB |
| 10 | 258,752 | 846,606   | 18,745,819  | 22.1 | 5.2 s  | 3.5 GB |
| 14 | 449,396 | 1,397,732 | 76,275,676  | 54.6 | 26.1 s | 9.1 GB |
| 20 | 850,502 | 2,527,002 | 176,550,811 | 69.9 | 57.5 s | **16.2 GB** |

**The binding constraint is LU fill-in, not variable count.** The Jacobian grows roughly `n^1.7`,
but its LU factors grow `~n^2.7`: the fill ratio itself climbs 15 → 22 → 55 → 70. At n=20 peak
memory already hits the machine's ~16 GB ceiling; extrapolating to n=34 gives ~6M Jacobian
nonzeros and a fill ratio well past 100, i.e. hundreds of millions of LU nonzeros. **34 regions is
not reachable on this hardware without condensation.** n=14 (9.1 GB) is the largest comfortable
size today; n=20 is the practical edge.

**Correction to the earlier plan.** The hypothesis recorded above — that sparse declaration alone
might suffice because "the thing making the model bigger is also making it emptier" — is
**refuted**. Two reasons:
1. Fill-in grows faster than the added sparsity helps (the table).
2. More fundamentally, **sparse declaration cannot reduce fill-in at all**: zero-flow cells are
   already `fix`ed by the guards in `build_equations.jl`, so they are already excluded from the
   Jacobian's free columns. Declaring them sparsely saves JuMP *build* memory only — it does not
   remove a single nonzero from the matrix being factorized.

So the priority flips: **Excerpt 49 `Substitute`/`Backsolve` condensation is the necessary lever**,
because it deletes whole variable *and equation* blocks (`xtradmar`, `xsuppmar`, `ppur`, …) from the
core system, genuinely shrinking what gets factorized. Worth pairing with an ordering investigation:
fill-in this severe usually means a few near-dense rows (the aggregation constraints that sum over
all regions/commodities grow longer with `nr`) are poisoning the elimination tree.

## Solver decision: Ipopt (not PATH) — reasoning

> **⚠️ SUPERSEDED 2026-07-25 — see "Gate #5c RESOLVED" above.** The reasoning below (why the model is
> a square system of equalities, and why PATH's complementarity machinery is unnecessary) remains
> **valid and worth reading**. What is now falsified is the *conclusion* that **Ipopt** should solve
> it. Empirically Ipopt hangs at `iter 0` on the exactly-square system, because zero degrees of
> freedom is degenerate for an interior-point optimizer. The very "square system of nonlinear
> equalities with no objective" property argued for here is precisely what makes a **Newton solver**,
> not an optimizer, the right choice. Current solver: sparse Newton on `F(x)=0` via MOI's nonlinear
> evaluator (~30 s at 25×6). PATH remains a reasonable future fallback if robustness under large
> shocks becomes the binding constraint.

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
is worth vendoring). **Status: done.** — one loose end: `HeaderArrayFile` is still declared in
`Project.toml` and is referenced by no `src/*.jl` file (only named in `read_data.jl`'s docstring as
the package that *couldn't* read three headers). Drop it before release.

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
Formula-to-Julia-array-computation translation WayangJulia used in its own Step 2. **Status: done
and verified** (2026-07-21) — `prepare_parameters.jl` computes **85** derived parameters, tabulated
Formula-block by Formula-block in `## Step 4 implementation details` below, and they feed the
benchmark solve that reproduces the base year at ‖F‖∞ ≈ 1.4e-9.

> This line read `**Status: not started** — depends on Steps 2-3 landing first` until 2026-07-30,
> long after the step was finished, contradicting both the milestone table and the implementation
> section below it. Steps 5-8 below had the same stale line for the same reason: the per-step Status
> lines were maintained through Step 3 and then abandoned in favour of the §N running log. If you
> are reading a Status line here, cross-check it against the milestone table before quoting it.

### Step 5 — Core equations (`build_model!.jl`)
Full translation of `TERM.TAB`'s Excerpts 1-49ish (production nests, Armington sourcing, regional
trade/margins via `ORG`/`DST`/`PRD`, CET export supply, government/household/investment demand,
market clearing, GDP both sides, Excerpt 38 labour-market closure), using the same three derivation
principles validated in WayangJulia (Principle A: common-price aggregation → plain sum; Principle
B: zero-pure-profits %-change → exact levels value identity; Principle C: genuine CES/CET nests via
a `ces()` helper), plus the log-differentiate-and-match technique for any %-change-linear equation
whose exact nonlinear levels form isn't immediately obvious. **Status: done and verified** —
`build_model!.jl` + `build_equations.jl` + `build_macros!.jl` translate Excerpts 1-49; the system is
exactly square (83,515 x 83,515 at 25 sectors x 6 regions) and the benchmark solve reproduces the
base year at ‖F‖∞ ≈ 1.4e-9. Reporting excerpts 30/33/34/36/37/42/43/45/46/48 are deliberately not
translated (diagnostics, not model content) — see the deferred list in §N.

### Step 6 — Dynamic + district extensions (`build_dynamics!.jl`, `build_district!.jl`)
Excerpts 50-54 (investment rule, real-wage adjustment) and Excerpt 55 (top-down district
extension — likely translatable as a WAYANG-style non-feedback satellite, given its name).
**Status: dynamics done; district extension N/A.** `build_dynamics!.jl` translates Excerpts 50-54
(investment rule, capital accumulation, `ELASTWAGE` real-wage adjustment — see §N.23). There is no
`build_district!.jl` and there should not be: Excerpt 55 needs district-level data that the shipped
database does not contain, so the step has nothing to run on. Not deferred — inapplicable.

### Step 7 — Closure & solve wiring (`initialize_model!.jl`, `solve_model!.jl`, `run_model!.jl`)
Encode one of the existing `.cmf` closures (`TERM.CMF` as the default, following WayangJulia's
`default_closure()` pattern) as an exogenous/endogenous variable split; Ipopt feasibility solve, no
objective, mirroring `WayangJulia/src/solve_model!.jl` verbatim in structure (same `set_attribute`
tuning: `max_iter`, `tol`, `constr_viol_tol`, adaptive `mu_strategy`); a `run_model!` homotopy driver
for shocks too large for a single direct solve, per WayangJulia's own conditioning findings (see its
`PLAN_SOLVER_CONDITIONING.md` — raw-vs-log CES residual scaling, positivity floors on prices).
**Status: done, but NOT as specified above.** The closure split lives in `closures.jl` +
`initialize_model!.jl` as planned. The solver does not: **Ipopt was abandoned** — it hangs on this
zero-DOF square system regardless of tuning — and replaced by a sparse Newton solve
(`solve_newton!.jl`) with natural-parameter homotopy (`continuation.jl`), pseudo-arclength
continuation past limit points (`arclength.jl`), and the `run_model!.jl` driver taking named
`Scenario`s (`scenarios.jl`, §N.22). `solve_model!.jl` retains the Ipopt path for reference only.

### Step 8 — Reporting (`calculate_gdp.jl`)
GDP income/expenditure decomposition and the GDP-both-sides invariant check, per this model's own
`CHKMOD.TAB` diagnostic conventions. **Status: done** — `calculate_gdp.jl`. Caveat recorded in §N:
the *benchmark database itself* fails the regional GDP identity by up to 36.1% (MalukuPapua). That
is a property of the shipped data, not of the translation, and no code change can fix it.

## Verification

> **See `VV_PLAN.md` for the full V&V plan (gates V1-V10).** The list below is the original
> Step-0 sketch and it is *incomplete as a V&V protocol*: every check in it is an
> internal-consistency check, so passing all of them would establish that the Julia model
> agrees with itself while establishing nothing about whether it agrees with `TERM.TAB`.
> `VV_PLAN.md` adds the Walras identity, numeraire invariance, path independence,
> aggregation consistency, analytic multipliers, closure ordering, a real regression suite,
> elasticity sensitivity, and the cross-check against the original GEMPACK model — the last
> being the single most valuable outstanding item in the project and, as of 2026-07-30,
> never attempted.

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
4. Bind into a solvable system (fix shock variables, set up solve)
5. Run a small shock (e.g. −5% export shift) and verify against expected signs
6. Benchmark replication (all shocks at 0 → ~0% change)
7. Price homogeneity test

## Immediate next steps (2026-07-28) — ⚠️ SUPERSEDED 2026-07-29

> Item 1 below is **cancelled**. See "Session findings (2026-07-29)" at the end of this file: the
> stall is a rank/consistency defect, not a step-quality defect, and further solver tuning treats a
> symptom. Follow plan **C0–C6** there instead. Items 3–4 survive as C6.

1. Make the second Newton step productive again: try adaptive LM μ (reduce as residual drops),
   Jacobi/diagonal preconditioning for CGNR, or a trust-region dogleg step.
2. Re-test Euler continuation once Newton can converge at a single small shock; use full Newton
   correction at each continuation step if affordable, or adaptive step size.
3. Apply the numeraire swap (`phi` free, `GDPPI` fixed) and run the price-homogeneity test.
4. Once a small shock converges, run a smoke-test scenario (e.g. −5% export shift) and compare
   signs/magnitudes against economic intuition.

## Current active blockages (2026-07-28)

### Active
1. **CGNR/Newton solver gives one good step on a shock, then stalls.** The primary linear
   solver is now regularized CGNR (Levenberg-Marquardt/Tikhonov, μ=0.1) on the augmented system
   `[J; sqrt(μ)I]·d ≈ [−F; 0]`. This fixed the earlier blow-up from numerically empty Jacobian
   columns and unregularized null-space drift. Result on `blabnat=0.9997` (0.03% shock):
   - First Newton step: scaled ||F||_inf drops 0.000300 → 0.0001937 (real progress).
   - Second Newton step: line search finds **no α** that reduces the residual below 0.0001937;
     CGNR itself converges to relative residuals ~3e-5 / 2.6e-4 in 200 iters, so the linear
     solve is healthy — the nonlinear model is simply a poor predictor past the first step.
   - **Underlying cause:** extreme CES/CET curvature and the very small shock domain; the
     linearized model at the post-shock point no longer predicts descent.
   - **Next options:** adaptive μ (start larger, shrink as residual drops), a true trust-region
     step, or Euler/homotopy continuation with smaller sub-steps and full correction at each.

2. **Euler continuation residual accumulates.** 50 steps with 3 inner corrections/step
   gives ||F||≈0.012 at step 20 (trending to ~0.07 by step 50); Newton cleanup stalls
   because backtracking alpha collapses to ~1e-8 (same fraction-to-boundary issue). Needs either:
   - Full Newton convergence at each Euler step (regularized CGNR may now make this feasible), OR
   - Adaptive step sizing with more inner corrections per step

3. ~~**Numeraire swap not applied.**~~ **STALE — RESOLVED.** The swap *is* applied, at
   [initialize_model!.jl:207-219](src/initialize_model!.jl:207): `phi` is unfixed and started at 1.0,
   `NatMacro("GDPPI")` is fixed at 1.0, with a fallback to fixing `phi` if `NatMacro` is absent.
   The price-homogeneity test is still unrun, but nothing is blocking it except the shock stall.

4. ~~**`build_pstras!.jl:173` `DISTANCE = reg1["MAKE"]`** — self-flagged as incorrect; needs proper
   distance data from `reg2`.~~ **NOT A DEFECT — RESOLVED (2026-07-29).** The assignment was dead:
   the very next line overwrote it with the correct `DISTANCE = ras["DIST"]` (ORG×DST, NR×NR) before
   anything read it, so the wrong value never reached the `DISTGONE` loop. Margin calculations were
   never affected. Removed the dead line and the stale "not right" comment that made it read like an
   open bug.

### Resolved
- **`ras_balance!` performance regression (was #1).** **FIXED (2026-07-26).** Root cause was
  not generator loops (those were already explicit) but **type-unstable array bindings**:
  `BASIC_U = parent(reg2["BSCU"])` etc. inferred as `Any` because `reg2` is `Dict{String,Any}`.
  Adding `::Array{Float64,N}` assertions on the five input arrays reduced Phase 1 from
  467s → 0.14s and Phase 2 from 177s → 0.22s. Full pipeline now completes in ~60s.

- **Newton step explosion from tiny-flow columns.** **FIXED (2026-07-28).** Original direct
  Newton/LU produced |d| ~ 100–10^5 because frozen columns and near-singular J made the linear
  solve unstable. Switched to **regularized CGNR** (μ=0.1) on `J'J + μI`. It now converges in
  ~200 iterations (vs. 0.0121 residual after 1000 unregularized iterations) and the first step
  is bounded and productive. Remaining issue is the second-iteration flat direction, not step
  explosion.

### Not started
5. **Gate 7 — closure & solve wiring beyond base.** `run_model!.jl` homotopy driver for shocks too
   large for a direct solve. Blocked on Newton shock convergence.

7. **Excerpt 49 condensation (34-region prerequisite).** 34-region Jacobian LU fill-in exceeds
   ~16 GB RAM (20 regions = 16.2 GB). Requires TABLO `Substitute`/`Backsolve` condensation to
   delete whole variable/equation blocks from the core system.

8. **RAS outer loop.** `converged=false` expected with single-pass RAS; multi-pass loop
   (feed `DGON` back into reg2) deferred until speed issue is resolved.

### Done
10. **Benchmark solve (gate 5c).** ✅ Solver switched Ipopt→Newton; converges 1 iter, 9.8s,
    scaled ‖F‖_inf=1.11e-9.
11. **Rank deficiency (gate 5c-ii).** ✅ 637→0 via zero-flow guards on `ppur_s`/`psuppmar_p`/`plnd`.
12. **Dynamic extension (gate 6).** ✅ 89,911 vars, converges 1 iter, benchmark replication unchanged.
13. **Macro reporting (gate 6c).** ✅ 90,460 vars, includes `NatMacro("GDPPI")` for numeraire swap.
14. **GDP reporting (gate 8).** ✅ Income/expenditure decomposition, model-side wgdpdiff=1.2e-12.
15. **Elasticities wired.** ✅ `SLAB`, `P028`, `SMAR`, `PO01`, `SCET`, `P018` stored in `prepare_parameters.jl`
    output dict (lines 198-204) and read by `build_model!.jl` via `haskey(params, ...)` fallback.

## Status summary (2026-07-28)

| Step | Description | Files | Status |
|------|-------------|-------|--------|
| 0 | Scaffold | `Project.toml`, module | ✅ done |
| 1 | Sets & data ingestion | `prepare_sets.jl`, `read_data.jl` | ✅ done & verified |
| 2 | Regionalization + RAS pipeline | `build_reg0!.jl`..`build_premod!.jl` (6 files) | ✅ done & verified |
| 3 | Aggregation 185→25 × 34 | `aggregation_data.jl`, `aggregate_model!.jl` | ✅ done & verified |
| 4 | Derived parameters | `prepare_parameters.jl` | ✅ done & verified (85 params, elasticities wired) |
| 5 | Core equations | `build_model!.jl`, `build_equations.jl` | ✅ 2,436,062 vars / 1,370,212 cons (levels form) |
| 5b | Base static closure | `initialize_model!.jl` | ✅ fix/free matches TERM.CMF |
| 5c | Benchmark solve | `solve_newton!.jl` | ✅ 9.8s, 1 iter, ‖F‖_inf=1.11e-9 |
| 5c-i | Self-executing solve | `solve_newton!.jl` | ✅ moved from scratchpad to src/ |
| 5c-ii | Rank deficiency removed | — | ✅ 637→0, full rank equilibrated Jacobian |
| 6 | Dynamic extension | `build_dynamics!.jl` | ✅ 89,911 vars, converges 1 iter |
| 6b | District extension | `build_district!.jl` | ⛔ N/A (data headers missing) |
| 6c | Macro reporting | `build_macros!.jl` | ✅ 90,460 vars, includes GDPPI |
| 7 | Shock solve wiring | `run_model!.jl` | ❌ blocked (CGNR step 1 works, step 2 stalls) |
| 8 | Reporting | `calculate_gdp.jl` | ✅ wgdpdiff=1.2e-12 |

* RAS pipeline runs and explicit-loop rewrite **speed-verified** (~60s full pipeline)

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

> **Superseded 2026-07-31 — and the reasoning above is where it went wrong.** The
> justification was sound *when written* (no `NatMacro` layer existed) and became stale the
> moment `build_macros!` was added, but nothing re-opened it, because a pinned orphan and an
> implemented equation give identical answers for every closure that leaves the variable
> endogenous. `origin/coalprice.CMF:83` swaps `houslack = natfhou`, making natfhou
> **exogenous**: the orphan disappeared, `solve_newton!` had nothing to pin, and the model came
> out 84,136 free variables against 84,135 equations. `E_natfhou!` now implements TERM.TAB:2053
> as `log(NatMacro["NomHou"]) == natfhou + log(NatMacro["NomGDPexp"])`, called after
> `build_macros!`. There are now **zero orphans and zero dead rows** at 25×6.
>
> The transferable rule: **an orphan variable is a missing equation until proven otherwise.**
> "It only feeds reporting" is not a reason to leave it unwritten — a later closure can make
> the variable exogenous, and then the missing equation is the difference between square and
> unsolvable.

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

## Course correction (2026-07-22): the entire Step 5 "executing" section above is off-plan

The user flagged, correctly, that the model had drifted from a **levels formulation** to a **Johansen
%-change formulation** — a deviation from this file's own original decision. Investigation confirmed
it in full:

- **Line 14** (original solver-decision writeup, written before Step 5 started): "Chosen: JuMP +
  Ipopt, **levels formulation**, solved as a square system of nonlinear equalities... exactly the
  architecture validated by `WayangJulia`."
- **Lines 291-300** (original Step 5 plan): laid out a real conversion methodology — Principle A
  (common-price aggregation → plain sum), Principle B (zero-pure-profits %-change → exact levels
  value identity), Principle C (genuine CES/CET nests via a `ces()` helper), plus log-differentiate-
  and-match for anything else — marked "Status: not started."
- **Line 417** ("Step 5 — Core equations (executing)," written once Step 5 actually started, commit
  `e5381a9`): silently redefines the approach as "Linearized %-change form... following TERM.TAB's own
  %-change convention," with no cross-reference to line 14 or acknowledgment that Principle A/B/C was
  abandoned. Every one of the ~130 `E_*!` functions since (all of Excerpts 6-39, including everything
  added in Step 5b above) is a direct, literal port of TABLO's own %-change variables (`p*` = price
  %-change, `x*` = quantity %-change, `w*` = value-share-weighted change, `del*` = level-delta, `f*` =
  shift/closure) as linear JuMP constraints — never converted to levels.
- **Independent confirmation**: `test/diagnose_gap.jl`'s constraint introspection filters on
  `F <: GenericAffExpr` and catches the *entire* model — every constraint in the ~2.4M-variable system
  is purely affine. A genuine levels CGE cannot be 100% affine (CES/CET nests, zero-profit price×
  quantity identities, and value aggregations are all nonlinear in levels) — this is hard evidence the
  levels conversion was never done.
- **The precedent this file cites didn't take the shortcut.** WayangJulia's own `plan.md` and its
  `src/` actually execute Principle A/B/C: `build_model!.jl` declares TABLO's own lowercase variable
  names as genuine dollar-levels JuMP variables (benchmark convention: price=1, quantity=benchmark
  value, away from benchmark both are real levels — no separate levels-naming scheme, just a semantic
  reinterpretation of the same names), computes CES/CET share-and-scale calibration parameters once via
  a closed-form `_ces_calibrate(quantities, sigma, output)` routine (`prepare_parameters.jl:28-52`),
  and expresses every nest via a shared `ces(y, p, α, σ, γ)` closed-form demand function (imported from
  `ComputableGeneralEquilibriumHelpers`) inside `log(...) .== log.(ces(...))` constraints. Zero-profit
  conditions become literal `log(aggregate_qty*aggregate_price) == log(sum(component_qty*component_price))`
  identities (Euler's theorem in levels, valid at any elasticity — not just at the benchmark cost
  shares). Regional/auxiliary "(change)"-type variables with no other economic pinning (`r1cap`,
  `omega`, and this project's own `zcon_reg`/`ztot_reg` analogues) are handled as log-differential,
  benchmark-0 variables built from `log()` of already-defined benchmark-1 price variables. None of this
  is exotic or unavailable — it's a fully worked, validated template sitting in the same Julia
  environment (`~/.julia/packages/GlobalTradeAnalysisProjectModelV7/GTAPJulia/WayangJulia/`).

**Decision, per explicit user direction: convert `build_equations.jl` to true levels, following
WayangJulia's methodology directly rather than re-deriving one from scratch.** This supersedes the
"Approach" bullets at line 417 and the "Status: not started" markers at lines 288/298 — those are now
historical record of the deviation, not the live plan. The %-change equations already written are the
correct *reference* (they encode the right economics, just linearized) — the conversion work is
translating each into its exact nonlinear levels counterpart using WayangJulia's four techniques, not
re-deriving the underlying CGE theory from `TERM.TAB` a second time.

### Levels-conversion plan (Step 5c)

1. **Infrastructure**: port `ces()`/`_ces_calibrate()` into IndotermJulia (new small helper file,
   adapted to the multi-region 25×34 scale — the same functions, no change to the closed-form math).
2. **Wire real elasticities**: `sigmadomimp`/`sigmalab`/`sigmaprim`/`sigmaout`/`exp_elast` are
   currently hardcoded placeholder constants in `build_model!.jl:283-287,390-391` — `SLAB`, `SCET`,
   `P018` are loaded into `prepare_parameters.jl` (lines 24, 26, 37) but never used again. These need
   to actually feed the CES/CET nests during the rewrite, closing a gap that predates this course
   correction.
3. **Convert excerpt-by-excerpt** (matching the catalog order already in this file): 6-12 (basic/
   purchaser prices, Armington CES, factor CES, zero-profit costs) → 13-17 (household LES, investment,
   government, export CET) → 19-23 (margins, regional-sourcing CES, MAKE/CET, market clearing) →
   24-29 (final-demand aggregates, tax revenue, GDP income/expenditure) → 27/38-39 (labour-market and
   household closure). Each block: identify which of the four techniques applies per equation
   (Principle A/B/C, or log-differential auxiliary), write the levels form, verify by taking the
   candidate nonlinear equation's differential at the benchmark point and confirming it reproduces the
   original %-change equation's coefficients exactly (WayangJulia's own verification technique).
4. **Update the base closure and solve wiring**: `initialize_model!.jl` currently sets every start
   value to `0.0` and fixes closure variables at `0.0` — correct for a %-change-at-benchmark
   formulation, wrong for levels (start values/fixes need to be the actual benchmark levels, matching
   WayangJulia's `generate_starting_values.jl` pattern). `solve_model!.jl` itself (Ipopt feasibility,
   zero objective) needs no change — that part of the architecture was never the problem.
5. **Re-verify**: benchmark replication becomes a real test again (a true levels model reproducing
   benchmark values exactly, not the trivially-true "all-zero %-change" check a linear homogeneous
   system gives for free), plus price homogeneity and GDP-both-sides on the converted model.

This is a large rewrite — WayangJulia's own equation-writing core (`build_model!.jl`, `build_fiscal!.jl`,
`build_regional!.jl` combined) is 1,267 lines / 151 named constraints for a national-only model;
IndotermJulia's current %-change `build_equations.jl` is already ~1,200 lines / ~130 functions *with*
the added regional dimension, so the levels rewrite is expected to be comparably large or larger.

---

## Session continuation (2026-07-26)

**Primary achievement:** Fixed the `ras_balance!` performance regression (467s → 3.4s) by adding concrete type assertions (`::Array{T,N}`) to all array bindings. The root cause was type instability: `reg2` is a `Dict{String,Any}`, so `parent(reg2[...])` returns `Any`-typed bindings. The compiler couldn't specialize inner loops, causing ~50× slowdown. Explicit type assertions restored specialization and brought phase 1 to 0.14s, phase 2 to 0.22s. Full pipeline (Steps 0-4) now runs in ~60s.

**Also done:** Updated PLAN.md status table and active blockages section. RAS pipeline speed is now **verified**, not just claimed.

**Next session priority order (from HANDOFF.md §7):**
1. Print elasticity values (30s) — closes the "wrong by construction under shock, invisible at benchmark" risk
2. Numeraire swap: free `phi`, fix `NatMacro("GDPPI")` per TERM.CMF; run price-homogeneity test
3. First real shock: `test/shock_blabnat_6reg.jl`, `blabnat = -3%` via 50-step Euler continuation
4. Fix `build_pstras!.jl:173` `DISTANCE = reg1["MAKE"]` defect (may explain 20% commodity imbalance)
5. Deferred: Excerpt 49 condensation (prereq for 34 regions — 20 regions already needs 16 GB), gate 7 `run_model!.jl`

**Environment note:** One idle bare `julia.exe` (PID 5640, started 07:57) still running — harmless but makes "is a run in flight?" checks ambiguous. Kill it for a clean signal.
Proceeding incrementally, excerpt by excerpt, committing and verifying each block rather than
attempting it in one pass.

---

## Session findings (2026-07-29) — the shock stall: evidence, verdict, and debug plan

This section supersedes "Current active blockages (2026-07-28) → Active #1" and item 1 of
"Immediate next steps (2026-07-28)". Both diagnosed the stall as *nonlinearity* and prescribed
*more solver tuning*. The evidence below says that diagnosis is wrong.

### A. What was measured

**A1. Every shock stalls, at a shock-specific floor — not just `blabnat`.**

| shock | closure | initial scaled ‖F‖∞ | stalled at | terminal state |
|---|---|---|---|---|
| `blabnat = 0.9997` (0.03%) | default | 3.00e-4 | **1.6897924105683306e-4** | `no_progress`, α→0 |
| `blabnat = 0.9997` | + 3 static swaps, `colnrm` freeze relaxed to 1e-12 | 3.00e-4 | **1.6897924105683306e-4** (bit-identical) | `no_progress` |
| `fgovgen = 0.01` | default | 9.95e-3 | **3.3927647286203875e-3** | `no_progress`, α=0 at iter 8 |
| `pfimp = 1.05` | default | 4.879e-2 | 4.807e-2 after 20 iters | α ∈ [2e-4, 8e-3], no progress |

Three economically unrelated shocks — a labour-productivity shifter, a government-spending
shifter, and an import price — all exhibit the same signature. A defect that is invariant to
*which* shock is applied is not a property of any one closure. **The defect is in the solver or in
the equation system's structure at a perturbed point, not in the closure.**

**A2. The stall floor is reproducible to the last bit under an unrelated code change.**
Relaxing the weak-column freeze `colnrm[i] < 1e-3` → `1e-12` at
[solve_newton!.jl:435](src/solve_newton!.jl:435) *and* applying three static closure swaps left the
`blabnat` stall at `0.00016897924105683306` — identical in every mantissa bit. A genuine
ill-conditioning / line-search / floating-point-path problem would perturb the floor. A
**hard structural obstruction** would not. This is the single most diagnostic observation of
the session.

**A3. CGNR reports a non-zero least-squares residual — F is not in the range of J.**
`_diag_rank3.jl` at the post-first-step point: CGNR relative residual **0.0169** at μ=0.01 after
200 iterations. For a *square, full-rank* system, CGNR must drive `‖J·d + F‖/‖F‖` to ~0. It does
not. There is a component of `F` in the **left null space of `J`** that no Newton step of any
length can remove — which is precisely a residual floor with α→0.

**A4. The μ (Levenberg-Marquardt) parameter is a pure tradeoff with no good setting.**

| μ | `cgnr_rel` | line-search α | ‖F‖ reduction |
|---|---|---|---|
| 0.01 | 0.008 – 0.02 (under-converged at maxit=200) | — | stalls |
| 0.1 | ~6e-6 (converged) | collapses to ~0.05 | stalls |
| 1.0 – 10.0 | ~1e-10 (exact) | — | step is steepest-descent-biased, barely moves ‖F‖ |

At every μ the outcome is the same floor. Regularization strength is *not* the free parameter that
unlocks this. Consistent with A3: raising μ makes the *linear* problem well-posed by changing the
problem, but the unreachable component of F is still unreachable.

**A5. Ruled out, individually, by direct experiment.**
- **Weak-column freeze** — relaxed 1e-3 → 1e-12; benchmark still passes (‖F‖=6.98e-10, solved) and
  the shock floor is bit-identical (A2). *Not the cause.*
- **CGNR iteration budget** — `maxit` 200 → 500 changed nothing (the apparent "timeout" was
  pipeline overhead, see A6). *Not the cause.*
- **Diagonal Jacobi preconditioning of CGNR** — implemented, appeared to hang; the hang reproduced
  with the preconditioner **reverted**, proving it innocent (again A6). Reverted and abandoned.
- **Closure choice** — three unrelated shocks, two different closures. *Not the cause* (A1).

**A6. ⚠️ Pipeline overhead is 429 seconds. Every previous "the solver hung" conclusion is void.**
Measured end-to-end with a timestamped harness:

```
T+0.00s   load data          T+228.7s  prepare_parameters      (17.7s)
T+37.6s   build_reg0         (34.3s)   T+246.4s  build_model_full!  (177.8s)  ← dominant
T+72.0s   build_reg1         (59.7s)   T+424.2s  benchmark_levels   (2.3s)
T+131.7s  build_reg2         (26.3s)   T+426.5s  initialize_model!  (1.7s)
T+158.0s  ras_balance        (9.0s)    T+428.2s  apply_static_swaps (0.9s)
T+167.0s  build_pstras       (17.0s)   T+429.1s  → solve_newton! ENTERED
T+184.0s  build_premod       (9.1s)
T+193.1s  aggregate_regions  (35.6s)
```

Two consequences: (i) the `ras_balance!` regression really is gone (9.0 s, as claimed at
2026-07-26); (ii) **`build_model_full!` at 177.8 s is 41% of the wall clock and every solver
experiment to date has paid 7 minutes of setup before the first Newton iteration.** Any test
budgeted under ~10 minutes was measuring the pipeline, not the solver.

**A7. The build is consistent with a square post-closure system (latest `eqcount.log`).**
At 25×6: `vars=122,001`, `cons=83,541`, `9,217` cells pinned by zero-flow guards at build time,
leaving `112,784` free pre-closure and a pre-closure gap of `29,243` — closed by the
`BASE_CLOSURE_*` fixes in `initialize_model!`, giving the square `83,541 = 83,541` the solver
reports. The `−6,919` equation delta vs. gate 6c's 90,460 is **expected, not a regression**: it is
the guards removing one equation *and* one unknown per zero-flow cell (`xmake` 96% pinned is
correct — MAKE is near-diagonal, 25×6 = 150 live cells of 3,750).

**A8. Closure facts established (keep — these are not the bug, but they were unclear before).**
- **The correct static-shock recipe** is: default closure **+ the 3 static `TERM.CMF` swaps applied
  *before* the benchmark solve.** Applied *after* → benchmark residual jumps to 0.944. Applied
  *before*, then re-solved → benchmark converges to 1.6e-6 and the shock starts at only 3.0e-4.
- **The 4 dynamic swaps need a multi-year driver.** `xcap=faccum` + `finv1=finv4` +
  `delfwage=flabsup_id` + `delUnity=1` gives benchmark ‖F‖=0.7958 and a shock initial residual of
  **42,618** (Excerpt 50's `CAPSTOK_OLDP·log(xcap/CAP) = CAPADD·delUnity` forces a one-shot capital
  jump). This confirms the design note in `build_dynamics!.jl` and `initialize_model!.jl` — the
  dynamic block is a *passive satellite* under a static closure, by construction.
- **`TERM.CMF`'s real shock is `blabnat = −3` (3%, levels 0.97).** Every test this session used
  `0.9997` — **100× smaller**. Even a fully working solver would still owe a continuation path to
  the real shock size.
- **`swap A = B` means A becomes endogenous, B becomes exogenous** (`TERM.CMF`'s own header comment:
  `! old exog   new exog`). The `delfwage = flabsup_id` swap was implemented backwards and has been
  corrected.

### B. Verdict — ❌ **REFUTED 2026-07-29 by measurement. See §E below for the corrected verdict.**

> Kept for the record because the reasoning is instructive about how it went wrong, not because any
> of it is actionable. Both premises turned out to be false: the "bit-identical floor" (A2) was a
> no-op code change, and the CGNR residual (A3) was non-convergence, not a null space. Skip to §E.

**The Newton solver is not failing to converge; it is converging to the least-squares solution of an
inconsistent system.** The evidence is A2 + A3 taken together: a residual floor that is invariant to
bit level under solver perturbation, plus a CGNR least-squares residual that will not go to zero.
That combination has one economical explanation — at a shocked point the Jacobian is
**rank-deficient**, and the shock pushes `F` into the resulting left-null direction.

Note the earlier finding that the *benchmark* Jacobian is full rank (gate 5c-ii, 637 → 0) is **not
in conflict**. The rank test was only ever run at the benchmark point. The zero-flow guards were
tuned so that every numerically-empty column disappears *there*. Nothing has ever verified rank at
a perturbed point.

**The leading mechanism, and why it hides at the benchmark:** a guard that pins the *variable* but
retains the *equation* leaves an over-determined row. At the benchmark that row reads `0 == 0` and
is invisible — the benchmark solves to 7e-10, as observed. Under any shock the row's terms move,
it becomes inconsistent, and it contributes an irreducible residual: exactly a shock-specific
floor that no step length can reduce. The guards touch **9,217 cells across 12 families**
(`xmake`, `xtradmar`, `delTAXint`, `ppur_s`, `xsuppmar`, `delTAXgov`, `delTAXexp`, `delTAXinv`,
`plnd`, `psuppmar_p`, `xtrad`, `delTAXhou`) — a large surface for exactly this error, and the
`if/else` pin-and-skip discipline was established retroactively, family by family, not uniformly.

Competing hypotheses, and why they now rank below it:
- *Severe nonlinearity / poor linear model* (the 2026-07-28 diagnosis) — **contradicted by A2**:
  nonlinearity would not produce a bit-identical floor across a changed code path, and by A4:
  a trust-region-like μ sweep finds no setting that helps.
- *Step caps / the second un-relaxed `colnrm < 1e-3` freeze at `solve_newton!.jl:~461`* — cheap to
  test, still worth eliminating, but A2 shows the *first* freeze is irrelevant, and a cap can only
  slow a step, never create a hard floor.

**Consequence for the plan:** step 1 of "Immediate next steps (2026-07-28)" — adaptive μ, Jacobi
preconditioning, trust-region dogleg — is **cancelled**. All three treat a symptom. Do not resume
solver tuning until the rank question is settled.

### C. Debug plan

Ordered so each step is cheap and each result changes what the next step is.

**C0 — Cache the built model (prerequisite; unblocks everything else).**
Serialize the pipeline output (`agg6` + `params`) to `data/cache_6reg.jls` and load it when present.
This removes ~250 s of the 429 s; add a `build_model_full!` cache (or accept 178 s) for the rest.
Rationale: every diagnostic below needs a built model, and at 7 minutes per attempt the debug loop
is unaffordable. *Pay this cost once.* Success: a solver experiment starts in <60 s.

**C1 — Fix the residual diagnostic (it has been silently lying).**
The per-equation residual reports returned `worst |resid| = 0.0` with 0 groups over 1e-9 — false.
`JuMP.value(con)` reads the **stale optimizer solution cache**; the Newton solver keeps its own `x`
and never calls `optimize!`, so the numbers were meaningless (accompanied by a flood of *"model has
been modified since the last call to optimize!"* warnings that should have been the tell).
The correct pattern already exists in [test/check_residual_6reg.jl](test/check_residual_6reg.jl):
pass an explicit point function, `JuMP.value(point, con)` where
`point(v) = JuMP.is_fixed(v) ? JuMP.fix_value(v) : JuMP.start_value(v)`. Port it into a reusable
`src/` diagnostic. Success: it reproduces the known 1.8e-4 worst-case at the benchmark.

**C2 — Name the offending equations (the decisive test).**
At the *stalled* `blabnat` point, write back the solver's `x` into the model's start values, then
run C1 and rank equations by `|residual|`, grouped by family via the `canon()` index-stripping
helper. The floor is 1.69e-4; the rows carrying it will stand out by orders of magnitude.
- **Expected outcome if the verdict is right:** the residual concentrates in one or two families,
  and each is a family that appears in the 12-family guard list of A7.
- **If instead it is smeared uniformly across thousands of equations,** the verdict is wrong and
  the problem is genuine ill-conditioning — go to C5.

**C3 — Confirm by rank, and identify the left-null direction.**
At the same stalled point, compute the equilibrated Jacobian `J_s` (reuse `compute_rowscale!` +
`ruiz_scales`) and get its numerical rank — sparse QR (`SuiteSparse.SPQR`) for a definitive count,
or a shifted `svds` probe on `J_sᵀ` for the smallest left singular vector `u` when a full
factorization is too costly. Then rank equations by `|u_i|`. This is the independent confirmation
of C2: **the equations with large `|u_i|` should be the same families C2 flags.** Two methods
agreeing is the standard the gate 5c-ii work already set.

**C4 — Repair, then re-verify.**
For each family C2/C3 names, read its builder in `build_equations.jl` and check the guard obeys the
invariant: *a guarded zero-flow cell must remove exactly one equation AND exactly one unknown.*
The canonical correct form is `E_plnd!`:
```julia
if LND[i,d] > 1e-10
    @constraint(m, xlnd[i,d] == alnd[i,d] * ces(...)[3])
else
    fix(plnd[i,d], 1.0; force=true)
end
```
Failure modes to look for: (a) `@constraint` emitted unconditionally while the variable is pinned →
over-determined row (**the predicted defect**); (b) `|| continue` skipping the constraint without
pinning the variable → under-determined column (the gate 5c-ii defect, already fixed in three
places — check the other nine). Then: re-run the benchmark (must stay ≤1e-9), re-run
`blabnat=0.9997` and confirm the floor is gone.

**C5 — Only if C2/C3 exonerate the guards.** Then the obstruction is genuine ill-conditioning at
the perturbed point, and the right lever is **Excerpt 49 condensation** (already required
independently for 34 regions) — eliminating `xtradmar`/`xsuppmar`/`ppur` removes the near-dense
aggregation rows that dominate both fill-in and conditioning. Not solver tuning.

**C6 — Then, and only then, resume the shock programme.**
Price-homogeneity test (the numeraire swap is already in place, see Active #3) → `blabnat = −3`
(the *real* `TERM.CMF` shock, 100× the test size) via Euler continuation with a full Newton
correction per step → a smoke-test scenario checked for sign/magnitude plausibility.

### E. ~~Corrected verdict~~ — **SUPERSEDED BY §F**

> The conclusion "the direction is dominated by a near-singular mode" is **wrong**, and the
> `‖d‖∞/‖b‖₂ ≈ 4e5` amplification that motivated it was measured on a **corrupted step**: at the
> time, `solve_newton!` converted the direction back to original units with `dy .* Dc_ruiz`, but
> `ruiz_scales` returns `Dr`/`Dc` as **divisors**, so every component was off by `Dc²`. With the
> correct `dy ./ Dc` the step satisfies `‖J·d + g‖₂/‖g‖₂ = 3.2e-10` and moves every variable by
> `O(shock)` except three, all explainable. §E's *measurements* stand; its *diagnosis* does not.
> Read §F instead. The rest of §E is kept as the record of what was believed and why.

### E (superseded). Corrected verdict (2026-07-29, measured) — the system is CONSISTENT; the step is the problem

§B was wrong. Three measurements settle it, all at the *stalled* `blabnat = 0.9997` point
(`test/scratch/diag_stall_nullspace.jl`, `test/scratch/diag_stall_direct.jl`):

| measurement | value | meaning |
|---|---|---|
| `lu(J_s)` | succeeds, 2.3 s, `nnz(L)+nnz(U)` = 3,391,666 | no structural singularity |
| `‖J·d − b‖₂ / ‖b‖₂` (direct LU) | **3.5e-10** | `F` **is** in `range(J)`; the system is **consistent** |
| columns with `‖·‖ < 1e-12` / `< 1e-3` | **0 / 0** | no empty columns; the zero-flow guards are clean |
| `‖J·d − b‖₂ / ‖b‖₂` (CGNR, μ=0, 4000 it) | 0.161 | **non-convergence**, not a left null space |

**Why §B's two premises were false.**

1. *A2, the "bit-identical floor across a changed code path."* The change was `solve_newton!.jl:435`
   `colnrm < 1e-3` → `< 1e-12`. But `colnrm` is measured on the **Ruiz-equilibrated** `J_s`, where
   every non-empty column has infinity-norm exactly 1. Both thresholds therefore select the identical
   set — the empty columns, of which there are **zero**. The edit was a **no-op by construction**, so
   bit-identical output is the expected result and carries no information about structure. The same
   argument retires the "second un-relaxed freeze at `:461`" concern: it is equally inert.
2. *A3, the CGNR least-squares residual.* CGNR's convergence rate is governed by κ(J)², and the LU
   step reveals an amplification `‖d‖∞ / ‖b‖₂ ≈ 1221 / 0.0033 ≈ 4e5`. 4000 iterations were nowhere
   near enough. A large CGNR leftover on an ill-conditioned matrix says nothing about rank — that is
   exactly why `diag_stall_direct.jl` exists as the control.

**What is actually wrong.** The exact Newton step is enormous and useless:

```
‖b‖₂ = 0.0033   →   ‖d‖∞ = 1221        (amplification ≈ 4e5)
```

and the line search finds **no descent at any α**, uncapped or capped. `‖F‖∞` falls off perfectly
linearly in α (6.79 → 3.39 → 1.70 → … , exactly halving) until α ≈ 8.5e-6, where it flattens onto the
floor 2.126e-4 — i.e. every α large enough to make first-order progress is already far outside the
region where the linearization holds, and every α small enough to be safe does nothing. The direction
is dominated by a **near-singular mode**, not a null one.

Corroborating this from the residual side: the dominant raw-residual family at the stalled point is
`alab_o` (46.5% of squared magnitude) at 2.13e-4, against a shock that demands a move of
`|log(0.9997)| = 3.0e-4`. The model has travelled ~30% of the way to the answer and stopped. This is
a step-taking failure, not an obstruction.

**Consequence for the plan.** C4 (repair zero-flow guards) is **moot** — there are no empty columns
and no inconsistency to repair. C5 (condensation) is not the fallback either; it addresses scale, not
conditioning. The live question is narrower: *which* near-singular mode, and why. `σ_min`, `v_min`
(the nearly-undetermined variables) and `u_min` (the nearly-dependent equations) are extracted by
inverse iteration in `test/scratch/diag_stall_nearnull.jl`; the fix follows from what that block names.

Note the §A5 line "weak-column freeze ruled out" was right for the wrong reason, and §A4's μ sweep is
now readable as the correct instrument pointed at the correct problem: LM regularization *is* the
right family of treatment for a near-singular Jacobian, and no μ helped because μ was being chosen
globally against a mode that needs to be removed structurally rather than damped.

### F. Current verdict (2026-07-29, later) — **the direction is right; the ACCEPTANCE TEST was wrong**

Four defects were found. The first three are fixed and verified; the fourth is the actual stall.

| # | defect | evidence it was real | status |
|---|---|---|---|
| 1 | `dy .* Dc_ruiz` instead of `dy ./ Dc_ruiz` — `ruiz_scales` returns divisors | round-trip `‖J·d+g‖₂/‖g‖₂` went from O(1) to **3.2e-10** | fixed |
| 2 | `ppur_s` zero-flow guard at `1e-10`, but the smallest positive `PUR_S` is 2.7e-8, so it only ever pinned exact zeros. `PUR_S[5,16,6]=9.95e-8` gave a price index amplified by ~1e7 | fraction-to-boundary `a0` went **0.0012 → 1.0**; non-finite rows → **0** at every α | fixed (`_MIN_PRICED_FLOW = 1e5·_TINY`) |
| 3 | line search evaluated a direction CGNR had not converged to (κ≈2.6e10 ⇒ rate ~κ²) | CGNR `max|dy| = 0.025` against a true **1221**; `cgnr_rel = 0.101` | fixed (direct sparse LU first) |
| 4 | **acceptance test = strict descent on `‖F‖∞`, with row scaling recomputed every iteration** | see below | **fixed, under verification** |

**Why #4 is the stall.** With `xcap`/`xlnd` exogenous in the base closure, the sector-specific fixed
factors of Livestock/BaliNusa absorb the whole adjustment in their *price* — so a 0.03% labour shock
moves `pcap[2,5] = plnd[2,5]` by −4.65%. (They move bit-identically because inverting two CES demands
with both quantities pinned gives `(pcap/plnd)^σ = const`. Expected, not a bug.) Row 38476,
`gret[i,d] = log(pcap+τ) − log(pinvitot+τ)`, therefore picks up ½(0.046)² ≈ **1.1e-3** of pure
log-curvature at α=1 — larger than the entire starting residual of 3.0e-4. It is the argmax of `‖F‖∞`
at every α ≥ 0.5. One row out of 83,515 vetoed the full step; the search settled on α=0.25, and α
shrank every iteration thereafter. That is linear convergence with a decaying rate — the stall.

Compounding it, `compute_rowscale!` was called **every iteration** (`:404`) while `‖F‖∞` was compared
**across** iterations (`:579`, `:584`). The yardstick changed length between readings.

**The fix (both now in `solve_newton!.jl`).** Freeze the row scaling for the whole solve — Ruiz
re-equilibrates the Jacobian each iteration anyway, so conditioning is unaffected. Accept steps on
**½‖F‖₂²** against a **non-monotone (Grippo) reference** — the worst merit in a 5-iteration window —
so a curvature overshoot can be corrected at the *next* iteration, which is what Newton does. `‖F‖∞`
remains the convergence criterion. The best iterate is tracked and restored, since a non-monotone
search may end on an uphill excursion. The CGNR retry after a rejected *exact* LU direction is
removed: at this conditioning it cannot produce a better direction than the one just rejected.

**Two hypotheses killed by measurement here — record them so they are not re-run.**
- *"The tiny CES coefficient 3.8e-4 is a degenerate capital value share (data defect)."* **False.**
  Capital's value share is 0.112–0.818 (median 0.398) across all 150 cells at 6 regions and
  0.0999–… (median 0.393) across all 850 at 34; **zero** cells below 0.01 at either scale
  (`test/scratch/diag_primary_shares.jl`). `0.00038` is a calibrated CES `ALPHA`, not a share.
- *"Undamped Newton will converge quadratically, so the line search is the whole story."* **Also
  false, but for a reason that was a measurement error on our side.** `max|d/x|` appeared to explode
  to 574× at iteration 2 — but every "runaway" was a `delTAXint`/`delPTX`/`delTAXhou` *additive*
  shift variable legitimately sitting at 1e-7…1e-27, where a relative step is meaningless. Likewise
  `xlab_o` spanning 0.012→801,678 is just the benchmark labour value-flow spread. The genuine
  economics stayed sane throughout (`phi`=1.00031, `realwage`∈[0.998,1.002], `xtot`∈[0.987,1.006]).
  **Nothing diverges.** Never rank a step by `|d/x|` on a variable whose benchmark is 0.

### G. Gate result for §F, and the defect it exposed (2026-07-29, later still)

The §F fix was implemented and gated with `test/verify_shock_fix.jl` (`scratchpad/gate.log`).
**The gate FAILED**, but it failed informatively, and part of §F is confirmed.

| | result |
|---|---|
| benchmark | ✅ `converged`, ‖F‖∞ = 1.4260876923799515e-9 — **no regression** |
| `blabnat = 0.9997` | ❌ `no_progress`, 12 iters, ‖F‖∞ = 2.649e-4 |
| `blabnat = 0.99` | ❌ `no_progress`, 1 iter, ‖F‖∞ = 9.752e-3 |

**Confirmed by the log — the merit change did exactly what §F predicted.** Iterations 1 and 2
of the 0.9997 shock now accept **α = 1.0** (the old monotone ∞-norm test settled on α = 0.25),
and the merit falls `6.752e-6 → 1.709e-6 → 5.26e-7`. Iteration 1's `ginf` *rises*
(3.0e-4 → 1.1e-3) and iteration 2 then absorbs it — the curvature overshoot §F described,
behaving as Newton is supposed to. The `elseif exact` early break and the best-iterate restore
both fired correctly.

**Stated honestly: on the ∞-norm, 0.9997 is now WORSE.** Prior code reached 7.03e-5 at
maxit(30); this reaches 2.649e-4 at 12. The merit fell, the reported norm did not. Do not
present the α=1 acceptance as a net win.

**The new defect — the Jacobian goes near-singular at later iterates.** From iteration 3 the
exact step magnitude grows monotonically while the *linear solve residual degrades*:

```
it  1  max|d|=1221    lin_rel=3.19e-10  α=1.0
it  2  max|d|=556     lin_rel=2.01e-10  α=1.0
it  3  max|d|=1670    lin_rel=9.58e-10  α=0.0625
it  4  max|d|=10760   lin_rel=8.81e-9   α=9.8e-4
it  6  max|d|=289300  lin_rel=1.63e-7   α=2.0e-6
it  7  LU REJECTED → CGNR+LM(mu=1e-7): max|d|=1.171  α=0.55  merit 4.99e-7 → 1.036e-7
it  8  max|d|=415500  lin_rel=5.51e-7   α=1.4e-6
it 12  max|d|=916800  lin_rel=9.73e-7   → no merit decrease at any α
```

An LU whose **relative** residual grows on the same-sized problem is the signature of genuine
near-singularity, not of bad scaling. Note the 0.99 shock does **not** start from the benchmark:
its initial residual 0.0097520 = |log 0.99 − log 0.9997|, i.e. it resumes from the 0.9997
stall point, which is why its very first step is already max|d| = 7.768e7.

**§E is therefore partially rehabilitated.** §E's "near-singular mode" was refuted *as an
explanation of the benchmark start point* — correctly, since there `lin_rel = 3.2e-10` and
`max|d| = 1221` are clean. It is the right description of the **late iterates**. Both statements
hold; they are about different points.

**The decisive datum is iteration 7.** It is the only iteration that fell through to damped
CGNR, and it produced the smallest step and the largest merit drop of the entire run. Near a
singular point a damped step beats the exact one *because* it is damped — the line search can
only rescale a direction, never redirect it, so an exact step aligned with the near-null mode
is unusable at every α.

**Three fixes applied in response** (`src/solve_newton!.jl`):
1. **Removed the `elseif exact` break.** A rejected exact ray now falls through to damped CGNR
   instead of stopping. The removed reasoning ("CGNR is garbage at κ ≈ 2.6e10") was measured at
   the *benchmark start point* and does not transfer to a near-singular iterate — iteration 7
   is the counterexample.
2. **μ adaptation now keyed on step QUALITY, not acceptance.** Old rule decayed μ whenever the
   merit fell at all; accepting α ≈ 1e-6 for a 1e-4 relative change drove μ to 1e-8 — no damping
   at all — exactly as J was going singular. New rule: decay only on `α > 0.5 && m < 0.9 m₀`;
   *increase* μ on `α < 0.1`.
3. **Trust region on ‖d‖∞.** An accurately-solved direction exceeding the radius is demoted to
   the damped path before the line search spends 20 residual evaluations rediscovering it.
   Radius follows the accepted move: doubled after a full step, pulled back to `α‖d‖∞` otherwise.

**Rule to carry forward:** *a small `lin_rel` certifies the linear solve, not the direction.* An
exact step can be useless when J is near-singular. Judge a direction by the move it produces,
not by the accuracy with which it was computed.

#### G2. Gate run 2 result (`scratchpad/gate2.log`) — the fixes work; a NEW floor at 1.7e-5

| | gate 1 | gate 2 |
|---|---|---|
| benchmark | 1.4260876923799515e-9 ✅ | **identical** ✅ |
| `blabnat = 0.9997` | `no_progress`, 12 it, 2.649e-4 | `maxit`, 30 it, **1.687e-5** |

**15× better than gate 1, and 4× better than the pre-fix best of 7.03e-5.** Status is `maxit`
(still moving), not `no_progress` (stuck). Iterations 1–4 are near-quadratic Newton —
merit `1.709e-6 → 5.26e-7 → 9.487e-8 → 3.269e-9`, ‖F‖∞ `1.1e-3 → 1.85e-5`. The near-singular
explosion is gone: the exact step magnitude is now **stable at ~1650** for all 30 iterations
instead of running away to 9.2e5. The drift diagnosis in §G was right.

**New defect, mine: the trust radius could contract but not expand.** Δ was grown only on
`α ≥ 0.99`, and the damped path essentially never returns α = 1, so Δ ratcheted monotonically
down — 1112 → 0.67 → … → 0.012 — against an LU step that stayed at ~1650. The exact path was
locked out for 28 consecutive iterations and the merit went flat at 2.5e-9. A trust region that
cannot expand after a good step is just a decaying step cap. Fixed: Δ now keys on the achieved
merit reduction (`red > 0.1` → `Δ ← max(2·moved, 2Δ)`; `red > 0.01` → hold; else contract).
`maxit` in the gate script raised 30 → 60, since iterations 1–4 were still productive.

**Open question for the next measurement, not to be guessed at.** From iteration 16 the merit
still falls (2.66e-9 → 2.53e-9) while `‖F‖∞` slowly *rises* (1.69e-5 → 1.82e-5) — the L2 merit
is shaving many small rows at the cost of the largest one. Either 1.7e-5 is a genuine floor in
one specific row, or it is the trust-region lockout above. Gate 3 separates these; if the floor
survives, run `test/scratch/diag_full_newton.jl` (task #19) to NAME the row carrying it. Do not theorise
about which row before measuring — three prior diagnoses in this section were refuted by
measurement.

#### G3. Gate run 3 (`scratchpad/gate3.log`) — ANSWER: the 1.7e-5 floor is REAL

`blabnat = 0.9997`: `maxit`, 60 iterations, ‖F‖∞ = **1.641919082126898e-5** (1160s).
Against 1.687e-5 at 30 iterations in gate 2: **52 extra iterations bought 2.7%.**

The Δ-expansion fix works exactly as designed — the radius now doubles on good steps
(10.51 → 21.02 → 42.04 at iterations 5–7) and iterations 5–8 take α = 1 and drop ‖F‖∞
from 5.98e-4 to 1.77e-5 in four steps. The residual then goes flat regardless. So the floor is
**not** a trust-region artifact, and it is not a line-search artifact either.

**Stop tuning the solver.** Three successive solver fixes (merit function, LM fallback +
μ-adaptation, trust region) moved the floor 2.6e-4 → 1.7e-5 and each was justified by
measurement, but the residual is now pinned by something that is not the step rule. Note also
that after iteration 8 the merit keeps falling (2.88e-9 → 2.29e-9) while ‖F‖∞ *rises*
(1.77e-5 → 1.83e-5): the L2 merit is shaving many small rows at the largest row's expense,
which is what a merit minimum that is **not a zero of F** looks like.

**Next measurement (task #19, no theorising first):** run `test/scratch/diag_full_newton.jl` with
`STALL_MAXIT=25` to reach the same floor under undamped Newton and print, per iteration, the
argmax row *by name and index* plus the residual mass by equation family. That names the
equation carrying the 1.7e-5. Only after it is named is it worth asking whether it is a
translation defect in that block or a genuine degeneracy of the closure.

### H. LOCALISED (2026-07-29) — 97% of the residual is ONE cell: Livestock × BaliNusa

`test/scratch/diag_full_newton.jl`, `STALL_MAXIT=25` (`scratchpad/rows2.log`). Every top-5 residual
row, at every iteration, is index **[2,5]**:

| share of ‖F‖₂² | family | variables |
|---|---|---|
| 51.7% | `alab_o\|pcap\|plab_o\|plnd\|xprim` | `plab_o[2,5] pcap[2,5] plnd[2,5]` |
| 36.0% | `gret\|pcap\|pinvitot` (row 38476) | `pcap[2,5] pinvitot[2,5] gret[2,5]` |
| 4.71% | `ggro\|xinvitot` | `xinvitot[2,5] ggro[2,5]` |
| 4.71% | `finv2\|xgdpexp\|xinvitot` | `xinvitot[2,5] finv2[2,5] xgdpexp[5]` |

**97% in sector 2 × region 5 = Livestock / BaliNusa**, stable across all 25 iterations. The
common variable in every family is `pcap[2,5]`. Undamped Newton does not converge here at all —
it oscillates (3e-4 → 5.5e-3 → 1.1e-2 → 8.4e-4) with `lin_rel ≈ 1e-10` throughout, so the
linear algebra is exact and the oscillation is genuine nonlinearity in this one cell.
`pcap` bottoms at 0.914 — an **8.6% price swing off a 0.03% shock**. This promotes the
`gret[2,5]` row from §F (then only "the row that vetoes the step") to the cell carrying
essentially all remaining error.

**Why this cell (`test/scratch/diag_cell25.jl`, measured).** It is the joint extreme of two ratios:

- **Fixed-factor share 0.172, rank 6 of 150** (median 0.457, max 0.818). CAP=1543.5, LND=686.0,
  LAB=10752 — the cell is 83% labour. With `xcap`/`xlnd` exogenous in the base static closure,
  the fixed factors' *prices* must absorb the entire adjustment, and a small fixed-factor share
  is large price leverage: ≈1/0.17 ≈ 5.8× versus ≈2.2× for a median cell. That is the 8.6% swing.
- **`CAP_v/INVEST_C` = 0.142, rank 0 of 150** (median 1.476). Investment is 7× the capital
  rental. Notably **all six lowest cells are sector 2** — [2,5], [2,6], [2,1], [2,4], [2,3],
  [2,2] occupy ranks 0–5 — so this is a Livestock-wide pattern, with BaliNusa its worst cell.
  Whether that is genuine or a pipeline allocation defect is NOT yet established; flagged, not
  concluded. (Related known defect in the same area: `build_pstras!.jl:173 DISTANCE = reg1["MAKE"]`.)

**Hypothesis killed on the way (the 4th in this section).** `E_pinvitot!` skips cells with
`INVEST_C[i,d] <= 1e-10` (`build_equations.jl:572`), which would leave `pinvitot` with no
defining equation while `E_gret!` still writes a row for every (i,d). Measured: **it skips 0 of
150 cells**, and `INVEST_C[2,5]` = 10,880 ranks 86/150 at 1.39× the median. No structural hole.

**Open, and the test now running (`test/verify_continuation.jl`).** Extreme data makes the solve
*hard*; it does not by itself make it *wrong*. A graduated ladder 0.999999 → 0.99999 → 0.9999 →
0.9997, each step starting from the previous solution, separates the readings: every step
converging to `tol` means the equations are right and the remedy is continuation (C6b, already
planned) rather than more solver tuning; a 1e-6 shock still stalling at ~1e-5 means the floor is
not step-size and the [2,5] block holds a real defect. **Do not act on the extreme-data story
until this returns** — it is currently correlation (rank 0 and rank 6 of 150 in exactly the
guilty cell), and the mechanism has not been confirmed.

### I. ✅ RESOLVED (2026-07-29) — the equations are RIGHT; shocks converge under continuation

`test/verify_continuation.jl` returned, and it settles §H. **The project's first ever converged
shocks:**

| `blabnat` | status | iters | ‖F‖∞ | `pcap[2,5]` | secs |
|---|---|---|---|---|---|
| 0.999999 | **converged** | 2 | 4.19e-9 | 0.99984518 | 62.6 |
| 0.99999  | **converged** | 2 | 3.26e-9 | 0.99843706 | 59.5 |
| 0.9999   | **converged** | 3 | 5.01e-9 | 0.98267596 | 77.8 |
| 0.9997   | maxit | 40 | 1.38e-5 | 0.9263198 | 881.8 |

Two or three Newton iterations to 1e-9 is not the behaviour of a mis-translated equation. **The
[2,5] equation-block audit is cancelled** — `E_gret!` / `E_xinvitot!` / `E_ggro!` / `E_finv2!` are
not to be edited on the strength of §H's data ranking. §H's "extreme data" story is confirmed as
*correlation with difficulty*, never as a defect. The remedy is continuation, which was already
planned as C6b.

**The script's own printed verdict ("the floor is not step-size") is OVERSTATED — ignore it.** It
fires on any incomplete ladder. Its own data does not support it: only the coarsest rung failed,
and that rung was a 3× jump (1e-4 → 3e-4) after three clean passes. What the data *does* show is
`pcap[2,5]` deviation per unit shock running **155× → 156× → 173× → 246×**, i.e. accelerating
superlinearly. That, combined with a Jacobian going near-singular *at the solution* (§E/§F), is
the signature of a **fold (limit point)** on the solution branch. A fold cannot be crossed by
natural-parameter continuation at any step size; it needs pseudo-arclength (Keller) continuation,
which parameterises by arclength and stays nonsingular through the turn.

**Discriminating measurement, running now:** a fine ladder 0.9999 → 0.99985 → 0.9998 → 0.99975 →
0.9997 (`LADDER` env var added to `verify_continuation.jl`). Walking all the way through ⇒ the
coarse break was only step size and plain continuation is enough. Converging and then stopping
abruptly with the response still accelerating ⇒ fold, and the driver needs arclength. Five
hypotheses in this investigation have already died to measurement — this one is bracketed, not
believed.

**Built in anticipation (needed under either outcome):**
- `src/continuation.jl` — `continuation_solve!(m, vars, shock_vr, target; h0, hmin, hmax, tol,
  maxit, grow_iters, report)` → `ContinuationResult`. Natural-parameter continuation with an
  adaptive step: grow after an easy step (≤ `grow_iters`), halve and roll back after a failure.
  Rollback is real — it snapshots and restores every start value, because `solve_newton!`
  overwrites them in place and a rejected step would otherwise strand the model on a
  non-solution. Exported from `IndotermJulia.jl`. This supersedes the legacy
  `euler_continuation!` in `solve_newton!.jl:825`, which carries its own inferior inner solver
  (10 % step cap, monotone ∞-norm acceptance, `colnrm` regularisation) predating every fix in
  §F/§G.
- `test/verify_continuation_driver.jl` — the gate. Benchmark must hold at ~1.4e-9, then
  `continuation_solve!` must reach `TARGET` (default 0.9997). On failure `s_reached` is the
  furthest genuinely *solved* point (residual ≤ tol, not a stalled iterate), which is exactly
  what locates a fold.

### J. FOLD CONFIRMED (2026-07-29) — and it is very likely a CLOSURE artifact, not economics

**Fine ladder** (`LADDER` env var added to `verify_continuation.jl`):

| `blabnat` | status | iters | ‖F‖∞ | `pcap[2,5]` |
|---|---|---|---|---|
| 0.99985 | converged | 4 | 2.68e-9 | 0.97186606 |
| 0.9998  | converged | 3 | 2.21e-9 | 0.95785452 |
| 0.99975 | maxit | 40 | 7.66e-7 | 0.92623464 |

Two more rungs converged in 3–4 iterations, so §I's conclusion stands and hardens. The
obstruction is a **fold (limit point)**, on three independent grounds:

1. **Square-root signature.** Local slope `d(pcap dev)/d(shock)` is **216** over 1e-4 → 1.5e-4
   and **280** over 1.5e-4 → 2e-4. Near a fold the response is `A − B·√(δ*−δ)`, so slope
   `∝ 1/√(δ*−δ)`. Two slopes fix both unknowns: `280/216 = √((δ*−1.25e-4)/(δ*−1.75e-4))` ⇒
   **δ* = 2.485e-4, `blabnat*` = 0.999752**.
2. **Out-of-sample hit.** The ladder converged at 0.9998 and broke at 0.99975. The break
   location was *not* used in the fit.
3. **Same stall point regardless of target.** The fit also gives `pcap[2,5]*` = 0.9249; the
   0.99975 run parks at **0.92623** and the 0.9997 run at **0.92632** — 0.14 % apart, and
   independent of where they were aiming. Being unable to pass a fixed point no matter what you
   aim at is what a turning point does. (Consistent with §E/§F's near-singular Jacobian *at the
   solution*, which a fold requires.)

Natural-parameter continuation cannot cross a fold at any step size. `continuation_solve!` will
therefore halve to `hmin` there — by design, that is how it *reports* the fold.

**But a fold at a 0.025 % labour shock is not economics — it is this closure.** The decisive
comparison is with the source's own scenario:

> `TERM.CMF:112` — `Shock blabnat = -3;   ! 3% increase labour productivity`

**The reference scenario is −3 %. We fold at −0.025 %, a factor of ~120 short.** No plausible
reading of the Indonesian economy has the equilibrium set ending a quarter of a basis point from
the benchmark. So the branch geometry is being imposed by the configuration, and reading
`TERM.CMF:69–113` shows exactly how — **the port applies only ONE of its seven active swaps:**

| `TERM.CMF` swap | line | in port? |
|---|---|---|
| `phi = NatMacro("GDPPI")` | 72 | ✅ applied |
| `xhouhtot = fhou` | 79 | ❌ |
| `houslack = shrBoTnom2` | 82 | ❌ |
| `fgovtot = fgovtot3` | 86 | ❌ |
| **`xcap = faccum`** (switch on capital accumulation) | **88** | ❌ |
| `finv1 = finv4` (dynamic investment rule) | 89 | ❌ |
| `delfwage = flabsup_id` (national labour adjustment) | 110 | ❌ |

Line 88 is the one that bites: **TERM.CMF runs the `blabnat` shock with `xcap` ENDOGENOUS.** Our
`BASE_CLOSURE_ARRAYS` fixes `xcap` *and* `xlnd`, so every fixed factor is frozen in all 25
sectors and the entire adjustment has to land in fixed-factor *prices* — with ~1/0.172 ≈ 5.8×
leverage in the [2,5] cell, measured at ~200× amplification. That is precisely the mechanism §H
identified, and it is switched on by a closure choice the source does not make for this shock.

Note also that `xlnd` has **no** swap in any source file, so freezing land IS faithful; freezing
capital is not. And `termsr.cmf` — the genuine short-run file — applies four further swaps of its
own (`invslack = NatMacro("RealInv")`, `xhouhtot = fhou`, `flab_i = flabsupA`,
`labslack = flabsup_id`), none of them in the port either. In other words the port's "base static
closure" is the raw TABmate automatic `Exogenous` list plus the numeraire swap — the state
*before* any closure choice is made, which is not a configuration anyone actually runs.

All seven swap partners (`faccum`, `finv4`, `fhou`, `shrBoTnom2`, `fgovtot3`, `delfwage`,
`flabsup_id`) are already declared in the port, so the swaps are implementable without new
equations.

**Consequence for the plan.** Do **not** write pseudo-arclength continuation yet. Arclength can
trace a branch around a fold; it cannot manufacture an equilibrium beyond one, and if the fold is
an artifact it would be sophisticated machinery aimed at the wrong target. The order is:
(1) confirm directionality (`test/scratch/diag_fold.jl`, running: `blabnat` → 0.9997 vs → 1.0003 —
one-sided fold vs symmetric obstruction); (2) apply TERM.CMF's active swaps and re-measure the
fold location; (3) only if a fold survives a faithful closure does arclength become the answer.

### K. §J's step (2) is DONE — the swap is right, and it moves the problem, not solves it

`test/verify_swapped_closure.jl`, 2026-07-29. Two results that must be held together:

**The swap implementation is correct.** All six of TERM.CMF's missing active swaps applied with
the expected shapes, and the benchmark still solves to
`‖F‖∞ = 1.4260876923799515e-9` — *bit-identical* to the base closure. That is the regression gate
every closure change has to pass, and it passed exactly. The reversed pair printed as
`swap flabsup_id → endogenous, delfwage → exogenous`, confirming that reading TERM.CMF:110
left-to-right would have produced the wrong closure and that the symmetric implementation was
necessary, not defensive.

| swap | shapes | exogenous side in the port |
|---|---|---|
| `xhouhtot = fhou` | 6 = 6 | `xhouhtot` |
| `houslack = shrBoTnom2` | 1 = 1 | `houslack` |
| `fgovtot = fgovtot3` | 6 = 6 | `fgovtot` |
| `xcap = faccum` | 150 = 150 | `xcap` |
| `finv1 = finv4` | 150 = 150 | `finv1` |
| `delfwage = flabsup_id` | 4 = 4 | **`flabsup_id`** ← reversed |

**But the continuation walk accepted ZERO steps** (`s_reached = 1.0`, `steps=0 rejects=6`, 796s):

```
h=0.5     ‖F‖∞ = 3.15e-5
h=0.25           6.56e-5   ← WORSE at a smaller step
h=0.125          6.70e-6
h=0.0625         2.22e-6
h=0.031          9.25e-7
h=0.016          2.71e-7
```

**This is NOT a fold, and the script's own printed verdict ("the fold survives a faithful
closure") was wrong** — it has since been corrected to detect `nsteps == 0` and say so. A fold is a
property of a point *on* the branch; s=1.0 is the benchmark, which solves to 1e-9. Two signatures
rule it out independently: the leftover residual falls only like `h^1.5` instead of collapsing to
tolerance, and it is **non-monotone in `h`**. A shock of ~1e-6 — numerically indistinguishable from
the benchmark — failing in 15 Newton iterations means the **Jacobian is near-singular at the base
year** under this closure. The solution exists; Newton cannot find the direction to it.

So §J's question ("does the fold survive a faithful closure?") is **not yet answerable**, because
the faithful closure cannot take even an infinitesimal step. The fold measurement at
`blabnat* = 0.9997525` stands as a property of the *base* closure and nothing more.

**Next: `test/scratch/diag_swap_bisect.jl`** applies each swap **alone** (not incrementally — an
incremental ladder confounds the culprit with its position in the ordering) and probes with a 1e-5
shock, which the base closure clears in 2 iterations. Prior suspicion, to be confirmed or killed
rather than assumed: `xcap = faccum` and `finv1 = finv4` each free 150 variables against equations
declared in the **dynamic** block (`build_dynamics!.jl:64-65`); in a single-period static solve
those accumulation relations may pin `xcap` only weakly, which is what a near-zero pivot looks
like. If every singleton is benign the damage is an interaction and the ladder has to be run
incrementally after all.

Arclength remains **not** the answer, now for a second and stronger reason: the obstruction under
the faithful closure is at the base year, where arclength has nothing to trace.

### L. LOCALISED — ONE scalar swap costs the conditioning, and my prior suspicion was WRONG

`test/scratch/diag_swap_bisect.jl`, 2026-07-29. Each swap applied **alone** to a fresh model (not
incrementally — an incremental ladder confounds the culprit with its position in the ordering),
probed with `blabnat = 0.99999`, which the base closure clears in 2 Newton iterations. Every case
kept the benchmark at exactly `1.4260876923799515e-9`, so all eight rows below are comparable.

| case | status | iters | ‖F‖∞ |
|---|---|---|---|
| none (control) | ✅ converged | 2 | 8.50e-9 |
| `xhouhtot = fhou` | ✅ converged | 2 | 1.63e-9 |
| **`houslack = shrBoTnom2`** | ❌ **maxit** | **20** | **3.24e-6** |
| `fgovtot = fgovtot3` | ✅ converged | 2 | 4.66e-9 |
| `xcap = faccum` | ✅ converged | 3 | 2.10e-9 |
| `finv1 = finv4` | ✅ converged | 2 | 9.20e-9 |
| `delfwage = flabsup_id` | ✅ converged | 2 | 8.50e-9 |
| ALL SIX | ❌ maxit | 20 | 1.03e-6 |

**§K's stated suspicion is refuted.** `xcap = faccum` and `finv1 = finv4` — the two that free 150
variables each against dynamic-block equations, and the two §K named as prior suspects — converge
in 3 and 2 iterations. The 150-variable swaps are harmless; a *single scalar* is not. Recorded
because the prediction was explicit and wrong, and the wrong prediction was the plausible one.

`ALL SIX` fails with the same signature as the singleton, so there is no interaction term to hunt:
one scalar accounts for the whole §K failure. **Therefore the other five swaps — including the
`xcap = faccum` that §J predicted would remove the fold — are usable now**, and the fold
re-measurement is running under `EXCLUDE=houslack` (`test/verify_swapped_closure.jl` gained an
`EXCLUDE` env knob for exactly this).

**Mechanism (predicted, being measured by `test/scratch/diag_houslack_null.jl`).** `TERM.TAB:2047-2050`:

```
E_fhou  (all,h,HOU)(all,d,DST)  whouhtot(h,d) = wlab_io(d) + fhou(h,d)  + houslack;
E_fhou2 (all,h,HOU)(all,d,DST)  whouhtot(h,d) = wgdpexp(d) + fhou2(h,d) + houslack;
```

`houslack` enters both blocks additively with coefficient 1 in every row, alongside the per-region
shifters. So whenever `houslack`, `fhou` and `fhou2` are all free,

    col(houslack) = Σ_d col(fhou[d]) + Σ_d col(fhou2[d])

**exactly** — the entries are literal 1s from an additive term, so the cancellation is exact in
floating point, not approximate. While `houslack` is exogenous this is inert. Free it and the
Jacobian is singular: the benchmark stays a solution (it always was), but Newton has no direction.
That is the observed signature precisely — bit-identical benchmark, zero accepted steps, residual
falling like h^1.5 and non-monotonically in h.

Note this is a property of the TABLO source, not a porting slip: both equations are in `TERM.TAB`
with `houslack` in each.

**The theory does not yet cover `ALL SIX`**, and that is stated rather than papered over: there
`xhouhtot = fhou` pins `fhou`, which should break the exact dependency — yet ALL SIX still failed.
`diag_houslack_null.jl` tests both cases by direct null check (‖J·v‖∞ against ‖ |J|·|v| ‖∞, the
magnitude before cancellation), so case B either exposes a second mechanism or shows the remaining
dependency is near-exact rather than exact.

**Measured (`test/scratch/diag_houslack_null.jl`, same day). Case A: PROVEN EXACTLY.**

```
houslack = shrBoTnom2 ALONE     ‖J·v‖∞ = 0.0        ‖ |J|·|v| ‖∞ = 2.0   ratio = 0.0
ALL SIX                         ‖J·v‖∞ = 1.0        ‖ |J|·|v| ‖∞ = 2.0   ratio = 0.5
```

`‖J·v‖∞` is not small in case A, it is **identically zero**. The singleton closure is
structurally singular, exactly as derived — no conditioning argument needed and none available.

Case B came out as predicted (not a null direction), and the way it failed is a precise lead
rather than a dead end: the residual is exactly `1.0` in six consecutive rows, indices
63122–63127, which is the `E_fhou` block carrying `houslack`'s uncompensated unit coefficient. So
**ALL SIX stalls for a different reason**, and three options remain open, in decreasing
plausibility: (a) it is merely HARD — the bisect's `maxit=20` was too small and the residual was
still falling at 1.03e-6 when it stopped; (b) singular along a different direction; (c) genuinely
near-singular. `test/scratch/diag_allsix_deep.jl` runs 60 iterations at a 1e-6 shock with the per-iteration
trajectory visible, which separates (a) from (b)/(c) without writing a null-space solver on spec.

**Bookkeeping caveat — RESOLVED, nothing to see.** The null run reported `83515 rows × 83516 free
columns`, one wider than the familiar square figure. The control (no swaps at all) reports
**83516 as well**, so the extra column is the pre-squaring orphan that `solve_newton!` pins, and it
is not caused by any swap. No under-determination.

### M. ALL SIX is NEAR-singular, not singular — and the closure question is now economic

`test/scratch/diag_allsix_deep.jl`, 2026-07-29. ALL SIX at a **1e-6** shock:
**converged, 2 iterations, ‖F‖∞ = 3.84e-9.** So option (b) is dead — there is no second exact null
direction — and the bisect's `maxit=20` failure at a 1e-5 shock was a *radius-of-convergence*
failure, not a rank failure.

But the script's own printed verdict, "(a) HARD, NOT singular", is **too generous**, and the trace
says why:

```
it 1 lmtry=1 mu=0.1 : max|d| = 36.61   a0=1.0   lin_rel=4.46e-8
it 2 lmtry=1 mu=0.01: max|d| =  2.70   a0=1.0   lin_rel=2.4e-9
```

My first reading of that line — "a step of 36.6 against a shock of 1e-6 is ~4×10⁷ amplification,
therefore option (c)" — was **wrong, and the movers run refuted it**. The ratio is meaningless:
`max|d|` is in *variable units* (mostly Rp billion) and the shock is dimensionless. This is the
same class of error as the phantom "574× runaway" in §A — a ratio taken across incommensurable
quantities — and it is worth naming twice because it has now cost two wrong verdicts.

**What actually moves (`test/scratch/diag_allsix_movers.jl`).** The prediction recorded here — that the
movers would be `houslack` and the `fhou2[d]` block, carrying the §L direction — is **refuted**.
`houslack` moved 9.2e-6 and `fhou2[1]` 1.9e-5, i.e. nothing. The top movers are:

| variable | benchmark | \|move\| | relative |
|---|---|---|---|
| `delVGDPEXP[2,7]`, `delGDPINC[1,2]`, `delXGDPEXP[2,7]` … | 0 | 33–39 | (additive) |
| `wlux[1]` | 903,646 | 33.40 | 3.7e-5 |
| `whouhtot[1]` | 1,644,636 | 31.51 | 1.9e-5 |
| `xuse[19,1,1]` | 617,391 | 16.20 | 2.6e-5 |
| `xinvi[2,2,1]`, `xinv_s[2,1]`, `xinv[2,1,1]` | ~25,000 | 15–17 | **6.2e-4** |

`delVGDPEXP` / `delGDPINC` / `delXGDPEXP` are declared at `TERM.TAB:1414`, `1461`, `1508` as
*"Ordinary change in … GDP component"* — additive **reporting** variables in money units, not
log-gaps. So the "a displacement of 39 means e³⁹" alarm is void: **no top mover is a log-gap**, and
575 free variables exceeding 1.0 in absolute terms just means 575 of them are large money numbers.

**But the amplification does not fully vanish under the honest metric.** A 1e-6 shock producing a
~2e-5 GDP response is ~20×, and ~600× in the investment block, whereas labour's ~0.4 cost share
predicts ~0.4e-6. Two candidates remain, and exactly one cheap test separates them:

  (i) a **genuine but stiff** derivative — plausible, since `xcap = faccum` and `finv1 = finv4` put
      the dynamic capital-accumulation block in charge of investment, the most elastic block in any
      CGE, and investment is precisely where the largest relative moves appear;
  (ii) **drift** along a weak direction, where the residual test is satisfied but the point is not
      a well-defined response.

A derivative is homogeneous: halve the shock and every response halves exactly. Drift has no reason
to. `test/scratch/diag_linearity.jl` solves the identical closure at δ and δ/2 from **fresh models** (since
`solve_newton!` mutates, re-shocking one model would measure a path, not two independent solves)
and prints the per-mover ratio. Ratios at 2.0 ⇒ (i), closure usable, continuation just needs small
steps. Scattered ratios ⇒ (ii), and the converged residual certifies nothing.

### N. The fold is real, it is a LIVESTOCK mode, and it is not any of the obvious defects

Four measurements, 2026-07-29. Together they answer §J step (2) and retire the `[2,5]` story.

**N.1 — The fold SURVIVES a faithful closure** (`EXCLUDE=houslack test/verify_swapped_closure.jl`).
Unlike ALL SIX (zero steps accepted), the five good swaps walk a real branch and then obstruct:

| s | δ | `pcap[2,5]` | Δ/δ |
|---|---|---|---|
| 0.99985 | 1.500e-4 | 1.0081684 | 54.5 |
| 0.9998125 | 1.875e-4 | 1.0122202 | 65.2 |
| 0.9998078 | 1.922e-4 | 1.0132692 | 68.3 |

`steps=3 rejects=7`, step size collapsing below `hmin` with an **accelerating** response — the
`A − B√(δ*−δ)` limit-point signature. So the answer to §J is: **the fold is NOT a closure
artifact.** It is also *earlier* than the base-closure fold (0.99981 vs 0.99975), so the `xcap =
faccum` swap that §J predicted would remove it does not.

**N.2 — The response is one coherent mode** (`test/scratch/diag_linearity.jl`). Doubling the shock:
`ratios: min=2.1909 median=2.1969 max=2.2095` over 20 movers — agreement to **0.9%**. That is not
drift (drift does not reproduce a ratio twenty times); it is a smooth response with strong
curvature, `bδ/a ≈ 0.10` at δ=1e-6. The script's own printed verdict said "NOT a clean derivative"
because it scored `|r − 2| < 0.05`; the discriminator is **uniformity**, not proximity to 2.0, and
it has been fixed. *(That was the fourth over-confident script verdict in this investigation — the
pattern is worth naming: a hard-coded threshold in a diagnostic is a hypothesis, and it gets tested
by the data like any other.)*

**N.3 — It is a SECTOR-2 mode, not a `[2,5]` cell** (`test/scratch/diag_pcap_rank.jl`). Ranking every
`pcap[i,d]` by relative response to a 1e-6 shock:

```
[2,3] 1187×  [2,1] 946×  [2,4] 313×  [2,5] 302×  [2,6] 212×  [2,2] 168×
[3,3]   86×  [3,6]  72×  ...  every other cell ≤ 86×
```

All six regions of sector 2, uniformly — which is exactly why N.2's ratios were uniform.
**`[2,5]` ranks 4th and `[2,3]` is four times larger.** `[2,5]` was never special; it is simply the
cell `verify_swapped_closure.jl` prints, inherited from older residual history. Reading "the fold is
at [2,5]" off a script that only ever looks at [2,5] was circular, and tasks #23/#24 were scoped on
that premise.

**N.4 — Four candidate causes, all refuted by measurement:**

| hypothesis | test | result |
|---|---|---|
| degenerate capital at `[2,5]` | `diag_cell25_scale.jl` | **No** — `VCAP[2,5]=1543.5`, rank 41/150, capital share 0.1189 vs 0.1117–0.1141 for other Livestock cells |
| dust capital cells drive it | `diag_pcap_rank.jl` | **No** — 0 of top 20 movers have `VCAP<100`; the real dust cells `VCAP[5,5]=0.040`, `VCAP[5,6]=0.033` never appear |
| near-Leontief `SIGMAPRIM[2]` | `diag_sector2.jl` | **Partial only** — `SIGMAPRIM[2]=0.239`, but sectors 3–7 are *lower* at 0.200 and respond 20× less |
| intermediate self-loop (feed) | `BSMR` diagonal | **No** — sector 2 self-share 0.038, near the bottom; sectors 18 (0.547) and 16 (0.494) dominate and are not stiff |

Also ruled out: multi-production (`MAKE` is exactly diagonal for all 25 sectors), and `SCET`,
`SGDD`, `SLAB`, `DPRC`, `QRAT`, `RADJ`, `ALFA`, `TFRO`, `XPEL`, all in-family for sector 2.

*Script defect found and fixed while re-checking this table:* `diag_sector2.jl`'s H2 branch probed
for a header named `"USE"`, which does not exist — the intermediate flow is `BSMR`
(25×2×29×6, COM×SRC×USR×DST, `build_reg1!.jl:348`). It found nothing and **skipped silently**, so
H2 had never actually been run. Now rewritten to read `BSMR` with the correct axis order — and to
restrict the denominator to the 25 industry columns, since `nu = na + 4` and the four final-user
columns are not industries. Re-run confirms the row above on real evidence: sector 2 self-share
0.0380, **rank 21 of 25 from the top** against a median of 0.2478. Livestock has one of the weakest
self-loops in the grid, not the strongest. H1 likewise measured properly: `SIGMAPRIM[2]=0.2391`
ranks 6 of 25 *from the bottom*, with sector 4 lower still at 0.200.

**N.5 — The one genuine sector-2 outlier: `TARG` (= `RNORMAL`, `TERM.TAB:2719`).**

```
TARG[2] = 0.030321   — the MINIMUM of all 25 sectors
           5.1× below the median (0.15506), 2.1× below the next lowest (sector 1, 0.064358)
```

`DPRC` (depreciation) is a uniform 0.05, so **sector 2 is the only sector of 25 whose normal rate of
return is below its depreciation rate** (net −0.0197; every other sector runs +0.014 to +0.156).
`build_premod!.jl:139` sets `CAPSTOK = CAP/RNORMAL`, giving Livestock a capital stock **46.6× its
rental** against 5–8× everywhere else, i.e. a gross return of 0.022–0.030 against 0.05 depreciation.
Capital in Livestock loses value faster than it earns.

One framing of this did **not** survive its own test: the steady-state check
`INVEST vs DPRC·CAPSTOK` gives `INV/need = 2.0880` for *every* sector (sector 1 alone differs at
1.8816), so investment is calibrated proportionally and does not single sector 2 out. The anomaly is
in the return, not in the investment ratio.

**Status and what it means.** The limit point is a genuine property of the calibrated model, not a
solver defect, a closure artifact, or a dust cell — so **pseudo-arclength is now unblocked** by the
standing rule (a fold has been shown to survive a faithful closure). But its *location* remains
economically implausible: a fold at a 0.02% labour shock, when published TERM models take ±5%
routinely, and an amplification of ~1000× on Livestock rentals. `TARG[2]` is the only defect-shaped
thing found, it is uniquely anomalous, and it sits in the block the active swaps govern — but the
causal chain from `RNORMAL[2] < DPRC[2]` to a 1187× `pcap` elasticity has **not** been demonstrated,
and sectors 4/7 with lower elasticities responding 20× less shows single-parameter reasoning has
already failed once here. Next: verify `TARG` against the source `INFILE` header and TERM.CMF, then
re-run `diag_pcap_rank.jl` with `TARG[2]` set to the grid median. If the Livestock mode collapses,
the fold is a data defect and arclength is the wrong build; if it survives, build the continuator.

**N.6 — `TARG` verified against the source. It is genuine INDOTERM data, faithfully carried.**

The first half of that next step is done, and it clears the translation:

- `TARG` is a **real header** in `national_data.csv`, not one of the hard-coded placeholders.
  `build_reg0!.jl:315` forwards it explicitly — the comment there records that these Step-6 dynamic
  parameters *were* previously dropped, so `build_premod!.jl:115` fell back to `_natvec("TARG", 0.10)`.
  The fallback is **not** what is in play now; the real values are carried through.
- The 34→6 / 185→25 aggregation weights it by `MAKE` output (`aggregate_model!.jl:154`,
  `_wavg_nd`), which is the correct treatment for a **rate** — summing it would have been the bug,
  and it is not what the code does.
- At the raw 185-sector level, **`Livestock` is itself the single lowest `TARG` of all 185 sectors**:

  ```
  Livestock    0.012162   ← minimum of 185
  Rubber       0.014870
  Cocoa        0.019972
  Coffee       0.022781
  Clove        0.024601
  ...          median over the 185 = 0.15773
  ```

So the aggregated 0.030321 is **not an aggregation artifact and not a translation defect** — it is
what the source database says, only diluted upward by the other Livestock-group members. Whether a
gross return of 0.0122 against a 0.05 depreciation rate is *economically* sensible is a question
about the INDOTERM database, not about this port. Note also that the four next-lowest sectors
(Rubber, Cocoa, Coffee, Clove) are likewise below `DPRC`; at the 25-sector grid they are averaged up
into their groups, which is why sector 2 emerges as the only aggregate with a negative net return.

**N.7 — the counterfactual** (`test/scratch/diag_targ_counterfactual.jl`): sets `RNORMAL[2]` to the grid
median and re-measures the ranking. `RNORMAL` reaches the model by two channels and a test moving
only one would be worthless — directly as a `GRETEXP` multiplier (`prepare_parameters.jl:843`), and
through `CAPSTOK = CAP/RNORMAL` (`build_premod!.jl:131–139`). **Trap avoided:**
`prepare_parameters.jl:814` *reads* `CAPSTOK` from `agg["STOC"]` and does not re-derive it, so
patching `TARG` alone would leave the stock stale and silently make the test partial. The script
therefore patches both source headers in `agg` and re-runs the real `prepare_parameters!`, so every
derived coefficient is recomputed by the production path; control and counterfactual each start from
an independent `deepcopy`; and a benchmark gate (≤1e-8) runs before either ranking is trusted, since
changing the capital calibration can invalidate the base year.

**RESULT — `TARG[2]` CAUSES the Livestock mode. The collapse is total.**

| | control | counterfactual |
|---|---|---|
| benchmark | converged, 1.42609e-09 | converged, **1.42609e-09** |
| shock 1e-6 | converged, 2 iters | converged, **1 iter** |
| worst sector-2 `pcap` | **1186.5×** | **2.2×** (×0.002) |
| worst non-sector-2 | 85.9× | 4.9× |
| sector-2 dominance | 13.8× | **0.4×** |

Sector 2 does not merely shrink — it leaves the top 12 entirely, which now reads `[3,6] 4.9`,
`[3,3] 4.8`, `[3,5] 4.4`, … The grid-wide worst amplification falls from 85.9 to 4.9, i.e. the whole
model becomes ordinary. The benchmark gate passed **bit-identically** in both cases, which is itself
the expected result and not a sign the patch failed to take: the calibration is self-consistent by
construction (`CAPSTOK = CAP/RNORMAL` preserves `GROSSRET = RNORMAL` at *any* level of `RNORMAL`),
so the base year holds regardless. That the patch took is proved by the 1-vs-2 iteration count and
the completely different ranking.

**One confound ruled out by the design.** The counterfactual scaled `CAPSTOK[2,:]` by `oldT/newT`,
which *shrinks* Livestock capital — and a small capital stock was the earlier "vertical supply
curve" story for a *large* response. It predicted the opposite sign of what happened, so it is not
what is driving the collapse.

**N.8 — A genuine aggregation defect, found while auditing N.7 — and it does NOT rescue the model.**

`build_premod!.jl:131–139` builds `CAPSTOK = CAP/RNORMAL` so that `GROSSRET = CAP/CAPSTOK` equals
`RNORMAL` exactly. Aggregation then treated the two sides differently: `STOC` is a value flow and is
**summed**, while `TARG` is a rate and was averaged by **MAKE output** (`aggregate_model!.jl:154`),
with `aggregate_regions!.jl` taking an **unweighted mean** over merged provinces. Neither weight is
the right one. `test/scratch/diag_rnormal_consistency.jl` measured the result:

```
GROSSRET/RNORMAL over 25 sectors:  min 0.7081, median 1.0056, max 1.1945
identity holds within 1% for only 10 of 25 sectors;  worst offender = sector 2
  sector 2:  CAP 22533.1  CAPSTOK 1.04951e6  GROSSRET 0.021470  RNORMAL 0.030321  ratio 0.7081
  sector 1:  ratio 0.8347      sector 13: ratio 1.1945      sector 25: ratio 0.9390
```

For a rate `r_j = flow_j / K_j`, the only aggregate consistent with summed components is
`Σflow_j / ΣK_j`, i.e. the arithmetic mean weighted by the **denominator** `K_j`. Fixed in both
stages: `DPRC`, `TARG`, `TFRO` now weight by `STOC` per (industry, region) cell, and `QRAT` — a
ratio of two investment/capital ratios — by trend investment `TFRO·STOC`. `ALFA`/`RADJ`/`REXP` are
behavioural parameters, not flow ratios, so they were left on output weighting rather than changed
on taste.

**Verified after a `refresh=true` rebuild: the identity is now EXACT.**

```
per-CELL  GROSSRET[i,d] == RNORMAL[i,d]
  150 cells: min 1.000000, median 1.000000, max 1.000000
  within 1e-8 of exact: 150 of 150     (was: 10 of 25 sectors within 1%)
```

*The diagnostic had the same bug it was auditing.* Its first version summed `CAP` and `CAPSTOK`
across regions while taking an **unweighted mean** of `RNORMAL` over that same axis — precisely the
mistake it was written to detect — which manufactured a discrepancy out of nothing and made the fix
look only partial (min 0.7081 → 0.9576, but sectors 9 and 12 apparently getting *worse*, at 1.2702
and 1.3370). The identity asserted by `build_premod!` is **per cell**, so the script now checks that
first and treats it as the verdict, with the by-sector table as a capital-weighted readout rather
than the test. Lesson, same family as the two earlier unit errors: a summary statistic that
aggregates a rate must use the same weights as the quantities it is being compared against.

**I predicted this fix would make the fold WORSE. It did the opposite, and that overturns N.7.**
Output weighting is biased *upward* whenever the rate varies within a group, because the low-rate
members carry disproportionately large stocks (`CAPSTOK` divides by the small rate). The consistent
Livestock rate is **0.021470**, *below* the carried 0.030321 — the same direction N.7 had just shown
made the mode explode. Re-running `test/scratch/diag_pcap_rank.jl` after the rebuild
(`scratchpad/rank_after_fix.log`):

| | before the fix | after the fix |
|---|---|---|
| benchmark | 1.4260876923799515e-09 | **1.4260876923799515e-09** (bit-identical) |
| shock 1e-6 | converged, 2 iters | converged, **1 iter** |
| worst cell | `[2,3]` **1186×** | `[2,4]` **11.27×** |
| sector-2 cells | 1187, 946, 313, 302, 212, 168 | 11.3, 7.4, 7.4, 5.1, and below |
| `[2,5]` | rank 4 of 150 | **rank 23 of 150**, 3.6× |
| grid spread | 1186× → 86× → … | 11.3× → 3.7× across the whole top 20 |

The model has become ordinary: nothing on the grid is more than ~3× the twentieth-place cell, and
sector 3 now interleaves with sector 2 instead of sector 2 standing alone at 14× the field.

**N.9 — Corrected verdict: the fold was OUR bug, and it is fixed.** N.7 and N.8 point in opposite
directions and only one of them is a clean experiment.

1. **N.8 moved exactly one thing.** `STOC` and `CAP` are untouched by the fix; the only quantity
   that changed is `TARG[2,d]`, and it went **down**. Under N.7's "the level is too low" story that
   is the worsening direction, yet the mode collapsed 105×. A single-variable experiment that
   contradicts the prediction beats the prediction.
2. **N.7 moved two things and cannot separate them.** It scaled `CAPSTOK[2,:]` by `oldT/newT`
   (shrinking Livestock capital ~5×) *and* raised the rate, deliberately holding the N.8
   inconsistency ratio fixed at 0.7081. Both a large increase (N.7) and a decrease (N.8) collapse
   the mode, so "the level" cannot be the operative variable in either. What N.8 changes and N.7
   does not is the **consistency** between the carried rate `RNORMAL` and the realised
   `GROSSRET = CAP/CAPSTOK`. That is the cause, and it was a defect in this port's aggregation.
3. **Therefore there is no data question to put to the user, and `TARG` must not be patched.**
   `TARG` is faithful source data (N.6); the sub-depreciation Livestock return is real and the model
   now handles it. The §N.9 "modelling decision" drafted before this run is **withdrawn** — it rested
   on N.7 alone.

Remaining check: this measures the fold *precursor* at the linearisation. Whether large shocks now
go through end-to-end is `test/verify_shock_fix.jl` (`blabnat = 0.9997` and `0.99`).

**N.10 — The direct-shock floor is a SEPARATE phenomenon, and the fix does not touch it.**
`scratchpad/shockgate_after_fix.log`:

```
benchmark          converged   1.4260876923799515e-9
blabnat = 0.9997   maxit 60    1.642115382120283e-5   (416s)   ❌
blabnat = 0.99     maxit 60    0.004737814103183395   (472s)   ❌
```

The 0.9997 failure is expected — §I established that the direct path stalls while the *ladder*
converges — and the 0.99 rung is not a clean reading at all, since `verify_shock_fix.jl` reuses one
model and so started 0.99 from the stalled 0.9997 iterate. What matters is the floor itself:

| | ‖F‖∞ at maxit, `blabnat = 0.9997` |
|---|---|
| §G3, **before** the aggregation fix | 1.641919082126898e-5 |
| now, **after** the fix | 1.642115382120283e-5 |

**Unchanged to four significant figures**, while the quantity the fix targeted moved 105×. Two
distinct phenomena have been conflated throughout §F–§N: the Livestock fold precursor (now fixed)
and the 1.7e-5 direct-shock floor (untouched). §H's localisation of 97% of that floor to `[2,5]`
was measured on the *old* database and must be re-measured before it is cited again.

The stall mechanism is also now legible, and it is a solver ratchet rather than a model property.
From it ≈15 the trust radius collapses monotonically 42 → 0.001 while the true LU step grows to
|d| ≈ 70,000, so every iteration the exact step is demoted to a damped CGNR step that barely moves;
merit then falls 0.08% per iteration while ‖F‖∞ *rises*. The §G3 fix made Δ grow on `red > 0.1`, but
the `red < 0.01` branch still contracts to `max(moved, 1e-3)` — and `moved` is tiny precisely
*because* the step was damped, which is self-reinforcing. Note this is the same signature §G3
already called out ("stop tuning the solver"): a merit minimum that is not a zero of F.

**N.11 — The N.8 fix does not touch the BASE closure at all, and that is the point.**
`scratchpad/ladder_after_fix.log`, re-running `test/verify_continuation.jl` on the fixed database
against the fine ladder already recorded in §I:

| rung | §I, before the fix | now, after the fix |
|---|---|---|
| 0.99985 | converged, 4 it, 2.68e-9, `pcap[2,5]` 0.97186606 | converged, 4 it, 2.67755e-9, **0.97186606** |
| 0.9998  | converged, 3 it, 2.21e-9, 0.95785452 | converged, 3 it, 2.21189e-9, **0.95785452** |
| 0.99975 | maxit 40, 7.66e-7, 0.92623464 | maxit 40, 7.67e-7, **0.92623464** |

**Bit-identical on every rung** — iteration counts, residuals, and `pcap[2,5]` to eight digits. So
the base-closure fold at `blabnat* = 0.999752` is entirely unaffected, and the 1.7e-5 direct floor
(N.10) with it.

*Correction to avoid repeating:* on first reading I compared this run against §I's **coarse** ladder
and reported "two extra rungs and an 18× better floor". The fine ladder had already been run and
recorded at §I; there is no improvement. Compare against the finest baseline that exists, not the
first table in the section.

**Why it is inert here, and why that is consistent rather than damning.** In the base static
closure `xcap` is exogenous, so the capital-accumulation block never enters the residual and
`CAPSTOK`/`RNORMAL`/`GRETEXP` are dead coefficients. `TARG[2]` moving 0.030321 → 0.021470 therefore
*cannot* change this solve — the bit-identical result is the expected one, and is a check on the fix
rather than a strike against it. The `xcap → endogenous, faccum → exogenous` swap is precisely what
makes those coefficients live, which is why `diag_pcap_rank.jl` (which calls
`apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))` at line 54) sees the 105× collapse and this ladder
sees nothing.

*This also resolves an apparent 43× contradiction:* the ladder's 1e-6 rung moves `pcap[2,5]` by
1.5e-4 (155×) while `diag_pcap_rank.jl` reports 3.6× at the same shock size. Not a bug in either —
**two different closures**. Any amplification number quoted from here on must name its closure.

**So there are two distinct folds and the N.8 result speaks only to the second:**

| | base static closure | TERM.CMF swapped closure |
|---|---|---|
| capital block | inert (`xcap` exogenous) | **live** (`xcap` endogenous) |
| fold at | `blabnat* = 0.999752` (§I) | 0.99981 (§N, *earlier*) |
| effect of the N.8 fix | none, bit-identical | 1186× → 11.27× at δ=1e-6 |

**N.12 — ✅ THE FOLD IS GONE IN THE LIVE-CAPITAL CLOSURE.** `scratchpad/ladder_swapped_after_fix.log`,
`SWAPS=1`, ten rungs:

```
0.999999  converged  1 it  2.68e-9      0.99981   converged  2 it  4.54e-9
0.99999   converged  1 it  7.33e-9      0.9998    converged  1 it  6.41e-9
0.9999    converged  2 it  4.19e-9      0.99975   converged  2 it  6.05e-9
0.99985   converged  2 it  9.78e-9      0.9997    converged  2 it  4.77e-9
                                        0.9995    converged  2 it  6.05e-9
                                        0.999     converged  2 it  6.75e-9
```

It walks **straight through 0.99981** (the swapped closure's own previously-recorded fold) and
through **0.999752** (the base closure's fold), out to 0.999 — 1000× the shock size that used to
break — every rung in **1–2 Newton iterations** to ~1e-9.

The branch is now linear, which is the strong form of the claim:

| δ | `pcap[2,5]` deviation | amplification |
|---|---|---|
| 1e-6 | 3.61e-6 | 3.610 |
| 1e-5 | 3.612e-5 | 3.612 |
| 1e-4 | 3.6123e-4 | 3.612 |
| 3e-4 | 1.08394e-3 | 3.613 |
| 1e-3 | 3.6157e-3 | 3.616 |

Constant to three significant figures across three orders of magnitude of shock size. Against §I's
base-closure trajectory accelerating 155× → 156× → 173× → 246× into its limit point, this is not a
weakened fold — there is no fold on this branch. **Pseudo-arclength continuation is off the critical
path** for the swapped closure. (The base-closure fold at 0.999752 is untouched and still real, but
the base closure is not the model TERM.CMF specifies.)

This also closes the §N.1 gate from the other direction: the six TERM.CMF swaps were the *faithful*
closure, and with the N.8 aggregation bug fixed they are also the *better-behaved* one — 1–2
iterations where the base closure needs 3–4 and then folds.

**Next:** ladder extended to the C6b target (`scratchpad/ladder_swapped_C6b.log`,
0.999 → 0.995 → 0.99 → 0.98 → 0.97), i.e. the real `blabnat = −3` scenario.

### N.13 Gate 7 — `run_model!` scenario driver (written)

`continuation_solve!` walks **one** `VariableRef`. That is the right primitive for the fold
experiments above, but a scenario generally moves several exogenous variables at once, and walking
them one after another is a *different experiment*: it drags the model through intermediate states
nobody asked about, and those states can be harder than the target.

`src/run_model!.jl` therefore walks every shock **together** along a scalar homotopy `t ∈ [0,1]`,
each shocked variable held at `base + t·(target − base)`, with the same adaptive step + rollback as
`continuation_solve!`. Three design decisions worth recording:

1. **The benchmark gate is unskippable.** The driver applies the closure, then re-solves the
   benchmark, and *raises* if it does not reach `tol`. §N.1–N.3 cost several sessions precisely
   because a closure that cannot reproduce its base year still produces numbers, and those numbers
   measure the closure bug rather than the scenario.
2. **The default closure is `TERM_CMF_SWAPS`, not the bare automatic one.** §N.12 settled that it
   is both the faithful closure and the better-behaved one. Passing `swaps = ()` is allowed but
   logs a warning that capital accumulation is inert in that closure — the trap that made the same
   cell read 155× and 3.6× depending on which closure was running.
3. **A shock on an endogenous variable raises.** Under the swapped closure `xcap`, `xhouhtot`,
   `fgovtot`, `finv1`, `delfwage` are all *endogenous*; shocking one silently would leave the
   system non-square.

`pct(-3)` converts a GEMPACK percentage shock on a bmk=1 variable to its levels target (0.97), so a
scenario file reads in the units the source states. `pct_shocks` handles variables whose benchmark
level is not 1 by scaling each variable's own base.

`test/run_scenario_6reg.jl` is the gate. It runs TERM.CMF's own scenario (`Shock blabnat = -3`)
under TERM.CMF's own closure and asserts four things in order: benchmark reproduced → homotopy
reaches `t = 1` (the FULL shock, not a fraction) → `wgdpdiff` still zero *under the shock* →
real GDP moved the right way and by a plausible multiple. That last check is the only one that is
not residual-based, and it is the only one that can catch a sign error: **every residual gate is
satisfied by any solution of the system, including one that moved the wrong way.** The levels
income-vs-expenditure gap is reported but not asserted — the benchmark database itself is off by up
to 36% in the worst region, so asserting it would be asserting the input data is perfect.

### N.14 The swapped ladder out to the C6b target — 0.995 reached, 0.99 is a *step-size* failure

`scratchpad/ladder_swapped_C6b.log` (`SWAPS=1`, `LADDER=0.999,0.995,0.99,0.98,0.97`), warm-starting
from the §N.12 point:

| `blabnat` | status | iters | ‖F‖∞ | `pcap[2,5]` | deviation | amplification |
|---|---|---|---|---|---|---|
| 0.999 | converged | 3 | 2.21e-9 | 0.9963843 | 3.6157e-3 | **3.616** |
| 0.995 | converged | 4 | 9.55e-9 | 0.98189557 | 1.810443e-2 | **3.621** |
| 0.99  | maxit 40 | — | 6.47e-5 | (0.96465129, not a solution) | — | — |

**−0.5% of national labour supply now solves.** That is 20× past §N.12's 0.999 and 500× past where
the base closure folds.

**The 0.99 rung is not evidence of a fold, and the script's own verdict line ("the [2,5] block needs
inspection") is wrong here.** `verify_continuation.jl` walks a *fixed* ladder and its pass/fail text
was written for a fixed ladder; it cannot distinguish a fold from a step that was simply too long.
Two facts say this one is too long:

1. ~~**Amplification is still linear at the last good point.** 3.616 at δ=1e-3 and 3.621 at δ=5e-3 —
   the same constant that held from δ=1e-6 in §N.12, now across *four* orders of magnitude. Near a
   fold the response goes as `A − B·√(δ*−δ)`; it does not stay linear right up to the wall.~~
   **OVERSTATED — see §N.15.** The sentence is true of a point *close* to a fold, and δ=5e-3 was
   not close: the obstruction sits at δ≈9.19e-3, so the last measurement was barely past the
   halfway mark, where `√(δ*−δ)` curvature is not yet visible. Linearity at 0.54·δ* was never
   evidence against a fold at δ*. The *conclusion* of §N.14 — "a fixed ladder cannot tell a fold
   from an overlong step, so re-measure adaptively" — stands and was the right call; this
   particular supporting argument does not.
2. **The step that failed was the largest yet, from the hardest point yet.** 0.995 → 0.99 is δ=5e-3,
   larger than the whole 0.999 → 0.995 leg, launched from a point that had itself needed 4
   iterations. A fixed ladder has no way to back off; it just fails.

This is exactly the failure mode adaptive continuation exists for, so the next measurement is
`run_model!`'s homotopy (§N.13) rather than another hand-written rung list. **Do not quote "the
swapped closure folds at 0.99" from this log** — the experiment cannot support that claim. If the
adaptive walk also collapses near 0.99 with `h` driven to `hmin`, *that* would be the fold evidence,
and it is a different measurement.

### N.15 The adaptive walk stops at `t = 0.30625` — this one IS a wall, not a long step

`scratchpad/gate7_blabnat3.log` — `run_model!`, TERM.CMF closure, `blabnat` 1.0 → 0.97, `h0 = 0.05`:

```
✅ t=0.05     3 it   3.49e-9        ✗ 0.3125    maxit  1.286e-5   h → 3.12e-3
✅ t=0.15     3 it   6.75e-9        ✗ 0.309375  maxit  6.365e-6   h → 1.56e-3
✅ t=0.25     4 it   1.63e-9        ✗ 0.307812  maxit  2.883e-6   h → 7.81e-4
✅ t=0.30     5 it   5.47e-9        ✗ 0.307031  maxit  1.148e-6   h → 3.91e-4
✅ t=0.30625  5 it   3.26e-9        ✗ 0.306641  maxit  5.656e-7   h → 1.95e-4
                                    ✗ 0.306445  maxit  1.269e-7   h → 9.77e-5  ⛔
```

Five clean steps, 5 accepted / 12 rejected, then **six successive halvings against the same point**.
Reached `t = 0.30625`, i.e. `blabnat = 0.9908125` — a **−0.92%** national labour-supply shock, ~30%
of the way to TERM.CMF's −3%. Unlike §N.14 this is not an overlong step: `h` was driven to `hmin` at
one parameter value while every earlier step converged, which is what `continuation.jl`'s own header
defines as the fold signature.

**Two readings remain, and the residual pattern cannot separate them.** The failed-step residual is
≈3e-3·h across all six halvings — Newton reached the same *relative* accuracy every time and
stopped. That is what a fold looks like, and it is equally what the §N.10 trust-region ratchet
(`solve_newton!.jl:636-638`, diagnosed, never fixed) looks like: it plateaus at a fixed fraction of
the step because `moved` is small only *because* the step was damped.

**The `pcap[2,5]` probe DECELERATES into the wall** — local slope `d(pcap)/dt` runs
0.1087 → 0.1072 → 0.0971 → **0.0410** over the last four points, and amplification falls
3.622 → 3.602 → 3.541 → 3.497. Folds accelerate. That probe is a leftover from the refuted "cell
[2,5]" framing (see [[indoterm-shock-residual-cell25]]) and is demonstrably **not** the singular
direction, so it says nothing about whether a fold exists — only that it is not there. This is the
fourth time a single-cell probe has misled; the rule stands: **rank all cells, never read a
localisation off a script that prints one.**

**The decisive measurement** is `test/scratch/diag_fold_tangent.jl`: at each converged point solve
`J·v = −∂F/∂λ` by **direct sparse LU** (never CGNR — §N.10) for the tangent `dx/dλ`, where `∂F/∂λ`
is the full-Jacobian column of the fixed shock variable.

* A fold makes `J` singular *at* the limit point, so `‖v‖∞ → ∞` like `1/(λ*−λ)`. No solver behaviour
  can suppress that — it is a property of the equations. `1/‖v‖∞` then falls ~linearly and
  extrapolates to the fold location, and the **argmax component names the direction**, which would be
  the first honest localisation of this obstruction.
* Bounded, flat `‖v‖∞` ⇒ not a fold; the branch continues and the defect is in step acceptance. Fix
  the ratchet first, and do not touch the model.

### §N.16 — ✅ DECIDED: it is a genuine FOLD, and the direction is region 4's household block

`test/scratch/diag_fold_tangent.jl` ran the §N.15 measurement. Result (`scratchpad/fold_tangent.log`),
TERM.CMF swapped closure, `blabnat` walking 1.0 → 0.97:

| t | blabnat | ‖dx/dλ‖∞ | 1/‖dx/dλ‖∞ | σ_min(J_scaled) | argmax |
|---|---|---|---|---|---|
| 0.05    | 0.9985    | 1.95338e6 | 5.11933e-7 | 3.854e-9  | `whouhtot[2]` |
| 0.15    | 0.9955    | 1.76704e6 | 5.65919e-7 | 3.195e-9  | `whouhtot[4]` |
| 0.25    | 0.9925    | 2.46976e6 | 4.04898e-7 | 2.089e-9  | `whouhtot[4]` |
| 0.30    | 0.9910    | 5.84110e6 | 1.71201e-7 | 7.61e-10  | `whouhtot[4]` |
| 0.30625 | 0.9908125 | 3.50683e7 | 2.85158e-8 | 1.159e-10 | `whouhtot[4]` |

**Verdict: FOLD.** Overall growth 17.95×, but **14× over the last stretch alone**, while `σ_min`
collapses 3.85e-9 → 1.159e-10. That is `‖v‖ ~ 1/(λ*−λ)`, not a plateau — and it is a property of
the equations, immune to anything `solve_newton!` does. The §N.10 trust-region ratchet is therefore
*not* what stops the walk at `t = 0.30625`; it may still be real, but it is not this.

**Location.** Linear extrapolation of `1/‖v‖∞ → 0` from the last two points (slope −2.283e-5) gives
**t\* ≈ 0.3075, blabnat\* ≈ 0.990775** (i.e. −0.92% labour supply). The `σ_min` column extrapolates
to t\* ≈ 0.3175. Trust the `1/‖v‖∞` figure: `σ_min` is an 8-step inverse power estimate on the
row-scaled matrix. So **t\* ≈ 0.307–0.318**, essentially where the walk actually stopped.

**Direction.** Top-5 tangent components at t=0.30625 — all region 4, all household/expenditure:

```
whouhtot[4]      -3.50683e7      delVGDPEXP[4,1]  -3.45144e7      wlux[4]  -3.2481e7
delGDPINC[4,2]   -2.99267e7      delXGDPEXP[4,1]  -2.39542e7
```

`delVGDPEXP`/`delXGDPEXP`/`delGDPINC` are *reporting* aggregates that read off `whouhtot`, so the
only two structural members of that list are `whouhtot[4]` (region 4 nominal household expenditure)
and `wlux[4]` (region 4 **supernumerary** expenditure). This is the first honest localisation of
this obstruction — and note it is **not** sector 2 / region 5, which §N.15 had already shown
decelerates into the wall.

**Hypothesis (UNVERIFIED — labelled as such until measured).** This is the Klein-Rubin/ELES
subsistence floor. `E_xlux!` is `xlux[c,d]·phou[c,d] == wlux[d]·alux[c,d]` and `E_xhouh_s_agg!` is
`xhou_s[c,d] == xlux[c,d] + xsub[c,d]`, with `xsub` pinned to `XSUB0·nhou·asub` — i.e. subsistence
demand is *inelastic*. If region 4's income falls far enough that `wlux[4] → 0`, all adjustment
margin is gone and the block goes singular. That would be an economically meaningful limit point,
not a translation bug. `test/scratch/diag_wlux_floor.jl` measures it.
**❌ REFUTED — see §N.17.**

**Remedy either way: pseudo-arclength continuation.** Natural-parameter continuation provably cannot
pass a fold; arclength parameterisation stays nonsingular through the turn. §N.12 had deprioritised
it; it is now back on the critical path.

### §N.17 — ❌ the subsistence-floor hypothesis is REFUTED; `whouhtot[4]` is running UP, not down

`test/scratch/diag_wlux_floor.jl` (`scratchpad/wlux_floor.log`). The §N.16 hypothesis was that region 4's
supernumerary expenditure `wlux[4]` runs to zero and exhausts the ELES adjustment margin. It does
the opposite.

**Premise fails.** Benchmark supernumerary share `WLUX0/HOUPUR_C` is **0.5495 in every region** —
`BLUX` is uniform, so no region is structurally thin and there is no "region 4 is different"
benchmark story to tell. (Worth knowing independently: the ELES budget split carries no regional
variation at all in this database.)

**Conclusion fails, with the sign reversed.** `wlux[4]/WLUX0[4]` along the branch:

| t | 0 | 0.05 | 0.15 | 0.25 | 0.30 | 0.30625 |
|---|---|---|---|---|---|---|
| `wlux[4]` | 1.000000 | 1.007270 | 1.023246 | 1.043151 | 1.059426 | **1.065307** |

It **rises**, and it **accelerates**: local slope `d(wlux[4]/WLUX0)/dt` goes 0.326 over
`t = 0.25→0.30` and then **0.941** over `0.30→0.30625` — nearly 3× in the final step, in exactly
the stretch where `‖dx/dλ‖∞` grew 14×. No other region does this (1,2,3,5 drift within ±1.1% and are
falling; 6 flips sign but stays inside 1%). The minimum `xlux/XLUX0` over all 150 cells never goes
below **0.9749**, and every `xsub` sits exactly at its benchmark, so nothing anywhere is near a
subsistence floor. The fold is a **runaway, not a floor**.

**What this leaves.** Under the TERM.CMF swaps `fhou` is exogenous, so `E_fhou!`
(`build_equations.jl:1777`) is what *determines* `whouhtot`:

```
log(whouhtot[d]/HOUPUR_C[d]) == log(wlab_io[d]) + fhou[d] + houslack
```

That equation is log-linear and cannot be singular on its own. The singularity must therefore be in
the **loop it closes**: `whouhtot[4] → regional household demand → regional output → labour demand →
wlab_io[4] → whouhtot[4]`. A fold is precisely where such a multiplier loop's gain reaches 1, and an
accelerating rise in the loop's own state variable is the expected signature. This is a **third
hypothesis and is NOT measured** — do not repeat §N.16's mistake of writing it up as if it were.

**Decision: stop mechanism-hunting here and build the remedy.** Three candidate mechanisms have now
been proposed and the measurement has killed two of them, but the remedy has been unchanged
throughout: **pseudo-arclength continuation**, which parameterises by arclength and stays
nonsingular through the turn. Knowing *why* the wall is there changes nothing about how to pass it —
it only tells us afterwards whether −3% labour supply is economically reachable in this closure at
all. Measure the loop gain *after* the branch can be continued, when the far side is available to
compare against.

### §N.18 — WHERE TO RESUME (stopped 2026-07-29 by request, mid-run)

**SUPERSEDED BY §N.19** — the run described here was attempted and failed for two reasons unrelated
to the fold. The env defaults quoted below are stale (`DS0` is now 0.005, `MAXSTEPS` 60).

`src/arclength.jl` (pseudo-arclength continuation) and `test/verify_arclength.jl` are **written,
load-tested, and NOT yet validated** — the first run was stopped before it produced output, so
nothing is known about whether the augmented assembly is correct. Resume with:

```
julia --project=. test/verify_arclength.jl
```

Env: `TARGET` (default 0.97), `DS0` (0.05), `MAXSTEPS` (120). ~30–45 s per corrector iteration at
25×6, so budget an hour. Julia block-buffers stdout when redirected — an empty log means running,
not stalled.

Three distinguishable outcomes, all handled by the script: `reached_target` (fold was
parameterisation-only, −3% becomes runnable); `turned` without reaching (the branch genuinely folds,
so `blabnat = 0.97` does not exist in this closure and `lam_fold` is the largest admissible shock);
neither (the only outcome that indicates a defect, since a fold cannot defeat arclength).

Nothing else is running. Open items, ranked: this validation → 34-region condensation (Excerpt 49,
the LU fill-in wall) → task #11 price homogeneity (written, never run) → §N.10 ratchet (diagnosed,
unfixed, exonerated for this fold but still misleading) → task #36 (fold mechanism, interpretive).

### §N.19 — the first arclength run failed, for TWO reasons, both mine and both now fixed

Run 1 of `verify_arclength.jl` (`scratchpad/arclength.log`) died after four accepted steps with
`OpenBLAS: malloc failed in gemm_driver`. **Note the shell still reported exit 0** — the "task
completed" notification was meaningless, and no assertion in the gate ever executed. Neither of the
three outcomes in §N.18 applies; the run never got far enough to test anything.

**Defect 1 — the metric. This is the `del*`-at-zero units trap for the SEVENTH time, reproduced
inside the fix for it.** Measured from the log:

| quantity | value |
|---|---|
| total `Δs` over 4 steps | 0.75 |
| total `Δλ` | −1.155e-7 |
| `τ_λ` (observed) | 1.54e-7 |
| implied `‖dz/dλ‖₂` (scaled) | 6.49e6 |
| largest *relative* response (`whouhtot[4]`) | **2.69** (−1.39e6 on a base of 515,617) |
| `Δs` required to reach λ=0.97 | **194,805** ⇒ ~2 months of wall-clock |

`τ_λ = 1/‖v‖₂`, and the 6.49e6 is **not economics**. `cscale = max(|x_bench|,1)` is deliberately
*absolute* for the additive `del*` variables that legitimately sit at zero — correct for step control
and for avoiding the phantom-runaway trap of §F. But those variables carry Rp-billion derivatives, so
on a unit-scaled column each contributes ~1e6 to the arclength metric. §N.16's own top-5 tangent
components name the culprits: after `whouhtot[4]` and `wlux[4]` come `delVGDPEXP`, `delXGDPEXP`,
`delGDPINC`. **Three reporting aggregates were consuming the entire arclength budget and leaving λ
1.5e-7 of it.** The fold was never the binding constraint on run 1.

Fix: arclength may be measured in *any* fixed positive-definite metric — the choice affects
efficiency, not validity. The metric is now built from the branch's own sensitivity at the start
point, `wt = 1/max(|vᵢ|, 1e-6‖v‖∞)` with `v = dz/dλ` (the bootstrap tangent, free — it was already
being computed and thrown away), normalised in the **∞-norm** so the measure carries no `√N` penalty.
`wt` is kept strictly separate from `cscale`; conflating a *metric* with a *state scaling* is the
whole defect.

**Defect 2 — the dense closing row.** With `τᵀ` as the last row the augmented matrix has one fully
dense row of 83,516 entries, and the factorisation exhausted memory in UMFPACK's dense frontal
kernels on its fifth call. This file's own notes had predicted exactly that and prescribed "scale the
tangent row down", which would not have helped — scaling does not make a dense row sparse.

Fix: close the system by **Rheinboldt local parameterisation** instead — advance the single
coordinate `k = argmax|τ|` and leave λ free, so the border row is `e_kᵀ`, one nonzero. Nonsingular at
the fold on the same grounds as pseudo-arclength (at a fold the tangent *is* the null direction, so
its argmax has `τ_k ≠ 0` by construction), and the factorisation stays genuinely sparse.

**A trap inside the fix, caught before running:** `k` must be the argmax of the **unweighted** tangent.
`wt` normalises every significant component to exactly 1.0, so an argmax over the *weighted* tangent
is a mass tie broken by index order — it would have selected the first variable above the floor,
which means nothing. The two choices control different things: `wt` sets the step size, `argmax|τ|`
maximises the nonsingularity margin. The largest mover among genuine value flows (`cscale > 1`) is
reported separately as `driver=`, since that — not the local coordinate — is the economically
meaningful direction, and it is how §N.16 identified `whouhtot[4]`.

Defaults rescaled to the shock, not to 1: `ds0=0.005`, `dsmax=0.02`, `dsmin=1e-7`, `MAXSTEPS=60`.
The routine now prints `dλ/ds` before stepping, as a tripwire: **if it reads ~1e-7 again, stop the
run — the metric is still being eaten and no amount of patience will help.**

**Both fixes verified by run 2** (`scratchpad/arclength2.log`), independently of anything downstream:

* the code computed "an unweighted metric would give dλ/ds≈1.54e-7" from `dz/dλ` at the start point,
  matching the 1.54e-7 derived from run 1's step log by hand — two independent routes to the same
  number, so the metric *was* the cause and the benchmark point is *not* near-singular;
* weighted `dλ/ds = −1.0`, a gain of 6.5e6;
* LU **1.9 s** against run 1's **28.4 s** — the single-nonzero border row, worth ~15× on its own.

### §N.19b — defect 3: a zero-length step was being accepted as corrector progress

Run 2 then rejected every step with the residual **bit-identical across all 20 corrector iterations**
(`‖F‖∞=0.001632`, `n=0.0`, ~2 s of LU each). Two clues fixed the cause without further measurement:
`n = 0.0` exactly (the closing row was satisfied from the predictor onward), and the rejected-step
residuals scaling as `0.001632 / 0.0004 = 4.08 ≈ 2²` between `ds=0.005` and `ds=0.0025` — i.e. pure
`O(ds²)` Euler-predictor truncation error. **The corrector was never moving the point at all.**

Cause: the fraction-to-boundary rule. `arclength.jl` computed `α = min(α, 0.95·gap/d)` with no special
case for `gap = 0`, so a variable sitting exactly ON its lower bound forced `α = 0` — and with 83,515
free variables at least one always does. Then the Armijo test `mt ≤ (1 − 1e-4·α)·base` degenerates at
α=0 to `base ≤ base`, which **passes**: a zero-length step was recorded as a successful corrector
iteration. One pinned variable vetoed the step for all the others, twenty times per attempt.

`solve_newton!.jl:566-572` — the solver that works — already had the missing branch: when
`gap <= 1e-8` it zeroes *that component* of the direction rather than letting it zero the global step.
`arclength.jl` now mirrors it exactly. Plus two guards so a no-op can never again masquerade as
progress: the line search requires `α ≥ αmin = 1e-10`, and an `α` that collapses below 1e-6 is
reported *with the binding variable's name and gap* and rejected rather than backtracked from.

**Lesson worth generalising: this is the third defect in a row where the new code diverged from
`solve_newton!` on a detail that file had already got right and documented** (merit function on
`‖F‖₂²` not `‖F‖∞`; frozen row scaling; at-bound component zeroing). When writing a second solver
against this model, diff against `solve_newton!` deliberately rather than reasoning each choice out
from scratch — its comments record measurements, not preferences.

### §N.19c — VERDICT: the fold is real, measured, and confirmed independently

> **🟢 SUPERSEDED BY §N.19e — read that first.** The fold measured below is real *for the scenario
> it was measured in*, and every number here stands. But that scenario was **incomplete**: it omitted
> `TERM.CMF`'s second shock, `delUnity = 1`. With the full scenario the branch passes straight
> through λ = 0.9908083 without turning and reaches `blabnat = 0.97`. **The fold is an artifact of
> the incomplete scenario, not a property of the model.** Do not cite λ* = 0.9908083 as a limit on
> productivity improvement.

> **⚠️ TERMINOLOGY CORRECTION — read before the rest of this section.** Everything below was first
> written calling `blabnat` a **labour-supply** shock. **That is wrong.** `blabnat` drives
> labour-augmenting **technical change** (`TERM.TAB:460` `alab_o # Labor-augmenting technical change #`;
> `:468` `blabnat # Driver for alab_o #`; `:496` equation `E_alab_o # Total labour productivity growth #`),
> and `TERM.CMF:112` states it outright: `Shock blabnat = -3; ! 3% increase labour productivity`.
> So λ → 0.97 is a **3% productivity IMPROVEMENT** — expansionary, not contractionary. The numerics
> are unaffected; the economic reading inverts. Read every "−0.92% labour supply" below as
> **"0.92% labour-productivity improvement"**. See §N.19d.

Run 3 (`scratchpad/arclength3.log`, TERM.CMF swapped closure, `ds0=0.005`, `dsmax=0.02`,
`MAXSTEPS=60`, tol 1e-8) traversed the limit point cleanly. **60 accepted steps, 0 rejections,
4.4 minutes.** Every accepted point satisfies `‖F‖∞ ≤ 1e-8` — the gate asserts this over the whole
path, so "passed the fold" cannot be an artifact of a loosened tolerance. On the descending arm and
through the turn the corrector converged in 3 iterations to `‖F‖∞ ≈ 2e-12`, four orders inside the gate.

**`λ* = 0.990831647`. Minimum λ on the branch = `0.9908083427`, i.e. a −0.9192% labour-supply shock.**

`dλ/ds` decays monotonically into the turn, changes sign, and re-grows on the far side — the textbook
saddle-node signature:

| descending | λ | dλ/ds | | ascending | λ | dλ/ds |
|---|---|---|---|---|---|---|
| | 0.9939150 | −0.0698 | | | 0.9908986 | +0.00219 |
| | 0.9928711 | −0.0556 | | | 0.9911147 | +0.00568 |
| | 0.9916764 | −0.0267 | | | 0.9916348 | +0.00955 |
| | 0.9910809 | −0.0137 | | | 0.9926186 | +0.01347 |
| | 0.9908359 | −0.00569 | | | 0.9961066 | +0.02001 |
| min λ | **0.9908083** | −0.00264 | | end | 1.0102314 | +0.02636 |

**Independent confirmation of §N.16.** §N.16 extrapolated `blabnat* ≈ 0.990775` from the collapse of
`1/‖dx/dλ‖∞`, and the natural-parameter walk (§N.15) stalled at `t = 0.30625 → 0.9908125`. Arclength
*traverses* the point and measures `0.9908083`. Three routes, two of them methodologically unrelated,
agreeing to **3.4e-5**. The fold is not an artifact of any one method.

**CONSEQUENCE — `blabnat = 0.97` does not exist in this closure.** Past λ* there is no solution at any
smaller `blabnat` on this branch, at any step size, with any solver. The −3% labour-supply scenario is
not a solve-harder problem. The largest labour-supply contraction the model admits in the TERM.CMF
swapped closure is **−0.92%**. This closes task #34 in the negative and definitively exonerates the
§N.10 trust-region ratchet for this obstruction.

**The far arm is a distinct branch, not a retrace.** Checked, because a sign-flip bug in the tangent
update would look identical in λ alone: at λ ≈ 0.99287 the descending arm has `dλ/ds = −0.0556`, the
ascending arm `+0.0143` — a factor of 3.9 apart. The local coordinate also flips `delVGDPEXP[4,1]` →
`delVGDPEXP[2,1]` across the turn. So **two solutions exist for every `blabnat` in (0.99081, 1)**, and
the ascending arm crosses λ = 1 (at 1.000100001) at a point that is *not* the benchmark — the benchmark
equilibrium is not unique in this closure. Worth knowing before any policy run: a scenario near this
region has a second equilibrium that continuation from the benchmark will never see.

Corroborating detail: the local coordinate approaching the fold is `delVGDPEXP[4,1]` — **region 4**,
the same region §N.16 localised via `whouhtot[4]`. `driver` (largest *relative* mover among genuine
value flows) stays in the region-4/sector-2 investment block throughout. This is consistent with the
task-#36 hypothesis that a region-4 income multiplier loop is what turns the branch.

**Measurement trap #8 — a metric normalised at a point cannot measure curvature at that point.**
Mid-run I extrapolated `dλ/ds ∝ √(λ*−λ)` from the first two steps and got λ* ≈ 0.99391, then flagged
it as contradicting §N.16. It was wrong, and the branch refuted it within minutes by sailing past
0.99391. Cause: the response-weighted metric is defined as `wt = 1/|dz/dλ|` **at the start point**, so
`dλ/ds = −1.0` there *by construction*. The 14× drop on step 1 is calibration, not curvature. Only
from step 2 onward is `dλ/ds` a curvature signal — and from there the decay is clean and the
extrapolation would have been sound. Do not feed a definitionally-fixed number into a curvature law.

### §N.19d — we have been running HALF the authors' scenario (`delUnity` is missing)

Found by finally reading the shipped GEMPACK control files instead of only `TERM.TAB`. `TERM.CMF`
defines the reference simulation with **two** shocks:

```
Shock blabnat =   -3; ! 3% increase labour productivity
Shock delunity=1;     ! move forward 1 year
```

**We have only ever applied the first.** `delUnity` is not decoration — in `TERM.TAB` it gates three
adjustment mechanisms simultaneously:

| TERM.TAB | mechanism |
|---|---|
| 2700 | capital accumulation — `0.01*CAPSTOK_OLDP*xcap = CAPADD*delUnity + faccum` |
| 2803 | rate-of-return adjustment — `RORADJ*[{GROSSRET0-GRETEXP0}*delUnity + delgret]` |
| 2909 | wage / employment-ratio adjustment — `ELASTWAGE*{[EMPRAT0-1]*delUnity + delempratio}` |

With `delUnity = 0` and `faccum` exogenous at zero, line 2700 forces the capital change to **exactly
zero**: `xcap` is endogenous under the swapped closure but *pinned* regardless. So every run so far
handed the economy a productivity boost while forbidding capital accumulation. TERM.TAB:2705 says this
in as many words — *"to switch off, put faccum endogenous; delUnity is exogenous but unshocked!"*

> **⚠️ CORRECTION to the rows above — I first overstated what `delUnity` gates**, in this section and
> in two chat messages, and it needs saying plainly because it weakens the hypothesis.
> * Rows 2803 and 2909 are **not** "switched off" at `delUnity = 0`. Read
>   `build_dynamics!.jl:116` and `:163-164`: the *endogenous* response terms `delgret` and
>   `delempratio` sit outside the `delUnity` factor and are present either way. Only the **constant**
>   catch-up term `(GROSSRET0-GRETEXP0)` and the **constant** momentum term `(EMPRAT0-1)` are gated.
>   Investment still sees rates of return, and wages still see employment. My earlier claim that
>   "investors get no signal from higher returns" is simply wrong.
> * `delUnity = 1` does not make capital *responsive* either. It sets it to a **predetermined** level:
>   `log(xcap/CAP) = CAPADD/CAPSTOK_OLDP`, last year's investment arriving. That is a different
>   exogenous number, not a new degree of freedom.
> * What survives intact is the **hard pin** at `delUnity = 0` (`build_dynamics!.jl:104-108`) — `xcap`
>   welded to the benchmark. That is real and is the one unambiguous difference between the two runs.
>
> So the run below is a genuine test, not a confirmation of an expected result. Predicting it would be
> guessing.

**This is the leading explanation for the fold** — with the caveat above materially reducing how
strongly it is held. It is neither a translation defect nor a solver defect but an incomplete
scenario. Note `TERM.CMF` is the *test simulation shipped with the
data package* (`readme.txt`: "concluding with a test simulation producing file Term.sl4"), so the
authors intended it to solve.

Translation status: `delUnity` exists and is correctly wired (`build_dynamics!.jl:74,107,116,164`) but
sits in `BASE_CLOSURE_SCALARS` at zero — `initialize_model!.jl:22` documents the base static closure
as leaving the dynamic mechanism switched off.

**NEXT ACTION, outranking task #36:** run the authors' actual scenario, `blabnat = -3` **with**
`delUnity = 1`. Decisive either way — if it converges the fold was an artifact of an incomplete
scenario and #36 is unnecessary; if it still folds we have a documented reference case to compare
against. **Caveat:** `delUnity` is a GEMPACK *change* variable shocked from zero, so check how
`build_dynamics!.jl` handles the change-variable convention in levels before trusting the first run.

**On the literature (searched 2026-07-30):** INDOTERM is real and published (Monash/CoPS with
Universitas Padjadjaran, built on Horridge's TERM; applications to West Java energy prices, rice
import policy, Industry 4.0, Nusantara relocation). **Nothing found documenting folds or limit points
in TERM** — but that is expected and is *not* evidence against one: GEMPACK solves by multi-step Euler
integration of the *linearised %-change* equations with Richardson extrapolation, never solving the
nonlinear levels system to a residual. That procedure has no mechanism for detecting a limit point.
A levels translation can see something the original solution method structurally cannot. Hold this
loosely — the missing `delUnity` remains the more likely story.

### §N.19e — RESOLVED: the authors' reference simulation runs end to end

`test/verify_delunity.jl`, `scratchpad/delunity3.log`. TERM.CMF swapped closure, both shocks.

**Stage A — `delUnity` 0 → 1** (one year forward, no productivity shock). Converged in **3 steps,
0 rejections**, ‖F‖∞ = 2.3e-9. Capital injection `CAPADD/CAPSTOK_OLDP` spans 3.69%–5.47% over the
150 cells (median 5.44%); Java's capital stock +5.56%, `realwage_id` +1.15%.

**Stage B — `blabnat` 1 → 0.97 by arclength, from there.** **Reached the target: 22 steps,
0 rejections**, ‖F‖∞ = 6.73e-10, λ range 0.99495 → 0.97000. Total 2.1 min.

| | `delUnity = 0` (§N.19c) | `delUnity = 1` (this run) |
|---|---|---|
| at λ = 0.9908083 | **turns** — limit point | passes through, no sign change in dλ/ds |
| λ attained | 0.9908083 (extremum) | **0.97 — target reached** |
| rejections | 0 / 60 | **0 / 22** |
| Newton iters/step | 3 | 2 (4 on the final landing step) |
| local coordinate | `delVGDPEXP[2,1]` | `delTAXhou[9,1,6]` → `delPGDPEXP[4,3]` |

`dλ/ds` decays smoothly as `C/s` (−0.0509 → −0.00119) with no sign change and no collapse toward
zero — a long gently-curving branch, not an approach to a limit point. **`blabnat = -3` exists.**

**Why it works — narrower than §N.19d first claimed.** The operative change is the **hard capital
pin**, and only that. At `delUnity = 0` the accumulation equation forces `log(xcap/CAP) = 0`
exactly, welding all 150 capital cells to the benchmark even though `xcap` is endogenous under the
swapped closure. Stage A measures what that pin suppresses: a 3.7–5.5% capital adjustment margin.
Deny it and the economy runs out of room at 0.92%; grant it and 3% is unremarkable. The ROR and
wage-employment mechanisms were **never** switched off (see the correction banner in §N.19d) and
played no part in the rescue. Right prediction, wrong mechanism — the mechanism is what gets
written down.

**Consequences.**
- Task #36 (region-4 income multiplier loop gain) is **unnecessary**: it existed only to explain a
  fold that is not a property of the model.
- `delUnity = 1` is the default for any year-forward scenario. A static/comparative scenario that
  deliberately holds capital fixed must say so explicitly, because `delUnity = 0` is not "neutral",
  it is a hard pin on 150 variables.
- **Process fix owed** — *partially paid 2026-07-30, see §N.22:* the scenario should be specified in
  one place mirroring `TERM.CMF` verbatim,
  so a missing shock is a visible absence rather than an unread line. Natural home is `run_model!`
  when arclength gets wired in as the homotopy fallback. Root cause of the miss: `Shock blabnat`
  and `Shock delunity` are *adjacent* lines (`TERM.CMF:112-113`); task #27 transcribed the closure
  down to the sixth swap at `:110`, took `:112` as "the shock", and stopped one line short. The
  one-scalar-λ shape of `continuation_solve!`/`arclength_solve!` made a two-shock scenario
  structurally invisible — the instrument determined what could be seen in the source.

### §N.20 — price homogeneity (task #11): the first two "violations" were both in the TEST

`test/verify_homogeneity.jl`, logs `scratchpad/homog1.log` … `homog3.log`. Base static closure
(this script does **not** call `apply_swaps!` — say so before comparing any number here against a
§N.19 number). Method: scale the numeraire `NatMacro("GDPPI")` by λ and require prices to scale by
λ while real quantities stay put.

**Run 1 — "6,858 additive violations". Every one was my classifier.** The script had been
hand-classifying by variable *name*: anything matching `del*` was assumed additive and therefore
required not to move. But `delVGDPEXP`, `delPGDPEXP`, `delTAXint`, `delPTX` are **nominal change**
variables — a change measured in rupiah. Under a 10% price scaling a nominal change *must* move by
10%; holding still would be the violation. Independent arithmetic settles it without any classifier
at all: Java's `delVGDPEXP[2,1] = 424,117.8` is exactly **10.000%** of Java household consumption
(4,241,178 Rp bn), and `delPGDPEXP[2,1]` equals it to nine digits, so `delXGDPEXP ≈ 0` — the entire
nominal change is price and the real change is zero. That *is* homogeneity holding.

**The fix — a two-λ proportionality discriminator, replacing the name-based rule.** For any variable
whose benchmark value is zero, solve at two scalings λ₁ = 1.05 and λ₂ = 1.10 and test

  v(λ₂) / v(λ₁) = (λ₂ − 1) / (λ₁ − 1)

It is parameter-free, scale-free, and strictly stronger than "does not move": it also catches a
variable that scales with the *wrong degree*. This is what distinguishes a genuine additive shifter
from a nominal change variable **without anyone deciding by name what a variable means.**

**Run 2 — the additive class went to 0 off, and exposed an ordering bug in the new classifier.**

| class | n | off | verdict |
|---|---|---|---|
| additive-invariant | 39,892 | **0** | ✅ |
| quantity | 33,653 | 267 | 1 defect propagated (below) |
| price | 22,195 | 977 | 801 closure-fixed, 36 zero-weight, **140 real** |
| nominal-change | 6,747 | 115 | |
| zero-base-unexplained | 111 | 111 | **all 111 were the ordering bug** |

The bug: the zero-benchmark branch tested *invariance* before *proportionality*, which gates a
scale-free test behind an **absolute** magnitude cutoff (`ZTOL = 1e-6`) — meaningless in a model
whose variables span rupiah-billions down to tax-rate dust. `delTAXint[4,2,4,5]` came in at
λ₁ = 9.92668e-07, λ₂ = 1.98534e-06 — a ratio of **exactly 2.000**, a perfect pass — and was filed as
a failure because λ₁ landed 0.7% under the cutoff. Fixed by asking proportionality first, in both
the scoring loop and the totals loop. *(Measurement trap #10: an absolute cutoff gating a scale-free
test. Ten of these are now on the record; the recurring shape is a metric whose units were never
checked against the quantity it is judging.)*

**The 176 flat `pdelivrd` prices are structurally indeterminate — measured, not asserted.**
`E_pdelivrd!` (`build_equations.jl:803`) is a share-weighted log average; where the flow is zero,
`BASSHR` and every `MARSHR` are zero too, the equation degenerates to `log(pdelivrd) = 0`, and the
price is pinned at 1.0 forever. Weight sum = 0 and base `DELIVRD` = 0 for 36 of the 48 failing free
cells. Nothing in the model multiplies these prices, so they cannot affect any result — but they are
excluded from the verdict on a **printed weight sum**, never on the strength of their name.

**All 267 quantity failures are ONE defect propagated.** `pcap[5,5]` and `pcap[5,6]` scale by
1.10041153 / 1.10044950 instead of 1.10 — a 4e-4 overshoot — while `pinvitot` in the same cells is
exact. The chain is `gret → ggro → xinvitot → xinvi`, and `E_xinvi!` (`:559`) is
`xinvi[c,i,d] = INVEST_v[c,i,d] * xinvitot[i,d]`, so one bad `xinvitot` fans out to 150 `xinvi`
cells. The independent "other" ratio bucket maxes at 1.10037515 — same magnitude, consistent origin.

### §N.21 — `pcap`: the structural hypothesis is REFUTED (task #38)

Hypothesis: a pinned `plnd` enters the primary-factor CES aggregator with a nonzero weight, giving
`pcap` a term that cannot scale — a genuine break of degree-zero homogeneity.
`scratchpad/probe_cells.jl` (parameters only, no model build) answers it directly:

```
cells with LND ≈ 0 but ALPHA[LND] > 0  (the structural-break signature):
   none — every zero-land cell also carries zero land share ✓
```

So no pinned price carries a live weight, and `E_pcap!` (`:296`) is a CES factor demand —
**homogeneous of degree 0 in prices by construction.** There is no structural route to the error.

What the two failing cells *do* have is size:

```
smallest CAP cells outright:
   Coal MalukuPapua  CAP = 0.0327748     ← fails
   Coal BaliNusa     CAP = 0.040271      ← fails
   OilGas BaliNusa   CAP = 8.97          ← 220× larger, passes
```

The failing set is exactly the two smallest capital stocks in the model, separated from the
third-smallest by a factor of 220. Their capital *share* is healthy (aCAP = 0.9753), so this is not
a small-share story — it is a row carrying ~0.03 Rp bn in a system whose rows run to 1e6, i.e. a
cell whose relative accuracy is limited by a scaled convergence test. That predicts the error
follows the tolerance down.

**Decisive test (running):** rerun at `NTOL = 1e-9` against the `NTOL = 1e-8` baseline in
`homog2.log`. Error shrinks ⇒ numerical, #11 closes as a pass. Error holds at 4e-4 ⇒ real defect in
the capital block. Read `pcap[5,5]`/`pcap[5,6]` against 1.10041153 / 1.10044950 and the "other"
bucket max against 1.10037515.

**A tolerance floor worth knowing.** `NTOL = 1e-11` is *below* the achievable residual floor at 25×6
(~6e-10). `continuation_solve!` accepts a step only when `residual <= tol`, so every step gets
rejected and the script reports "homogeneity is UNTESTED" — which reads exactly like a model problem
and is not one. `verify_homogeneity.jl` now guards this explicitly and exits with a diagnostic
instead. Measured floor: benchmark `status=converged ‖F‖∞ = 6.98e-10`.

### §N.22 — scenarios get one declarative home (task #39, partial)

`src/scenarios.jl`. Pays the process debt §N.19e recorded as owed. `TERM_CMF_REFERENCE` carries
**both** of TERM.CMF's `Shock` lines with their source line numbers attached:

```julia
shocks = [
    "blabnat"  => pct(-3),   # TERM.CMF:112  ! 3% increase labour productivity
    "delUnity" => 1.0,       # TERM.CMF:113  ! move forward 1 year
],
```

Verified against `origin/TERM.CMF` that these are the **only** two `Shock` statements in the file, so
"the list is complete" is a checkable claim rather than an assumption. `describe(sc)` prints closure
and shocks into the run log. The incomplete variant is kept as `TERM_CMF_NO_DELUNITY`, labelled
`INCOMPLETE`, so the §N.19c comparison is reproducible and cannot be rebuilt by accident.

**Wired (2026-07-30).** `scenarios.jl` is in the module include list (after `run_model!.jl`, which
defines `pct`, and `closures.jl`, which defines `_side_name` — both are called at constant-definition
time). `run_model!(agg, params, sc::Scenario)` is the form scripts should use; it prints `describe(sc)`
into the log and **rejects** `name`/`swaps`/`shocks`/`pct_shocks` as keyword overrides, because a
scenario editable at the call site is not a single declarative place.

**The completeness claim is now mechanical, not documentary.**
`test/verify_scenario_complete.jl` parses `origin/TERM.CMF` — stripping `!` comments, matching
`shock` case-insensitively, since the file writes `Shock blabnat` and `Shock delunity` with different
capitalisation — and fails if the source and `TERM_CMF_REFERENCE` disagree on shock count or on which
variables are shocked. It reports `:112 blabnat`, `:113 delunity`, and passes. A **negative control**
was run: fed `TERM_CMF_NO_DELUNITY`, the same comparison reports `detected missing: ["delunity"]`, so
the check can go red — a green test that cannot fail would have been worth nothing here.

**Arclength fallback wired, with its limitation encoded rather than hidden.** On a step collapse
below `hmin`, `run_model!` now retries with `arclength_solve!` — but **only when exactly one shock is
live**, because arclength traces one `VariableRef` and `t` is a bookkeeping scalar, not a model
variable. With several shocks the driver says so and stops. The alternative — promoting `t` to a real
variable tied to each shocked variable by a linear constraint (square, solution-set-identical) —
would compose in general but changes the built model's structure, so it is **not** done here.

For the multi-shock case the working recipe is the staged run §N.19e actually performed: apply the
switch-type shocks first (`delUnity` 0→1, 3 steps), then arclength the continuous one from there
(`blabnat` → 0.97, 22 steps, 0 rejections). Staging is not an approximation — each stage ends on a
genuine solution. Note this means the fallback does **not** fire for `TERM_CMF_REFERENCE` itself,
which has two shocks and does not need it.

On the turned-without-reaching path the driver prints the §N.19e warning inline: check the scenario
is complete before reporting a fold as a property of the model.

**Verified so far:** package precompiles, both scenarios print, the `Scenario` method resolves, the
`arclength` kwarg is on the primitive, and the completeness test passes with a working negative
control.

**Exercised 2026-07-31, and it found a real bug — task #39 is open again, not done.**
`test/verify_arclength_fallback.jl` ran `TERM_CMF_NO_DELUNITY` end to end (35.4 min, 204 steps, 10
rejects, final ‖F‖∞ = 1.6e-9 — the solver itself converges cleanly at every accepted point). It
correctly detected the collapse and correctly fired the fallback, but the branch it then traced is
wrong: `arclength_solve!` reports `turned = true` at λ ≈ 1.00749 — a fold **above** the benchmark
value of 1.0, on the opposite side from the target (`blabnat = 0.97`) — and then, rather than
folding back toward the target as the branch's far side, keeps climbing further away, ending at
`blabnat = 1.0194`. See VV_PLAN.md "Open items" for the standing hypothesis (the fallback restarts
`arclength_solve!`'s bootstrap tangent *at* the natural-homotopy collapse point, which is itself
already at/near the documented fold λ*≈0.9908 — exactly where the bootstrap's own precondition,
"valid HERE because the start point is far from the fold" (`arclength.jl:342-344`), is violated).
Not yet fixed. This is a solver/numerical-method bug, in scope to fix without asking — it does not
change what the model *is*, only whether this one fallback path finds the right branch.

### §N.23 — labour market: `ELASTWAGE`, and a correction

The reference simulation's employment/real-wage split is a **constant 0.4985 across all four
occupations** (0.502, 0.498, 0.499, 0.499), including occupation 1, which moves the *opposite*
direction from the other three. That constant is `ELASTWAGE`, read straight off `TERM.TAB:2899`:
`ELASTWAGE(o) # Elasticity of wage to employment: i.e. 0.5 #`. Recovering a **named input
parameter** from simulated output is a much stronger check on Excerpt 54 than an inferred elasticity.

**Correction to what was reported in chat.** This was first described as "the national labour supply
curve with elasticity 2" — wrong equation, and the reading inverted. There are two candidates and
only one is live under TERM.CMF's closure:

- `E_labslack` (`TERM.TAB:1984`), `realwage_id = [1/0.5]*xlab_id + flabsup_id`, i.e. wage = 2 ×
  employment, is **INERT**: TERM.CMF:110 swaps `delfwage = flabsup_id`, making `flabsup_id`
  endogenous, and `TERM.TAB:1987` says so outright — *"nb above equation does nothing if flabsup_id
  endogenous"*. It absorbs a residual and nothing else.
- `E_delfwage` (`TERM.TAB:2907`, Excerpt 54) is live, translated at `build_dynamics!.jl:158` and
  checked line-by-line against the TABLO. Its coefficient gives wage = **0.5** × employment.

Note the two differ by inversion (2 vs 0.5), which is exactly why the misattribution was easy: the
measured number is consistent with a *misread* of the inert equation. **If a coefficient looks
inverted, check which of the pair the closure left active before concluding the translation is
wrong.**

### D. Housekeeping carried forward
- ~~`solve_newton!.jl:435` … the **second** freeze at `~461` is still `1e-3` — decide both together.~~
  **RESOLVED (§E):** both freezes are inert. `colnrm` is measured on the Ruiz-equilibrated `J_s`,
  where every non-empty column has infinity-norm 1, and the empty-column count is 0. Either
  threshold selects the empty set. Harmless to leave, but neither is a lever.
- `MOI.initialize(evaluator, [... :Constraints ...])` is invalid — `:Constraints` is not an MOI
  feature for `ReverseAD.NLPEvaluator`. This was a diagnostic-script bug, not a model bug.
- ~~The real `NLPEvaluator` construction is **not** in `solve_newton!.jl`; locate it before C3.~~
  **RESOLVED:** it is at `solve_newton!.jl:248`, built from a `MOI.Nonlinear.Model` populated by
  iterating `list_of_constraint_types(m)` → `all_constraints(m, Ftype, S)`, skipping
  `Ftype <: VariableRef`. The three `test/scratch/diag_stall_*.jl` scripts reproduce that iteration order
  exactly so Jacobian row *i* maps to a known equation signature.
- ~~`build_pstras!.jl:173` `DISTANCE = reg1["MAKE"]` remains self-flagged as wrong (Active #4).~~
  **RESOLVED: not a defect.** It was dead code — overwritten by the correct `DISTANCE = ras["DIST"]`
  on the very next line, before the `DISTGONE` loop ever read it, so margins were never affected.
  Removed, with a comment recording why the self-flagged line was harmless.

### E. Publication-readiness verdict, 2026-07-31

Overall fit of the model for publication work, based on progress as of this date (see
`VV_PLAN.md` for the full gate-by-gate detail this verdict is drawn from):

**What blocks a publication claim right now:**
- **No sensitivity analysis (V8)** — every result is a point estimate off 2016-vintage
  elasticities. Any reviewer will ask for this; it doesn't exist yet.
- **No regression suite (V7)** — no `runtests.jl`, no pinned baselines. Nothing currently
  stops a refactor from silently moving a published number.
- **V6 (closure ordering)** was investigated at length (three solver attempts plus one
  closure-matching fix) and ended time-boxed, not confirmed — see `VV_PLAN.md` §V6. Not a
  pass.
- **V1 (Walras identity)** is blocked on methodology, not just unimplemented — the substitute
  test surfaced a genuine near-singularity that hasn't been explained.
- **V4 (34-region aggregation)** is blocked on a memory wall (~16GB LU fill-in); everything
  validated so far is at the aggregated 6-region resolution only.
- **Only one external validation scenario exists** (coalprice, comparative-static, short-run
  closure). V9 licenses confidence in that scenario family, not in every closure/shock
  combination the model might be asked about.
- A known **data defect** (regional GDP identity off up to 36% in the shipped benchmark) is
  documented but not a translation bug — still needs a caveat in any manuscript using regional
  breakdowns.

**Practical read.** This is publishable as a *methods/translation* paper today ("we ported
TERM to Julia in levels form and validated against the authors' own published results to
within measurement-method tolerance") with V9 as the centerpiece. It is **not** ready to
support new policy counterfactuals until V6 is resolved (or the model is shown to genuinely
lack a solution on that branch, which is itself a publishable finding but a different claim),
V1 is resolved or explicitly scoped out, and V8 exists in at least a minimal one-at-a-time
elasticity sweep — those three are what a referee would ask about first if the paper claims
anything beyond "we reproduced a known result."
  Removed, with a comment recording why the self-flagged line was harmless.
