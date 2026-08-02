# HANDOFF — IndotermJulia, 2026-07-26

Written at the end of a session whose only instruction was *"run the model."* Read
`AGENTS.md` (repo root) first for the hard constraints and the hard-won lessons; read
`PLAN.md`'s gate table for project state. **This file covers only what changed today and
what the next person should do first** — it is deliberately narrow and will go stale fast.

---

## 1. One-paragraph state

**The 25×6 gate is green.** `test/solve_benchmark_6reg.jl` passes end-to-end in ~54 s
(43.9 s build + 9.8 s solve), and the pipeline's RAS performance regression is fixed. Getting
there took two things today. First, the run failed the benchmark-replication assertion (drift
2.5e-3 against a 1e-4 threshold) — caused by two defects in `src/solve_newton!.jl` introduced
since the last verified run, **not** in the model equations; both fixed (§3). Second, the
equation count had dropped 90,460 → 83,541 with no explanation on record; that is now
accounted for exactly, and is correct behaviour rather than a regression (§5).

What that does **not** mean: no shock has ever been run to completion, the numeraire swap is
not applied, and the model has only ever been solved at the point it was calibrated to. Every
result so far is a *consistency* result. §7 is ordered to attack that.

---

## 2. What I actually ran

```bash
julia --project=IndotermJulia IndotermJulia/test/solve_benchmark_6reg.jl
```

Verbatim output of the first (failing) run:

```
--- Running data pipeline (Steps 0-4) ---
  max DIFFIND_sc = 1.0302519982092393e-7, max DIFFCOM_sc = 0.20427646684500125
    COM[134] (Aircraft): MAKE_I=4932.66, SALES=4095.95, DIFFCOM_sc=-0.2043
    COM[13]  (Tobacco):  MAKE_I=14628.9, SALES=18010.5, DIFFCOM_sc= 0.1878
    COM[133] (RailwayEqp)                              DIFFCOM_sc= 0.1479
    COM[84]  (LeatherPrd)                              DIFFCOM_sc= 0.1369
    COM[76]  (TobaccoPrd)                              DIFFCOM_sc= 0.1191
--- aggregate_regions! (34→6) ---
--- build_model_full! @ 25×6 ---
Build: 67.7s, vars=122001, cons=83541
--- solve_newton! (square levels system F(x)=0) ---
  squaring: pinned 1 orphan var(s) [natfhou] , deleted 0 dead row(s)
  system: 83541 equations, 83541 unknowns (square)
  initial scaled ||F||_inf = 1.4260876923799515e-9
  it 1: scaled ||F||_inf = 3.725290298461914e-9  alpha = 1.0
Solve: 23.0s  status=converged  iters=1  scaled ||F||_inf=3.725290298461914e-9

Max |Δ| vs benchmark (absolute): 0.004582019988447428  at wlux[2]
Max |Δ| vs benchmark (relative): 0.0024731657976853625  at delVGDPEXP[8]
Raw ||F||_inf: 9.918212890625e-5

ERROR: AssertionError: solution deviates from benchmark by 0.00247 (relative)
       — not replicated                          [solve_benchmark_6reg.jl:79]
```

Logs are committed in the repo (the original scratchpad paths were session-scoped and are
not findable from another session):

- `IndotermJulia/logs/2026-07-26_bench_FAIL.log` — the failing run above
- `IndotermJulia/logs/2026-07-26_bench_PASS.log` — the post-fix re-run

**Post-fix re-run: PASSED.**

```
  initial scaled ||F||_inf = 1.4260876923799515e-9
  start point already satisfies F(x)=0 to tol — no step taken
Solve: 9.8s  status=converged  iters=0
Max |Δ| vs benchmark (relative): 0.0
PHASE 1 GATE PASSED: 25×6 levels model replicates its benchmark ✓
```

Build also dropped 67.7s → 43.9s and solve 23.0s → 9.8s.

> **Read that `0.0` correctly.** It is exact *because the solver now takes no step at all*,
> so the replication assertion is trivially satisfied and carries no information. The
> number that means something is the **initial** residual, **1.426e-9** — evaluated at the
> benchmark point with no solver involvement, and identical in both runs. That is the real
> evidence the model was never at fault. Do not quote the `0.0` as a quality metric.

---

## 3. What I changed — two edits, both in `src/solve_newton!.jl`

Nothing else was touched. No equation file, no pipeline file, no test file.

### Edit 1 — the solver stepped when it was already converged

The start point residual was **1.426e-9**, comfortably under `tol = 1e-8`. The loop had no
pre-iteration convergence test, so it took a full Newton step anyway. That step moved
variables by up to 4.6e-3 and made the inf-norm residual **worse** (1.43e-9 → 3.73e-9).

The step got *accepted* because of a norm mismatch that is worth internalising:

| | norm used |
|---|---|
| line-search acceptance test (`fnew < f0*(1 - 1e-4·alpha)`) | **2-norm** (`norm(gs)`) |
| convergence test (`norm(gs, Inf) < tol`) | **inf-norm** |

At a residual of ~1e-9 you are solving float noise. A step can lower the 2-norm slightly
while blowing up one component — which is exactly what happened. This is why a benchmark
solve, which *starts at the solution by construction*, must check before it steps.

Fix: a guard before the loop, and `for it in 1:(status == :converged ? 0 : maxit)`.

### Edit 2 — `qr` had been promoted from fallback to primary linear solver

Line was `d_eq = -(qr(Jeq) \ gs)`. The documented design (module docstring, PLAN.md §7) is
**`lu` first, `qr` least-squares only as a fallback.**

This matters because SuiteSparseQR's `\` is *rank-revealing*: when it decides the matrix is
rank-deficient it silently returns a **basic** solution with columns zeroed — no error, no
warning, just a wrong step in exactly the near-null-space directions. And PLAN.md records
that SPQR's rank estimate on this Jacobian is unreliable (it reported **71,101** against a
real deficiency of **637**, because its tolerance keys off the largest column norm).

Fix: `lu` first, `qr` in a nested `catch` with a log line when it triggers.

> Either edit alone probably clears the gate. Edit 1 is the proximate cause; edit 2 is why
> the bad step was so large. Both are correct on their own merits — keep both.

---

## 4. Dead code found (left in place, flag for review)

`cgnr` and `cgnr_resid` (`solve_newton!.jl:110–139`) are **defined but never called**. Their
docstring says CGNR is *"critical for the ~637 rank-deficient singular values in this model's
Jacobian."*

That comment is **stale and misleading**: gate 5c-ii resolved the deficiency 637 → 0, and the
run above confirms `lu`-with-full-rank works. My read is that another agent hit the drift
symptom, wrote CGNR to suppress it via minimum-norm least squares, then found the real cause
elsewhere and never wired it in — but left the docstring asserting a deficiency that no
longer exists.

**Recommendation:** delete both functions, or at minimum rewrite the docstring. Someone will
otherwise read it and conclude the rank problem is still open. I did not delete it because it
wasn't in scope for "run the model" and it is harmless at runtime.

---

## 5. RESOLVED — the ~6,900-equation drop is four new zero-flow guards

*(Settled 2026-07-26 by `test/diag_eqcount_6reg.jl`. Kept in full because the reasoning is
reusable, and because the "why you cannot skip this" argument at the end still applies to
any future count change.)*

| when | equations |
|---|---|
| gate 5c-i (core) | 88,599 |
| + dynamics (gate 6) | 89,911 |
| + macros (gate 6c) | 90,460 |
| **today** | **83,541** |

That is **−6,919** against the last verified figure. The system stayed square, so ~6,900 free
variables disappeared *too*. `build_model_full!` does still call `build_dynamics!`
(`build_model!.jl:627`) and `build_macros!` (`:633`) — I checked; they are not simply
unwired.

**Why symmetric loss is the important clue.** The zero-flow guards added in gate 5c-ii have
the shape:

```julia
if FLOW > 1e-10
    @constraint(m, ...)        # one equation
else
    fix(var, 1.0; force=true)  # removes one free variable
end
```

Every extra cell that falls into the `else` branch removes exactly one equation *and* one
unknown. So a symmetric drop of ~6,900 is the fingerprint of **~6,900 more benchmark cells
testing as zero than before**. Either the upstream data changed, or a guard's threshold or
predicate changed.

**Timeline evidence pointing at code, not data:** `build_equations.jl` (Jul 25 19:47) and
`build_model!.jl` (Jul 25 18:19) were both edited *after* the 15:19/15:32 run that recorded
90,460. `ras_balance!.jl` was also rewritten (Jul 26 10:59), but I diffed it and the rewrite
is arithmetically faithful — `sum(generator)` → explicit `@inbounds for` with identical
accumulation order. So I lean toward the evening-of-Jul-25 equation edits.

### The answer

`test/diag_eqcount_6reg.jl` counts variables fixed at *build* time — i.e. guard pins, before
`initialize_model!` applies any closure-exogenous fixes. 9,217 total. Four of the families
are new since the 90,460 measurement (all four guards were added to `build_equations.jl` on
Jul 25 evening), and they account for the delta exactly:

| family | site | pinned | of total | |
|---|---|---:|---:|---|
| `xmake` | `build_equations.jl:922` | 3,600 | 3,750 | 96.0% |
| `xtradmar` | `:773` | 2,992 | 16,200 | 18.5% |
| `xsuppmar` | `:879` | 288 | 1,944 | 14.8% |
| `xtrad` | `:829` | 39 | 1,800 | 2.2% |
| | **sum** | **6,919** | | |

**6,919 == 90,460 − 83,541.** Not approximately — exactly. The remaining 2,298 pins
(`ppur_s` 594, `plnd` 108, `psuppmar_p` 48, `delTAX*` 1,548) are the gate 5c-ii guards that
were already in place when 90,460 was recorded.

So this is a **continuation of the gate 5c-ii zero-flow work by another agent, not a
regression** — and it is economically correct: a cell with zero benchmark flow has zero CES
share and can never become positive under any shock, which is exactly GEMPACK's own
semantics. The 96% figure for `xmake` is not alarming either — `xmake[c,i,d]` is
commodity × industry × region and the make matrix is nearly diagonal by construction.

Regenerate any time with:

```bash
julia --project=IndotermJulia IndotermJulia/test/diag_eqcount_6reg.jl
```

(It caches the pipeline output to `agg6_cache.jls`, so repeat runs skip the ~3-minute data
stage. Delete that file if upstream data changes.)

**Why you could not have skipped this.** A passing benchmark replication will *not* catch it. Pinning a
variable that should be free is invisible at the benchmark — the pinned value **is** the
benchmark value, so the residual stays at 1e-9 and every assertion passes. It only shows up
under a shock, as a variable that refuses to move. Given that the next milestone is precisely
"run the first shock," an unexplained 8% shrinkage of the system is a trap laid directly
across the path.

---

## 6. Corrected status of PLAN.md's four "active blockages"

| # | PLAN.md says | Reality after today |
|---|---|---|
| 1 | `ras_balance!` 124 s regression, fix unconfirmed | ✅ **CONFIRMED FIXED** — pipeline completes; update PLAN.md |
| 2 | Elasticities not wired | ✅ **STALE — they ARE wired**, `prepare_parameters.jl:198-204` writes `SLAB`/`P028`/`P015`/`SMAR`/`PO01`/`SCET`/`P018`. See caveat below |
| 3 | Numeraire swap not applied | ❌ still true — unchanged |
| 4 | Precompile 14–40 s | ❌ still true — unchanged |

PLAN.md's "Current active blockages (2026-07-26)" section should have #1 and #2 moved to Done.

**#2 verified numerically (2026-07-26).** `prepare_parameters.jl:24-28` falls back to
`zeros(T, na)` — *not* to the old placeholder — when a key is missing from `agg`, and σ = 0
is Leontief, a worse silent failure than the σ = 0.5 placeholder it replaced because nothing
downstream complains. So the values were printed:

| key | n | min | max | zeros |
|---|---:|---:|---:|---:|
| `SLAB` (labour CES) | 25 | 0.35 | 0.35 | 0 |
| `P028` (primary factor CES) | 25 | 0.2 | 1.68 | 0 |
| `P015` (Armington dom/imp) | 25 | 0.9 | 9.33 | 0 |
| `SMAR` (margin substitution) | 9 | 0.2 | 0.2 | 0 |
| `PO01` (population) | 6 | 7.00e6 | 1.46e8 | 0 |
| `SCET` (CET output mix) | 25 | 0.5 | 0.5 | 0 |
| `P018` (export demand) | 25 | 1.8 | 20.51 | 0 |

No zeros anywhere, and the dispersed ones (`P028`, `P015`, `P018`) are clearly real data.
`PO01` at 7.0e6–1.46e8 is the right order for Indonesian island-group populations, which
also confirms `aggregate_regions!` summed it rather than averaging it.

One residual ambiguity worth knowing: `SCET` is uniform 0.5 — numerically identical to the
old hardcoded fallback `fill(0.5, na)`. It reaches `build_model!` through the `haskey` branch,
so it *is* coming from the parameter dict, but a uniform value cannot itself distinguish
"GEMPACK's uniform default, faithfully carried" from "someone wrote the placeholder into the
dict". `SLAB` is uniform 0.35 ≠ the 0.5 fallback, which is mild evidence for the former.

---

## 7. Recommended order of work

Items 1 and 2 as originally written (read `bench2.log`; explain the 83,541) are **both done** —
see §2 and §5. What remains:

1. **Print the elasticity values** (§6 caveat). One command, 30 seconds, and it closes the
   last "wrong by construction under shock, invisible at benchmark" risk. Same trap shape as
   §5: benchmark replication cannot detect it, because at a balanced benchmark the elasticity
   is never exercised.
2. **The numeraire swap.** Free `phi`, fix `NatMacro("GDPPI")` (built in gate 6c), per
   `TERM.CMF`. Then a price-homogeneity test: scale the numeraire, confirm all real
   quantities are unchanged and all nominal prices move proportionally. This is the strongest
   correctness check available that does *not* require a real shock, and it exercises exactly
   the code paths benchmark replication leaves untested.
3. **The first shock** — `test/shock_blabnat_6reg.jl`, `blabnat = -3` via 50-step Euler
   continuation. Never run to completion by anyone. Expect this to be where real problems
   surface, which is why 1 and 2 come first.
4. **`build_pstras!.jl:173`** carries a self-flagged defect:
   `DISTANCE = reg1["MAKE"] # not right — need DISTANCE from reg2`. Weigh it against today's
   measured asymmetry — `DIFFIND_sc = 1.03e-7` (essentially exact) vs `DIFFCOM_sc = 0.204`.
   Industry balance is perfect while commodity balance is off 20% on a handful of commodities
   (Aircraft, Tobacco, RailwayEqp, LeatherPrd, TobaccoPrd — all small, all import-heavy or
   excise-heavy). That pattern is suggestive but not proof; do not assume the two are related
   without checking.
5. Deferred, unchanged: Excerpt 49 condensation (prerequisite for 34 regions — 20 regions
   already needs 16.2 GB), gate 7 `run_model!.jl`, reporting excerpts
   30/33/34/36/37/42/43/45/46/48.
6. **Housekeeping:** everything since `c02f2a6` is uncommitted (13 modified + 20+ untracked).
   The tree currently has no restore point.

---

## 8. Two claims in the docs I would not repeat without checking

Recorded because both are mine from earlier sessions and both are weaker than they read.

- **"The equations are verified correct."** True in a specific, limited sense: they are
  internally consistent and correctly calibrated. Benchmark replication **cannot** detect an
  equation that is economically wrong but benchmark-consistent, because at a balanced SAM
  almost any plausible market-clearing form holds exactly. No shock has ever run. Treat
  "verified" as "verified against itself."
- **"The GDP income/expenditure imbalance is a data property, not a translation defect."**
  I asserted this and it was an overclaim — the database is produced by this project's own
  Julia pipeline, and I never separated "present in the original source data" from
  "introduced by the port." What *is* established: `CHKMOD.TAB` computes the same quantity and
  writes it to a DIAG header as a **diagnostic**, not an `Assertion` — so the original authors
  instrumented it deliberately. The **magnitude** (36.1% for MalukuPapua, 1.07% nationally)
  has never been cross-checked against the original's output, which needs a GEMPACK licence
  the project does not have.

Related and still open: `src/build_pstras!.jl:173` carries a self-flagged defect —
`DISTANCE = reg1["MAKE"]  # not right — need DISTANCE from reg2`. `datnotes.txt` documents
`pstras.tab`'s job as *"make tiny changes to factor payments, MAKE and IMPORTS to make
supply = demand."* Today's run shows `DIFFIND_sc = 1.03e-7` (industry balance essentially
exact) against `DIFFCOM_sc = 0.204` (commodity balance off 20%). A known-wrong line in the
stage whose documented purpose is balancing, with a 20% commodity imbalance downstream of it,
is worth chasing before anyone trusts regional GDP numbers.

---

## 9. Environment notes specific to today

- One **idle bare `julia.exe` REPL** (PID 5640, started 07:57, no command-line args) was
  running throughout and is not mine. It is ~187 MB and idle; harmless, but it makes
  "is a run in flight?" checks by process-name ambiguous. I used `count ≥ 2` as the
  in-flight test. Kill it if you want a clean signal.
- Redirected stdout is **block-buffered** — a log file stays completely empty until the
  process exits. Do not read an empty log as a crash. (This is in `AGENTS.md`; it caught me
  again anyway in a previous session.)
- Uncommitted work is extensive: 13 modified files and 20+ untracked ones, including every
  file added since `c02f2a6`. Nothing in this project since Step 5b is in a commit. That is a
  standing risk independent of everything above.
