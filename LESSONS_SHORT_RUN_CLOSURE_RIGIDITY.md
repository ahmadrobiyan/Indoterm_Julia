# LESSONS — short-run closure rigidity, and the solver failure it causes

**Audience:** any agent (or human) about to run an **economy-wide** shock in this
Julia port, or about to reach for a short-run / rigid-factor closure because it
"looks like the right economic story".

**Status:** empirical, reproduced, with the log evidence quoted below. Not theory.

---

## 1. TL;DR

> **An economy-wide shock needs an economy-wide adjustment margin. A rigid-factor
> closure removes the margins and forces the whole adjustment into *factor prices*,
> where the cells with the smallest fixed-factor shares amplify it by ~200×. That
> produces a near-singular Jacobian, and the Newton solver grinds instead of
> failing — it burns hours at full CPU with no progress and no error.**

Two corollaries that cost the most time:

1. **A grinding solve is not a bug report about your model.** The benchmark still
   solves to 1e-9 and the shock is applied correctly. Nothing is *wrong*; the
   closure simply has nowhere to put the shock.
2. **If the grind is identical at 1/5th the shock size, the closure is the
   problem, not the magnitude.** Do not "fix" it by shrinking the shock. Change
   the closure.

---

## 2. The failure, as it actually appeared

Scenario: government demand **+10%** (`fgovgen => logpct(10)`) under
`COALPRICE_SWAPS` — the project's short-run comparative static.

### 2a. The benchmark is perfect — so nothing is *invalid*

```
system: 84136 equations, 84136 unknowns (square)
initial scaled ||F||_inf = 1.4260876923799515e-9
start point already satisfies F(x)=0 to tol — no step taken
benchmark: status=converged  ‖F‖∞=1.4260876923799515e-9
shocks:
  fgovgen              0.0 → 0.09531018
```

The closure is valid, the shock variable is exogenous, and the shock target is
correct. **Every downstream conclusion must start from here** — if you diagnose
this as "the shock variable must be endogenous" or "the closure is wrong", you
will chase a non-existent bug.

### 2b. Every *shocked* point grinds

```
homotopy: t = 0 → 1  (h0=0.25, tol=1.0e-8)
  initial scaled ||F||_inf = 13887.12318505987          ← huge for a 25%-of-shock step
  it 1 lmtry=1 mu=0.1: LU 0.5s  max|d|=25650.0  ... lin_rel=1.31e-15  Dc_spread=1.04e7
      ls 1 a=1.0 merit=0.0006794 ...
      ls 2 a=0.5 merit=3.554e7 ...       ← line search walks alpha down to 1.9e-6
      ... (ls 2..20)
  it 1: scaled ||F||_inf = 0.017425367972488004
  it 2 ... it 5: LU step |d|=50210.0 exceeds trust radius 33.32 — using damped CGNR
  it 5 lmtry=1 mu=0.01: CGNR 0.6s  max|d|=15.34 ... cgnr_it=200/200     ← maxed
  it 6 ... it 30: same pattern
  ✗ t=0.25 failed (maxit, ‖F‖∞=0.0003559) — h → 0.125
```

The signature, which you should learn to recognise on sight:

| log signature | what it means |
|---|---|
| `initial scaled ‖F‖∞` in the **thousands** for a small t-step | the model is violently sensitive to this instrument |
| `LU step |d|=… exceeds trust radius … using damped CGNR` | the Newton step is unusable; trust region collapsed |
| `cgnr_it=200/200` | the fallback solver is **maxed** — it cannot make progress |
| `Delta=… on red=… (poor-step contract)` repeatedly | steps accepted only after shrinking to nothing |
| `✗ t=… failed (maxit, …)` then `h → h/2`, never clearing | genuine near-singularity, **not** a bad predictor |
| `Dc_spread = 1.04e7` | the documented bad Jacobian scaling (rows span 1e-10…1e12) |

### 2c. The decisive test: same failure at 1/5th the shock

| shock | homotopy step | initial scaled ‖F‖∞ |
|---|---|---|
| gov **+10%** | t = 0.25 | **13 887** |
| gov **+2%** | t = 0.25 | **2 858** |

Ratio of residuals ≈ 4.86; ratio of shocks = 5. **Residual is linear in shock
size ⇒ the obstruction is structural, not a magnitude problem.** At +2% the
corrector still stalls around ‖F‖∞ ≈ 1e-5 with CGNR maxed and never reaches tol.

### 2d. The same shock, same model, different closure — clean

```
closure swaps:
  swap xhouhtot → endogenous, fhou → exogenous (6)
  swap houslack → endogenous, shrBoTnom2 → exogenous (1)
  swap fgovtot  → endogenous, fgovtot3 → exogenous (6)
  swap xcap     → endogenous, faccum → exogenous (150)
  swap finv1    → endogenous, finv4 → exogenous (150)
  swap flabsup_id → endogenous, xlab_id → exogenous (4)
...
  ✅ t=0.25  it=4  ‖F‖∞=6.054e-9  h=0.25  43.3s
  ✅ t=0.625 it=4  ‖F‖∞=9.546e-9  h=0.375 50.9s
  ✅ t=1.0   it=4  ‖F‖∞=3.260e-9  h=0.562 83.3s
reached t = 1 in 3 step(s), 0 rejection(s)
```

**Nothing about the shock changed.** Only the closure. The short-run closure was
the entire problem.

---

## 3. Why it happens — the rigidity mechanism

From `closures.jl`'s own documentation of the bare/rigid closure:

> with `xcap` AND `xlnd` frozen in all 25 sectors, fixed-factor *prices* carry the
> whole adjustment, at **~200× amplification in the cell with the smallest
> fixed-factor share**.

That is the mechanism, stated precisely:

1. An economy-wide demand shock (government spending, import prices, productivity)
   requires **resources to move** between sectors and uses.
2. A rigid-factor closure freezes the two quantity margins that would normally
   move: capital (`xcap`) and land (`xlnd`) — and fixes real wages too
   (`flab_i = flabsupA`, `labslack = flabsup_id` in `COALPRICE_SWAPS`).
3. The market-clearing burden therefore falls on **factor prices**.
4. In sectors/cells where the fixed-factor share is tiny, the price must move
   ~200× to clear — a near-singular Jacobian in exactly those rows.
5. The solver's linear solve degenerates, the trust region collapses, and CGNR is
   maxed out. It does not diverge and it does not error. **It grinds.**

### Important nuance — the closure is not "broken"

`COALPRICE_SWAPS` is correct and well-behaved **for what it was built for**: a
single-commodity world price shock (`COALPRICE_REFERENCE`, externally validated
under V9) and targeted interventions (`HILIRISASI_*`). The rigidity only bites
when the shock is **economy-wide**. Choose the closure to match the *shock's
scope*, not just the story you want to tell.

---

## 4. Diagnostic decision tree

```
Solve misbehaves.
│
├─ Does the BENCHMARK (zero shock) solve to tol?
│    NO  → CLOSURE BUG. Under/over-determined system, or a swap mis-applied.
│          Fix the closure. Stop. (run_model! errors out here by design.)
│    YES → the closure is VALID. Continue.
│
├─ Is the shocked var actually exogenous under this closure?
│    NO  → run_model! raises immediately ("shock X is ENDOGENOUS under this
│          closure"). Pick a closure in which your instrument is fixed.
│    YES → continue.
│
├─ Does the failure residual scale LINEARLY with shock size?
│    YES → STRUCTURAL / CLOSURE RIGIDITY. Shrinking the shock will not help.
│          Change the closure.  ← this session's case
│    NO  → possibly a genuine limit point / fold. Use the arclength fallback
│          (single-shock scenarios only) to test whether the branch turns.
│
└─ Did the corrector hit `cgnr_it=200/200` and `exceeds trust radius`?
     YES → trust-region collapse ⇒ near-singular Jacobian. See §5 fixes.
```

**Remember the two failure regimes are different things:**

- *Benchmark fails* → a **modelling/closure bug**. Legitimately yours to fix.
- *Benchmark exact, shocks grind* → **branch hardness / conditioning**. The model
  is fine; your closure choice is wrong for this shock.

---

## 5. What to avoid (anti-patterns)

### 5.1 Do not reach for the short-run closure for economy-wide shocks

`COALPRICE_SWAPS` (and `TERM_SR_SWAPS`) are the wrong tool for a national
demand/supply shock. The project already recorded this twice before this session —
**read `closures.jl` and V6 before choosing**:

- `verify_closure_ordering.jl` (V6): *"it ran for 35+ minutes burning full CPU with
  zero step output, i.e. stuck deep inside a single Newton solve at one homotopy
  point"* — under `TERM_SR_SWAPS`.
- `demo_scenario_development.jl`: *"do NOT use `TERM_SR_SWAPS` for an economy-wide
  labour shock — its branch folds at λ≈0.9998."*

### 5.2 Do not shrink the shock to make a rigid closure "work"

The linear-in-magnitude result above shows this is futile, and it silently changes
the experiment (a 2% spending shock is a different policy question from a 10% one).
If you have to shrink the shock to get a solve, the closure is the problem.

### 5.3 Do not transcribe a source `.CMF` closure blindly

`origin/LONGRUN.CMF` opens its long-run block with `swap flabsup = flab_io;`.
**There is no variable `flabsup`** — `TERM.TAB:1983–1996` (Excerpt 38) declares
`flabsupA`, `flabsupB` and `flabsup_id`; the plain name is a pre-rename leftover.
`apply_swaps!` raises on undeclared names, so that closure list is **not
implementable** in this port.

The `.CMF` files in this repo are partly recycled from other countries (region
names like `Busan`, `KAINU`, `EKalimantan`; commodity `"construction"`), and they
also carry **cancelling duplicate swaps** — `LONGRUN.CMF` and `NOSHOCK.CMF` both
state `NatMacro("RealInv") = invslack` at lines 39 **and** 52, the same pair both
ways, so they cancel and the net closure does **not** fix national investment.

Read a `.CMF` as a *hint*, verify every name against `TERM.TAB`, and prefer the
tested `closures.jl` constants.

### 5.4 Do not invent a closure and label it as a source file's

If `LONGRUN.CMF`'s closure cannot be transcribed, use the port's own tested
long-run closure and **say so in `notes`**. Never produce something that *claims*
to be a source's closure while being a different experiment.

### 5.5 Do not report an unreachable variant as a result

The short-run government-spending variant does not solve in this port. That is a
finding to *state*, not a number to report. The project's own culture (see
`TERM_CMF_NO_DELUNITY`, labelled `INCOMPLETE`) is to name unreachable/deliberate
variants for what they are.

### 5.6 Do not confuse `pct` and `logpct`

Both solve silently. The declaration is the test:

| declaration | benchmark | use |
|---|---|---|
| `@variable(m, v >= 1e-6)` | 1.0 | `pct(p)` = `1 + p/100` |
| unbounded / additive in a `log(...)` or `exp(...)` | 0.0 | `logpct(p)` = `log1p(p/100)` |

Concretely: `xgov(c,s,d) = XGOV0(c,s,d) · exp(fgovtot + fgov + fgov_s + fgovgen)`
(`build_equations.jl:684`) ⇒ `fgovgen` is **log-space** ⇒ `logpct`. Passing
`pct(10)` = 1.10 would apply `exp(1.10) − 1` ≈ **+200%** instead of +10%, and the
solve would still look perfectly healthy.

### 5.7 Do not forget `delUnity = 1` when accumulation is switched on

Any closure containing `xcap = faccum` needs `delUnity = 1`. At `delUnity = 0` the
accumulation equations **weld all 150 `xcap` cells to their benchmark values** and
the capital-mobility margin is inert — reproducing the pinned-closure failure that
produced the spurious fold at λ\* = 0.9908083.

### 5.8 Do not assume a bare array variable can take a "uniform" shock

`_shock_ref` requires an explicit index for any array with >1 element:

```julia
shocks = ["pfimp" => pct(-10)]                                # ✗ ERROR
shocks = [("pfimp", c) => pct(-10) for c in eachindex(AGGCOM)] # ✓
```

### 5.9 Do not expect the arclength fallback for a multi-shock scenario

`arclength_solve!` traces **one** `VariableRef`. With several shocks live, the
fallback is skipped (and `run_model!` says so). The working recipe for a
multi-shock scenario that folds is to **stage** it: apply the discrete/switch
shocks first, then arclength the continuous one (PLAN.md §N.19e).

---

## 6. The solutions applied here

1. **Match the closure's scope to the shock's scope.** Use
   `TERM_LR_SWAPS_MATCHED` for economy-wide shocks — it restores the
   economy-wide margins: capital mobile (`xcap = faccum`), employment fixed with a
   *flexible real wage* (`flabsup_id = xlab_id`), plus `delUnity = 1`.

   | | rigid (symptom) | long run (solution) |
   |---|---|---|
   | capital `xcap` | exogenous (frozen in place) | **endogenous** (`xcap = faccum`) |
   | adjustment margin | factor prices only, ~200× amplified | capital reallocation + real wage |
   | result | grinds, CGNR maxed, never reaches tol | t = 1 in 3 steps, 0 rejections |

2. **Always subtract a baseline.** Every long-run scenario carries
   `delUnity = 1`, which is a large event on its own (capital stock +5.5646%, real
   GDP +2.4277% with **no policy**). Report the **policy-attributable** effect as
   the difference from a `delUnity = 1`-only run:

   | scenario | raw RealGDP | minus baseline | policy effect |
   |---|---|---|---|
   | baseline (`delUnity=1` only) | +2.4277% | — | — |
   | gov demand +10% | +2.3245% | −2.4277% | **−0.10 pp** |
   | import price −10% | +2.3286% | −2.4277% | **−0.10 pp** |
   | labour productivity +3% | +4.0597% | −2.4277% | **+1.63 pp** |

   Note `CapStock` is **identical** (+5.5646%) in all three: it is the
   accumulation baseline, not a policy response. Any table showing it varying by
   scenario is wrong.

3. **Expect supply-side results from a supply-side closure.** With employment
   fixed and capital set by accumulation, real GDP is supply-determined — demand
   shocks move **prices and composition, not output**. The signal lives in the CPI
   (import price −1.85 pp; government −0.67 pp), not the GDP column. Reporting
   only GDP would badly misrepresent these scenarios.

4. **Verify the instrument is exogenous under the chosen closure.** `fgovtot`, for
   example, is exogenous under the automatic closure and `COALPRICE_SWAPS` but
   **endogenous** under `TERM_CMF_SWAPS` / `TERM_LR_SWAPS_MATCHED` (because of
   `fgovtot = fgovtot3`). That is why the government scenario moves the scalar
   `fgovgen` instead.

5. **Run the internal consistency check.** The national income-vs-expenditure
   identity held at 0.86% / 2.28% / 0.85% across the three scenarios — inside the
   project's 5% sanity bound. Do this before believing any headline number.

---

## 7. Checklist before running any scenario

- [ ] Is the shock **economy-wide** or **targeted**? → choose closure by scope, not
      by story.
- [ ] Did you read `closures.jl`'s docstring for the closure you picked?
- [ ] Is your shock variable **exogenous** under that closure? (`run_model!` will
      tell you immediately.)
- [ ] `pct` or `logpct`? Check the **declaration** (`>= 1e-6` ⇒ ratio ⇒ `pct`;
      unbounded/additive ⇒ `logpct`).
- [ ] Array instrument with >1 element ⇒ **indexed** shock entries.
- [ ] Accumulation closure ⇒ **`delUnity = 1`** included.
- [ ] Multi-shock ⇒ **arclength fallback unavailable**; plan to stage.
- [ ] Got a **baseline** run to subtract?
- [ ] Numeraire recorded (`:gdppi` vs `:exrate`) — real quantities are immune, the
      CPI column is not.
- [ ] If it grinds: check whether residual scales linearly with shock size
      **before** touching the magnitude.

---

## 8. Procedure / tooling lessons (cost real time this session)

1. **Never trust an empty log as "no output".** Julia block-buffers stdout when
   redirected, and in this environment **background-task log capture stays empty
   until process exit** — even for a `println` + `flush(stdout)`, and even with
   shell `>` redirection. Consequences: a `run_in_background` solve gives you *no*
   live signal, and a completed background task can hand back a 0-byte log.

   **Workaround that works:** have the Julia script write its **own** results file
   with an explicit `open(...)` / `println` / `flush` after each stage, then poll
   *that* file. Foreground runs do stream output normally (subject to the tool's
   10-minute cap), and a redirect test (`println` → file) confirms redirection
   works in the foreground.

2. **Distinguish "grinding" from "deadlocked" by CPU time**, not by silence:

   ```powershell
   Get-Process -Name julia | Select-Object Id,CPU,WorkingSet64,StartTime
   ```
   CPU time ≈ wall time ⇒ genuinely computing (this session: 1960 s CPU over
   ~32 min wall). Flat CPU ⇒ actually stuck.

3. **Know the cost model.** `cached_pipeline(6)` is ~0.5 s (1.1 MB cache);
   `build_model_full!` has run 40–180 s; each homotopy step ~40–85 s. Budget
   accordingly and prefer few, well-chosen trials. Three concurrent runs are fine
   on 8 cores / ~2 GB each.

4. **Prefer self-labeling scenarios.** `describe(sc)` prints the closure, every
   shock and the notes before every run, so a log always records which closure
   produced a number — the cheapest guard against the "quoted without naming the
   closure" failure this project has already suffered once.

---

## 9. Where this is already documented in the repo

| source | point |
|---|---|
| `closures.jl` (module docstring) | the ~200× fixed-factor price amplification; the base closure folds at `blabnat ≈ 0.99975` |
| `closures.jl` (`TERM_SR_SWAPS`) | written for one targeted investment shock, **not** an economy-wide shock |
| `test/verify_closure_ordering.jl` (V6) | 35+ min full-CPU grind with zero step output under `TERM_SR_SWAPS` |
| `test/demo_scenario_development.jl` | "do NOT use `TERM_SR_SWAPS` for an economy-wide labour shock"; `blabnat = -5%` "not cheaply reachable" |
| `test/verify_arclength.jl` / PLAN.md §N.19 | a plain homotopy walk gives no signal about *why* a point is hard — it just grinds |
| `TERM_CMF_REFERENCE` notes | `delUnity = 1` is load-bearing; at 0 it hard-pins all 150 `xcap` cells |
| `test/reference/draftreport_coalprice.md` | numeraire choice rebases nominal results; do not compare across numeraires |
| `VV_PLAN.md` V4 | regional near-zero cells are resolution-sensitive — screen before quoting regional signs |

**The single durable rule:** *the closure is the experiment.* A closure that
cannot absorb your shock is not a modelling error to debug — it is a closure to
change.
