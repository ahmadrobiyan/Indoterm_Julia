# AGENTS.md — working guide for INDOTERM (Julia port)

Orientation for anyone (human or agent) picking this project up. Read this, then
`IndotermJulia/PLAN.md` (the gate table at the top is the single source of truth for project state).

## What this project is

A Julia translation of **INDOTERM** — a TERM-family, ORANI-G-derived multi-region CGE model of
Indonesia — from its GEMPACK/TABLO source. Authoritative sources live in `IndotermJulia/origin/`:

- `TERM.TAB` — the model itself (equations, `Omit`/`Substitute`/`Backsolve` declarations)
- `TERM.CMF` — the closure (`Exogenous` list, swaps, shocks)

When Julia behaviour and these files disagree, **the TAB/CMF files win**. Several bugs found so far
were the port silently departing from TABLO semantics.

## Hard constraints (user-set — do not violate)

1. **"I need this model translation using level as long as possible."** The **levels** formulation is
   mandatory. Phase 4 (Johansen / %-change linearisation) is **rejected**. Do not propose it.
2. Region count for the validation model is **6 island groups** (34 → 6), locked.
3. Phase 1 (validate at reduced scale) comes **before** scaling back up.

## Current state (2026-07-26)

- **The equations are verified correct.** Benchmark relative residual **1.86e-9**. Do not re-open
  equation correctness without new evidence — run the residual check first (below).
- **Solver is a sparse Newton, not Ipopt.** Ipopt hangs at `iter 0` on the square system.
- **Jacobian rank deficiency: RESOLVED (637→0).** All 637 null directions were numerically empty
  columns from zero-flow variables (`ppur_s`, `psuppmar_p`, `plnd`). Now guarded in `build_equations.jl`
  — missing equation blocks converted to `if/else` with `fix(...; force=true)` at benchmark value.
  Equilibrated Jacobian is **full rank (88,599)**; Newton converges in **1 iter / 2.2 s**.
- Squaring logic + Newton solver are in **`src/solve_newton!.jl`** (moved from scratchpad, gate #5c-i).
- **Active performance issue:** `ras_balance!` slowed from <1s to ~124s after `Pkg.update()`. Root
  cause: `sum(generator)` overhead in Julia. Generator-based loops being rewritten to explicit
  `@inbounds for` — see PLAN.md "Current active blockages".
- **Elasticities not wired:** `SLAB`, `P028`, `SMAR`, `PO01`, `SCET`, `P018` are read from data but
  never stored in `prepare_parameters.jl`'s output dict. `build_model!.jl` silently uses hardcoded
  placeholders. Shock magnitudes will be wrong until fixed — benchmark replication is unaffected.
- **No real shock has been run yet.** The pipeline + benchmark solve works; shock drivers exist in
  `test/shock_blabnat_6reg.jl` but are blocked by the RAS performance regression.

## Hard-won lessons (each cost real time — don't rediscover them)

### Diagnose before solving
A CGE benchmark solve should start *at* the solution. So when a solve misbehaves, **evaluate the
residuals at the start point first** (`test/check_residual_6reg.jl`) — no solver call, fast, and it
settles "are the equations right?" immediately. Every solver pathology so far turned out to be
structural, never a wrong equation.

### Scaling is not optional
Jacobian row scales span **2.7e-10 to 1.7e12**. Value identities carry ~1e6 coefficients; price
equations are O(1). Consequences if you ignore it: absolute tolerances become meaningless, `lu`
throws `SingularException`, and SuiteSparseQR reports a **fake** rank deficiency (71,101 vs the real
637 — its tolerance keys off the largest column norm). Always equilibrate rows before drawing any
numerical conclusion.

### A square system needs a square solver
`#free variables` must equal `#equality constraints`. But note the twist: **zero degrees of freedom
hangs Ipopt** (it is an interior-point *optimizer*; it wants `n_x > n_c`), while a Newton solver on
`F(x)=0` *requires* exactly that. Signature of the hang: log stops after `iter 0`, plus
`Too few degrees of freedom (n_x, n_c)`. Don't wait it out — 56 min produced nothing.

### Two recurring bug classes from TABLO semantics
1. **`Omit` means DELETE.** An `Omit`-ed variable has no equation and sits at its neutral benchmark
   (ratio-type → 1.0, additive → 0.0). Leaving it endogenous under-determines the system. See
   `OMIT_NAMES` in `src/initialize_model!.jl`.
2. **Guarded zero-flow skips.** Blocks that emit `@constraint` only where a benchmark flow is
   non-zero (`PUR > 1e-10 || continue`) while the variable is declared **densely** leave every
   skipped cell an unknown with no equation. Fix: `if/else`, pinning the skipped cell at its
   benchmark. Fixed for `E_delTAX*!`; also fixed for `ppur_s`/`psuppmar_p`/`plnd` (gate 5c-ii).

### `sum(generator)` is 50× slower than an explicit `@inbounds for` loop in Julia
JuMP's generator syntax (`sum(TRADE[c,s,:,d])` or `sum(TRADE[c,s,r,d] for d in 1:NR)`) creates
an anonymous closure and allocates per-iteration — fast for <100 iterations, catastrophic at
4386+ iterations. `ras_balance!.jl`'s Phase 1 loops dropped from ~2s to <1s just by rewriting
5 nested-generator sums as explicit `@inbounds for` loops. Always use explicit loops for hot
accumulator code in Julia. Profile first with `@time` on a representative subsection.

### Localising under-determination
Use **Hopcroft-Karp bipartite maximum matching** (constraints × free variables). Unmatched free
variables name the under-determined families; unmatched constraints are redundant or dead rows.
Script: `scratchpad/diag_matching.jl`. Check the fingerprint of unmatched constraints — if it is
*empty*, they are all-fixed `0 == 0` dead rows, **not** Walras redundancy, and
`dependency_detector="mumps"` will not remove them (use `JuMP.delete`).

## Environment gotchas (Windows)

- **Julia build + run is ~3-4 min**; stdout is **block-buffered** when redirected, so a log stays
  empty until the process exits. Ipopt's own `output_file` flushes per iteration — use it for live
  progress.
- **Ipopt locks its `output_file`.** Reusing the name across runs gives
  `INVALID_OPTION: Error opening output file`, and the log reads as `Device or resource busy` until
  the process exits. Use a fresh filename per run. (Note: the default solver is now Newton/Ipopt;
  this only matters if you run Ipopt directly.)
- **Julia soft scope**: accumulator loops (`n += 1`) at top level in a script silently become local
  and throw `UndefVarError`. Wrap loop drivers in a function. This has bitten repeatedly.
- **Don't stack waiters.** A `Bash` wait that times out keeps running in the background. Stop the old
  one (`TaskStop`) before launching another, or you accumulate orphaned pollers watching log strings
  that will never appear.

## Verification commands

```bash
julia --project=IndotermJulia IndotermJulia/test/run_pipeline_6reg.jl
```

```bash
julia --project=IndotermJulia IndotermJulia/test/check_residual_6reg.jl
```

```bash
julia --project=IndotermJulia IndotermJulia/test/solve_benchmark_6reg.jl
```

```bash
julia --project=IndotermJulia IndotermJulia/test/shock_blabnat_6reg.jl
```

## Conventions

- Flow-vs-index distinction matters: genuine value **flows** are seeded at their benchmark flow via
  `benchmark_levels(params)`; **index/ratio** variables stay at 1.0; additive shifters/`(change)`
  deltas at 0.0. Getting this wrong is what originally made the benchmark infeasible.
- Ratio-type variables are identified by their `>= 1e-6` lower bound (`_is_ratio_type`).
- `u`-index order: `u_hou = na+1`, `u_inv = na+2`, `u_gov = na+3`, `u_exp = na+4`.
