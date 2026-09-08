# AGENTS.md — working guide for INDOTERM (Julia port)

Orientation for anyone (human or agent) picking this project up. Read this, then
`IndotermJulia/VV_PLAN.md` (the gate table at the top is the single source of truth for
project/validation state — this file covers orientation, hard constraints, and lessons learned;
it does not duplicate gate status).

## What this project is

A Julia translation of **INDOTERM** — a TERM-family, ORANI-G-derived multi-region CGE model of
Indonesia — from its GEMPACK/TABLO source. Authoritative sources live in `IndotermJulia/origin/`
(read-only reference — do not edit):

- `TERM.TAB` — the model itself (equations, `Omit`/`Substitute`/`Backsolve` declarations)
- `TERM.CMF` — the closure (`Exogenous` list, swaps, shocks)

When Julia behaviour and these files disagree, **the TAB/CMF files win**. Several bugs found so far
were the port silently departing from TABLO semantics.

## Hard constraints (user-set — do not violate)

1. **"I need this model translation using level as long as possible."** The **levels** formulation is
   mandatory. A %-change/Johansen linearisation rewrite is **rejected**. Do not propose it.
2. Region count for the validation model is **6 island groups** (34 → 6), locked. Sector count
   locked at **25**. Changing either is a modelling decision, not a test-methodology choice — it
   requires explicit user authorization, not an agent's unilateral call.
3. Phase 1 (validate at reduced scale) comes **before** scaling back up to finer region resolution
   (12, 34 regions remain test instruments only — see `VV_PLAN.md` V4 for how 12-region was used
   this way).

## Current state

Status is tracked in `VV_PLAN.md`'s gate table (V1-V9+), not duplicated here — read that file for
what currently passes, fails, or is a disclosed limitation. As of the last update:

- **Benchmark reproduces cleanly.** Sparse Newton solver, not Ipopt (Ipopt hangs at `iter 0` on
  this square, zero-degrees-of-freedom system — see "A square system needs a square solver" below).
- **Elasticities are wired.** `SLAB`, `P028`, `SMAR`, `PO01`, `SCET`, `P018` are read from data and
  stored in `prepare_parameters.jl`'s output dict (85 parameters total, elasticities included).
- **Real shocks have been run**, not just benchmark replication — e.g. `coalprice.CMF` (+12% coal
  export price) and `TERM_CMF_REFERENCE`/`TERM_CMF_NO_DELUNITY` scenarios (`src/scenarios.jl`),
  including full continuation/homotopy paths and multi-round regional-consistency diagnostics
  (`VV_PLAN.md` V4). `delUnity=1` (a full year of capital accumulation) and `ELASTWAGE=0.5`
  (labour/wage split) are both load-bearing, verified closure parameters — not placeholders.
- **A known, disclosed limitation exists (V4):** region-level *sign* results for a small set of
  near-zero flow/GDP variables (concentrated in BaliNusa) are resolution-sensitive — see
  `VV_PLAN.md`'s V4 section for the full mechanism, magnitude analysis, and root cause. National
  results and the large majority of regional results are unaffected. `src/regional_confidence.jl`
  provides `print_regional_confidence_report(vals, bmk, nr)` as a cheap, single-solve screen for
  flagging candidate near-zero/resolution-fragile region-level results in **any** scenario's
  output — run it after any `run_model!` call whose region-level results you intend to report by
  sign/direction.

## Hard-won lessons (each cost real time — don't rediscover them)

### Diagnose before solving
A CGE benchmark solve should start *at* the solution. So when a solve misbehaves, **evaluate the
residuals at the start point first** — no solver call, fast, and it settles "are the equations
right?" immediately. Every solver pathology found in this project turned out to be structural,
never a wrong equation.

### Scaling is not optional
Jacobian row scales span roughly **1e-10 to 1e12**. Value identities carry large coefficients;
price equations are O(1). Consequences if you ignore it: absolute tolerances become meaningless,
`lu` throws `SingularException`, and SuiteSparseQR can report a **fake** rank deficiency (its
tolerance keys off the largest column norm). Always equilibrate rows before drawing any numerical
conclusion.

### A square system needs a square solver
`#free variables` must equal `#equality constraints`. But note the twist: **zero degrees of
freedom hangs Ipopt** (it is an interior-point *optimizer*; it wants `n_x > n_c`), while a Newton
solver on `F(x)=0` *requires* exactly that. Signature of the hang: log stops after `iter 0`, plus
`Too few degrees of freedom (n_x, n_c)`. Don't wait it out.

### Two recurring bug classes from TABLO semantics
1. **`Omit` means DELETE.** An `Omit`-ed variable has no equation and sits at its neutral benchmark
   (ratio-type → 1.0, additive → 0.0). Leaving it endogenous under-determines the system. See
   `OMIT_NAMES` in `src/initialize_model!.jl`.
2. **Guarded zero-flow skips.** Blocks that emit `@constraint` only where a benchmark flow is
   non-zero (`PUR > 1e-10 || continue`) while the variable is declared **densely** leave every
   skipped cell an unknown with no equation. Fix: `if/else`, pinning the skipped cell at its
   benchmark.
3. **Never translate a TABLO equation suffix positionally when it indexes a set by NAME.**
   Positional translation silently misaligns the equation the moment set ordering differs from
   assumption — caught late in this project (Step 6/8, regional GDP identity work).

### `sum(generator)` is much slower than an explicit `@inbounds for` loop in Julia
JuMP's generator syntax (`sum(TRADE[c,s,:,d])` or `sum(TRADE[c,s,r,d] for d in 1:NR)`) creates an
anonymous closure and allocates per-iteration — fast for <100 iterations, very slow at thousands.
Always use explicit loops for hot accumulator code in Julia. Profile first with `@time` on a
representative subsection.

### Localising under-determination
Use **Hopcroft-Karp bipartite maximum matching** (constraints × free variables). Unmatched free
variables name the under-determined families; unmatched constraints are redundant or dead rows.
Check the fingerprint of unmatched constraints — if it is *empty*, they are all-fixed `0 == 0`
dead rows, **not** Walras redundancy, and a symbolic dependency detector will not remove them
(use `JuMP.delete`).

### Never pin a price at a constant
Pinning a price variable at a fixed value (rather than letting the model's own numeraire/closure
determine it) breaks the model's degree-zero homogeneity in prices — it silently freezes what
should be a live variable and produces plausible-looking but wrong results on downstream cells.
Anchor to a share/ratio quantity instead if a variable needs a numerical floor.

### Never infer rank or singularity from an iterative-solver leftover
A stalled or slow-converging iterative solve (e.g. CGNR) is not evidence the underlying system is
rank-deficient or inconsistent — it may simply be dominated by a near-singular mode that a direct
solve handles fine. Check with a direct factorization (LU) residual before concluding a defect.

### Region-scaling has a hard wall — LU fill-in, not variable count
Sparse LU fill-in ratio grows steeply with region count (measured: ~15× at 6 regions, ~70× at 20),
driving peak memory well past typical machine limits before regions get anywhere near the model's
native 34. Sparse-matrix declaration alone cannot fix this — it is a structural property of the
elimination ordering. Reaching the full 34-region resolution requires a separate condensation
effort, not just "more patience" or "more RAM" at the margin. Treat any resolution above 6 as a
test instrument (see `VV_PLAN.md` V4), not a target for the locked validation model, unless that
condensation work is explicitly undertaken.

## Environment gotchas (Windows)

- **Julia build + run is a few minutes**; stdout is **block-buffered** when redirected, so a log
  stays empty until the process exits. If a solver has its own `output_file` option, it flushes
  per iteration — use it for live progress instead.
- **Julia soft scope**: accumulator loops (`n += 1`) at top level in a script silently become local
  and throw `UndefVarError`. Wrap loop drivers in a function. This has bitten repeatedly.
- **Don't stack waiters.** A shell wait that times out keeps the underlying process running in the
  background. Stop the old one before launching another, or you accumulate orphaned pollers
  watching log strings that will never appear.

## Verification

`VV_PLAN.md`'s gate table is the authoritative status source — read it first. The gates themselves
are Julia scripts under `test/verify_*.jl` (one file per gate, e.g.
`verify_aggregation_consistency.jl` for V4, `verify_multiplier.jl` for V5, `verify_homogeneity.jl`,
`verify_walras.jl`, `verify_numeraire.jl`, `verify_closure_ordering.jl`, `verify_delunity.jl`,
`verify_coalprice_reference.jl`, and others — see `test/` for the full list). Run an individual
gate directly, e.g.:

```bash
julia --project=. test/verify_aggregation_consistency.jl
```

Or run a whole tier through the runner (`test/gates.tsv` assigns each gate to `fast`,
`integration`, `full` or `nongating`):

```bash
bash scripts/run_gates.sh fast
```

**Read this before trusting a green run.** Most `verify_*` scripts print their verdict and
**exit 0 even when they fail** — an exit status alone certifies nothing. `run_gates.sh`
compensates for gates marked `marker` by requiring exit 0 **and** a `✅` **and** no
`❌`/`FAILED`/`FAILS`. That is a stopgap, and it means a green tier reads as "no gate printed a
failure", which is weaker than "every gate asserted a pass". Do not simplify the runner back to
an exit-code check, and do not add a new gate without deciding its row in `gates.tsv`. The real
fix is converting those scripts to throwing `Test.jl` `@testset`s; until then the weaker
guarantee is the one you have.

`test/pipeline_cache.jl` provides `cached_pipeline(nr; rmap)` to avoid rebuilding the data pipeline
on every run — most gate scripts use it. For a quick region-level sanity check on any solved
scenario, use the early-warning screen described above:

```julia
using IndotermJulia
agg, params = cached_pipeline(6)
r = run_model!(agg, params, scenario)
print_regional_confidence_report(r.values, benchmark_levels(params), 6)
```

## Conventions

- Flow-vs-index distinction matters: genuine value **flows** are seeded at their benchmark flow via
  `benchmark_levels(params)`; **index/ratio** variables stay at 1.0; additive shifters/`(change)`
  deltas at 0.0. Getting this wrong is what originally made the benchmark infeasible.
- Ratio-type variables are identified by their `>= 1e-6` lower bound (`_is_ratio_type`).
- `u`-index order: `u_hou = na+1`, `u_inv = na+2`, `u_gov = na+3`, `u_exp = na+4`.
- Always name the model closure explicitly when reporting a result or discussing a discrepancy —
  the base static closure and `TERM_CMF_REFERENCE` differ in which blocks are inert (e.g. the
  capital block), and the same cell can read very differently under each (`VV_PLAN.md`, "two
  closures" finding). "The result under closure X" is a meaningfully different claim from "the
  result."
