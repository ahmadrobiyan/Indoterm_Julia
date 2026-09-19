# HANDOFF: INDOTERM Scaling to 34 Regions (Stage 2b)
**Date:** 2026-09-13 — **addendum 2026-09-17: Stage 2b.4/2b.5 CLOSED as `non-gating` (user decision). 6 island groups stays the publication model.**
**Status (2026-09-13):** Stage 2a (Seam) COMPLETE; Stage 2b (Scaling) IN PROGRESS — memory wall DISPROVED 2026-09-13 18:03, blocker is shock convergence, not memory.
**Status (2026-09-17):** Stage 2b.4/2b.5 closed `non-gating` — 34-province shocks not required for validation.

## CORRECTION 2026-09-13 18:03 — OOM claim withdrawn
`logs/memledger_34_b.log` (later than this handoff) proves 34-region Schur factorization SUCCEEDS:
- Square system 1,196,221 x 1,196,221; free-J nnz 5,709,864; Ktilde nnz 5,242,253
- `lu(Ktilde) nnz(L)+nnz(U)=234,197,687 peakRSS=15.66 GiB` on ~23.8 GiB machine
- Benchmark CONVERGED `||F||=8.67e-9`, no step needed
- The "structural memory wall" verdict below is therefore WRONG. Do not pursue condensation/Krylov/RAM for memory reasons.

## Convergence status (from memledger_34_b.log)
- Benchmark (zero shock): YES, converged 8.675669960211962e-9.
- Shock `blabnat 1.0 -> 0.99` with `maxit=1, linsolve=:schur`: NO — factorization OK (261.8s, lin_rel 5.51e-11) but line search stalls over 20 backtracks, merit stuck ~0.043, `status=maxit`.

## Plan — Stage 2b revised (convergence, not memory)
- 2b.1 Pipeline scaling: DONE (cache_34reg.jls HIT, 34 regions / 25 sectors verified, no flattening).
- 2b.2 Sparsity audit: DONE (numbers above; 234M < 250M rule, 15.66 GiB < 23.8 GiB).
- 2b.3 Benchmark solve: DONE (converged 8.6e-9).
- 2b.4 Shock convergence — **CLOSED `non-gating` 2026-09-17**: small-shock probe — `blabnat=0.999` (0.1%) stalls at `t=0.6` / `1%` stalls outright; +12% coal under proper `COALPRICE_SWAPS` closure stalls at `t=0.05` (`||F||=1.9e-5`). Line-search / conditioning wall, not memory. See `VV_PLAN.md` addendum.
- 2b.5 Real scenario — **CLOSED `non-gating` 2026-09-17**: `COALPRICE` at 34 regions not required for validation; 6-island stays the publication model. Probe kept as `test/scratch/_solve_34_coalprice.jl` (default `+12%` `@ tol=1e-4` — honest floor at this size).
- Constraint (AGENTS.md): 6-region stays the locked validation model; 34 remains a test instrument until 2b.4-2b.5 pass. No `src/` solver changes without authorization; scratch + logs only.

## 🎯 Goal
Transition the INDOTERM model from 6-region validation to full 34-region production, ensuring it can be solved on available hardware.

## 🛠️ What has been achieved (Stage 2a)
1. **Schur Complement Solver:** Implemented a high-performance linear solver (`src/schur_linsolve.jl`) that eliminates large families of variables to reduce LU fill-in.
2. **The Seam:** Integrated the Schur solver into `solve_newton!.jl` and `arclength.jl` via the `linsolve = :schur` flag.
3. **Verification:** Proven bit-identical equivalence to direct LU on the 6-region `coalprice` scenario and confirmed path-independence (continuity) through turning points (folds).
4. **Numerical Stability:** Implemented `SCHUR_DFLOOR = 1e-3` to prevent divergence on ill-conditioned matrices.

## 🚩 The "Wall" (Stage 2b Findings)
Attempts to run the model at **34 regions** resulted in repeated **Hard OOM (Out of Memory) crashes** during the first Jacobian factorization.

- **Dimension Verification:** Confirmed the model is correctly sized at 34 regions / 25 sectors (earlier "blow-up" was a log interpretation error).
- **Sparsity Analysis:** The "PRIMARY" elimination set, while reducing fill, is still insufficient to bring the 34-region Schur core under the hardware's RAM limit.
- **Verdict:** The system is hitting a structural memory wall. The current "linear algebra shortcut" is not aggressive enough for the 34-region resolution.

## 📋 Remaining Work & Recommendations
To reach the 34-region target, the following paths are available:

### 1. Aggressive Structural Condensation (High Effort, High Payoff)
Move beyond the current "PRIMARY" set. Perform a deep symbolic analysis of the Jacobian to identify a larger set of variables for elimination. This is the "Excerpt 49" approach mentioned in the project history.

### 2. Iterative Solver / Krylov Methods (Medium Effort)
Replace the direct LU factorization of the Schur core with an iterative method (e.g., CGNR or GMRES) combined with a block-preconditioner. This would avoid the memory-intensive factorization step entirely.

### 3. Hardware Expansion
Increase available system RAM or expand the page file/swap space to accommodate the current Schur core size.

## 📂 Critical Files
- `src/schur_linsolve.jl`: The Schur engine.
- `src/solve_newton!.jl` & `src/arclength.jl`: The solver interfaces.
- `src/prepare_parameters.jl`: Now contains hard-coded dimensions to prevent flattening errors.
- `SCALE34_STAGE1_NOTES.md`: contains the original scaling projections.
