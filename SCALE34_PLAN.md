# Reaching 34 regions: Krylov vs Excerpt 49 condensation

**Status: 🟡 ASSESSMENT COMPLETE 2026-09-12. Stage 1 (decisive measurement) proposed, not started.**

Companion to `VV_PLAN.md` and `PLAN.md` §"Region-scaling measurement". Everything the
model currently reports is 25 sectors × **6 island groups**; the data pipeline already produces
the 25 × **34 province** database (`data/cache_34reg.jls`) and the JuMP model builds at that
scale (2.4 M variables, 98 s). What stops a 34-region *solve* is the sparse LU of the Newton
Jacobian. This document assesses the two levers the 2026-09-08 handoff named — a
Krylov/block-structure linear-algebra change vs GEMPACK-style condensation — and lays out a
staged plan whose first stage settles the choice with one measurement.

---

## 1. The wall, measured

### 1.1 Fill-in growth (`test/scale_probe.jl`, 2026-07-25, default `lu`)

| regions | free vars | jac nnz | LU nnz | fill | LU s | peak RSS |
|---|---|---|---|---|---|---|
| 6  | 121,452 | 427,015   | 6.6 M   | 15.4 | 3.0  | 2.5 GB |
| 10 | 258,752 | 846,606   | 18.7 M  | 22.1 | 5.2  | 3.5 GB |
| 14 | 449,396 | 1,397,732 | 76.3 M  | 54.6 | 26.1 | 9.1 GB |
| 20 | 850,502 | 2,527,002 | 176.6 M | 69.9 | 57.5 | 16.2 GB |

LU nnz grows ≈ n^2.7. Fitting free-vars = a·n + b·n² to the 6- and 20-region rows gives
≈ 2.2 M free variables and ≈ 8 M Jacobian nonzeros at n = 34; extrapolating LU nnz at the same
exponent gives ≈ 0.7–0.8 **billion** LU nonzeros, ≈ 12 GB for the factors alone plus the JuMP
model. The machine has **23.8 GB** (not 16 as older notes assumed), so 34 regions is not
categorically impossible on this hardware — it is at the edge, with no headroom for the
arclength fallback or a second model in memory.

### 1.2 Is it the ordering? — No. (`test/scratch/_probe_lu_ordering.jl`, 2026-09-12)

Same Jacobian, every ordering × strategy UMFPACK offers (`logs/lu_ordering_2026-09-12.log`):

| regions | ordering | LU nnz | fill | s |
|---|---|---|---|---|
| 6  | AMD (default)        | 3.32 M  | **10.97** | 0.67 |
| 6  | METIS = CHOLMOD = BEST | 4.42 M | 14.59 | 1.3–2.9 |
| 10 | AMD                  | 15.63 M | **24.74** | 5.6 |
| 10 | METIS = CHOLMOD = BEST | 26.15 M | 41.39 | 7–13 |

Strategy (auto / unsymmetric / symmetric) changes nothing. METIS nested dissection is *worse*
than AMD by 30–70 % and the AMD fill still grows 11 → 25 from 6 to 10 regions. **Ordering is
not the lever.** The default is already the best available.

### 1.3 Is it a few dense rows? — No.

Degree profile of the free-column Jacobian:

| | 6 regions | 10 regions |
|---|---|---|
| rows × cols | 84,136² | 166,582² |
| avg nnz per row | 3.6 | 3.6 |
| rows with degree ≥ 1000 | 6 (all 1193) | 10 (all 1193) |
| cols with degree ≥ 200 | 1 (`phi` = 300) | 71 (`phi` = 500, `pdom` = 236 ×70) |
| densest 1 % of rows hold | 14.3 % of nnz | 15.1 % of nnz |

One degree-1193 row per region (a regional aggregate), no dense columns, and a matrix that is
*extremely* sparse on average. Ordering algorithms handle a handful of dense rows well; that is
not what is filling the factors. What fills them is the **inter-regional coupling structure**:
every `xtrad[a,s,r,d]`, `pdelivrd[a,s,r,d]`, `xtradmar[a,s,m,r,d]` links region `r` to region
`d` at the commodity level, so the region graph is complete. A complete graph has no small
vertex separator, which is exactly the condition under which nested dissection fails and
minimum-degree fill grows super-linearly. This is structural, not a tuning problem.

---

## 2. The two levers, assessed

### 2.1 Krylov / block-structure linear algebra

Replace `lu(J_s)` in `solve_newton!` with a preconditioned iterative solve (GMRES or BiCGSTAB)
whose preconditioner exploits the regional block structure — e.g. block-Jacobi with one LU per
region's intra-regional block, or an additive-Schwarz / Schur-complement preconditioner over
the trade coupling.

*For:* memory is O(nnz(J)) plus the preconditioner, so 34 regions fits trivially; the
per-region blocks are 6-region-sized problems the current code already factorizes in < 1 s.

*Against, and these are serious:*
- **κ(J_s) ≈ 2.6 × 10¹⁰ after Ruiz equilibration** (measured, `solve_newton!.jl:457`). CGNR
  was already tried and abandoned for exactly this reason — at 200 iterations it returned
  `max|dy| = 0.025` against a true `1221`. GMRES converges on singular-value *clusters*, not
  κ², so it is not doomed the way CGNR was, but a preconditioner strong enough for κ ~ 10¹⁰ is
  usually a near-complete factorization, which re-introduces the fill.
- **Inter-regional coupling is strong, not weak.** Block-Jacobi preconditioners work when the
  off-diagonal blocks are small perturbations. For Indonesian provinces, inter-regional trade
  is a large share of every region's use; the coupling block is comparable in weight to the
  diagonal. Expect many GMRES iterations, each a full Jacobian-vector product.
- **Inexact Newton changes the solver's contract.** The trust region, line search, and the
  `lin_rel` acceptance test in `solve_newton!` all assume an exact LU direction. An inexact
  direction near a fold (V6 fold at λ* ≈ 0.99984; the P028 near-singular stalls) will
  interact with those guards in ways nobody has characterised. The arclength fallback also
  factorizes a bordered Jacobian.
- **No precedent in this model family.** GEMPACK solves TERM with a direct solver (MA48).

*Cost:* days for a prototype on saved Jacobians; **weeks** to make it robust inside the
Newton/continuation machinery. *Probability it reaches 34 regions robustly:* low-to-medium.

### 2.2 Excerpt 49 condensation

`TERM.TAB` Excerpt 49 lists the variables GEMPACK eliminates before factorization
(`Substitute x using E_x` — removed from the system, recovered afterwards; `Backsolve` —
removed and recovered). The large-dimension ones:

| variable | dims | size at 34 regions | directive |
|---|---|---|---|
| `xtradmar[a,s,m,r,d]` | na·ns·nm·nr² | 520,200 | Substitute |
| `xtrad[a,s,r,d]`      | na·ns·nr²    | 57,800  | Substitute |
| `pdelivrd[a,s,r,d]`   | na·ns·nr²    | 57,800  | Backsolve |
| `xsuppmar[m,r,d,p]`   | nm·nr³       | 353,736 | Backsolve |
| `xint[a,s,i,r]`       | na·ns·na·nr  | 42,500  | Substitute |
| `ppur[a,s,u,r]`       | na·ns·nu·nr  | 49,300  | Substitute |
| `tuser`, `xmake`, `pmake`, `xint_s`, `ppur_s`, `xinvi`, `aint_s`, … | | | Substitute |

Every one of these is defined by a single explicit equation (`x = f(other variables)`), so
eliminating it is exact block Gaussian elimination — a Schur complement — not an
approximation. *Solving the condensed system gives the same Newton step.*

*For:*
- **Proven at exactly this scale.** GEMPACK solves 34-region TERM with this list on ordinary
  PCs. The condensed core is what MA48 factorizes, and it fits.
- It removes the nr² and nr³ families — the ones whose coupling makes the region graph
  complete — from the factorized system. The remaining core (`pdom`, `pbasic`, `puse`,
  `xtrad_d`, factor markets, macro) is O(na·nr) with coupling through `xtrad_d`/`puse` sums.
- The Newton direction is unchanged, so every guard in `solve_newton!`, the arclength
  fallback, and every validated gate keeps its semantics.

*Against:*
- **In principle a minimum-degree ordering could find the same elimination** — AMD is free to
  pivot the definitional rows first. The measurement above says it does not achieve it (fill
  still explodes), but the *reason* — heuristic ordering vs. a problem-aware block order — is
  exactly what Stage 1 below has to demonstrate. Condensation's advantage is not guaranteed;
  it is likely, and cheap to test.
- Done the GEMPACK way — symbolically substituting `E_xtrad` into every equation that uses
  `xtrad` in `build_equations.jl` — it is **weeks** of work, touches the file the standing
  instruction says not to refactor, and re-opens every V1–V9 gate.

### 2.3 A third route: Jacobian-level condensation (recommended for Stage 1)

The condensation does not have to be symbolic. At each Newton iteration the Jacobian is
assembled as a sparse matrix; partition its columns into **S** (the Excerpt 49 substitute
set) and **C** (the core), and its rows into the defining equations of S and the rest:

```
J = [ D  B ]   D = ∂E_S/∂x_S  (block-triangular: each substituted variable has one
    [ A  K ]                   defining equation, and the definitions form a DAG)
Schur complement   K̃ = K − A·D⁻¹·B      (same size as the core)
solve   K̃ Δx_C = r_C − A·D⁻¹·r_S ,   then   Δx_S = D⁻¹(r_S − B·Δx_C)
```

`D⁻¹` is a sparse triangular solve, never a factorization. `K̃` is factorized with AMD. This
is precisely what GEMPACK's condensation computes, obtained from the assembled Jacobian
instead of from rewritten equations. It needs **no change to `build_equations.jl`**, is a
drop-in replacement for the single `lu(J_s)` call, and can be measured on the existing
6/10/14/20-region Jacobians before any solver code is touched.

Risk: `K̃` may be denser than `K` (the substitution spreads `xtrad`'s CES terms into the
balance rows). The question is whether `nnz(LU(K̃)) ≪ nnz(LU(J))` — an empirical number,
available in a day.

---

## 3. Recommendation

**Do not start with Krylov. Do not start with symbolic condensation.** Start with the
Jacobian-level condensation measurement (§2.3): it costs a day, changes no model code, and its
result is decisive either way:

- If LU fill of the condensed core at 14–20 regions is a small multiple of the core's nnz and
  extrapolates to well under ~8 GB at 34 → implement it as the linear solver in
  `solve_newton!` (one injection point) and go to 34 regions with the direct solver GEMPACK
  itself uses. Krylov is never needed.
- If the condensed core *still* fills catastrophically → the coupling is intrinsic to the core
  and no direct method reaches 34 on this machine. Then Krylov is the only route, and the
  condensed core `K̃` is the operator to precondition (much better conditioned than `J_s`,
  because the identity-like S block is gone). Prototype GMRES on saved `K̃` matrices offline
  before touching the Newton loop.

Either way the symbolic Excerpt 49 rewrite of `build_equations.jl` is **not** on the critical
path and should stay off it.

---

## 4. Staged plan

### Stage 0 — ordering and dense-row check ✅ done 2026-09-12
`test/scratch/_probe_lu_ordering.jl`. Verdict above: ordering is not the lever; fill is
structural. (Runs at 14/20 regions would complete the table but change no conclusion; run
only if a reviewer asks.)

### Stage 1 — Jacobian-level condensation measurement (1–2 days, no `src/` changes)
`test/scratch/_probe_condense.jl`:
1. Build the model at n ∈ {6, 10, 14}, assemble the benchmark Jacobian as
   `_probe_lu_ordering.jl` does, keep the column → variable-family map and a row →
   equation-family map (the latter needs the constraint names `build_equations.jl` already
   attaches, or a positional record of which `E_*!` produced each row — check first).
2. Define S = the Excerpt 49 Substitute + Backsolve set that exists in the Julia model
   (`PLAN.md` line ~985 cross-check lists the name correspondences; `xtradmar` →
   `E_xtradmar_na!`).
3. Verify `D` is nonsingular and block-triangular under a topological order of the
   definitional DAG (if it is not, the substitute set must shrink until it is — record which
   variables drop out and why).
4. Form `K̃`, factorize with AMD, report: size of core, nnz(K̃), LU nnz, fill, time, and the
   error of a full-system solve reconstructed through the Schur path vs direct `lu(J)`
   (must agree to ~1e-10 relative — this is the correctness gate for the method).
5. Repeat with S = only the nr²/nr³ families (`xtradmar`, `xtrad`, `pdelivrd`, `xsuppmar`)
   to learn which eliminations carry the benefit.
6. Fit growth at 6/10/14; if the extrapolation to 34 is under 8 GB, run n = 20 to confirm.

**Decision rule:** proceed to Stage 2a if projected 34-region LU nnz of `K̃` ≤ 250 M
(≈ 4 GB factors); otherwise Stage 2b.

### Stage 2a — condensed direct solver in `solve_newton!` (3–5 days)
- Add a linear-solve seam at the single `F_lu = lu(J_s)` site: `linsolve(J_s, rhs) → dy`,
  default = current LU (bit-identical behaviour), option = Schur path from Stage 1.
- The seam is the **only** change to `solve_newton!.jl`; the standing "do not refactor"
  instruction must be explicitly lifted for this one edit before it is made.
- Arclength fallback: `src/arclength.jl` factorizes the bordered `[J b; cᵀ d]`; route it
  through the same seam (bordered Schur is one extra rank-1 update).
- Regression: every gate in `test/gates.tsv` at 6 regions under `linsolve=:schur` must match
  `linsolve=:lu` to the gate tolerances (V3 path-independence is the sharpest check).

### Stage 2b — Krylov on the condensed core (only if 2a's decision rule fails; 2–4 weeks)
- Offline first: save `K̃` at 10/14/20 regions to disk; prototype GMRES with (i) ILU(τ) on
  `K̃`, (ii) regional block-Jacobi with per-region LU, (iii) block-Jacobi + a dense
  Schur complement over the national/aggregate variables. Report iterations to 1e-10 relative
  residual and memory. Abandon any preconditioner that needs > 50 iterations at 20 regions.
- Only if one converges: wire it behind the Stage 2a seam as `linsolve=:gmres`, then
  characterise inexact-Newton behaviour against the trust-region guards on the 6-region
  gates before any 34-region attempt.

### Stage 3 — first 34-region solve
- Benchmark replication (zero shock) at 25 × 34 with the chosen solver; record RSS, LU time,
  Newton iterations. Requires the 34-region `params` from `cached_pipeline(34)` — already
  cached.
- Then `COALPRICE_REFERENCE` at 34 regions. The V9 reference report is at 6 island groups, so
  the 34-region result is validated by **aggregation consistency**: 34-region results
  aggregated to the 6 groups vs the validated 6-region run (this is the V4 apparatus,
  `test/verify_aggregation_consistency.jl`, pointed at a new pair).

### Stage 4 — re-gate at 34
V1 (Walras), V2 (numeraire), V3 (path independence), V5 (multiplier) re-run at 34 regions; V8
sweep is *not* repeated at 34 (13 solves × 34-region cost) — the 6-region ranges stand.

---

## 5. Budget and hardware

- RAM 23.8 GB. Peak RSS at 20 regions with plain LU was 16.2 GB; the JuMP model itself is a
  significant share and does not shrink under either lever. Stage 1 must report RSS.
- No new packages are required for Stages 0–2a (SparseArrays/UMFPACK only). Stage 2b needs
  GMRES: either a ~100-line in-house implementation or `Krylov.jl`/`IterativeSolvers.jl`,
  neither of which is in the depot — adding one requires registry access and is a decision
  for that stage.

## 6. Not authorized by this document

- Any edit to `src/solve_newton!.jl`, `src/build_equations.jl`, or `src/arclength.jl`.
- Any 34-region solve attempt.
- Adding packages.

Stage 1 is a scratch probe and logs only.
