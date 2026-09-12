# Stage 1 progress notes (Jacobian-level condensation probe)

Running record kept by the opencode agent for easy handover. Newest entries at the bottom.
Related: `HANDOFF-2026-09-12.md` (authorization + traps), `SCALE34_PLAN.md` (assessment + plan).
Probe: `test/scratch/_probe_condense.jl` (untracked scratch). Logs: `logs/` (git-ignored).

## 2026-09-12 ~07:10 — Took charge, updated handoff

- User authorized Stage 1 and put opencode in charge, sandboxed (scratch + logs only, no `src/` edits).
- Recorded in `HANDOFF-2026-09-12.md`: header status note, §1 "Authorized ... in a sandbox",
  §7 renamed to "Decisions and Status" with Stage 1 marked AUTHORIZED.
- Deleted stray `test/scratch/_probe_condense.jl` + backtick file (a broken duplicate with
  `\"`-escaped quotes; the good draft `test/scratch/_probe_condense.jl` was kept).
- Sandbox used: in-repo scratch isolation (`test/scratch/` + `logs/` only). No branch created;
  `master` still holds the 3 unpushed commits plus the handoff edit (uncommitted).

## 2026-09-12 ~07:15 — Diagnosed the earlier failed run

- `logs/probe_condense_2026-09-12.log` (12,206 lines, ended mid-run, process gone) shows the old
  draft's greedy row matcher failed: 12,123 "Could not find unique defining row" warnings.
- Measured shortfall: n=6 PRIMARY S=18,425 cols → 17,711 rows; FULL S=64,766 → 61,641 rows
  (n=10 PRIMARY 56,133 → 53,956). Rectangular D then crashed `D \ B` with
  `MethodError: no method matching ldiv!(::SparseArrays.SPQ...)` (SPQR on non-square).
  The "singularity" label in that draft was wrong — it was a shape error.
- Root cause hypothesis: greedy max-|entry| with row dedup steals shared aggregation rows
  (e.g. `E_xsuppmar_p[m,r,d] == sum(xtradmar[...])`, one row shared by na·ns cells),
  cascading into ~4% unmatched. Structural check of `src/build_equations.jl`
  (`E_xtradmar_na!` :812, `E_pdelivrd!` :826, `E_xtrad!` :864, `E_xsuppmar!` :915) shows every
  free S-cell HAS its own defining equation (fixed cells are `JuMP.fix`ed out of the free
  Jacobian), so a perfect matching should exist — the loss is algorithmic, not structural.
- Within-S dependencies (PRIMARY): xtradmar→xtrad, xtrad→pdelivrd; xsuppmar is a leaf.
  So D should be block-triangular under a topological order.

## 2026-09-12 ~07:20 — Rewrote `test/scratch/_probe_condense.jl`

- Kept `region_map` / `free_jacobian` verbatim from `_probe_lu_ordering.jl`; sizes from ARGS.
- New: per-family free-column census; Hopcroft-Karp max matching on the stored pattern with
  rows explored strongest-first; per-family matched report with shrink-to-matched S;
  topological order of the S-dependency graph (cycles → skip set with family report);
  strict-upper dust check; sparse triangular `D⁻¹B` via hand-written forward substitution
  (dense workspace + touched-list accumulator, output stays sparse — no SuiteSparse API risk);
  1e-14 dust drop on Y; `K̃ = K − A·Y` sparse; default-AMD `lu(K̃)`; one `lu(J)` baseline per
  size; Schur-vs-direct verify (gate ≤1e-10); `Sys.maxrss()` after each heavy step.
- Syntax verified via `Meta.parseall` (scratch `parsecheck.jl` in Temp dir). Scratch files
  are LF (no CRLF trap for these).

## 2026-09-12 07:23 — Launched n=6 smoke run

- Detached: `Start-Process julia --project=. test/scratch/_probe_condense.jl 6`
  (PID 21924; my bash tool has no run_in_background flag, Start-Process is the equivalent).
- Logs: `logs/condense_2026-09-12_n6.log` + `.err.log` (0 bytes so far = still in
  build/Jacobian phase; `say()` flushes every line once emitting starts).
- Expected ~8–15 min. Next: read log tail; if PRIMARY matches perfectly and verifies,
  launch n=10, then n=14 (one size per process).

## 2026-09-12 ~07:30 — n=6 smoke run DONE (4 min). Matching solved, cycles found

- `logs/condense_2026-09-12_n6.log` (116 lines): build+init 76 s, J 84,136² nnz=302,815,
  LU(J) baseline nnz=3,320,701 fill=10.97 (reproduces ordering probe exactly).
- Hopcroft-Karp matched 100% in 0.04–0.09 s: PRIMARY 18,425/18,425, FULL 64,766/64,766.
  The old 4% shortfall was purely the greedy algorithm. Census (n=6 free cols):
  PRIMARY xtradmar 13,208 / xtrad 1,761 / pdelivrd 1,800 / xsuppmar 1,656.
- BUT both topo sorts cycled, so no K̃ was measured yet:
  PRIMARY: 4,811 cyclic nodes (families xtrad, xtradmar);
  FULL: 1,775 cyclic nodes (delPTX, pcst, plab, plab_o, pprim, realwage, xlab, xlab_o,
  xprim, xsuppmar, xsuppmar_d).
- Diagnosis: unweighted HK pairs cells with SHARED aggregation rows (e.g. an xtradmar
  cell matched to `E_xsuppmar_p[m,r,d] == sum(xtradmar)`, a row with 625 S-entries),
  which poisons the S-dependency DAG. Evidence: FULL (different matching) is acyclic on
  xtrad/xtradmar — matching quality is luck-of-the-draw.
- Fix decided: tiered S-degree-filtered matching. True defining rows have O(1) S-entries
  (E_xtradmar_na 2, E_xtrad 2, E_pdelivrd 1, E_xsuppmar 1); shared rows have 625+ (sums
  over c,s) or nr (sums over one region index). Tier 1 keeps rows with S-degree ≤ 6,
  tier 2 (for leftovers) ≤ 4n. No src replay needed; per-family report exposes any
  family whose true rows are legitimately wide. Also switching log separators/symbols
  to ASCII (the UTF-8 box chars mojibake when viewed as ANSI).

## 2026-09-12 ~07:35 — Probe v2 written, n=6 relaunched

- `hopcroft_karp` now takes an optional initial matching; `run_set` gained `nreg` and does
  HK tier1 (rows with S-degree <= 6) then tier2 for leftovers (S-degree <= 4n), plus a
  matched-row S-degree diagnostic (max, count > K1).
- Log symbols ASCII-only (`-` separators, `DinvB`, `Ktilde`, `LU(Ktilde)`).
- Relaunched n=6 detached (PID 26248, 07:31:45) -> `logs/condense_2026-09-12_n6b.log`.
- Launcher note: `Start-Process` prints `Unknown: ChildProcess.kill` but the julia child
  survives (verified: PID live, CPU accumulating, log files created).

## 2026-09-12 ~23:00 — Hermes took over (sandbox agent cold 15 h), fixed the cyclic veto
- n6b result: tiered matching works (FULL 64,766 → 64,718, core |C|=19,418) but a
  1,225-node cyclic cluster (delPTX, pcst, plab_o, pprim, xlab, xlab_o, xprim)
  vetoed the whole set — no K̃ measured yet.
- Fix (scratch-only): cyclic nodes drop back to the core, acyclic remainder measured;
  only a fully-cyclic S skips. Rationale: the cyclic cluster is small; the decision-rule
  number (fill of LU(K̃) for the nr²/nr³ families) doesn't need it. Positional rows
  (§3.2a) remain the follow-up for the cluster itself.
- Relaunched n=6 (`proc_0ac3a43e72b2`) -> `logs/condense_2026-09-12_n6c.log`.
  Next on green n=6: n=10, then n=14, one size per process.

## 2026-09-12 ~23:05 — n6c: acyclic path runs, K̃ singular; D diag min 9.74e-14
- Cardinality-maximal matching pairs columns with dust-coefficient rows. One 1e-14
  pivot poisoned D⁻¹B (verify 4.56e+64 — usefully diagnostic, not a number).
- Fix: pivot-dominance floor (matched |J| < 1e-9 drops to core, per-family report).

## 2026-09-12 ~23:10 — n6e: PRIMARY VERIFIES STRUCTURALLY, fill 4.18 vs 10.97
- PRIMARY (16,869 acyclic cols): D exactly triangular (strict-upper 0), D diag min
  3.82e-07, K̃/K=1.74, **`lu(K̃)` direct success, fill 4.18 (baseline 10.97)**,
  but verify 9.23e-03 — misses the 1e-10 gate along the weak-pivot modes.
- FULL still singular (D diag min 1.10e-09 — weak pivot in an unidentified family;
  floor dropped only 2 pdelivrd pairs, so the offender is elsewhere).
- Next: raise floor toward 1e-6 with a D-diagonal spectrum report, so the level is
  set empirically instead of guessed.

## 2026-09-12 ~23:15 — census: weak pivots are share-composites, not wrong rows
- 3,652 weak (<1e-3) pivots concentrate in delTAXint (16% of family), xint_s (23%),
  xint (9%), xinvi, delTAXhou, xsub — families whose TRUE defining equations carry
  small calibrated shares. Symbolic substitution is indifferent; numerical D⁻¹
  divides by them. Positional rows (§3.2a) would return the SAME equations, so they
  cannot fix these families — only a dominance floor (drop to core) can.
- `xtradmar` has ZERO weak pivots (13,208 cols, sdeg max 1, D diag min 1.29e-02,
  exactly triangular D, Y with zero densification). Added TRADENEST set.

## 2026-09-12 ~23:20 — n6h TRADENEST: clean elimination, fill GETS WORSE
- 13,208/13,208 matched, 0 dominance drops, D exactly triangular, verify 2.98e-07
  (not gate-clean but reasonable).
- **`lu(K̃)` fill 16.06 vs 10.97 baseline — absolute LU nnz 4.44M vs 3.32M.**
  Single-family elimination REFUTED as a fill lever at n=6. Mechanism: the
  eliminated diagonal block was separator material for AMD; deleting it merges the
  graph, and the trade coupling re-enters the core through A·Y (+5% edges, but the
  strategically worst 5%). K̃/K=1.05 in count, far worse in structure.
- Consequence for the decision rule: the question is no longer "eliminate the
  trade nest" but whether the FULL substitute set (all coupling leaves the
  system together) restores separability — which needs the weak-pivot families
  handled (floor keeps them in core, so FULL-as-measurable excludes them;
  their coupling stays). Open: does partial elimination ever win, or must the
  coupling leave entirely? n=10 next for the scaling slope.
