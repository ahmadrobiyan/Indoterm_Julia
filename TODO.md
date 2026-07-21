# IndotermJulia — TODO

> Working task list companion to `PLAN.md` (which holds the full design rationale, learnings, and
> status history). This file is the short, checkable list of what's next. Update it as work lands;
> keep the narrative/rationale in `PLAN.md`.

## Housekeeping (do first)

- [x] ~~**Commit or stash Step 4/5 work.**~~ — committed 2026-07-22 as `e5381a9` ("Add Step 4/5 core
      translation: derived parameters and CGE equations"), covering `build_equations.jl`,
      `build_model!.jl`, `prepare_parameters.jl`, the `PLAN.md`/`TODO.md` updates, the P021/FRISCH
      fix, and `test/run_full_model.jl`. Not yet pushed to `origin/master` — ask before pushing.

## Step 5 completion (core equations) — current focus

- [x] ~~Confirm `test/run_pipeline.jl` (Steps 0-4) still passes end-to-end on this machine~~ —
      **it failed 2026-07-22** with `MethodError: no method matching Float64(::Vector{Float64})` in
      `prepare_parameters.jl:29`. Root cause: `P021` (the Frisch/marginal-budget-share parameter) is
      genuinely region-specific in the data (`build_premod!.jl` stores it as
      `NamedArray(FRISCH_d, REG, (:DST,))`, length `nr`=34), but `aggregate_model!.jl`'s pass-through
      block (line 216) unwraps it to a plain `Vector{Float64}` before `prepare_parameters!` sees it.
      `prepare_parameters.jl` assumed a scalar national Frisch parameter — even its "working" branch
      (`vec(...)[1]`) would have silently kept only region 1's value and discarded the other 33.
      **Fixed**: `P021_v` is now carried as an `nr`-length vector throughout, and `FRISCH` in the
      `BLUX`/`SLUX` loop is indexed by region (`FRISCH[d]`) instead of being a single scalar. Re-ran
      full pipeline (Steps 0-4) after the fix — passes end-to-end, "Derived 83 parameters ✓".
- [x] ~~**Step 5a — build at full 25×34 scale.**~~ — **done 2026-07-22.** First attempt hit
      `KeyError: key "INVEST_C" not found` in `build_model!.jl:49`. Root cause: `2PUR` (the
      investment-by-commodity-by-industry matrix) passes through unchanged from `reg1` → `reg2` →
      `ras_balance!` (confirmed present in all three), but `build_pstras!.jl`'s output dict simply
      never included a `"2PUR"` entry — the only stage in the chain that dropped it. Downstream,
      `build_premod!.jl`'s `haskey(pstras,"2PUR") ? ... : nothing` silently degraded to `nothing`,
      `aggregate_model!.jl` never populated `agg["2PUR"]`, and `prepare_parameters.jl`'s
      `INVEST_I`/`INVEST_C` block (guarded on `V2PUR !== nothing`) silently skipped — so the params
      dict was missing two keys that `build_model!.jl` reads unconditionally. **Fixed** by adding
      `"2PUR" => haskey(ras,"2PUR") ? ras["2PUR"] : nothing` to `build_pstras!.jl`'s output dict,
      matching the existing `BSMR`/`UTAX` pass-through pattern. Re-ran: full pipeline + full model
      build now succeeds — **98.0s build, 2,376,562 variables, 1,368,512 constraints, 169 vars-dict
      entries.** The ~1.01M vars-minus-constraints gap is expected at this stage (not a bug): GEMPACK
      CGE closures always have more variables than core equations, with the closure (Step 5b) fixing
      exactly that many variables exogenously to square the system — this hasn't been wired yet.
- [ ] Audit variable/equation coverage: `build_model!.jl` declares ~168 variable arrays; confirm
      every variable either has a defining equation in `build_equations.jl` or is intentionally left
      for closure (exogenous fix). Gaps here will show up as an under-determined system at solve time
      — better to catch them by inspection first.
- [x] ~~Check whether Excerpts 30–48 are needed for a basic square closure~~ — **read directly,
      2026-07-22: they are not.** Excerpts 30–37 and 40–48 are entirely regional→national
      aggregation, contribution/decomposition reporting (`MainMacro`, `NatMacro`, Keller
      decomposition, terms-of-trade contributions, inter-regional trade reporting) and accounting
      `Assertion`/`Write ... to file SUMMARY` diagnostics (Excerpt 41's `CHECKA`-`CHECKE`, of which
      `CHECKA`/`CHECKB` are already ported to `prepare_parameters!`; `CHECKC`/`D`/`E` are not, but are
      the same pattern and low priority). None of these define a variable that any core Excerpt
      1–29/38–39 equation depends on for its own solution — they're all one-way reads *from* the core
      solution. **Excerpt 49 ("Condensation actions")** is GEMPACK's own `Substitute`/`Backsolve`
      variable list — a useful cross-check: every name it lists should already have a corresponding
      `E_*!` function in `build_equations.jl` (`Substitute` names are algebraically-eliminated in
      GEMPACK for solver efficiency only; JuMP/Ipopt has no need to eliminate them, keeping them as
      ordinary variables+constraints is correct and simpler — same call PLAN.md's solver-decision
      section already made for the overall MCP-vs-Ipopt question). Remaining excerpts 50–55 are Step
      6 (dynamic 50–54, district 55), out of Step 5's scope.
- [x] ~~**Cross-check Excerpt 49's variable list against `build_equations.jl`.**~~ — **done
      2026-07-22, no real gaps found.** All 11 initial grep "MISS" hits were false positives:
      - `xtradmar`: implemented, just under the function name `E_xtradmar_na!` (matches the historical
        naming-mismatch bug already logged in `PLAN.md`'s bug-fix list), called from
        `build_model!.jl:468`.
      - `xhou_s`, `xhouh_s`: both defined by `E_xhouh_s_agg!` (`build_equations.jl:304-308`), which
        constrains `xhou_s[c,d]` directly — the single-household simplification collapsed the
        household-type-disaggregated equation into this one, it didn't drop the constraint.
      - `contCPI, continccom, contincind_d, contMainMacro, contnatxtot, contxprim_i, xrowdem,
        xrowdem_d`: none of these are declared as variables anywhere in `build_model!.jl` or
        `build_equations.jl` — confirms they're genuinely Excerpt 30-40 reporting/decomposition-only
        and out of scope for the core square system, per the point above.
      All remaining Excerpt 49 names were already "OK" in the first pass. **Step 5's equation coverage
      against Excerpt 49 is complete** — no missing core equations identified.

## Step 5b — Ipopt solve wiring

- [ ] `initialize_model!.jl` — encode one `.cmf` closure (start with `TERM.CMF`, per `PLAN.md`'s
      `default_closure()` convention) as the exogenous/endogenous variable split; set starting values
      (0 for %-changes at benchmark).
- [ ] `solve_model!.jl` — Ipopt feasibility solve (no objective), mirroring WayangJulia's
      `set_attribute` tuning (`max_iter`, `tol`, `constr_viol_tol`, adaptive `mu_strategy`).
- [ ] `run_model!.jl` — homotopy/warm-start driver, only if a direct solve fails to converge on a
      real shock.

## Verification (per PLAN.md's protocol — do in this order)

- [ ] **Benchmark replication**: solve with all shocks at 0 → expect ~0% change everywhere.
- [ ] **Price homogeneity test**: shock only the numéraire → nominal moves, real variables don't.
- [ ] **GDP-both-sides check**: after any shock, income-side GDP == expenditure-side GDP
      (`calculate_gdp.jl`, Step 8 — can be stubbed early just for this check).
- [ ] **Smoke-test scenario**: run one named `.cmf` shock (`TERM.CMF` or `sim1.cmf`) and sanity-check
      signs/magnitudes.

## Step 6 — Dynamic + district extensions (after Step 5 verified)

- [ ] `build_dynamics!.jl` — Excerpts 50-54 (investment rule, real-wage adjustment).
- [ ] `build_district!.jl` — Excerpt 55 (top-down district extension; likely a WAYANG-style
      non-feedback satellite based on its name — confirm by reading that excerpt).

## Step 8 — Reporting

- [ ] `calculate_gdp.jl` — GDP income/expenditure decomposition, benchmark-check reporting, per
      `CHKMOD.TAB`'s own diagnostic conventions.

## Smaller loose ends (from PLAN.md's "Immediate follow-up items")

- [ ] Write aggregated 25×34 output to a GEMPACK `.har` file for cross-checking against the
      original TABLO run (if a licensed GEMPACK install ever becomes available).
- [ ] RAS outer loop: current `build_pstras!` is single-pass (`converged=false` expected); a
      multi-pass loop feeding `DGON` back into reg2's distance weighting would tighten
      DOMERR/IMPERR discrepancies. Not urgent — data-quality refinement, not a blocker.
- [ ] Add `P015` (`ARMSIGMA`, the Armington elasticity) to `build_premod!` output — currently read
      but not stored; `build_model!.jl` currently hardcodes `sigmadomimp = fill(5.0, na)` instead.
- [ ] Once `HeaderArrayFile.jl` is confirmed unused elsewhere, drop it from `Project.toml` (currently
      kept only because Step 1 originally depended on it before the harpy-CSV fallback).
