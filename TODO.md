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
- [x] ~~Audit variable/equation coverage~~ — **done 2026-07-22 via the Exogenous-list cross-check
      below; found and fixed one real gap (`srctwist`/`avesrctwist`).**
- [x] ~~**Cross-check `TERM.CMF`'s default-closure `Exogenous` list (48 names) against the 169
      declared JuMP variables.**~~ — **done 2026-07-22.** 7 names came up unmatched:
      `delfwage_o, delUnity, emptrend, frnorm, frnorm_id, gtrend, srctwist`. Traced each by line
      number against `TERM.TAB`'s own `! Excerpt N of TABLO input file: !` markers:
      - `delfwage_o` (line 2906), `delUnity` (2696), `emptrend` (2878), `frnorm`/`frnorm_id`
        (2724-2725), `gtrend` (2732) all fall inside Excerpts 50-54 (lines 2664-2920) — confirmed
        Step 6 dynamic scope (investment rule / real-wage adjustment), correctly out of Step 5.
      - `srctwist` — genuine gap, see below.
- [x] ~~**Fix `E_xtrad!`: missing `srctwist`/`avesrctwist`/`SIGMADOMDOM` terms.**~~ — **done
      2026-07-22.** `TERM.TAB` lines 985-1005 contain *two* candidate `E_xtrad` formulations back
      to back. Read byte-by-byte for TABLO's `!...!` comment delimiters: the first (lines 992-996,
      using `srctwist`/`avesrctwist`/`SIGMADOMDOM(c)`) is live code, closed by the real
      `Substitute xtrad using E_xtrad;` directive. The second (lines 998-1005, the `twistsrc(i,s,k)`
      "alternative form") is entirely inside one unclosed TABLO comment — it opens at line 998's
      lone `!` and doesn't close until the `!` at the very end of line 1005 — so `twistsrc` and the
      second `E_xtrad`/`Variable` declarations are dead documentation, never compiled. Confirms the
      **first** formulation is canonical, and `build_equations.jl`'s old `E_xtrad!`
      (`xtrad-atrad == xuse-(pdelivrd+atrad-puse)`) was missing three things: the `srctwist`/
      `avesrctwist` regional-sourcing-preference shift terms (and their defining equation,
      `E_avesrctwist`, was absent entirely — its variables aren't referenced anywhere else in the
      Julia code), and the `SIGMADOMDOM(c)` CES elasticity coefficient (implicitly using 1.0
      instead). The elasticity itself (`SGDD`, aggregated correctly in `build_premod!.jl`/
      `aggregate_model!.jl`) was *already computed* by `prepare_parameters.jl:26` but never added to
      its output dict — a silent drop of the same kind as the P021/2PUR bugs. **Fixed**: added
      `p["SGDD"] = SGDD` to `prepare_parameters.jl`; added `srctwist[na,ns,nr,nr]` and
      `avesrctwist[na,ns,nr]` `@variable`s to `build_model!.jl`; added `E_avesrctwist!` to
      `build_equations.jl` (mirrors `E_puse!`'s `ID01(DELIVRD_R)*lhs = sum_r DELIVRD*rhs` pattern);
      rewrote `E_xtrad!` to include all three missing terms. `srctwist` is left as a free variable
      with no defining equation — correct, since `TERM.CMF` lists it `Exogenous` in the base closure
      (to be fixed at 0 by Step 5b's `initialize_model!.jl`, same as any other closure-exogenous
      variable). Re-ran full pipeline + full model build: **61.6s, 2,436,062 variables (+59,500 —
      exactly `na*ns*nr*nr` + `na*ns*nr` = 57,800+1,700), 1,370,212 constraints (+1,700, all
      `E_avesrctwist` instances active), 171 vars-dict entries (+2)** — matches hand-calculated
      expectations exactly.
- [ ] **Systemic issue found while fixing the above, not yet addressed**: `prepare_parameters.jl`
      reads `SLAB, P028, SMAR (as SMAR_v), PO01, SCET, P018` from `agg` (lines 24-38) but — like
      `SGDD` before the fix above — never adds any of them to its output dict `p`. They are
      genuinely dead reads (grep confirms zero other uses of these five bindings anywhere in the
      file). `build_model!.jl` compensates with hardcoded placeholder elasticities instead
      (`sigmalab = fill(0.5, na)`, `sigmaprim = fill(0.5, na)`, `sigmaout = fill(0.5, na)`, and
      `sigmadomimp = fill(5.0, na)` per the pre-existing `P015`/`ARMSIGMA` item below). This means
      `E_xlab_o!`/`E_xprim!`/`E_xmake!` (or whichever equations these feed) are running on
      placeholder elasticities even though the real, data-derived values (`SLAB`=labour CES,
      `P028`=primary-factor CES, `SCET`=CET output transformation) are computed correctly upstream
      and simply never wired in. Lower priority than `srctwist` was (this doesn't change the
      variable/constraint *count*, only solved magnitudes on a real shock — benchmark replication
      with all-zero shocks won't catch it either), but should be fixed before trusting shock
      magnitudes. Same root-cause pattern each time: fix by adding `p["SLAB"]=SLAB` etc. to
      `prepare_parameters.jl` and threading them into `build_model!.jl` in place of the hardcoded
      fills — bundle with the `P015`/`ARMSIGMA` fix below since it's the same fix shape.
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

- [x] ~~`initialize_model!.jl` — encode one `.cmf` closure~~ — **done.** Encodes `TERM.CMF`'s
      48-name default `Exogenous` list (42 modeled as `BASE_CLOSURE_SCALARS`/`BASE_CLOSURE_ARRAYS`,
      6 dynamic-only names in `DYNAMIC_ONLY_CLOSURE_NAMES` documented as Step-6-scope and skipped);
      fixes each at a shock value (default 0.0), sets a 0.0 start value on every declared variable,
      leaves everything else free. Deliberately does **not** apply `TERM.CMF`'s numéraire swap
      (`phi = Natmacro("GDPPI")`, needs Excerpt 30-40 reporting, out of Step 5 scope) or dynamic-closure
      swaps (need Step 6) — `phi` stays fixed as the default numéraire instead.
- [x] ~~Verify `initialize_model!` against predicted gap~~ — **done, and this is what surfaced the
      much bigger systemic bug below.** Built a touched/untouched constraint-term introspection
      diagnostic (`test/diagnose_gap.jl`): walks every affine constraint, collects every `VariableRef`
      it references, and reports any *free* (non-fixed) variable that never appears in any constraint
      at all — a variable with literally no defining or using equation. First run found **13,877
      untouched free variables**, breaking down by base name as `wlab_o=>850, xsuppmar_d=>10404,
      psuppmar_p=>1258, xsuppmar_rd=>306, xtrad_d=>162, fgret=>850, fhou2=>34, plab_id=>4, wlab_id=>4,
      rlab_id=>4, natfhou=>1` — none of this was expected from the vars-minus-constraints gap alone,
      so each entry needed root-causing individually (see the two bugs below).
- [x] ~~Fix `LAB_O`=0 bug~~ — **done, root cause was a leftover stub, not a genuine data gap.**
      `wlab_o=>850` in the untouched list (`E_wlab_o!`/`E_plab_o!`/`E_wprim!` all guard on
      `LAB_O[i,d] > 1e-10`) traced to `params["LAB_O"]` being identically zero. Traced upstream to
      `build_pstras!.jl`'s `V1LAB_iod` computation, which had literally been left as
      `for o in 1:NO; V1LAB_iod[i,o,d] = 0.0; end` — a "We'll recompute below" placeholder that was
      never filled in. The correct value needs `reg0`'s national labour-occupation shares (`OSHR`,
      IND×OCC) applied to `reg1`'s regional labour factor payment (`FAC_r[i,g_lab,d]`), but `OSHR`
      wasn't in `reg1`'s output dict at all. **Fixed**: added `"OSHR" => reg0["OSHR"]` to
      `build_reg1!.jl`'s output dict; replaced the zero-stub in `build_pstras!.jl` with
      `V1LAB_iod[i,o,d] = FAC_r[i,g_lab,d] * OSHR[i,o]`. Re-ran `test/diagnose_gap.jl`: untouched count
      dropped to 13,027 (−850, `wlab_o` fully gone from the breakdown), all other counts unchanged —
      confirms no side effects.
- [x] ~~Investigate remaining untouched-variable gap after the LAB_O fix~~ — **done, found the single
      most severe bug in the project so far: 10 of 11 lookup-dict `*_setup!` functions in
      `build_equations.jl` were defined and exported but never called anywhere.** Many equations
      depend on a module-level `const XXX_idx = Dict(...)` populated by a dedicated `*_setup!()`
      function that must run before the equation reads it via `get(dict, key, 0.0)`. Confirmed by grep
      that `build_model!.jl` (the only place any equation-builder is invoked) called just one of these
      11 setup functions (`E_delXGDPEXPa_setup!`, which is differently-shaped and doesn't count) —
      `E_plab_o_setup!`, `INVEST_setup!`, `USE_IS_setup!`, `USE_usc_setup!`, `TRADMAR_setup!`,
      `SUPPMAR_setup!`, `SUPPMAR_D_setup!`, `TRADE_setup!`, `PUR_src_setup!`, `TAX_PUR_setup!`,
      `STOCKS_setup!` were all dead code, meaning every equation reading their dicts
      (`E_plab_o!`, `E_wlab_o!`, `E_wprim!`, `E_pinvitot!`, `E_xint_i!`, `E_xuse!`, `E_xsuppmar_p!`,
      `E_psuppmar_p!`, `E_xsuppmar_d!`, `E_xsuppmar_rd!`, `E_xtrad_d!`, `E_xtrad_r!`, `E_xfind!`,
      `E_delTAXint!`, `E_delTAXhou!`, `E_delTAXinv!`, `E_delXGDPEXPb!`, `E_delPGDPEXPb!`) had been
      silently summing over an empty dict (`get(...,0.0)` always falling back to zero) — spanning
      labour prices, margin supply, market clearing, investment allocation, tax revenue, and GDP
      decomposition. Worse, this interacted directly with the LAB_O fix just above: making `LAB_O`
      nonzero let `E_plab_o!`/`E_wlab_o!`/`E_wprim!`'s guards start passing, but since `V1LAB_idx` was
      still empty, those equations began wrongly constraining `plab_o`/`wlab_o`/`wprim` to exactly 0
      instead of just being absent — "touched but wrongly valued," invisible to the touched/untouched
      diagnostic. Also found, while wiring the fix: `E_xinv_s!` (defines `xinv_s`, consumed by
      `E_xfinb!`) was a wholesale missing equation-block call, a different bug shape (no call at all,
      vs. a missing setup call). **Fixed**: read `prepare_parameters.jl` in full to identify the
      correct pre-computed source array for each dict (`params["USE"]`, `params["TAX"]`,
      `params["PUR"]`, `params["INVEST"]`, `params["STOCKS"]`, plus `TMAR`/`MARS` newly extracted from
      `agg`), added all 8 remaining setup!() calls (`E_plab_o_setup!` was one of these 8) to
      `build_model_full!` in `build_model!.jl` right before the Excerpt 6 block, guarded with
      `haskey(params, ...)` where the source is itself conditionally populated; added the missing
      `E_xinv_s!(m, vars, na, nr, params)` call to the Excerpt 14 block.
- [x] ~~Fix `STOK` (stock-change) pass-through gap, found during the same investigation~~ — **done.**
      `params["STOCKS"]` (feeds `STOCKS_idx`, `GDPEXPSUM`, `CHECKA`/`CKRATA`,
      `E_delXGDPEXPb!`/`E_delPGDPEXPb!`) was silently defaulting to `zeros(T,na,nr)` via
      `prepare_parameters.jl`'s `agg["STOK"] !== nothing` guard — even though `reg1["STOK"]`/
      `reg2["STOK"]` (regional industry stock-change data) existed correctly upstream. The one gap:
      `build_pstras!.jl`'s output dict never included a `"STOK"` key at all — same silent-drop shape
      as the earlier `2PUR`/`UTAX` bugs. **Fixed** by adding
      `"STOK" => haskey(reg1, "STOK") ? reg1["STOK"] : nothing,` to `build_pstras!.jl`'s output dict,
      matching the existing `UTAX`/`2PUR` pass-through lines.
- [x] ~~Verify all of the above together~~ — **done 2026-07-22.** Re-ran `test/diagnose_gap.jl`:
      untouched free vars dropped from **13,027 → 3,575** (−9,452). `xsuppmar_rd` (was 306) is
      completely gone; `xsuppmar_d` fell from 10,404 → 1,258 (now matching `psuppmar_p`'s count
      exactly, consistent with both sharing the same sparse margin-flow index set). Remaining
      breakdown: `xsuppmar_d=>1258, psuppmar_p=>1258, fgret=>850, xtrad_d=>162, fhou2=>34, plab_id=>4,
      wlab_id=>4, rlab_id=>4, natfhou=>1` — not yet individually root-caused, but each looks like
      genuine data sparsity (zero-valued margin/labour/region combinations in the real 25×34 dataset)
      rather than a new missing-wiring bug, since none of these dropped or changed shape when the
      setup calls were added. Also wrote `test/check_setup_dicts.jl` to numerically confirm every one
      of the 12 lookup dicts is now non-empty with a real nonzero sum after `build_model_full!` (e.g.
      `V1LAB_idx` sum 6.43e6, previously 0.0 for all 3,400 entries), and that `agg["STOK"]` and
      `params["STOCKS"]` now carry matching nonzero sums (33,956.96) end-to-end — confirms `plab_o`/
      `wlab_o`/`wprim` are no longer wrongly forced to zero (the touched/untouched diagnostic alone
      can't detect "touched but wrongly valued," so this numeric check was the necessary complement).
- [ ] Root-cause the remaining 3,575 untouched vars (lower priority — likely genuine data sparsity,
      not a wiring bug, but not yet confirmed row-by-row).
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
      **Same bug now confirmed for 5 siblings** (see Step 5 section above,
      2026-07-22): `SLAB, P028, SMAR, PO01, SCET, P018` are all read into local bindings in
      `prepare_parameters.jl` and then never added to its output dict, so `build_model!.jl` falls
      back to hardcoded placeholder elasticities for all of them. Fix all six in one pass.
- [ ] Once `HeaderArrayFile.jl` is confirmed unused elsewhere, drop it from `Project.toml` (currently
      kept only because Step 1 originally depended on it before the harpy-CSV fallback).
