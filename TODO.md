# IndotermJulia — TODO

> Working task list companion to `PLAN.md` (which holds the full design rationale, learnings, and
> status history). This file is the short, checkable list of what's next. Update it as work lands;
> keep the narrative/rationale in `PLAN.md`.

## Housekeeping (do first)

- [x] ~~**Commit or stash Step 4/5 work.**~~ — committed 2026-07-22 as `e5381a9` ("Add Step 4/5 core
      translation: derived parameters and CGE equations"), covering `build_equations.jl`,
      `build_model!.jl`, `prepare_parameters.jl`, the `PLAN.md`/`TODO.md` updates, the P021/FRISCH
      fix, and `test/run_full_model.jl`. Not yet pushed to `origin/master` — ask before pushing.

## CURRENT FOCUS (2026-07-25) — read this first

Everything below this section is historical unless linked from here. Full rationale in `PLAN.md`
("Gate #5c RESOLVED").

**Where we are: the levels model is VERIFIED CORRECT and solves at 25×6.**

- Benchmark relative residual **1.86e-9**; Newton converges in **1 iteration / 2.2 s**.
- Equilibrated Jacobian is **full rank** (was deficient by 637).
- Solver changed **Ipopt → sparse Newton** (Ipopt hangs at `iter 0` on a square system — zero
  degrees of freedom is degenerate for an interior-point optimizer).
- The model is now **naturally square from `src/` alone**: squaring collapsed from
  "pin 157 + delete 108 dead rows" to **pin 1, delete 0** — and that one pin is `natfhou`, the
  stub already documented under gate 5b.

**Next, in order:**

- [x] ~~**Gate 5c-i — productionise the solve.**~~ — **done 2026-07-25.** New
      `src/solve_newton!.jl` (`solve_newton!` → `NewtonResult`), exported; `Project.toml` gains
      `SparseArrays`/`LinearAlgebra`; `test/solve_benchmark_6reg.jl` rewritten to use it and now
      passes unaided: square 88,599=88,599, converged in 1 iter (42.8 s), scaled `‖F‖∞` 1.11e-9,
      max relative Δ vs benchmark 5.4e-7. Squaring is defensive and automatic (pins orphan vars —
      today just `natfhou`; deletes dead rows — today none) and errors out with a diagnostic rather
      than solving a non-square system silently.
- [ ] **First real economic test** — run `TERM.CMF`'s `Shock blabnat = -3` and sanity-check signs and
      magnitudes. Unblocked now that the null space is gone (shock results would previously have been
      non-unique). **This is the active task.** Note `initialize_model!` already accepts `shocks`.
- [ ] **Price homogeneity test** — shock only the numéraire; nominal variables move, real ones must
      not. This is the natural check that the remaining closure/numéraire wiring is right.
- [ ] **Re-test at full 34 regions.** Needs Phase 3 condensation first: 2.4M variables, and the
      sparse QR fallback will not scale as-is. Profile whether `lu` alone suffices now that the
      Jacobian is full rank — it may, which would make Phase 3 cheaper than planned.
- [ ] Fold the `OMIT_NAMES` / zero-flow-guard patterns into a short "porting TABLO faithfully" note
      in `AGENTS.md` if more instances turn up (two classes found so far, both now fixed).

**Do not re-litigate:** the *solver* choice (Ipopt hung 56 min at `iter 0`, measured), or the
correctness of the equations *that are present*.

### Excerpt coverage audit (2026-07-25) — what is and is not translated

`TERM.TAB` has **55 excerpts**. Audited each against `src/`:

| group | excerpts | status |
|---|---|---|
| Sets, data reads, sign checks, flow updates | 1, 2, 3, 5 | ✅ handled by the pipeline (`prepare_sets`/`read_data`/`build_*`) |
| Core model | 4, 6-29 | ✅ translated (levels) |
| Regional macro reporting — **partial** | 31 | ⚠️ parameters only (`TRADE_CR`, `IMPUSED_C`, `IMPLANDED_C`); `ximps`/`ximpused`/`pimpused`/`ximplanded`/`pimplanded` and `MainMacro`/`SelMacro` were **missing** |
| National macro aggregation | 32 | ⚠️ **missing** (`NatMacro`, `contMainMacro`, `shrBoT`, `shrBoTnom`, `shrBoTnom2`, `finvgdp`) |
| Labour/household closure, GDP | 35, 38-41, 44 | ✅ translated |
| Condensation actions | 49 | n/a — a Phase 3 *lever*, not equations |
| Dynamic extension | 50, 51, 53, 54 | ✅ translated (`build_dynamics!.jl`) |
| Dynamic diagnostics summary tables | 52 | ⏭️ diagnostics only, deferred |
| District extension | 55 | ⛔ N/A — no data (see Step 6 section) |
| **Contribution/decomposition reporting** | 30, 33, 34, 36, 37, 42, 43, 45, 46, 48 | ⏭️ **not translated** (~55 equations) |

- [x] ~~**Excerpts 31-32 — `MainMacro`/`NatMacro`.**~~ — `src/build_macros!.jl`, **verified
      2026-07-25**: square 90,460 = 90,460, converges in 1 Newton iteration, benchmark replication
      unchanged (max relative Δ 5.42e-7). Note the *raw* residual rose 9.2e-5 → 0.146 while the
      *scaled* one stayed 1.106e-9 — the macro equations carry ~1e7 weights, which is exactly why
      convergence is judged on the equilibrated residual.
      This is *not* cosmetic: `NatMacro("GDPPI")` is the target of `TERM.CMF`'s **active numeraire
      swap** (`swap phi = Natmacro("GDPPI")`, `TERM.CMF:72`), which the port has been unable to apply
      — `initialize_model!.jl` says so in its own docstring. Implementing it unblocks using the
      model's real numeraire instead of leaving `phi` fixed.
- [ ] **Apply the real numeraire swap** now that `NatMacro("GDPPI")` exists: free `phi`, fix
      `NatMacro("GDPPI")`. Re-verify benchmark replication and then run the price-homogeneity test,
      which is the check this swap actually matters for.
- [ ] Excerpts 30/33/34/36/37/42/43/45/46/48 — contribution & decomposition reporting (~55 equations:
      commodity contributions to national results, endowment/tech/tax contributions to real GDP,
      COM/REG contributions to trade price indices, inter-regional trade reporting). Pure reporting —
      none feeds back into the core solve, so they are safe to add incrementally. Excerpt 46 is
      weights-only (0 equations).

### Limits of the benchmark-replication proof — keep in mind for shock results

Benchmark replication proves the equations **that exist** are internally consistent and correctly
calibrated. It **cannot** detect an equation that is *economically wrong but benchmark-consistent*,
because at a balanced benchmark SAM almost any plausible market-clearing form holds exactly. Such a
bug only surfaces **under a shock**. Keep this in mind when the first shock results land.

Two instances were flagged in earlier sessions. **Both were re-verified against `TERM.TAB` on
2026-07-25 and are ALREADY FIXED** — the notes claiming otherwise (further down this file) are
stale and are retained only as history:

- [x] ~~`E_pdomB` missing~~ — **implemented** in `E_pdomA_sum!` (`build_equations.jl`). Margin
      commodities use the levels form `xcom == xtrad_d + Σ xsuppmar_rd + K` with
      `K = MAKE_I − TRADE_D − Σ SUPPMAR_RD` carrying the benchmark imbalance. Verified: at the
      benchmark RHS collapses to `MAKE_I` (replicates), and differentiating gives
      `MAKE_I·x̂com = TRADE_D·x̂trad_d + SUPPMAR_RD·x̂suppmar_rd` — exactly `TERM.TAB:1200-1203`.
      Non-margins use the ratio form `xcom/MAKE_I == xtrad_d/TRADE_D`, which linearizes to
      `x̂com = x̂trad_d` (`TERM.TAB:1196-1198`). The MAR→aggregated-COM mapping lands the 9 margins
      on commodities 20/21/22 (Trade / Transport / InfoComm).
- [x] ~~`xuse` purchaser-vs-basic wedge~~ — **resolved** in `E_xuse!`. Each user's
      purchaser-valued quantity is converted to basic values by a benchmark ratio
      (`β_u = USE_u / PUR_u`) before summing, and `xuse` enters scaled by
      `α = USE_U / TRADE_R`: `α·xuse == Σ_u β_u·x_u`. Replicates because
      `USE_U = USE_I + Σ USE_final`, and `xuse` keeps its `TRADE_R` (basic-value) benchmark, which is
      what the `E_xtrad` sourcing nest expects.

## Step 5 completion (core equations) — historical

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

## Step 5b — closure + solve wiring (historical; solver since changed Ipopt → Newton)

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
- [x] ~~Fix missing `E_ggro!`/`E_fgret!` equations (Excerpt 15 investment rule).~~ — **done
      2026-07-22.** `fgret=>850` in the untouched breakdown exactly equaled `na*nr` (25×34), a red
      flag for a wholesale missing equation rather than data sparsity. Grep confirmed `fgret` had
      zero defining equations anywhere. Read `TERM.TAB` lines 738-777: `ggro` also lacked its own
      equation — it only appeared as a term inside `E_xinvitot!`'s constraint, so it was "touched"
      but underdetermined by `na*nr` equations (same masking pattern as the `xlab_id`/`realwage_id`
      case below). **Fixed**: added `E_ggro!` (`ggro == finv1 + 0.33*(2*gret - invslack)`) and
      `E_fgret!` (`gret == fgret + capslack`) to `build_equations.jl`, wired both into
      `build_model_full!`'s Excerpt 15 block in `build_model!.jl` right after `E_xinvitot!`. All
      referenced variables (`gret, ggro, finv1, invslack, fgret, capslack`) already existed as JuMP
      variables — no new declarations needed. Re-ran `test/diagnose_gap.jl`: untouched dropped
      3,575 → 2,725 (−850, exactly `na*nr`), `fgret` completely gone from the breakdown, no other
      counts changed.
- [x] ~~Fix missing "_id" labour-aggregate family (Excerpt 27: `E_plab_id!`, `E_realwage_id!`,
      `E_xlab_id!`, `E_wlab_id!`, `E_rlab_id!`).~~ — **done 2026-07-22.** `plab_id=>4, wlab_id=>4,
      rlab_id=>4` in the untouched breakdown each exactly equaled `no` (4 labour occupations) — same
      100%-of-array red flag as `fgret` above. Grep confirmed all 5 equations in the family were
      missing; `xlab_id`/`realwage_id` didn't show up as untouched only because they're *consumed*
      (not defined) by an unrelated existing constraint at `build_equations.jl:1173-1174` — the same
      "touched but not well-defined" masking as `ggro`. Also found the required coefficient
      (`SLAB_ID`/`LAB_ID`, the region-summed labour aggregate and its share) had never been computed
      in `prepare_parameters.jl` at all — only the sibling `LAB_I`/`SLAB_I` (per-industry) and
      `LAB_IO` (per-region) existed. **Fixed**: added `p["LAB_ID"]`/`p["SLAB_ID"]` computation to
      `prepare_parameters.jl` right after the existing `SLAB_I_arr` block (`LAB_ID[o] = sum_d
      LAB_I[o,d]`, `SLAB_ID[o,d] = LAB_I[o,d]/LAB_ID[o]`); added the 5 equation functions to
      `build_equations.jl` mirroring the sibling `_i`-family functions but summing over region `d`
      with the `SLAB_ID` share weight instead of summing over industry `i`; wired all 5 calls into
      `build_model_full!` between the `_i` and `_io` family calls in `build_model!.jl`. Re-ran
      `test/diagnose_gap.jl`: untouched dropped 2,725 → 2,713 (−12, exactly `3*no`), `plab_id`/
      `wlab_id`/`rlab_id` completely gone from the breakdown.
- [x] ~~Root-cause the remaining 2,713 untouched vars.~~ — **done 2026-07-22.** Checked each
      remaining entry individually:
      - `fhou2=>34` matched `nr` exactly (100%-of-array red flag, same heuristic as `fgret`/
        `plab_id`). Grep confirmed `build_equations.jl` had `E_fhou!` and `E_natfhou!` but no
        `E_fhou2!` at all, even though `build_model!.jl:270` declares `fhou2[1:nr]` as a variable.
        `TERM.TAB` lines 2038-2053 (Excerpt 39) show `E_fhou`/`E_fhou2` as a genuine TABLO chain: both
        equations share `whouhtot(h,d)` on the LHS, but `whouhtot` is already pinned by the
        pre-existing `E_whouhtot!` (`whouhtot == phouhtot + xhoutot`), so in practice `E_fhou!` solves
        for `fhou` and `E_fhou2!` solves for `fhou2` — a "propensity to consume from regional GDP"
        ratio, mirroring `fhou`'s "propensity to consume from labour income." **Fixed**: added
        `E_fhou2!` (`whouhtot[d] == wgdpexp[d] + fhou2[d] + houslack`, mirroring `E_fhou!`'s exact
        shape with `wgdpexp` in place of `wlab_io`) to `build_equations.jl`, wired the call into
        `build_model_full!` in `build_model!.jl` right after `E_fhou!`. Re-ran `diagnose_gap.jl`:
        untouched dropped 2,713 → 2,679 (−34, exactly `nr`), `fhou2` completely gone.
      - `natfhou=>1` is a **documented, intentional gap, not a bug**: `E_natfhou!` exists but its
        body is empty (`function E_natfhou!(...) end`) and is still called — a deliberate stub.
        `TERM.TAB` line 2053 requires `NatMacro("NomHou")`/`NatMacro("NomGDPexp")`, national aggregates
        computed by the Excerpt 30-40 `MainMacro`/`NatMacro` reporting layer already established as
        out of Step 5's core scope (same category as the un-ported `CHECKC`/`D`/`E` diagnostics —
        "one-way reads *from* the core solution," nothing downstream depends on `natfhou`). Left as
        is; revisit only if/when the Excerpt 30-40 reporting layer is built.
      - `xtrad_d=>162` — confirmed genuine data sparsity, not a bug. `E_xtrad_d!` guards on
        `TRADE_D[c,s,r] > 1e-10`; total array size is `na*ns*nr`=1,700, so 162 zero-valued
        commodity/source/region combinations (real gaps in the 25×34 trade data) is a small fraction,
        not the 100%-of-array pattern that flagged every real bug this session.
      - `xsuppmar_d=>1258, psuppmar_p=>1258` — left as genuine data sparsity per the prior segment's
        conclusion (both share the same sparse margin-flow index set, guarded the same way as
        `xtrad_d`); not re-verified row-by-row this round since neither changed shape across three
        consecutive unrelated fixes, which would be surprising for a live wiring bug.
      Final state: **2,679 untouched free vars**, of which only `natfhou=>1` is a known,
      documented, deliberately-out-of-scope gap; everything else is confirmed real data sparsity.
- [x] ~~`solve_model!.jl` — Ipopt feasibility solve (no objective), mirroring WayangJulia's
      `set_attribute` tuning (`max_iter`, `tol`, `constr_viol_tol`, adaptive `mu_strategy`).~~ —
      written 2026-07-22, wired into `IndotermJulia.jl`; not yet committed. Superseded in priority by
      the levels-conversion course correction below — revisit once `build_equations.jl` is in levels
      form (the %-change-era `solve_model!.jl` itself needs no change, just a levels model to call it on).
- [x] ~~Fix `initialize_model!.jl` closure-benchmark bug (all closure vars fixed at 0.0, wrong for
      ratio-type vars).~~ — **done 2026-07-23.** `test/check_start_point.jl` (new diagnostic: evaluates
      every constraint at the fixed/start initial point, no solver call) found 124,100/1,396,876
      constraints evaluating to `Inf`, concentrated in `log(...)` reads of `tuser`/`ppur`/`puse`/
      `aint_s`/`bint_scd`/`bint_s`. Root cause: 17 of 41 `TERM.CMF` closure names are ratio-type
      (`>=1e-6`-bounded, correct benchmark `1.0`) but `initialize_model!.jl` fixed all 41 at `0.0`
      uniformly — right for shift-type names, wrong for ratio-type ones (`log(0)=-Inf`). **Fixed**:
      added `_is_ratio_type` helper (same `>=1e-6` convention `solve_benchmark.jl` already used for
      free vars); closure fix now defaults to `1.0` for ratio-type, `0.0` for shift-type. Confirmed via
      before/after `check_start_point.jl` runs: 124,100 → 0 non-finite constraints. Real Ipopt solve no
      longer terminates `INVALID_MODEL` either.
- [ ] **New blocker, found 2026-07-23, paused per user direction — do not resume without being asked**:
      with the `Inf`-constraint bug fixed, `test/solve_benchmark.jl` now runs real Newton iterations but
      fails with `MUMPS returned INFO(1) =-13 - out of memory` (5,288–38,784 MB depending on
      `mumps_mem_percent`), ending `OTHER_ERROR`/`INFEASIBLE_POINT`. Confirmed NOT a genuine memory
      shortage (host has ~16.68GB free; a raw Julia allocation test succeeded to 20GB cleanly).
      Diagnosed as most likely a MUMPS-internal scaling limit (probably 32-bit integer indexing) for a
      KKT system this large (2.4M vars / 1.4M cons). Only `Ipopt v1.15.0`+bundled MUMPS is installed —
      no `HSL_jll`/Pardiso. Full writeup in `PLAN.md`'s "Gate #5c status detail" section. Options once
      resumed: (a) install HSL (MA57), needs external STFC academic registration; (b) investigate/reduce
      KKT fill-in; (c) stay paused. **User chose to pause on 2026-07-23** — awaiting further direction.
- [ ] Recalibrate `solve_model!.jl`'s `print_level` default — documented as avoiding Ipopt's
      full-vector per-iteration dumps at level 8, but empirically level 8 still dumps them
      (`new vars[i]`/`curr_c[i]`/`curr_x[1][i]`/`final y_c/y_d/z_L/z_U` for all 2.4M variables),
      producing gigabyte-scale log files. Only `print_level>=11` was assumed to do this; not true at
      this model's scale. Low urgency, fold in whenever solver work resumes.
- [ ] `run_model!.jl` — homotopy/warm-start driver, only if a direct solve fails to converge on a
      real shock.

## Course correction (2026-07-22): levels conversion — current top priority

User caught that `build_equations.jl` had silently drifted from the documented **levels formulation**
(`PLAN.md` line 14) to a direct port of TABLO's own **%-change (Johansen)** equations (`PLAN.md` line
417) — confirmed via the touched/untouched diagnostic (100% of the model's constraints are
`GenericAffExpr`, impossible for a genuine nonlinear levels CGE) and by comparing against WayangJulia's
actual source, which *did* execute the levels conversion the plan called for. Full writeup in
`PLAN.md`'s "Course correction" section. **Decision: convert to true levels, following WayangJulia's
methodology directly** (`ces()`/`_ces_calibrate()` helper, Principle A/B/C, log-differential auxiliary
variables). This is now the primary line of work; everything below in this section works toward it.

> **✅ THIS WHOLE SECTION IS COMPLETE (verified 2026-07-25).** The model is genuinely in levels form:
> `src/ces_helper.jl` provides `ces`/`ces_calibrate`, the equations use nonlinear `ces(...)`/`log(...)`
> constructions throughout, and the converted model reproduces its own benchmark to **1.86e-9**
> relative — the real (non-trivial) replication test this conversion was undertaken to make possible.
> Boxes below ticked retrospectively; kept for the rationale.

- [x] ~~Port `ces(y, p, α, σ, γ)` and `_ces_calibrate(quantities, sigma, output)` from WayangJulia~~
      — done, `src/ces_helper.jl:14` (`ces`) and `:38` (`ces_calibrate`).
- [x] ~~Fold in the pre-existing elasticity-wiring gap~~ (`SLAB, P028, SMAR, PO01, SCET,
      P018`/`P015`/`ARMSIGMA`) — done; these are wired as the `σ` inputs to the CES calibration
      (e.g. `ALPHA_FAC`/`GAMMA_FAC` from `P028` in `prepare_parameters.jl`).
- [x] ~~Convert Excerpts 6-12 (basic/purchaser prices, Armington CES, factor CES, zero-profit
      costs) to levels form.~~
- [x] ~~Convert Excerpts 13-17 (household LES demand, investment, government, export CET) to levels
      form.~~
- [x] ~~Convert Excerpts 19-23 (margins, regional-sourcing CES, MAKE/CET, market clearing) to levels
      form.~~ — **done 2026-07-23.** All 18 functions rewritten against `origin/TERM.TAB` lines
      895-1219 directly (not just the old %-change Julia port). Highlights: `E_xsuppmar_p!`,
      `E_xsuppmar_d!`, `E_xsuppmar_rd!`, `E_xtrad_d!`, `E_xtrad_r!`, `E_xcomA_B!` all collapse to
      plain physical-quantity sums (their fixed weights equal the summed variable's own benchmark
      exactly — confirmed against TABLO, **not** `MARS*DIST` as the old `SUPPMAR_idx` Dict assumed);
      the now-dead `TRADMAR_idx`/`SUPPMAR_idx`/`SUPPMAR_D_idx`/`TRADE_idx` Dict-caching machinery and
      their `*_setup!` calls/exports were deleted. `E_xmake!` became a genuine `ces()`-based CET nest
      (new `ALPHA_MAKE`/`GAMMA_MAKE` calibration in `prepare_parameters.jl`, negative-`SCET` CET dual
      convention); `E_xtotA_B!` is its exact Shephard's-lemma revenue-identity dual. **Bug fix along
      the way**: `E_xsuppmar!` had completely dropped TERM.TAB's `SIGMAMAR` CES elasticity term (was a
      pure Leontief `xsuppmar==xsuppmar_p+asuppmar`) — restored per the real TABLO equation. Added
      `xtradmar >= 0` bound (missed in the original bound audit — needed now that
      `E_xtradmar_na!` reads `log(xtradmar)`), plus the other 9 bounds already flagged
      (`xtrad`/`pdelivrd`/`xcom`/`xmake`/`pmake`/`xsuppmar`/`xsuppmar_p`/`psuppmar_p`/`xsuppmar_d`/
      `xsuppmar_rd`). **Known gap intentionally NOT fixed (pre-existing, out of scope for this
      conversion pass)**: TABLO splits market clearing into `E_pdomA` (NONMAR commodities only) and
      `E_pdomB` (MAR/margin commodities: `MAKE_I*xcom == TRADE_D*xtrad_d + SUPPMAR_RD*xsuppmar_rd`,
      TERM.TAB lines 1196-1203) — the Julia port only ever had `E_pdomA_sum!`, applied uniformly to
      *all* commodities including margin ones. `E_pdomB` is still entirely missing. [**STALE — this
      was FIXED in a later session; `E_pdomA_sum!` now implements both split forms. Verified against
      `TERM.TAB:1196-1203` on 2026-07-25. See "Limits of the benchmark-replication proof" at the top
      of this file.**] This predates the
      levels conversion (it's a pre-existing gap in the original %-change port, not something this
      pass rewrote from scratch), so left as-is per the "big-bang, one pass" mandate; needs a real fix
      before trusting margin-commodity market-clearing results.
- [x] ~~Convert Excerpts 24-29 (final-demand aggregates, tax revenue, GDP income/expenditure) to
      levels form.~~
- [x] ~~Convert Excerpts 27/38-39 (labour-market closure, household closure) to levels form.~~
- [x] ~~Update `initialize_model!.jl` (start values / closure fixes need to be benchmark levels, not
      `0.0`)~~ — done via `src/benchmark_levels.jl` + the `bmk_levels` argument to
      `initialize_model!`, which seeds each genuine flow variable at its benchmark value-flow
      (index/ratio variables stay 1.0, additive shifters 0.0).
- [x] ~~Re-run the verification protocol on the converted model~~ — done 2026-07-25: benchmark
      replication holds to **1.86e-9** relative. See "CURRENT FOCUS" at the top.

## Verification (per PLAN.md's protocol — do in this order)

- [x] ~~**Benchmark replication**: solve with all shocks at 0 → expect ~0% change everywhere.~~ —
      **passed 2026-07-25 at 25×6.** Relative residual 1.86e-9; Newton converges in 1 iter / 2.2 s;
      max drift off benchmark 5.4e-7. Still to repeat at full 34 regions. See the ⚠️ caveat under
      "CURRENT FOCUS" — this test cannot detect a benchmark-consistent but economically wrong
      equation (`E_pdomB`).
- [ ] **Price homogeneity test**: shock only the numéraire → nominal moves, real variables don't.
- [ ] **GDP-both-sides check**: after any shock, income-side GDP == expenditure-side GDP
      (`calculate_gdp.jl`, Step 8 — can be stubbed early just for this check).
- [ ] **Smoke-test scenario**: run one named `.cmf` shock (`TERM.CMF` or `sim1.cmf`) and sanity-check
      signs/magnitudes.

## Step 6 — Dynamic + district extensions (after Step 5 verified)

- [x] ~~`build_dynamics!.jl` — Excerpts 50-54 (investment rule, real-wage adjustment).~~ —
      **done 2026-07-25.** `src/build_dynamics!.jl` translates Excerpt 50 (capital accumulation),
      51 (investment rule), 53 (national dynamic reporting) and 54 (real-wage adjustment) to levels,
      plus `update_dynamics!` for the inter-period `Update` statements a recursive multi-period
      driver needs. Verified: system stays square (89,911 = 89,911), converges in 1 iteration, and
      benchmark replication is unchanged (max relative Δ 5.42e-7, identical to the pre-dynamics run).
      Under the base static closure the block is a passive satellite — `faccum`/`finv4`/`delfwage`
      absorb — exactly as TERM.TAB's own "to switch off" note prescribes.
- [ ] **`build_district!.jl` — Excerpt 55: NOT APPLICABLE to this dataset.** The top-down district
      extension requires a sub-provincial layer that INDOTERM simply does not have. Excerpt 55 reads
      `RGN` (district set), `MRGN` (district→region mapping), `SREG` (split flag), `LOCI` (local-industry
      flag) and `MVTO` (district outputs/value-added); **none of those five headers exists** in any of
      `national.har`, `regsupp.har` or `DISTGONE.HAR` (66 headers total, checked 2026-07-25). The
      finest spatial unit in the data is the 34 provinces = `REG`. Writing the module now would give
      an untestable, unrunnable stub, so it is deliberately deferred rather than faked.
      **To enable later**: obtain district-level data providing those five headers, then implement it
      as a non-feedback satellite computed *after* the core solve (regional industries grow at the
      region rate; local/"municipal" industries follow local demand — see TERM.TAB:2921-3010).

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
- [x] ~~Once `HeaderArrayFile.jl` is confirmed unused elsewhere, drop it from `Project.toml`~~ — **done**: confirmed no `using`/`import`/qualified call in `src/`, `test/` or `analysis/` (only a docstring in `read_data.jl` explaining why it is NOT used); removed from `[deps]` and `[compat]`. (was: currently
      kept only because Step 1 originally depended on it before the harpy-CSV fallback).
