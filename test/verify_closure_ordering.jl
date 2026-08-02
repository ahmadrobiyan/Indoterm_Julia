"""
Closure ordering (VV_PLAN.md V6).

Standard CGE theory: the same shock produces a LARGER real output response under
a long-run closure (capital mobile across sectors) than under a short-run closure
(capital fixed in place). This project already maintains two closures and has
been bitten once by quoting a number without naming which one — `PLAN.md` §N
found the same cell reading 155× differently between the base static closure and
TERM.CMF's swaps, purely from which variables were free. This gate turns that
lesson into a standing check: run the same complete shock under both, and verify
the ordering theory demands actually holds.

**Method** (VV_PLAN.md §V6). Run `TERM_CMF_REFERENCE`'s shocks (`blabnat`,
`delUnity`) under two closures:
  * long-run — `TERM_LR_SWAPS_MATCHED` (`closures.jl`), NOT `TERM_CMF_SWAPS`.
    2026-07-31: an earlier version of this gate used `TERM_CMF_SWAPS` directly
    and got an INCONCLUSIVE verdict at every shock magnitude tried (VV_PLAN.md
    §V6 RESULT) — `TERM_SR_SWAPS`'s branch folded at λ*≈0.9998 regardless of
    target size. The official-documentation check traced this to a mismatched
    pair: `TERM_CMF_SWAPS` uses `TERM.CMF`'s own labour swap
    (`delfwage=flabsup_id`, a third variant that is neither textbook long-run
    nor short-run per Excerpt 38), while `TERM_SR_SWAPS` uses the textbook
    short-run labour swap. `TERM_LR_SWAPS_MATCHED` is identical to
    `TERM_CMF_SWAPS` except it swaps in `termlr.cmf`'s own textbook long-run
    labour recipe (`flabsup_id = xlab_id`) instead, so both legs are now
    matched on labour AND capital, as Excerpt 38 actually defines the
    short-run/long-run distinction. `xcap` is still swapped endogenous
    (`xcap = faccum`): capital can move.
  * short-run — `TERM_SR_SWAPS` (`closures.jl:59`), termsr.cmf's short-run swaps.
    `xcap`/`xlnd` stay exogenous: capital is welded to its benchmark allocation.

**Why the short-run run is STAGED, not a straight homotopy on both shocks
together.** `closures.jl`'s own docstring warns `TERM_SR_SWAPS` was written for
termsr.cmf's own scenario (a targeted investment shock), not an economy-wide
labour shock. A first attempt at this gate ran both shocks together through
plain `t∈[0,1]` homotopy under `TERM_SR_SWAPS` and it did not fail cleanly — it
ran for 35+ minutes burning full CPU with zero step output, i.e. stuck deep
inside a single Newton solve at one homotopy point, not obviously converging or
diverging. That is exactly the natural-parameter-continuation failure mode
`arclength.jl`'s header describes for the long-run branch's own fold (§N.19c):
a plain homotopy walk gives no signal about WHY a point is hard, it just grinds.
So this run follows the exact two-stage recipe §N.19e used for
`TERM_CMF_REFERENCE` itself:
  A. `delUnity` 0 → 1 alone, by natural continuation (`continuation_solve!`) —
     cheap, and worth seeing on its own.
  B. from there, `blabnat` 1 → 0.97 by pseudo-arclength (`arclength_solve!`),
     which can locate and pass a fold instead of stalling on one.
Staging does not change the scenario or the closure — `TERM_SR_SWAPS` and both
shock targets are identical to a direct run. It only changes HOW the solver gets
there, and only for the closure that struggled; the long-run run is unaffected.

Compare national RealGDP (`NatMacro(RealGDP)`, `build_macros!.jl`'s
`"RealGDP" => d -> xgdpexp[d]`) and sectoral output (`xtot[i,d]`, benchmark = 1,
`build_equations.jl:453-459`) between the two runs' final points.

**Pass.** |ΔRealGDP| is larger under the long-run closure than the short-run
closure (capital immobility damps the response), AND a majority of the 25×6
`xtot` cells show the same ordering. Sector-level exceptions are legitimate
(relative-price effects can push an individual sector's short-run response above
its long-run one even while the aggregate orders correctly) and are listed, not
suppressed or treated as failures on their own.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_closure_ordering.jl
Env:  NTOL (default 1e-8)   DSMAX (default 0.2)   MAXSTEPS (default 200)
      BLABNAT_PCT (default -3) — 2026-07-31: at BOTH -3 and -0.1 the short-run
      branch under `TERM_SR_SWAPS` folded at λ*≈0.9998 regardless of target
      size, when the long-run leg used `TERM_CMF_SWAPS`. That was a matching
      bug (see Method above), not a magnitude problem — this run now uses
      `TERM_LR_SWAPS_MATCHED` for the long-run leg instead. Retest at -3 first;
      only fall back to a smaller BLABNAT_PCT if the fold persists even with
      the matched pair.
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, Printf

const NTOL       = parse(Float64, get(ENV, "NTOL", "1e-8"))
const DSMAX      = parse(Float64, get(ENV, "DSMAX", "0.2"))
const MAXSTEPS   = parse(Int, get(ENV, "MAXSTEPS", "200"))
const BLABNAT_PCT = parse(Float64, get(ENV, "BLABNAT_PCT", "-3"))

agg6, params = cached_pipeline(6)

println("="^72)
println("V6 — closure ordering: TERM_CMF_REFERENCE's shocks, long-run vs short-run")
println("="^72)

# ── Run 1: long-run closure (TERM_LR_SWAPS_MATCHED) — straight homotopy, as usual ──
println("\n── run 1: long-run closure (TERM_LR_SWAPS_MATCHED — xcap endogenous, labour matched to TERM_SR_SWAPS) ──")
sc_lr = Scenario(name = "TERM_LR_SWAPS_MATCHED, blabnat=$(BLABNAT_PCT)% delUnity=1 (V6)",
                  source = "closures.jl TERM_LR_SWAPS_MATCHED — TERM_CMF_SWAPS with " *
                            "termlr.cmf's textbook long-run labour swap (flabsup_id=xlab_id) " *
                            "in place of TERM.CMF's own delfwage=flabsup_id",
                  swaps = TERM_LR_SWAPS_MATCHED,
                  shocks = ["blabnat" => pct(BLABNAT_PCT), "delUnity" => 1.0],
                  numeraire = TERM_CMF_REFERENCE.numeraire,
                  notes = "V6 — see verify_closure_ordering.jl header for why this replaced TERM_CMF_SWAPS.")
el1 = @elapsed r1 = run_model!(agg6, params, sc_lr; tol = NTOL, gdp = false)
@printf("solved=%s  t_reached=%.6g  steps/rejects=%d/%d  ‖F‖∞=%.3e  (%.1f min)\n",
        r1.solved, r1.t_reached, r1.nsteps, r1.nrejects, r1.residual, el1/60)
@assert r1.solved "run 1 (long-run) did not reach the target — nothing to compare"
@assert r1.residual <= max(NTOL, 1e-8) "run 1 (long-run) residual $(r1.residual) is not converged"

# ── Run 2: short-run closure (TERM_SR_SWAPS) — staged, per §N.19e's recipe ──
println("\n" * "─"^72)
println("run 2: short-run closure (TERM_SR_SWAPS — xcap exogenous), STAGED")
println("─"^72)

t = @elapsed ((m2, vars2) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t; digits=1))s")
initialize_model!(m2, vars2; bmk_levels=benchmark_levels(params), numeraire=sc_lr.numeraire)
println("closure: TERM_SR_SWAPS (4 swaps)")
apply_swaps!(m2, vars2, collect(TERM_SR_SWAPS); verbose=true)

r0 = solve_newton!(m2, vars2; maxit=5, tol=NTOL, verbose=false)
println("benchmark: status=$(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= NTOL "run 2 benchmark solve failed under TERM_SR_SWAPS (‖F‖∞=$(r0.residual)) — " *
                            "a closure that cannot reproduce its own base year makes every number past " *
                            "this point meaningless, before the ordering question is even reached"

dU = vars2["delUnity"]
@assert JuMP.is_fixed(dU) "delUnity is not exogenous under TERM_SR_SWAPS — cannot stage the shock"
println("\nSTAGE A — delUnity 0 → 1 (natural continuation)")
elA = @elapsed rA = continuation_solve!(m2, vars2, dU, 1.0;
                                        h0=0.25, hmin=1e-4, tol=NTOL, maxit=30)
@printf("STAGE A: %s  steps/rejects %d/%d  ‖F‖∞=%.3e  (%.1f min)\n",
        rA.reached_target ? "✅ reached" : "❌ STOPPED SHORT",
        rA.nsteps, rA.nrejects, rA.residual, elA/60)
@assert rA.reached_target "STAGE A failed — the one-year-forward baseline does not solve under " *
                          "TERM_SR_SWAPS with NO productivity shock applied at all. This is a " *
                          "separate, more serious problem than the shock itself, and V6 cannot " *
                          "proceed to a comparison."

const BLABNAT_TARGET = pct(BLABNAT_PCT)
lam = let v = vars2["blabnat"]; v isa AbstractArray ? first(v) : v end
println("\nSTAGE B — blabnat 1 → $(BLABNAT_TARGET) (pseudo-arclength)")
elB = @elapsed rB = arclength_solve!(m2, vars2, lam, BLABNAT_TARGET;
                                     ds0=0.005, dsmin=1e-7, dsmax=DSMAX, tol=NTOL,
                                     maxit=20, maxsteps=MAXSTEPS)
@printf("STAGE B: reached=%s  turned=%s  λ=%.8g  steps/rejects=%d/%d  ‖F‖∞=%.3e  (%.1f min)\n",
        rB.reached_target, rB.turned, rB.lam_reached, rB.nsteps, rB.nrejects, rB.residual, elB/60)

if !rB.reached_target
    if rB.turned
        println("\n↩ STAGE B TURNED at λ* ≈ $(round(rB.lam_fold; sigdigits=10)) without reaching $(BLABNAT_TARGET) — " *
                "the short-run closure genuinely cannot sustain this shock at full strength. " *
                "V6 cannot compare a converged short-run point that does not exist; report this " *
                "as the finding rather than forcing a comparison.")
    else
        println("\n❌ STAGE B stopped without reaching the target and without turning — not a fold, " *
                "points at a bifurcation, a domain boundary, or a defect in the staged solve itself.")
    end
    println("\nDone (V6 INCONCLUSIVE — short-run closure did not reach a comparable point).")
    exit(0)
end
@assert rB.residual <= NTOL "STAGE B \"reached\" but residual $(rB.residual) is not converged"

# ── Snapshot run 2's final point in the same shape as r1.values ─────────────
function snapshot_values(vars::Dict{String,Any})
    out = Dict{String,Any}()
    for (nm, v) in vars
        if v isa AbstractArray
            a = Array{Float64}(undef, size(v))
            for idx in eachindex(v)
                vr = v[idx]
                a[idx] = JuMP.is_fixed(vr) ? JuMP.fix_value(vr) : JuMP.start_value(vr)
            end
            out[nm] = a
        else
            out[nm] = JuMP.is_fixed(v) ? JuMP.fix_value(v) : JuMP.start_value(v)
        end
    end
    out
end
V2 = snapshot_values(vars2)

# ── Compare ───────────────────────────────────────────────────────────────
k_gdp = findfirst(==("RealGDP"), MAINMACROS)
@assert k_gdp !== nothing "MAINMACROS has no \"RealGDP\" entry — cannot compare aggregate response"

gdp_lr = r1.values["NatMacro"][k_gdp]
gdp_sr = V2["NatMacro"][k_gdp]
dgdp_lr = gdp_lr - 1.0   # benchmark = 1
dgdp_sr = gdp_sr - 1.0

println("\n" * "="^72)
println("NATIONAL REAL GDP RESPONSE (NatMacro(\"RealGDP\"), benchmark = 1)")
@printf("  long-run:  %.8f   (ΔRealGDP = %+.6f)\n", gdp_lr, dgdp_lr)
@printf("  short-run: %.8f   (ΔRealGDP = %+.6f)\n", gdp_sr, dgdp_sr)
gdp_ordered = abs(dgdp_lr) > abs(dgdp_sr)
println(gdp_ordered ?
    "  ✅ long-run response is larger in magnitude, as theory requires" :
    "  ❌ long-run response is NOT larger — ordering violated at the aggregate level")

xtot1, xtot2 = r1.values["xtot"], V2["xtot"]
na, nr = size(xtot1)
@assert size(xtot2) == (na, nr) "xtot shape mismatch between runs"

ok = 0; bad = 0
worst = Tuple{Float64,Int,Int,Float64,Float64}[]   # (margin, i, d, dlr, dsr)
for d in 1:nr, i in 1:na
    dlr = xtot1[i,d] - 1.0
    dsr = xtot2[i,d] - 1.0
    if abs(dlr) >= abs(dsr)
        global ok += 1
    else
        global bad += 1
    end
    push!(worst, (abs(dsr) - abs(dlr), i, d, dlr, dsr))
end

println("\n" * "="^72)
n = na * nr
@printf("SECTORAL OUTPUT (xtot[i,d], %d×%d = %d cells)\n", na, nr, n)
@printf("  %d/%d (%.1f%%) cells order long-run ≥ short-run in magnitude\n", ok, n, 100ok/n)
majority_ok = ok > n / 2
println(majority_ok ?
    "  ✅ majority of sectors order the same way as the aggregate" :
    "  ❌ majority of sectors do NOT order the same way — not just isolated exceptions")

println("\n  worst 10 exceptions (short-run response exceeds long-run, by margin):")
sort!(worst; by = x -> -x[1])
for (margin, i, d, dlr, dsr) in first(worst, 10)
    margin > 0 || break
    sname = i <= length(AGGCOM) ? AGGCOM[i] : "sector$i"
    @printf("    %-14s reg%-2d  Δlr=%+.6f  Δsr=%+.6f  margin=%.6f\n", sname, d, dlr, dsr, margin)
end

println("\n" * "="^72)
ok_all = gdp_ordered && majority_ok
println(ok_all ?
    "✅ CLOSURE ORDERING HOLDS.\n" *
    "   The long-run closure produces a larger-magnitude real output response than the\n" *
    "   short-run closure, both at the aggregate level and for a majority of sectors." :
    "❌ CLOSURE ORDERING VIOLATED — see the detail above.")
println("\nDone.")
