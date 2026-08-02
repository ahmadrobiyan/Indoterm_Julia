"""
Run the authors' ACTUAL reference scenario: `blabnat = -3` **with** `delUnity = 1`.

Every run before this one applied only half of `TERM.CMF`'s shipped test simulation:

    TERM.CMF:112   Shock blabnat =   -3; ! 3% increase labour productivity
    TERM.CMF:113   Shock delunity=1;     ! move forward 1 year

`delUnity` is a GEMPACK (change) variable, benchmark 0, additive in this port
(`build_dynamics!.jl:19-21`), so the shock is simply fixing it at 1.0. It gates three
constant terms:

  * `TERM.TAB:2700` capital accumulation — `CAPSTOK_OLDP*log(xcap/CAP) = CAPADD*delUnity + faccum`.
    With `delUnity = 0` AND `faccum` exogenous at zero the right side is **exactly zero**, so
    `xcap` is pinned at the benchmark even though the swapped closure makes it endogenous.
    This is a hard pin and is the one unambiguous difference.
  * `TERM.TAB:2803` the rate-of-return **catch-up** constant `(GROSSRET0-GRETEXP0)*delUnity`.
  * `TERM.TAB:2909` the wage/employment **momentum** constant `(EMPRAT0-1)*delUnity`.

Note carefully what `delUnity = 0` does NOT do: the endogenous response terms in those last
two equations (`delgret`, `delempratio`) are present either way, so investment and wages still
respond to the shock. Only the constant catch-up/momentum terms vanish. Do not describe
`delUnity = 0` as "adjustment switched off" — it is not.

**Two stages, because both shocks move away from the benchmark.**

  A. `delUnity` 0 → 1 alone. This is the model's own one-year-forward baseline with no
     productivity shock, and it is worth seeing on its own: it says how much capital arrives
     from last year's investment.
  B. from there, `blabnat` 1 → 0.97 by **arclength**, so that if this still folds we measure
     where instead of merely failing.

**Reading the outcome.** PLAN.md §N.19c measured a fold at λ* = 0.9908083 with `delUnity = 0`.
  1. Stage B reaches 0.97 ⇒ the fold was an artifact of running an incomplete scenario. The
     −3% productivity result exists, and task #36 is unnecessary.
  2. Stage B folds again ⇒ the fold survives the authors' own closure and is a real property.
     Compare `lam_fold` against 0.9908083: if it moved, `delUnity` matters but does not rescue;
     if it barely moved, `delUnity` is irrelevant to this obstruction.
  3. Stage A itself fails ⇒ the one-year-forward baseline does not solve. That is a separate
     and more serious problem than the productivity shock, and it must be reported as such
     rather than blamed on `blabnat`.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_delunity.jl
Env:  TARGET (default 0.97)   DSMAX (default 0.02)   MAXSTEPS (default 80)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, LinearAlgebra, Printf

const TARGET   = parse(Float64, get(ENV, "TARGET", "0.97"))
# Run 1 exhausted an 80-step budget at dsmax=0.02 while still descending smoothly
# (2 Newton iterations, 0 rejections, every step). dλ/ds decays like C/s along this
# branch — 0.0509 at s≈0.03 down to 0.0051 at s≈1.6 — so λ ≈ λ₀ − C·ln(s) and reaching
# 0.97 needs roughly a doubling of arclength, not more tiny steps. Hence dsmax up 10×.
# The branch is smooth enough to earn it, and a rejection just halves ds anyway.
const DSMAX    = parse(Float64, get(ENV, "DSMAX", "0.2"))
const MAXSTEPS = parse(Int, get(ENV, "MAXSTEPS", "200"))
const TOL      = 1e-8

println("="^72)
println("TERM.CMF reference scenario  —  blabnat → $TARGET  WITH  delUnity = 1")
println("="^72)
flush(stdout)

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))
apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))
println("closure: TERM.CMF all six swaps"); flush(stdout)

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("benchmark: $(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= TOL "benchmark regressed — fix that before reading anything below"
flush(stdout)

dU  = vars["delUnity"]
lam = let v = vars["blabnat"]; v isa AbstractArray ? first(v) : v end
@assert JuMP.is_fixed(dU) "delUnity is not exogenous — the closure is not what this test assumes"
println("delUnity starts at $(JuMP.fix_value(dU))  (TERM.CMF shocks it to 1)")

# How much capital does delUnity=1 actually inject? Worth knowing BEFORE the solve,
# because it bounds how much of an adjustment margin this can possibly provide.
if haskey(params, "CAPADD")
    CAPADD = parent(params["CAPADD"]); CSO = parent(params["CAPSTOK_OLDP"])
    g = [CSO[i] > 1e-10 ? CAPADD[i] / CSO[i] : 0.0 for i in eachindex(CAPADD)]
    nz = filter(!iszero, g)
    if !isempty(nz)
        @printf("capital injection log(xcap/CAP) = CAPADD/CAPSTOK_OLDP: min %.4g  median %.4g  max %.4g\n",
                minimum(nz), sort(nz)[cld(length(nz), 2)], maximum(nz))
    end
end
flush(stdout)

println("\n" * "─"^72)
println("STAGE A — delUnity 0 → 1  (one year forward, no productivity shock)")
println("─"^72); flush(stdout)
elA = @elapsed rA = continuation_solve!(m, vars, dU, 1.0;
                                        h0=0.25, hmin=1e-4, tol=TOL, maxit=30,
                                        report=[("natcapstok", 2), ("realwage_id", 1)])
@printf("\nSTAGE A: %s  delUnity=%.6g  steps/rejects %d/%d  ‖F‖∞=%.4g  (%.1f min)\n",
        rA.reached_target ? "✅ reached" : "❌ STOPPED SHORT",
        rA.s_reached, rA.nsteps, rA.nrejects, rA.residual, elA/60)   # ContinuationResult calls it s_reached, not lam_reached
flush(stdout)

if !rA.reached_target
    println("""
⛔ STAGE A FAILED — the one-year-forward baseline does not solve, with NO productivity
   shock applied at all. This is a separate and more serious problem than `blabnat`, and
   it must not be reported as "the -3% scenario failed". Outcome 3 of the three above.
   Stage B is not run, because starting it from a non-solution would be meaningless.""")
    println("\nDone.")
    exit(0)
end

println("\n" * "─"^72)
println("STAGE B — blabnat 1 → $TARGET  by arclength, with delUnity = 1")
println("─"^72); flush(stdout)
elB = @elapsed rB = arclength_solve!(m, vars, lam, TARGET;
                                     ds0=0.005, dsmin=1e-7, dsmax=DSMAX, tol=TOL,
                                     maxit=20, maxsteps=MAXSTEPS,
                                     report=[("whouhtot", 4), ("wlux", 4)])

println("\n" * "="^72)
println("RESULT   (stage A $(round(elA/60, digits=1)) min + stage B $(round(elB/60, digits=1)) min)")
println("="^72)
@printf("reached target : %s\n", rB.reached_target ? "✅ yes" : "no")
@printf("λ reached      : %.10g   (target %.10g)\n", rB.lam_reached, rB.target)
@printf("turned         : %s\n", rB.turned ? "yes, at λ* ≈ $(round(rB.lam_fold; sigdigits=10))" : "no")
@printf("steps/rejects  : %d / %d\n", rB.nsteps, rB.nrejects)
@printf("final ‖F‖∞     : %.4g\n", rB.residual)

# Same guard as verify_arclength.jl: a "passed the fold" claim must not be manufactured
# by a loosened tolerance.
bad = [p for p in rB.path if !(p[3] <= TOL)]
@assert isempty(bad) "accepted points that are NOT solutions: $bad"

if !isempty(rB.path)
    lams = [p[1] for p in rB.path]
    @printf("λ range traced : %.10g … %.10g  (extreme %.4g%%)\n",
            maximum(lams), minimum(lams), (minimum(lams) - 1) * 100)
end

const FOLD_WITHOUT = 0.9908083427   # PLAN.md §N.19c, delUnity = 0
println()
if rB.reached_target
    println("""✅ THE FOLD WAS AN INCOMPLETE-SCENARIO ARTIFACT. With `delUnity = 1` — the
   authors' actual reference simulation — `blabnat = $TARGET` solves. The fold measured at
   λ* = $FOLD_WITHOUT in §N.19c was a property of running half the scenario, not of the
   model. Task #36 (region-4 loop gain) is unnecessary. Update §N.19c/§N.19d and make
   `delUnity = 1` the default for any year-forward scenario.""")
elseif rB.turned
    shift = rB.lam_fold - FOLD_WITHOUT
    println("""↩ THE FOLD SURVIVES the authors' own closure, at λ* ≈ $(round(rB.lam_fold; sigdigits=10))
   versus $FOLD_WITHOUT without delUnity (shift $(round(shift; sigdigits=3))).
   $(abs(shift) < 1e-4 ?
     "That shift is negligible: `delUnity` is IRRELEVANT to this obstruction, and the missing
   shock was never the explanation. Task #36 becomes the live question again." :
     "`delUnity` moves the fold but does not remove it, so it is part of the story and not
   the whole of it. Task #36 stays live.")
   Either way the -3% productivity scenario does not exist on this branch.""")
elseif rB.nsteps >= MAXSTEPS
    # Distinguish "ran out of budget" from "stopped for an unexplained reason". Run 1
    # conflated them and printed the defect verdict for a branch that was descending
    # perfectly happily — 0 rejections, 2 Newton iterations per step, ‖F‖∞ ~7e-10.
    println("""⏱ STEP BUDGET EXHAUSTED at λ = $(round(rB.lam_reached; sigdigits=10)) — NOT a fold and NOT a
   defect. $(rB.nrejects) rejection(s) over $(rB.nsteps) step(s) says the branch was still tracking
   cleanly; it simply ran out of steps. Raise MAXSTEPS/DSMAX and continue.
   $(rB.lam_reached < FOLD_WITHOUT ?
     "NOTE: this already passed λ = $FOLD_WITHOUT, where the delUnity=0 branch folded, without
   turning. The fold is gone. That conclusion does not depend on reaching 0.97." : "")""")
else
    println("""❌ stopped without reaching the target, without turning, and with step budget to
   spare. A fold cannot defeat arclength, so this points at a bifurcation, a domain boundary,
   or a defect — NOT at a fold. Do not report this as confirming anything about delUnity.""")
end
println("\nDone.")
