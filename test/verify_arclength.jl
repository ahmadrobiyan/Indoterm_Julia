"""
Can pseudo-arclength continuation pass the fold that stopped everything else?

PLAN.md §N.16 measured a genuine fold on the `blabnat` branch at λ* ≈ 0.9908
(‖dx/dλ‖∞ 1.95e6 → 3.51e7, σ_min 3.85e-9 → 1.16e-10, direction `whouhtot[4]`).
Natural-parameter continuation — `continuation_solve!`, and the homotopy inside
`run_model!` — provably cannot pass it: past a limit point there is no solution
at that λ at all, so no step size helps. `arclength_solve!` parameterises by
arclength instead, and its augmented Jacobian stays nonsingular through the turn.

**Read the outcome carefully — three different things count as success here, and
only one of them is "the shock ran".**

  1. `reached_target` — the branch actually gets to `blabnat = 0.97`. That means
     the fold was a fold in the *parameterisation* only and the target was always
     on this branch; every earlier method just could not steer there.
  2. `turned` and NOT `reached_target` — the branch genuinely folds back. Then
     `blabnat = 0.97` **does not exist** in this closure, and `lam_fold` is the
     largest labour-supply shock the model admits. That is a real result about
     the model, not a failure of the code, and it is reported as such.
  3. Neither — arclength stopped for some other reason. This is the only outcome
     that indicates a defect, because a fold cannot defeat arclength. Suspect a
     bifurcation (where the augmented system is ALSO singular), a domain
     boundary, or a bug in the augmented assembly.

Guard against the obvious self-deception: assert the benchmark still solves
first, and assert every accepted point is a genuine solution (‖F‖∞ ≤ tol), so a
"passed the fold" claim cannot be manufactured by a loosened tolerance.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_arclength.jl
Env:  TARGET (default 0.97)   DS0 (default 0.05)   MAXSTEPS (default 120)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, LinearAlgebra, Printf

const TARGET   = parse(Float64, get(ENV, "TARGET", "0.97"))
const DS0      = parse(Float64, get(ENV, "DS0", "0.005"))
const DSMAX    = parse(Float64, get(ENV, "DSMAX", "0.02"))
const MAXSTEPS = parse(Int, get(ENV, "MAXSTEPS", "60"))
const TOL      = 1e-8

println("="^72)
println("pseudo-arclength continuation  —  blabnat → $TARGET")
println("="^72)

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))
apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))
println("closure: TERM.CMF all six swaps"); flush(stdout)

# solve_newton! is also what SQUARES the model (pins orphans, deletes dead rows),
# which arclength_solve! relies on — so this call is structural, not just a check.
r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("benchmark: $(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= TOL "benchmark regressed — fix that before reading anything below"
flush(stdout)

lam = let v = vars["blabnat"]; v isa AbstractArray ? first(v) : v end

el = @elapsed r = arclength_solve!(m, vars, lam, TARGET;
                                   ds0=DS0, dsmin=1e-7, dsmax=DSMAX, tol=TOL,
                                   maxit=20, maxsteps=MAXSTEPS,
                                   report=[("whouhtot", 4), ("wlux", 4)])

println("\n" * "="^72)
println("RESULT   ($(round(el/60, digits=1)) min)")
println("="^72)
@printf("reached target : %s\n", r.reached_target ? "✅ yes" : "no")
@printf("λ reached      : %.10g   (target %.10g)\n", r.lam_reached, r.target)
@printf("turned         : %s\n", r.turned ? "yes, at λ ≈ $(round(r.lam_fold; sigdigits=10))" : "no")
@printf("steps/rejects  : %d / %d\n", r.nsteps, r.nrejects)
@printf("final ‖F‖∞     : %.4g\n", r.residual)
@printf("local coord    : %s\n", r.local_coord)

println("\n── branch traced ──")
println("λ                  iters   ‖F‖∞")
for (l, it, res) in r.path
    @printf("%-19.10g %-7d %.5g\n", l, it, res)
end

# Every accepted point must be a real solution — otherwise "passed the fold" is
# just a loosened tolerance wearing a hat.
bad = [p for p in r.path if !(p[3] <= TOL)]
@assert isempty(bad) "accepted points that are NOT solutions: $bad"

if !isempty(r.path)
    lams = [p[1] for p in r.path]
    @printf("\nλ range traced: %.10g … %.10g  (extreme shock %.4g%%)\n",
            maximum(lams), minimum(lams), (minimum(lams) - 1) * 100)
end

println()
if r.reached_target
    println("""✅ THE FOLD WAS PASSABLE. blabnat = $TARGET is on this branch after all —
   natural-parameter continuation could not steer there, arclength could. The
   $(abs((TARGET-1)*100))% labour-PRODUCTIVITY scenario is now runnable end to end; wire
   arclength into run_model! as the fallback when the homotopy stalls.""")
elseif r.turned
    println("""↩ THE BRANCH GENUINELY FOLDS at λ ≈ $(round(r.lam_fold; sigdigits=10))
   (= $(round(-(r.lam_fold-1)*100; sigdigits=4))% labour PRODUCTIVITY). Past that point there is
   NO solution at any smaller blabnat on this branch, so `blabnat = $TARGET` is not a
   solve-harder problem — it does not exist in this closure.

   ⚠️ BEFORE reporting that as a property of the model, check the SCENARIO is
   complete. A fold measured at λ* = 0.9908083 was reported exactly that way and
   was wrong: the run had omitted TERM.CMF's second shock, `delUnity = 1`. With
   both shocks the branch passes straight through and reaches 0.97 in 22 steps
   (PLAN.md §N.19e). `delUnity = 0` is not neutral — it hard-pins all 150 capital
   cells to the benchmark, and suppressing that adjustment margin is what created
   the fold. Compare also against the base static closure, which folds elsewhere
   entirely (0.999752 vs 0.99981 — the two closures are NOT comparable).""")
else
    println("""❌ stopped without reaching the target AND without turning. A fold cannot
   defeat arclength, so this does NOT confirm §N.16 — it points at something else:
   a bifurcation (where the AUGMENTED system is singular too), a domain boundary
   hit by the predictor, or a defect in the augmented assembly in src/arclength.jl.
   Do not report this as "the fold is impassable".""")
end
println("\nDone.")
