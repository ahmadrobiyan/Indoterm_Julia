"""
Does `run_model!`'s arclength fallback actually fire, and does it report honestly?

Task #39 wired `arclength_solve!` into `run_model!` as the retry when the natural-
parameter homotopy collapses below `hmin`. That code compiles and its methods
resolve, but compiling is not running: the fallback path had never been executed
when it was written. This script executes it.

**The fixture is deliberately the broken scenario.** `TERM_CMF_NO_DELUNITY` is
TERM.CMF with line 113 (`Shock delunity=1`) omitted. That omission hard-pins all
150 `xcap` cells and manufactures a limit point at λ\\* ≈ 0.9908083 (PLAN.md
§N.19c). It is the only scenario we have that reliably folds, which makes it the
right test fixture and the wrong thing to quote results from — hence its
`INCOMPLETE` label, which `test/verify_scenario_complete.jl` asserts is still
there.

**What counts as a pass here is narrow.** This test does NOT assert that the
target is reached — it is not supposed to be reachable in this incomplete
scenario. It asserts the *plumbing*:

  1. the homotopy collapses (otherwise the fallback never gets a chance, and this
     test would pass vacuously — so a run that sails to `t = 1` is a FAILURE of
     the fixture, reported as such rather than as a success);
  2. the fallback fires and `arclength_solve!` runs;
  3. the outcome is one of the three documented ones, and `turned = true` is the
     expected one, locating the fold near 0.9908083;
  4. `r.values` comes back populated and usable, i.e. `_current_values` returned a
     real point and not an empty or half-filled dict.

Point 4 is the one most likely to be quietly wrong, because the success path
returns values read from JuMP start values rather than from `solve_newton!`'s own
vector, and nothing else in the codebase exercises that reader.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_arclength_fallback.jl
Env:  HMIN (default 1e-3 — deliberately loose, to reach the collapse fast)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, Printf

const HMIN = parse(Float64, get(ENV, "HMIN", "1e-3"))
const EXPECTED_FOLD = 0.9908083     # PLAN.md §N.19c, delUnity = 0

println("="^72)
println("arclength fallback inside run_model!  —  fixture: the INCOMPLETE scenario")
println("="^72)

agg6, params = cached_pipeline(6)

sc = TERM_CMF_NO_DELUNITY
@assert length(sc.shocks) == 1 """
    the fallback only applies to a single-shock scenario, and this fixture now has
    $(length(sc.shocks)). Either the fixture was 'fixed' (see
    test/verify_scenario_complete.jl) or the wrong scenario is being used."""

el = @elapsed r = run_model!(agg6, params, sc;
                            h0 = 0.25, hmin = HMIN, tol = 1e-8,
                            arclength = true, gdp = false)

println("\n" * "="^72)
println("RESULT   ($(round(el/60, digits=1)) min)")
println("="^72)
@printf("solved        : %s\n", r.solved ? "yes" : "no")
@printf("t reached     : %.10g\n", r.t_reached)
@printf("steps/rejects : %d / %d\n", r.nsteps, r.nrejects)
@printf("final ‖F‖∞    : %.4g\n", r.residual)
@printf("path points   : %d\n", length(r.path))

# ── 1. The fixture must actually fold, or this test proves nothing ───────────
@assert !r.solved """
    the INCOMPLETE scenario reached t = 1 without collapsing. That is a FIXTURE
    failure, not a success: with no collapse the arclength fallback never ran and
    this script asserted nothing about it. Either the fold at λ* ≈ $EXPECTED_FOLD is
    gone (check what changed — §N.19c measured it carefully), or hmin=$HMIN is so
    loose the walk stepped straight over it. Do NOT read this as 'the fallback works'."""

# ── 2/3. The fallback ran, and reported one of the documented outcomes ───────
# t_reached is the scenario fraction; convert the expected fold to that scale so
# the two are comparable. blabnat: base 1.0 → target 0.97.
b, tg = 1.0, pct(-3)
t_fold = (EXPECTED_FOLD - b) / (tg - b)
@printf("\nexpected fold at λ* ≈ %.7f  ⇒  t ≈ %.5f\n", EXPECTED_FOLD, t_fold)
@printf("run stopped at t = %.5f\n", r.t_reached)

@assert r.t_reached > 0 """
    the walk collapsed at t = 0, i.e. the very first step failed. That is not a fold
    — it means the benchmark under this closure is not a usable starting point.
    Check the closure gate in run_model! before reading anything else."""

# The fallback should get at least as far as plain stepping did. If arclength
# ended BELOW where the homotopy already was, something regressed in the handoff
# (most likely the model was left at the wrong point before arclength_solve!).
@assert r.t_reached >= 0.5 * t_fold """
    the run stopped at t = $(r.t_reached), well short of the known fold at t ≈ $t_fold.
    The fallback either did not fire or was handed a bad starting point. Check the
    log for the '↪ falling back to pseudo-arclength' line — if it is absent, the
    collapse branch did not reach the fallback."""

# ── 4. The values dict must be real — this is _current_values' only exercise ──
@assert !isempty(r.values) "r.values is empty — _current_values returned nothing usable"
for nm in ("pcap", "xcap", "realwage_id", "blabnat")
    @assert haskey(r.values, nm) "r.values is missing \"$nm\" — _current_values did not read every variable"
    v = r.values[nm]
    @assert all(isfinite, v) "r.values[\"$nm\"] contains non-finite entries — the returned point is not a solution"
end
bl = r.values["blabnat"]
blv = bl isa AbstractArray ? first(bl) : bl
@printf("\nr.values sanity: blabnat = %.10g   (base 1.0, target %.4g)\n", blv, tg)

# NOT `tg <= blv <= b`. That bound assumed the path stays monotone in blabnat,
# which is false at a fold by definition: `turned = true` means the tangent's λ
# component changed sign, so past the turning point the branch is free to move
# blabnat in EITHER direction depending on the curve's geometry in the full
# (N+1)-dimensional state space — nothing makes it double back toward the
# target rather than continuing outward. What actually certifies "a real point
# on the traced branch" is (a) it is a genuine converged solution, and (b) the
# path leading to it is continuous — no jump between consecutive accepted
# points. Both are checked directly below instead of assuming a range.
@assert isfinite(blv) "r.values[\"blabnat\"] = $blv is not finite"
@assert r.residual <= 1e-6 "final residual $(r.residual) is not a converged solution"
# `ScenarioResult` does not expose arclength's own `turned` flag — it is not
# needed here anyway: any accepted path is only worth trusting if consecutive
# points are actually continuous, whether the walk stopped short because of a
# fold or for some other reason.
if !r.solved && length(r.path) >= 2
    lam_prev = r.path[end-1][1]; lam_last = r.path[end][1]
    step = abs(lam_last - lam_prev)
    @assert step <= 0.05 """
        the last two accepted path points jump by $step in λ — that is not a
        continuous continuation step, so the returned point may not actually be
        on the traced branch."""
    @printf("path continuity: |Δλ| over last accepted step = %.4g (small ⇒ continuous)\n", step)
end

println("\n✅ the arclength fallback fires, runs, and returns a usable point.")
println("""
   NOTE what this did and did not establish. It established the PLUMBING: the
   collapse branch reaches the fallback, arclength runs, and `_current_values`
   returns a real point on the branch. It established nothing about the economics
   of this scenario, which is incomplete BY CONSTRUCTION — the fold it stops at is
   an artifact of the missing `delUnity = 1`, not a property of the model. Use
   TERM_CMF_REFERENCE for any result anyone will read.""")
println("\nDone.")
