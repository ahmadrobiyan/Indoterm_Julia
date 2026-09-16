"""
Do the 2026-09-10 arclength guards turn a bogus verdict into an honest one?

`test/verify_arclength_fallback.jl` already shows the guards do not obstruct a
traversal that works: on the `TERM_CMF_NO_DELUNITY` fixture the trust region
fires 26 times, the discontinuity guard never fires, and the run still lands the
documented fold at λ* = 0.9908086 (expected 0.9908083). That is the no-false-
positive half of the evidence.

This probe is the other half: the case the guards were written for. Before the
patch, every `P028` case reached an `hmin` collapse, handed over to arclength,
and came back with `TURNED ... so the target does not exist in this closure` —
reached through a first accepted step that moved λ by −0.269 against a predicted
1.22e-6 (`logs/p028.log`). The claim under test is that the patched routine now
either tracks the branch or says it cannot, and in neither case reports a
turning point off a discontinuous jump.

**`×0.75` is the fixture, not `×0.5`, because it is the cheapest.** All three of
`×0.5`, `×0.6`, `×0.75` fold; `×0.75` folds earliest (`t ≈ 0.1024` against
`0.6426` for `×0.5`), so it reaches the handover in the least wall-clock. The
guards are properties of `arclength_solve!`, not of the elasticity, so the
cheapest case that exercises them is the right one.

What a PASS looks like, in order of preference:

  1. arclength reaches λ = 0.4054651 — the earlier "does not exist" was an
     artifact outright, and the missing V8 point is recoverable.
  2. `⛔ collapsed below dsmin` with NO `TURNED` line — an honest "this method
     cannot get there from here", which is the usable negative result.
  3. `⚠ DISCONTINUOUS step rejected` lines appear and the run ends without a
     turning point — the guard caught the jump it was written for.

What a FAIL looks like: a `TURNED`/`does not exist` verdict again, or a λ
trajectory that walks off to ≈ −7.3 as before. Either means the guards missed.

Run:  julia --project=IndotermJulia IndotermJulia/test/scratch/_probe_arclength_guards.jl
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using Printf

const FACTOR = 0.75

function main()
    agg6, params = cached_pipeline(6)
    agg2 = copy(agg6)
    agg2["P028"] = vec(parent(agg6["P028"])) .* FACTOR
    p2 = prepare_parameters!(agg2)
    @assert vec(parent(p2["P028"])) != vec(parent(params["P028"])) "P028 was not perturbed"

    println("\n── arclength-guard validation: P028 ×$FACTOR, hmin=1e-6, h0=0.05, maxit=60")
    println("   pre-patch verdict for this case (logs/p028.log): ⛔ hmin at t=0.102435 → ")
    println("   arclength → ⛔ dsmin at λ=-0.3966 → TURNED at λ≈-0.0724 → INCOMPLETE, 8795s")
    # Both flushes are load-bearing; see the note in _probe_p028.jl. Redirected
    # stdout is block-buffered, and verbose=false would leave nothing to flush.
    flush(stdout)

    t0 = time()
    res = run_model!(agg2, p2, COALPRICE_REFERENCE;
                     tol = 1e-8, gdp = false, verbose = true,
                     h0 = 0.05, hmin = 1e-6, maxit = 60)
    dt = time() - t0

    @printf("\n── RESULT (%.1f min)  solved=%s  t=%.4f\n", dt / 60, res.solved, res.t_reached)
    flush(stdout)
end

main()
