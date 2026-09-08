"""
Does the `_MIN_PRICED_FLOW` guard explain the pcap homogeneity violation? (#38)

MECHANISM. `build_equations.jl:141-146` pins `ppur_s[c,i,d] = 1.0` whenever the
benchmark flow `PUR_S[c,i,d]` falls at or below `_MIN_PRICED_FLOW = 1e-4` Rp bn.
That is a PRICE held at a constant. Scale the numeraire by λ and every other
price in the model becomes λ×itself while these stay at 1.0. Where a pinned cell
still carries a nonzero flow, it enters the cost aggregations with a small but
real weight, and the cost chain for that (industry, region) is left short of λ by
roughly the pinned share of its purchases.

Note what is and is not wrong here. Pinning a price whose flow is EXACTLY zero is
harmless and correct — nothing multiplies it. The leak is the band
`0 < PUR_S <= 1e-4`: nonzero weight, frozen price.

THE PREDICTION IS ARITHMETIC, WHICH IS WHY THIS IS WORTH RUNNING. If the
mechanism is the dominant one, then for each (industry, region)

    relative homogeneity error  ≈  (λ − 1) × (pinned purchases / total purchases)

Measured at λ = 1.10: pcap[5,5] overshot by 3.741e-4 and pcap[5,6] by 4.086e-4.
So the pinned shares must come out near 3.74e-3 and 4.09e-3 respectively, and
near zero for the control cells. A match across two failing cells AND the
controls, on a number this script never sees, would convict. An order-of-
magnitude miss means the mechanism is present but not dominant — which is a
result too, and must be reported as such rather than rounded toward the story.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_pinned_price_share.jl
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using Printf

const MIN_PRICED_FLOW = 1e-4      # = 1e5 * _TINY, build_equations.jl:76
const LAMBDA = 1.10               # the λ the overshoots below were measured at

# (cell, measured pcap overshoot at λ=1.10). NaN where nothing was measured, so
# the control rows still print their share without inventing a target for it.
const CELLS = [((5, 5), 3.741e-4), ((5, 6), 4.086e-4),
               ((5, 2), NaN), ((5, 1), NaN), ((1, 2), NaN)]

agg6, params = cached_pipeline(6)
PUR_S = parent(params["PUR_S"])
na, nu, nr = size(PUR_S)
@printf("PUR_S dims: %d commodities × %d users × %d regions\n", na, nu, nr)

npin  = count(x -> 0 < x <= MIN_PRICED_FLOW, PUR_S)
nzero = count(==(0.0), PUR_S)
@printf("cells pinned with a LIVE flow (0 < PUR_S <= %g): %d of %d\n",
        MIN_PRICED_FLOW, npin, length(PUR_S))
@printf("cells pinned with a dead flow (PUR_S == 0):      %d  (harmless — nothing multiplies them)\n\n",
        nzero)

@printf("%-8s %14s %14s %12s   %12s %12s %s\n",
        "cell", "total PUR_S", "pinned PUR_S", "pinned frac",
        "predicted", "measured", "ratio")
println("-"^104)

results = NTuple{3,Float64}[]
for ((i, d), measured) in CELLS
    col = @view PUR_S[:, i, d]
    tot = sum(col)
    pin = sum(x for x in col if 0 < x <= MIN_PRICED_FLOW; init = 0.0)
    frac = tot > 0 ? pin / tot : NaN
    pred = (LAMBDA - 1) * frac
    ratio = isfinite(measured) && pred != 0 ? measured / pred : NaN
    @printf("%-8s %14.6g %14.6g %12.5g   %12.4g %12s %s\n",
            "($i,$d)", tot, pin, frac, pred,
            isfinite(measured) ? @sprintf("%.4g", measured) : "—",
            isfinite(ratio) ? @sprintf("%.3g", ratio) : "")
    isfinite(measured) && push!(results, (pred, measured, ratio))
end

println("\n" * "="^74)
println("VERDICT")
println("="^74)

if isempty(results)
    println("no measured cells supplied — nothing to convict or acquit.")
else
    ratios = [r[3] for r in results if isfinite(r[3])]
    if isempty(ratios)
        println("""
        ACQUITTED. The failing cells buy nothing through pinned-price cells, so this
        mechanism contributes exactly zero to them. Whatever breaks homogeneity is
        elsewhere.""")
    elseif all(r -> 0.5 <= r <= 2.0, ratios)
        println("""
        CONVICTED. Predicted and measured agree within a factor of 2 on every cell,
        from a share this script computed out of the benchmark database alone and
        never fitted. The _MIN_PRICED_FLOW guard is the cause of the homogeneity
        violation in task #38.

        This is a translation-side defect, not a property of TERM: GEMPACK has no
        such guard. It was added (build_equations.jl:60-76) to stop dust flows from
        throttling the Newton step, and it does that job — but it pays for it by
        freezing a price that should scale, which is exactly a degree-of-homogeneity
        error. The fix is to pin those prices to the NUMERAIRE rather than to the
        constant 1.0, so they scale with everything else while still contributing no
        Jacobian information.""")
    elseif all(r -> 0.05 <= r <= 20.0, ratios)
        println("""
        PRESENT BUT NOT PROVEN DOMINANT. Predicted and measured are within an order
        of magnitude but do not match cleanly. The mechanism is real and contributes,
        but something else contributes comparably. Do NOT close #38 on this — the
        honest statement is that a second cause is still unaccounted for.""")
    else
        println("""
        REJECTED AS THE DOMINANT CAUSE. The predicted and measured errors differ by
        more than an order of magnitude, so whatever the guard does to these cells it
        is not what verify_homogeneity_tol.jl measured. Report the miss and look
        elsewhere rather than adjusting the story to fit.""")
    end
end
println("\nDone.")
