"""
What is different about the two cells that break price homogeneity? (task #38)

`verify_homogeneity_tol.jl` established the overshoot is REAL: it is completely
independent of solver tolerance (identical to 4 significant figures while the
residual fell 525×) and it grows with λ. It is confined to `pcap[5,5]` and
`pcap[5,6]` — industry 5 in regions 5 and 6 — against a control cell that sits at
1e-9, six orders of magnitude smaller.

This script does no solving. It reads the benchmark parameters for those cells
and for controls, looking for the structural asymmetry.

THE LEADING HYPOTHESIS is `E_plnd!` (`build_equations.jl:305-325`): when
`LND[i,d] <= 1e-10` it drops the land demand row and pins `plnd[i,d] = 1.0`.
A price pinned at 1.0 does NOT scale with λ, while every other price in the same
CES aggregator does. If the land share `ALPHA_FAC[i,3,d]` is exactly 0 the pinned
price is harmless — `α^σ = 0` kills its term in the aggregator, which is why the
guard is written the way it is. But if `LND` is below the 1e-10 gate while
`ALPHA_FAC[i,3,d]` is merely *small and nonzero*, the pinned price survives into
the aggregator with a small weight and breaks homogeneity by roughly that weight.

That would predict all of what was measured: an error independent of tolerance
(it is structural), growing with λ (it is a fixed price against scaling ones),
and confined to a couple of cells (only those hit the guard inconsistently).

The test is therefore a joint one, and BOTH halves must hold to convict:
  (a) the failing cells are gated (`LND <= 1e-10`), and
  (b) their land share is nonzero anyway.
If (a) holds but (b) gives exactly 0, the hypothesis is dead and the printout
below is the evidence that killed it — read the control rows before concluding.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_pcap_homog_cells.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using Printf

const FAILING = [(5, 5), (5, 6)]
const CONTROL = [(1, 2), (5, 2), (5, 1)]   # (1,2) healthy; (5,·) same industry, big regions

agg6, params = cached_pipeline(6)

CAP   = parent(params["CAP"])
LND   = parent(params["LND"])
PRIM  = parent(params["PRIM"])
LAB_O = parent(params["LAB_O"])
ALPHA = params["ALPHA_FAC"]
GAMMA = params["GAMMA_FAC"]

const LND_GATE = 1e-10      # E_plnd!'s threshold, build_equations.jl:311

@printf("%-8s %12s %12s %12s %12s  %-6s %10s %10s %10s\n",
        "cell", "CAP", "LND", "LAB_O", "PRIM", "gated", "aLAB", "aCAP", "aLND")
println("-"^110)

function row(i, d, tag)
    gated = LND[i, d] <= LND_GATE
    @printf("%-8s %12.5g %12.5g %12.5g %12.5g  %-6s %10.4g %10.4g %10.4g   %s\n",
            "($i,$d)", CAP[i,d], LND[i,d], LAB_O[i,d], PRIM[i,d],
            gated ? "YES" : "no",
            ALPHA[i,1,d], ALPHA[i,2,d], ALPHA[i,3,d], tag)
    return (gated = gated, aLND = ALPHA[i,3,d])
end

fail = [row(i, d, "← FAILS homogeneity") for (i, d) in FAILING]
ctrl = [row(i, d, "control") for (i, d) in CONTROL]

println("\n" * "="^74)
println("VERDICT ON THE PINNED-plnd HYPOTHESIS")
println("="^74)

gated_fail = all(f.gated for f in fail)
live_share = [f.aLND for f in fail]

if !gated_fail
    println("""
    DEAD. The failing cells are NOT gated by E_plnd! (LND > $LND_GATE), so their
    plnd is a free, scaling price like any other and the pinned-price mechanism
    cannot be the cause. Look elsewhere — the next candidates are any other place
    a price is fixed at a constant rather than scaled, or a parameter with price
    units baked in.""")
elseif all(==(0.0), live_share)
    println("""
    DEAD, and cleanly so. The failing cells ARE gated, but their land share is
    EXACTLY zero, so α^σ = 0 removes the pinned plnd from the CES aggregator and
    it cannot influence pcap. The guard is behaving as its comment claims.
    The cause is elsewhere.""")
else
    @printf("""
    CONVICTED, pending an arithmetic check. The failing cells are gated by
    E_plnd! AND carry a nonzero land share (%.4g, %.4g), so a plnd pinned at 1.0
    sits inside the CES price aggregator alongside prices that scale by λ.

    Order-of-magnitude check before believing it: the induced relative error in
    pcap should be roughly the land share times (λ−1). At λ=1.10 the measured
    overshoots were 3.741e-4 and 4.086e-4, so the shares would need to be of
    order %.4g and %.4g. Compare against the actual shares above — if they are
    orders of magnitude apart, this mechanism is present but is NOT the dominant
    one, and saying so is the honest result.\n""",
            live_share[1], live_share[2], 3.741e-4 / 0.10, 4.086e-4 / 0.10)
end

# Whether or not the hypothesis survives, count how widespread the gate is. A
# guard that fires on two cells is a curiosity; one that fires on hundreds is a
# systematic issue and changes how much the 4e-4 matters.
ngated = count(<=(LND_GATE), LND)
nlive  = count(i -> LND[i] <= LND_GATE && ALPHA[CartesianIndex(Tuple(i)[1], 3, Tuple(i)[2])] != 0,
               CartesianIndices(LND))
@printf("\nE_plnd! gate fires on %d of %d (industry,region) cells; %d of those have a nonzero land share.\n",
        ngated, length(LND), nlive)
println("Done.")
