"""
Is the fold at δ ≈ 2e-4 an economic feature, or a degenerate-data artifact at [2,5]?

`test/verify_swapped_closure.jl` under the five faithful swaps walks three steps and
then obstructs, with `pcap[2,5]` (Livestock / BaliNusa) accelerating all the way:

    s=0.99985     δ=1.500e-4   pcap[2,5]=1.0081684   Δ/δ = 54.5
    s=0.9998125   δ=1.875e-4   pcap[2,5]=1.0122202   Δ/δ = 65.2
    s=0.9998078   δ=1.922e-4   pcap[2,5]=1.0132692   Δ/δ = 68.3

and `test/diag_linearity.jl` independently finds every mover sharing one doubling
ratio (2.191–2.210, spread 0.9%), i.e. a single dominant mode with ~10% quadratic
content at δ=1e-6. Both point at a limit point.

The fold is therefore real IN THIS SYSTEM. But its LOCATION is the problem: a fold
at a 0.02% labour-supply shock is not credible economics. Published TERM models
take ±5% shocks routinely. A limit point 250x smaller than the shocks the model
exists to run is far more likely to come from a near-degenerate cell than from
Indonesian macroeconomics — and it is worth settling BEFORE building a
pseudo-arclength continuator, because arclength would faithfully trace a branch
that should not be there.

The mechanism to test: primary-factor demand is CES over labour/capital/land, and
a sector whose capital stock is tiny has a nearly vertical supply curve for
capital — a small quantity change forces a large rental-price change. If
VCAP[2,5] is orders of magnitude below the other cells, the huge `pcap[2,5]`
elasticity is arithmetic, not economics.

This prints the capital, labour and land value flows for the whole 25x6 grid,
ranks the cells by capital share and by absolute size, and marks where [2,5]
sits. It is pure data inspection — no solve — so it is cheap and cannot be
confounded by the solver.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_cell25_scale.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using Printf

agg6, params = cached_pipeline(6)

# The factor value flows live in the aggregated DATABASE, not in `params` —
# `build_model!.jl:16` reads them as `parent(agg["1CAP"])`, and `parent` strips the
# NamedArray wrapper so plain 1-based indexing applies. Labour carries an extra
# occupation dimension (`build_premod!.jl:76` sums it out), so collapse it the
# same way rather than assuming a 2-D array.
VCAP = parent(agg6["1CAP"])
VLND = parent(agg6["1LND"])
_lab = parent(agg6["1LAB"])
VLAB = ndims(_lab) == 3 ? dropdims(sum(_lab; dims=2); dims=2) : _lab
na, nr = size(VCAP)
@printf("\ngrid: %d sectors x %d regions\n", na, nr)

println("\n══ capital value flow VCAP[i,d] ══")
@printf("  %-4s", "sec")
for d in 1:nr; @printf(" %12d", d); end
println()
for i in 1:na
    @printf("  %-4d", i)
    for d in 1:nr; @printf(" %12.4g", VCAP[i,d]); end
    println(i == 2 ? "   <== Livestock" : "")
end

# Rank by ABSOLUTE size: a vertical supply curve comes from a small STOCK, and
# the share can look ordinary while the level is negligible.
cells = [(i, d, VCAP[i,d]) for i in 1:na, d in 1:nr] |> vec
sort!(cells; by = x -> x[3])
println("\n══ 12 smallest capital cells ══")
for (i, d, v) in first(cells, 12)
    @printf("  VCAP[%2d,%d] = %-14.6g%s\n", i, d, v,
            (i == 2 && d == 5) ? "   <== the fold cell" : "")
end
k = findfirst(x -> x[1] == 2 && x[2] == 5, cells)
@printf("\n  VCAP[2,5] = %.6g  — rank %d of %d smallest\n", VCAP[2,5], k, length(cells))
@printf("  median capital cell = %.6g   ratio [2,5]/median = %.4g\n",
        sort([c[3] for c in cells])[cld(length(cells), 2)],
        VCAP[2,5] / sort([c[3] for c in cells])[cld(length(cells), 2)])

println("\n══ primary-factor mix in the fold cell and its neighbours ══")
@printf("  %-12s %14s %14s %14s %10s\n", "cell", "VLAB", "VCAP", "VLND", "cap share")
for (i, d) in ((2,5), (2,1), (2,2), (2,3), (2,4), (2,6), (1,5), (3,5))
    tot = VLAB[i,d] + VCAP[i,d] + VLND[i,d]
    @printf("  [%2d,%d]%-6s %14.6g %14.6g %14.6g %10.4f%s\n", i, d, "",
            VLAB[i,d], VCAP[i,d], VLND[i,d], tot > 0 ? VCAP[i,d]/tot : NaN,
            (i == 2 && d == 5) ? "   <==" : "")
end

tiny = count(x -> x[3] < 1.0, cells)
@printf("\n  %d of %d capital cells are below 1.0 in value-flow units\n", tiny, length(cells))
println("\nRead this as: if VCAP[2,5] is orders of magnitude below the median, the fold")
println("is a degenerate-cell artifact and the fix is a minimum-stock floor in the data")
println("(the same remedy as _MIN_PRICED_FLOW for dust trade flows), NOT a continuator.")
println("If VCAP[2,5] is unremarkable, the limit point is a genuine property of the")
println("calibrated model and pseudo-arclength is the right next build.")
println("\nDone.")
