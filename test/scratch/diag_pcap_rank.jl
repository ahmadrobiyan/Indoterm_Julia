"""
Which capital cell actually drives the fold? (Not necessarily [2,5].)

`test/verify_swapped_closure.jl` prints `pcap(2,5)` as its probe, and that number
accelerates as the walk approaches its obstruction. But it prints that cell because
an EARLIER residual investigation flagged Livestock/BaliNusa — it is a fixed probe,
not a measurement of what dominates. Reading "the fold is at [2,5]" off a script
that only ever looks at [2,5] would be circular.

`test/diag_cell25_scale.jl` then removed the obvious reason to blame that cell:

    VCAP[2,5] = 1543.5, rank 41 of 150, capital share 0.1189
    other Livestock cells:                capital share 0.1117 - 0.1141

so [2,5] is completely ordinary. What that scan DID find is real degeneracy
elsewhere, in the same corner of the grid:

    VCAP[5,5] = 0.0403   VCAP[5,6] = 0.0328   VCAP[6,5] = 8.97   VCAP[5,4] = 20.0

against a region-5 column of ~1,500-43,000 — four orders of magnitude down. A
sector with essentially no capital has a near-vertical supply curve for it: a
small quantity change forces a huge rental-price change. That is the classic
source of an early limit point, and it is the same pathology `_MIN_PRICED_FLOW`
already fixes for dust trade flows.

`test/diag_linearity.jl` showed every mover sharing one doubling ratio (2.191 -
2.210, spread 0.9%), i.e. the response is dominated by a SINGLE mode. So ranking
the `pcap` components of that mode names the cell that drives it — no fold-chasing
walk required, and a small converged shock is enough to expose the direction.

This solves ALL SIX at a small shock and ranks every pcap[i,d] by RELATIVE move
(relative is right here, and safe: pcap is a price with a benchmark near 1, not
an additive variable sitting at 0 — the trap that produced the phantom "574x
runaway" applies to `del*` variables, not to these). Each cell is printed beside
its VCAP so the correlation between tiny capital and huge price response is
directly visible.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_pcap_rank.jl
Env:  SHOCK (default 0.999999), TOPN (default 20)
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, Printf

SHOCK = parse(Float64, get(ENV, "SHOCK", "0.999999"))
TOPN  = parse(Int,     get(ENV, "TOPN", "20"))

agg6, params = cached_pipeline(6)
bmk = benchmark_levels(params)
VCAP = parent(agg6["1CAP"])

(m, vars) = build_model_full!(agg6, params)
initialize_model!(m, vars; bmk_levels=bmk)
apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("\nbenchmark: status=$(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= 1e-8 "benchmark regressed"

pcap = vars["pcap"]
na, nr = size(VCAP)
base = [JuMP.start_value(pcap[i,d]) for i in 1:na, d in 1:nr]

sv = vars["blabnat"]
svr = sv isa AbstractArray ? first(sv) : sv
JuMP.fix(svr, SHOCK; force=true)
r = solve_newton!(m, vars; maxit=30, verbose=false)
println("shock blabnat=$SHOCK: status=$(r.status)  iters=$(r.iters)  ‖F‖∞=$(r.residual)")
r.status == :converged || println("⚠️  not converged — the ranking below is not trustworthy.")

rows = Tuple{Int,Int,Float64,Float64,Float64}[]   # i, d, base, relmove, VCAP
for i in 1:na, d in 1:nr
    b = base[i,d]
    (b === nothing || !isfinite(b) || abs(b) < 1e-12) && continue
    now = JuMP.start_value(pcap[i,d])
    now === nothing && continue
    push!(rows, (i, d, b, abs(now - b) / abs(b), VCAP[i,d]))
end
sort!(rows; by = x -> -x[4])

δ = 1.0 - SHOCK
@printf("\n── top %d pcap cells by RELATIVE move (shock δ = %.3g) ──\n", TOPN, δ)
@printf("  %-10s %-14s %-14s %-14s %s\n", "cell", "pcap base", "rel move", "amplif (/δ)", "VCAP")
for (i, d, b, rel, vc) in first(rows, TOPN)
    @printf("  [%2d,%d]%-4s %-14.8g %-14.6g %-14.4g %.6g%s\n", i, d, "", b, rel, rel/δ, vc,
            (i == 2 && d == 5) ? "   <== the printed probe" : "")
end

k = findfirst(x -> x[1] == 2 && x[2] == 5, rows)
k === nothing ? println("\n  [2,5] absent from the ranking") :
    @printf("\n  [2,5] ranks %d of %d — rel move %.6g, amplification %.4g\n",
            k, length(rows), rows[k][4], rows[k][4]/δ)

# Does a small capital stock predict a large price response? If the top movers are
# the dust cells the scan found, that is the vertical-supply-curve mechanism and
# the remedy is a data floor, not a continuator.
top = first(rows, TOPN)
@printf("\n  median VCAP over the top %d movers = %.6g\n", TOPN,
        sort([x[5] for x in top])[cld(length(top), 2)])
@printf("  median VCAP over the whole grid     = %.6g\n",
        sort(vec(VCAP))[cld(length(VCAP), 2)])
dust = count(x -> x[5] < 100.0, top)
@printf("  %d of the top %d movers have VCAP < 100 (grid has %d such cells of %d)\n",
        dust, TOPN, count(<(100.0), vec(VCAP)), length(VCAP))

println("\n", dust >= cld(TOPN, 2) ?
    "➡ The fold is driven by DUST CAPITAL CELLS, not by [2,5]. A sector with almost no\n" *
    "  capital has a near-vertical rental supply curve, so an ordinary shock produces an\n" *
    "  enormous price move and an early limit point. The fix is a minimum-stock floor in\n" *
    "  the data — the same remedy as _MIN_PRICED_FLOW for dust trade flows — NOT a\n" *
    "  pseudo-arclength continuator, which would faithfully trace a branch that should\n" *
    "  not exist." :
    "➡ The top movers are NOT dust cells, so the vertical-supply-curve story does not\n" *
    "  explain this fold. The limit point looks like a genuine property of the calibrated\n" *
    "  model, and pseudo-arclength is the right next build.")
println("\nDone.")
