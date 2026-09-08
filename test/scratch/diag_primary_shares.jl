"""
Is the near-singular mode a DATA defect?

`diag_step_quality.jl` localised the α-cap to row 38476,
`gret[i,d] - (log(pcap[i,d]+τ) - log(pinvitot[i,d]+τ))`, for (i,d) = (2,5).
At the stalled point `pcap[2,5] = 0.94473` — a 5.5% fall from its benchmark of 1
for a 0.03% labour shock — and the Newton step wants a further -4.1%. The
curvature of `log` at that size is ½(0.041)² = 8.6e-4, which is exactly the
amount by which that row overshoots the Newton prediction at α=1. One row with a
runaway variable is throttling α for the whole system: the shock equation
`log(alab_o) = log(blabnat)+log(blab)+log(blab_d)` is linear-in-logs and holds
96.7% of the residual mass, yet only gets 25% of its step per iteration.

`pcap[2,5]` and `plnd[2,5]` move *identically* (both x=0.94473, d=-0.038885), and
their CES demand equations carry share coefficients of 3.8e-4 and 1.3e-5. This
script checks the primary-factor value shares in the aggregated database to see
whether those are genuine or a regionalisation/RAS artifact.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_primary_shares.jl
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))

agg, params = cached_pipeline(6)
for k in ("1LAB", "1CAP", "1LND")
    println("$k  size=$(size(agg[k]))")
end

LAB = agg["1LAB"]; CAP = agg["1CAP"]; LND = agg["1LND"]
na, nr = size(CAP)
vlab = [sum(@view LAB[i, :, d]) for i in 1:na, d in 1:nr]   # 1LAB is (IND, OCC, REG)

println("\nna=$na  nr=$nr")
println("totals:  LAB=$(sum(LAB))  CAP=$(sum(CAP))  LND=$(sum(LND))")

prim = vlab .+ CAP .+ LND
scap = [prim[i,d] > 0 ? CAP[i,d]/prim[i,d] : NaN for i in 1:na, d in 1:nr]
slnd = [prim[i,d] > 0 ? LND[i,d]/prim[i,d] : NaN for i in 1:na, d in 1:nr]

println("\n── capital value share of primary factors: distribution ──")
ok = filter(isfinite, vec(scap))
println("  cells=$(length(ok))  min=$(minimum(ok))  median=$(sort(ok)[cld(end,2)])  max=$(maximum(ok))")
for thr in (1e-5, 1e-4, 1e-3, 1e-2, 0.05)
    println("  cells with capital share < $thr : $(count(<(thr), ok))")
end

println("\n── the 20 lowest capital shares (i, d) ──")
idx = sortperm(vec(scap); by = v -> isfinite(v) ? v : Inf)
for k in idx[1:20]
    i, d = fldmod1(k, 1)[1], 0
    i = ((k - 1) % na) + 1; d = ((k - 1) ÷ na) + 1
    println("  i=$(rpad(i,3)) d=$(rpad(d,3))  scap=$(rpad(round(scap[i,d],sigdigits=4),12)) " *
            "slnd=$(rpad(round(slnd[i,d],sigdigits=4),12)) " *
            "VLAB=$(rpad(round(vlab[i,d],sigdigits=6),14)) VCAP=$(rpad(round(CAP[i,d],sigdigits=6),14)) " *
            "VLND=$(round(LND[i,d],sigdigits=6))")
end

println("\n── the suspect cell and its region-5 / industry-2 neighbours ──")
println("  scap[2,5]=$(scap[2,5])  slnd[2,5]=$(slnd[2,5])")
println("  VLAB[2,5]=$(vlab[2,5])  VCAP[2,5]=$(CAP[2,5])  VLND[2,5]=$(LND[2,5])")
println("\n  industry 2 across regions:")
for d in 1:nr
    println("    d=$d  scap=$(rpad(round(scap[2,d],sigdigits=4),12)) VLAB=$(rpad(round(vlab[2,d],sigdigits=6),14)) VCAP=$(round(CAP[2,d],sigdigits=6))")
end
println("\n  region 5 across industries (10 smallest capital shares):")
for i in sortperm([isfinite(scap[i,5]) ? scap[i,5] : Inf for i in 1:na])[1:10]
    println("    i=$i  scap=$(rpad(round(scap[i,5],sigdigits=4),12)) VLAB=$(rpad(round(vlab[i,5],sigdigits=6),14)) VCAP=$(round(CAP[i,5],sigdigits=6))")
end

# Does the same pathology exist in the 34-region source, or is it created by the
# 34->6 region aggregation / RAS?
println("\n── does the pathology predate region aggregation? (25x34) ──")
agg34, _ = cached_pipeline(34)
L34 = agg34["1LAB"]; C34 = agg34["1CAP"]; N34 = agg34["1LND"]
na2, nr2 = size(C34)
v34 = [sum(@view L34[i, :, d]) for i in 1:na2, d in 1:nr2]
p34 = v34 .+ C34 .+ N34
s34 = [p34[i,d] > 0 ? C34[i,d]/p34[i,d] : NaN for i in 1:na2, d in 1:nr2]
ok34 = filter(isfinite, vec(s34))
println("  cells=$(length(ok34))  min=$(minimum(ok34))  median=$(sort(ok34)[cld(end,2)])")
for thr in (1e-4, 1e-3, 1e-2)
    println("  cells with capital share < $thr : $(count(<(thr), ok34))")
end
println("\nDone.")
