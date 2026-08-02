"""
Why is 97% of the shock residual in ONE cell, sector 2 / region 5?

`test/diag_full_newton.jl` (rows2.log) showed every top-5 residual row at [2,5]:
gret[2,5]/pcap[2,5]/pinvitot[2,5] (36%), plab_o|pcap|plnd[2,5] (52%),
xinvitot[2,5]|ggro[2,5] (4.7%), xinvitot[2,5]|finv2[2,5] (4.7%).

Candidate cause: `E_pinvitot!` SKIPS cells with INVEST_C[i,d] <= 1e-10, leaving
pinvitot[i,d] with no defining equation, while E_gret! still writes
gret = log(pcap+tau) - log(pinvitot+tau) for EVERY (i,d). Print the raw flows so
the cell is characterised by data, not by hypothesis.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_cell25.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
const IJ = IndotermJulia

agg6, params = cached_pipeline(6)

INVEST_C = parent(params["INVEST_C"])
CAP_v    = parent(params["CAP"])
na, nr = size(INVEST_C)
println("dims: na=$na nr=$nr")

# How many cells does E_pinvitot! skip, and is [2,5] one of them?
skipped = [(i, d) for i in 1:na, d in 1:nr if INVEST_C[i, d] <= 1e-10]
println("\nE_pinvitot! skips $(length(skipped)) of $(na*nr) cells (INVEST_C <= 1e-10)")
isempty(skipped) || println("  skipped: ", skipped)

println("\n── INVEST_C by cell (investment value, industry x region) ──")
println("        ", join([lpad("r$d", 12) for d in 1:nr]))
for i in 1:na
    println(lpad("i$i", 5), "  ", join([lpad(round(INVEST_C[i, d]; sigdigits=4), 12) for d in 1:nr]))
end

println("\n── the flagged cell [2,5] against its row and column ──")
println("INVEST_C[2,5] = ", INVEST_C[2, 5])
println("CAP_v[2,5]    = ", CAP_v[2, 5])
srt = sort(vec(INVEST_C))
println("INVEST_C: min=", srt[1], "  p05=", srt[max(1, length(srt) ÷ 20)],
        "  median=", srt[length(srt) ÷ 2], "  max=", srt[end])
rk = count(<(INVEST_C[2, 5]), vec(INVEST_C))
println("INVEST_C[2,5] rank = $rk of $(na*nr) (0 = smallest)")
println("ratio INVEST_C[2,5] / median = ",
        round(INVEST_C[2, 5] / max(srt[length(srt) ÷ 2], eps()); sigdigits=4))

# gret = log(pcap+tau) - log(pinvitot+tau): sensitivity is 1/(pinvitot+tau).
# With pinvitot a benchmark=1 index this is harmless UNLESS the cell is degenerate.
println("\n── smallest 10 INVEST_C cells ──")
for k in first(sortperm(vec(INVEST_C)), 10)
    i, d = Tuple(CartesianIndices(INVEST_C)[k])
    println("  [$i,$d]  INVEST_C=", rpad(round(INVEST_C[i, d]; sigdigits=5), 14),
            "  CAP_v=", round(CAP_v[i, d]; sigdigits=5))
end

println("\n── CAP_v/INVEST_C ratio (drives gret curvature) — largest 10 ──")
rat = [INVEST_C[i, d] > 1e-10 ? CAP_v[i, d] / INVEST_C[i, d] : Inf for i in 1:na, d in 1:nr]
for k in first(sortperm(vec(rat); rev=true), 10)
    i, d = Tuple(CartesianIndices(rat)[k])
    println("  [$i,$d]  CAP_v/INVEST_C=", rpad(round(rat[i, d]; sigdigits=5), 14),
            "  CAP_v=", rpad(round(CAP_v[i, d]; sigdigits=5), 12),
            "  INVEST_C=", round(INVEST_C[i, d]; sigdigits=5))
end

println("\n── CAP_v/INVEST_C ratio — SMALLEST 10 (the norm is ~2.0) ──")
for k in first(sortperm(vec(rat)), 10)
    i, d = Tuple(CartesianIndices(rat)[k])
    println("  [$i,$d]  CAP_v/INVEST_C=", rpad(round(rat[i, d]; sigdigits=5), 12),
            "  CAP_v=", rpad(round(CAP_v[i, d]; sigdigits=5), 12),
            "  INVEST_C=", round(INVEST_C[i, d]; sigdigits=5))
end
let r25 = rat[2, 5], v = sort(vec(rat))
    println("rat[2,5] = ", round(r25; sigdigits=5),
            "   rank ", count(<(r25), vec(rat)), " of ", length(rat),
            "   median=", round(v[length(v) ÷ 2]; sigdigits=5))
end

# The 52%-of-mass family is alab_o|pcap|plab_o|plnd|xprim — the primary-factor
# nest. Its curvature is set by the CAPITAL/LAND value shares in the cell, since
# xcap and xlnd are exogenous in the base static closure and the whole adjustment
# has to land in pcap/plnd. Print the shares for [2,5] against the distribution.
println("\n── primary-factor value shares ──")
LND = parent(params["LND"])          # (na, nr)
labtot = parent(params["LAB_O"])     # (na, nr): labour summed over occupations
if LND !== nothing
    tot = [CAP_v[i, d] + LND[i, d] + labtot[i, d] for i in 1:na, d in 1:nr]
    capsh = [tot[i, d] > 1e-10 ? CAP_v[i, d] / tot[i, d] : NaN for i in 1:na, d in 1:nr]
    fixsh = [tot[i, d] > 1e-10 ? (CAP_v[i, d] + LND[i, d]) / tot[i, d] : NaN
             for i in 1:na, d in 1:nr]
    println("cell [2,5]: CAP=", round(CAP_v[2, 5]; sigdigits=5),
            "  LND=", round(LND[2, 5]; sigdigits=5),
            "  LAB=", round(labtot[2, 5]; sigdigits=5))
    println("  capital share = ", round(capsh[2, 5]; sigdigits=4),
            "   fixed-factor (cap+lnd) share = ", round(fixsh[2, 5]; sigdigits=4))
    fv = sort(filter(isfinite, vec(fixsh)))
    println("  fixed-factor share across cells: min=", round(fv[1]; sigdigits=4),
            " median=", round(fv[length(fv) ÷ 2]; sigdigits=4),
            " max=", round(fv[end]; sigdigits=4))
    println("  rank of [2,5] = ", count(<(fixsh[2, 5]), fv), " of ", length(fv))
    println("\n  cells with the HIGHEST fixed-factor share (least able to adjust):")
    for k in first(sortperm(vec(fixsh); rev=true), 8)
        i, d = Tuple(CartesianIndices(fixsh)[k])
        isfinite(fixsh[i, d]) || continue
        println("    [$i,$d]  fixed share=", rpad(round(fixsh[i, d]; sigdigits=4), 10),
                "  CAP=", rpad(round(CAP_v[i, d]; sigdigits=5), 12),
                "  LND=", rpad(round(LND[i, d]; sigdigits=5), 12),
                "  LAB=", round(labtot[i, d]; sigdigits=5))
    end
else
    println("  (LND key absent; available keys matching cap/lnd/lab:)")
    println("  ", filter(k -> occursin(r"(?i)cap|lnd|lab", k), sort(collect(keys(params)))))
end

println("\nDone.")
