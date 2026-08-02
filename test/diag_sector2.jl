"""
Why is sector 2 (Livestock) 10-20x stiffer than every other sector?

`test/diag_pcap_rank.jl` settled what the dominant mode actually is, and it was not
what either standing hypothesis predicted. Ranking every pcap[i,d] by relative
response to a 1e-6 shock under the full TERM.CMF closure:

    [2,3] 1187x   [2,1] 946x   [2,4] 313x   [2,5] 302x   [2,6] 212x   [2,2] 168x
    [3,3]   86x   [3,6]   72x   ... every other cell <= 86x

Two hypotheses died there:
  * It is NOT a [2,5] cell problem. [2,5] ranks 4th and [2,3] is four times larger;
    [2,5] was simply the cell `verify_swapped_closure.jl` happens to print, chosen
    from older residual history — circular evidence.
  * It is NOT dust capital. 0 of the top 20 movers have VCAP < 100, the dust cells
    VCAP[5,5]=0.040 / VCAP[5,6]=0.033 do not appear at all, and the top movers'
    median VCAP (5,012) is the same order as the grid median (10,016).

What is left is a SECTOR-WIDE property: all six regions of sector 2, uniformly.
That also explains why `test/diag_linearity.jl` found every mover sharing one
doubling ratio to within 0.9% — one sectoral mode, not a scattering of cells.

Sector 2's capital share is 0.1189, i.e. ordinary (`diag_cell25_scale.jl`), so the
factor mix is not the cause. Two candidates remain, both cheap to read straight
out of the calibrated data with no solve at all:

  (H1) NEAR-LEONTIEF primary-factor substitution. `TERM.TAB:136/159` declares
       SIGMAPRIM(i) from header "P028", used at TERM.TAB:475/479/483 and ported to
       `build_equations.jl:296` E_pcap!. As SIGMAPRIM -> 0 factor demand becomes
       Leontief, the rental supply curve turns vertical, and pcap must move a long
       way to clear a small quantity change. If SIGMAPRIM[2] is far below the other
       sectors, that alone explains the ranking.

  (H2) A STRONG SELF-LOOP in intermediate demand — livestock consuming feed, i.e.
       sector 2 buying its own output. A large diagonal USE share raises the
       Leontief multiplier (1-a)^-1 for that sector specifically, amplifying any
       shock that reaches it. If USE[2,2,d]/output is far above the other diagonals,
       that explains it instead.

They are not exclusive and they are distinguishable by inspection, so this prints
both: SIGMAPRIM for all 25 sectors with sector 2 marked, and the intermediate
self-loop share by sector. No model is built.

RESULT (both refuted):
  H1  SIGMAPRIM[2] = 0.2391, rank 6 of 25 from the BOTTOM, and sector 4 is lower
      still at 0.200 while responding ~20x less. Not near-Leontief, and not
      monotone in the elasticity anyway.
  H2  self share 0.0380, rank 21 of 25 from the TOP against a median of 0.2478.
      Livestock has one of the WEAKEST self-loops in the grid, not the strongest.
      (This branch previously probed a header named "USE", which does not exist,
      so it skipped silently and H2 was never actually tested. It is BSMR.)

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_sector2.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using Printf

agg6, params = cached_pipeline(6)

# `build_model!.jl:302` reads the elasticity as params["P028"], falling back to a
# flat 0.5 when the header is absent. Reproduce that exactly, and say which branch
# applied — a fallback would make the whole H1 question moot.
na = size(parent(agg6["1CAP"]), 1)
has028 = haskey(params, "P028")
SIGMAPRIM = has028 ? parent(params["P028"]) : fill(0.5, na)
println(has028 ? "\nSIGMAPRIM from header P028 (as build_model!.jl:302 reads it)" :
                 "\n⚠️  params has NO \"P028\" — build_model! falls back to a flat 0.5, so H1 is moot")

println("\n══ H1: primary-factor CES elasticity by sector ══")
for i in 1:na
    @printf("  SIGMAPRIM[%2d] = %-10.6g%s\n", i, SIGMAPRIM[i], i == 2 ? "   <== Livestock" : "")
end
srt = sort(collect(enumerate(SIGMAPRIM)); by = x -> x[2])
@printf("\n  smallest: sector %d at %.6g;  sector 2 at %.6g;  median %.6g\n",
        srt[1][1], srt[1][2], SIGMAPRIM[2], sort(collect(SIGMAPRIM))[cld(na, 2)])
rank2 = findfirst(x -> x[1] == 2, srt)
@printf("  sector 2 ranks %d of %d from the bottom\n", rank2, na)

# ── H2: intermediate self-loop ────────────────────────────────────────────────
# The intermediate-input flow is NOT called "USE" in this database — an earlier
# version of this script probed for that header, found nothing, and silently
# skipped H2 entirely, so H2 was never actually tested. The flow is `BSMR`, laid
# out COM×SRC×USR×DST (`build_reg1!.jl:348`), where USR is the na industries
# followed by 4 final users (`prepare_parameters.jl:256`: nu = na + 4). So the
# self-loop for industry j is the c=j, u=j slice summed over source and region,
# and its denominator is industry j's whole intermediate column — final-user
# columns are not industries and must stay out of both.
if !haskey(agg6, "BSMR")
    println("\n⚠️  no \"BSMR\" header in the aggregated database — H2 cannot be tested.")
else
    U = parent(agg6["BSMR"])
    @printf("\n══ H2: intermediate self-loop (BSMR is %s, COM×SRC×USR×DST) ══\n",
            join(size(U), "×"))
    ncom, _, nu, _ = size(U)
    n = min(ncom, na)
    selfflow = [sum(@view U[j, :, j, :]) for j in 1:n]
    colsum   = [sum(@view U[:, :, j, :]) for j in 1:n]
    diagshare = [colsum[j] > 0 ? selfflow[j] / colsum[j] : 0.0 for j in 1:n]
    @printf("  (%d of the %d USR columns are industries; the other %d are final users)\n",
            n, nu, nu - na)
    @printf("  %-8s %-16s %-16s %s\n", "sector", "BSMR[j,·,j,·]", "total inputs", "self share")
    for j in 1:n
        @printf("  %-8d %-16.6g %-16.6g %.4f%s\n", j, selfflow[j], colsum[j], diagshare[j],
                j == 2 ? "   <== Livestock" : "")
    end
    s = sort(collect(enumerate(diagshare)); by = x -> -x[2])
    @printf("\n  largest self-share: sector %d at %.4f;  sector 2 at %.4f;  median %.4f\n",
            s[1][1], s[1][2], diagshare[2], sort(diagshare)[cld(n, 2)])
    @printf("  sector 2 ranks %d of %d from the top\n",
            findfirst(x -> x[1] == 2, s), n)
end

println("""

Read this as:
  * SIGMAPRIM[2] far below the others  -> H1: near-Leontief factor demand. The
    rental supply curve is nearly vertical, so pcap[2,*] must swing hard to clear.
    That is a CALIBRATION input, so check it against TERM.CMF / the source header
    before treating the fold as real — a mistyped or defaulted elasticity would
    produce exactly this signature.
  * A large sector-2 self share      -> H2: Leontief multiplier amplification.
  * Neither unusual                  -> the stiffness is emergent, the limit point
    is a genuine property of the calibrated model, and pseudo-arclength is the
    right next build.
""")
println("Done.")
