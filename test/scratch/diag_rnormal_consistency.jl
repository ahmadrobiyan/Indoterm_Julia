"""
Is the steady-state identity GROSSRET == RNORMAL still true after aggregation?

`build_premod!.jl:131-139` constructs the capital stock as CAPSTOK = CAP/RNORMAL,
and says why in its own comment: it "makes the benchmark GROSSRET = CAP/CAPSTOK
equal RNORMAL exactly — the steady-state condition". So at the 185-sector level
the identity holds by construction.

Aggregation to 25 sectors then treats the two sides DIFFERENTLY:

  * CAPSTOK travels as the header "STOC", a value flow, and is SUMMED;
  * RNORMAL travels as the header "TARG", a rate, and is averaged by MAKE output
    (`aggregate_model!.jl:154`, `_wavg_nd`).

Averaging a rate rather than summing it is the right instinct, but MAKE output is
the wrong weight for THIS rate. The rate of return is a ratio of two flows, so the
only aggregate consistent with the summed components is

    RNORMAL_agg = sum(CAP_j) / sum(CAPSTOK_j)
                = sum(CAP_j) / sum(CAP_j / RNORMAL_j)

i.e. a CAP-weighted HARMONIC mean. An arithmetic mean weighted by output is a
different number whenever the rate varies within a group, and it is biased UPWARD
relative to the harmonic mean — the more so the wider the spread. Livestock is
exactly the bad case: its members run from 0.012162 (the minimum of all 185) up to
~0.066, and the low-rate member carries the largest stock by construction, since
CAPSTOK = CAP/RNORMAL divides by that small rate.

If the identity is broken, then in the aggregated model the capital block is
calibrated against a rate of return the capital stock does not actually earn, and
that inconsistency sits precisely in the block that `test/diag_targ_counterfactual.jl`
just showed drives the fold (removing the TARG anomaly collapsed the Livestock
pcap amplification from 1186.5x to 2.2x).

This is pure data inspection — no model is built, nothing is solved. It prints, for
every sector, the realised GROSSRET against the carried RNORMAL, so the question is
settled by reading the two columns rather than by argument.

NOTE ON DIRECTION: if the correct harmonic aggregate is LOWER than the arithmetic
one, then fixing the aggregation makes sector 2's rate of return smaller still, and
the fold WORSE — which would mean the stiffness is a genuine property of the source
database rather than an aggregation artifact. This test is therefore capable of
refuting the convenient answer, not just confirming it.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_rnormal_consistency.jl
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using Printf

agg6, params = cached_pipeline(6)

CAP     = parent(agg6["1CAP"])          # capital rental value flow, IND x REG
CAPSTOK = parent(agg6["STOC"])          # capital stock, IND x REG
RNORMAL = parent(agg6["TARG"])          # rate of return as carried by aggregation
DPRC    = parent(agg6["DPRC"])
na, nr  = size(CAP)

# ── The identity is per CELL, and the test must not repeat the sin it audits ──
# An earlier version of this script summed CAP and CAPSTOK over regions while
# taking an UNWEIGHTED mean of RNORMAL over the same axis — exactly the mistake it
# was written to detect. That manufactures a discrepancy out of nothing whenever a
# sector's rate varies across regions, so it cannot distinguish "the aggregation
# is wrong" from "my summary statistic is wrong". The identity actually asserted by
# `build_premod!.jl:131-139` is per (industry, region) cell, so check that first
# and treat it as the verdict; the by-sector table below is a readout, not the test.
println("\n══ per-CELL identity  GROSSRET[i,d] == RNORMAL[i,d] ══")
cellratios = Float64[]
wcell = Ref((0, 0, 0.0))
for i in 1:na, d in 1:nr
    s = CAPSTOK[i,d]; c = CAP[i,d]; rn = RNORMAL[i,d]
    (s > 1e-10 && rn > 1e-12) || continue
    q = (c / s) / rn
    push!(cellratios, q)
    q > 0 && abs(log(q)) > wcell[][3] && (wcell[] = (i, d, abs(log(q))))
end
let srt = sort(cellratios)
    @printf("  %d cells: min %.6f, median %.6f, max %.6f\n",
            length(srt), srt[1], srt[cld(length(srt), 2)], srt[end])
    @printf("  within 1e-8 of exact: %d;  within 1%%: %d;  worst cell [%d,%d]\n",
            count(x -> abs(x - 1) < 1e-8, cellratios),
            count(x -> abs(x - 1) < 0.01, cellratios),
            wcell[][1], wcell[][2])
end

println("\n══ steady-state identity, per sector (summed over the 6 regions) ══")
println("  (RNORMAL here is CAPSTOK-weighted across regions — an unweighted mean")
println("   would reintroduce the very inconsistency this script is testing for.)")
@printf("  %-4s %-14s %-14s %-12s %-12s %-9s %s\n",
        "sec", "CAP", "CAPSTOK", "GROSSRET", "RNORMAL", "ratio", "net of DPRC")
worst = Ref((0, 0.0))     # Ref, not a bare global: a top-level `for` opens a soft
                          # scope, so plain assignment would create a new local.
for i in 1:na
    c  = sum(@view CAP[i, :])
    s  = sum(@view CAPSTOK[i, :])
    gr = s > 1e-10 ? c / s : NaN
    kw = sum(@view CAPSTOK[i, :])
    rn = kw > 1e-10 ? sum(RNORMAL[i,d] * CAPSTOK[i,d] for d in 1:nr) / kw : 0.0
    ratio = (isfinite(gr) && rn > 1e-12) ? gr / rn : NaN
    isfinite(ratio) && ratio > 0 && abs(log(ratio)) > worst[][2] &&
        (worst[] = (i, abs(log(ratio))))
    @printf("  %-4d %-14.6g %-14.6g %-12.6g %-12.6g %-9.4f %+.4f%s\n",
            i, c, s, gr, rn, ratio, gr - DPRC[i,1],
            i == 2 ? "   <== Livestock" : "")
end

# How far off is the identity overall? If aggregation were consistent every ratio
# would be exactly 1.0000.
ratios = Float64[]
for i in 1:na
    c = sum(@view CAP[i, :]); s = sum(@view CAPSTOK[i, :])
    kw = sum(@view CAPSTOK[i, :])
    rn = kw > 1e-10 ? sum(RNORMAL[i,d] * CAPSTOK[i,d] for d in 1:nr) / kw : 0.0
    (s > 1e-10 && rn > 1e-12) && push!(ratios, (c / s) / rn)
end
srt = sort(ratios)
@printf("\n  GROSSRET/RNORMAL over %d sectors: min %.4f, median %.4f, max %.4f\n",
        length(ratios), srt[1], srt[cld(length(srt), 2)], srt[end])
@printf("  sectors within 1%% of the identity: %d of %d\n",
        count(x -> abs(x - 1) < 0.01, ratios), length(ratios))
@printf("  worst offender: sector %d\n", worst[][1])

broken = count(x -> abs(x - 1) > 0.01, ratios)
println("\n", broken == 0 ?
    "➡ The identity HOLDS. Aggregation preserves the steady-state condition, so the\n" *
    "  MAKE-output weighting of TARG is harmless here and the low Livestock rate of\n" *
    "  return is a faithful property of the source database — not an aggregation bug.\n" *
    "  The fold is then a DATA question for the modeller, not a translation defect." :
    "➡ The identity is BROKEN for $broken of $(length(ratios)) sectors. In the aggregated\n" *
    "  model the capital block is calibrated against a rate of return that the capital\n" *
    "  stock does not actually earn. Since CAPSTOK and CAP are both summed flows and\n" *
    "  are therefore mutually consistent, the inconsistent quantity is RNORMAL, and the\n" *
    "  fix is to DERIVE it after aggregation as sum(CAP)/sum(CAPSTOK) rather than\n" *
    "  carrying an output-weighted average of the disaggregated rates.\n" *
    "  Check the sign before concluding this helps: if the derived rate for sector 2 is\n" *
    "  LOWER than the carried 0.030321, the fix makes the fold worse, and the stiffness\n" *
    "  is a genuine property of the source data.")
println("\nDone.")
