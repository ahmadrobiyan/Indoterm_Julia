"""
Is the fold at t* ≈ 0.307 the Klein-Rubin SUBSISTENCE FLOOR in region 4?

PLAN.md §N.16 established — by direct tangent measurement, not inference — that the
obstruction is a genuine fold and that its direction is region 4's household block:
`whouhtot[4]` carries ‖dx/dλ‖∞, with `wlux[4]` third behind a reporting aggregate.

That names a block. It does not yet name a MECHANISM. The ELES/Stone-Geary
structure in build_equations.jl Excerpt 13 offers an obvious candidate:

    E_xsub!         xsub[c,d]  == XSUB0[c,d] * nhou[d] * asub[c,d]     (INELASTIC)
    E_xlux!         xlux[c,d] * phou[c,d] == wlux[d] * alux[c,d]
    E_xhouh_s_agg!  xhou_s[c,d] == xlux[c,d] + xsub[c,d]

Subsistence demand does not respond to price or income at all (ahou_s is Omitted,
so asub is constant). ALL household adjustment therefore has to come out of the
supernumerary part, whose scale is the single per-region scalar `wlux[d]`. If a
region's income falls far enough that `wlux[d]` approaches zero, that margin is
exhausted: xlux → 0, the household block loses its only degree of freedom, and the
system goes singular. That is a real economic limit point, not a translation bug.

PREDICTION, stated before the measurement so it can be wrong:
  * If the hypothesis holds, `wlux[4]/WLUX0[4]` falls steeply and heads for 0
    around t ≈ 0.307, and region 4 has the LOWEST supernumerary share of the six.
  * If the hypothesis fails, `wlux[4]` is still comfortably positive at the last
    converged point and something else in the block is singular — in which case
    read the `xlux` minimum and the per-commodity table below instead.

Both outcomes are informative, so print the full six-region table either way and
rank every region. Do NOT read a localisation off a single printed cell — that
mistake has now been made four times (§N.15).

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_wlux_floor.jl
Env:  TS=0.05,0.15,0.25,0.30,0.30625   (defaults to the §N.16 points)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, LinearAlgebra, Printf

const TS = let s = get(ENV, "TS", "")
    isempty(s) ? [0.0, 0.05, 0.15, 0.25, 0.30, 0.30625] : parse.(Float64, split(s, ','))
end
const SPAN = -0.03          # blabnat: 1.0 → 0.97

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))
apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("benchmark: $(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= 1e-8 "benchmark regressed"

const na, nr = size(vars["xhou_s"])
const WLUX0 = parent(params["WLUX0"])          # 1 × nr
const XSUB0 = parent(params["XSUB0"])          # na × 1 × nr
const XLUX0 = parent(params["XLUX0"])
const HOUPUR_C = parent(params["HOUPUR_C"])    # 1 × nr

val(v) = JuMP.is_fixed(v) ? JuMP.fix_value(v) :
         (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v))

# ── Benchmark structure: which region has the least supernumerary margin? ────
# This is a pure parameter fact — known before any solve — so it is the cheapest
# possible test of the hypothesis's premise.
println("\n" * "="^78)
println("BENCHMARK ELES STRUCTURE  (supernumerary share = WLUX0 / HOUPUR_C)")
println("="^78)
sub0 = [sum(XSUB0[c,1,d] for c in 1:na) for d in 1:nr]
lux0 = [WLUX0[1,d] for d in 1:nr]
shr0 = [lux0[d] / HOUPUR_C[1,d] for d in 1:nr]
println("region   HOUPUR_C        WLUX0(supernum)   subsistence     supernum share")
for d in 1:nr
    @printf("  %-6d %-15.6g %-17.6g %-15.6g %.4f%s\n", d, HOUPUR_C[1,d], lux0[d], sub0[d],
            shr0[d], d == argmin(shr0) ? "   ← LOWEST" : "")
end
println("\nregions ranked by supernumerary share (thinnest margin first): ",
        join(sortperm(shr0), " < "))
println(argmin(shr0) == 4 ?
    "  ⇒ region 4 IS the thinnest — consistent with the hypothesis's premise." :
    "  ⇒ region 4 is NOT the thinnest (region $(argmin(shr0)) is). The fold direction is\n" *
    "     region 4 anyway, so a pure benchmark-share story is already insufficient —\n" *
    "     what matters is how far each wlux MOVES, measured below.")

# ── Walk the branch ─────────────────────────────────────────────────────────
println("\n" * "="^78)
println("wlux(d) ALONG THE BRANCH   (reported as wlux/WLUX0, benchmark = 1)")
println("="^78)
rows = Tuple{Float64,Vector{Float64},Float64,Tuple{Int,Int}}[]
for tt in TS
    JuMP.fix(vars["blabnat"] isa AbstractArray ? first(vars["blabnat"]) : vars["blabnat"],
             1.0 + tt * SPAN; force=true)
    r = solve_newton!(m, vars; maxit=40, verbose=false)
    if r.residual > 1e-8
        println("\nt=$tt  ❌ did not converge (‖F‖∞=$(r.residual)) — stopping")
        break
    end
    w = [val(vars["wlux"][d]) / WLUX0[1,d] for d in 1:nr]
    xl = [val(vars["xlux"][c,d]) for c in 1:na, d in 1:nr]
    xs = [val(vars["xsub"][c,d]) for c in 1:na, d in 1:nr]
    # smallest luxury quantity relative to its own benchmark — the cell that runs
    # out of adjustment margin FIRST, wherever it is
    rel = [XLUX0[c,1,d] > 0 ? xl[c,d] / XLUX0[c,1,d] : Inf for c in 1:na, d in 1:nr]
    k = argmin(rel)
    push!(rows, (tt, w, rel[k], (k[1], k[2])))
    @printf("\n── t=%-8g blabnat=%-11.7g conv %d it, ‖F‖∞=%.3g\n", tt, 1.0 + tt*SPAN, r.iters, r.residual)
    print("   wlux/WLUX0 :")
    for d in 1:nr; @printf(" [%d]=%.6f", d, w[d]); end
    println()
    @printf("   min xlux/XLUX0 = %.6f at (c=%d, d=%d);  its xsub/XSUB0 = %.6f\n",
            rel[k], k[1], k[2], XSUB0[k[1],1,k[2]] > 0 ? xs[k[1],k[2]]/XSUB0[k[1],1,k[2]] : NaN)
end

# ── Verdict ─────────────────────────────────────────────────────────────────
println("\n" * "="^78)
println("VERDICT")
println("="^78)
if length(rows) < 2
    println("Not enough converged points to judge.")
else
    print(rpad("t", 11)); for d in 1:nr; print(rpad("wlux[$d]", 12)); end
    println(rpad("min xlux/XLUX0", 16), "cell")
    for (tt, w, mn, cell) in rows
        print(rpad(tt, 11))
        for d in 1:nr; print(rpad(round(w[d]; digits=6), 12)); end
        println(rpad(round(mn; digits=6), 16), cell)
    end

    w_last = rows[end][2]; w_prev = rows[end-1][2]
    dt = rows[end][1] - rows[end-1][1]
    println("\nlocal slope d(wlux/WLUX0)/dt over the last step (dt=$(round(dt;sigdigits=3))):")
    slopes = (w_last .- w_prev) ./ dt
    for d in 1:nr
        # linear extrapolation of THIS region's wlux to zero
        z = slopes[d] < 0 ? rows[end][1] + w_last[d] / (-slopes[d]) : Inf
        @printf("   region %d: %+10.5f   ⇒ hits zero at t ≈ %s\n", d, slopes[d],
                isfinite(z) ? string(round(z; sigdigits=5)) : "never (rising)")
    end

    println("""

    HOW TO READ THIS  (t* ≈ 0.307–0.318 from §N.16)
      * A region whose wlux extrapolates to zero AT t* — and region 4 in particular —
        confirms the subsistence-floor mechanism: households there have exhausted the
        only margin the ELES gives them, and the fold is economics, not a bug.
      * wlux[4] still well above zero with a mild slope REFUTES it. The fold is then
        somewhere else in the whouhtot[4] block (candidates: the houslack/shrBoTnom2
        swap, which pins ONE aggregate for all six regions, or the fhou/xhouhtot swap
        that makes household spending endogenous in the first place). Do not invent a
        replacement story from this table — measure the next one.
      * Either way the REMEDY is unchanged: pseudo-arclength continuation to turn the
        fold. Confirming the mechanism only tells us whether the turn is economically
        meaningful or whether a closure choice put the wall there.""")
end
println("\nDone.")
