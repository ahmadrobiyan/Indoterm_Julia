"""
Economic results of the TERM.CMF reference simulation — and whether they make sense.

`test/verify_delunity.jl` established that the scenario *solves*. That is a numerical
statement. This script asks the separate and more important question: are the answers
economically coherent?

**The decomposition is the whole point.** The scenario has two shocks and they do
different things, so reporting only the combined endpoint would confound them:

    benchmark ──(delUnity 0→1)──► YEAR-FORWARD ──(blabnat 1→0.97)──► FULL SCENARIO
                 capital arrives                  3% productivity gain

  * A = year-forward vs benchmark — last year's investment landing as capital.
    Nothing to do with productivity.
  * B = full scenario vs benchmark — what TERM.CMF reports.
  * **B − A = the pure productivity effect**, which is what to sense-check against
    economic logic. Reading B alone would credit productivity with capital's work.

**What theory says a 3% labour-productivity gain should do** (write the predictions
down BEFORE looking, so the check is a test and not a rationalisation):

  1. Real GDP up. Same labour, more output per unit.
  2. Real wage up — competition for now-more-productive labour bids it up. In a
     TERM/ORANI-G labour market with the national supply mechanism switched on
     (`delfwage = flabsup_id`), that is where the gain mostly lands.
  3. Employment roughly flat, or mildly up. Productivity cuts labour needed per unit
     of output but output rises; the two partly cancel. A large employment FALL would
     be suspicious; a large rise would be too.
  4. Prices down — unit labour costs fall, so GDP deflator and CPI should soften.
  5. Exports up — lower costs improve competitiveness.
  6. Capital UNCHANGED between A and B. This is a sharp, falsifiable prediction: with
     `delUnity = 1` the accumulation equation pins `log(xcap/CAP) = CAPADD/CAPSTOK_OLDP`,
     a predetermined number that does not depend on `blabnat`. If capital moves between
     A and B, the accumulation block is wired wrong.

Prediction 6 is the most valuable line in this report, because it can only come out
right for the right reason.

Run:  julia --project=IndotermJulia IndotermJulia/test/report_reference_sim.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, Printf

const TOL = 1e-8

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))
apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
@assert r0.residual <= TOL "benchmark regressed"
println("benchmark: $(r0.status)  ‖F‖∞=$(r0.residual)")

na, nr = length(AGGCOM), length(REG6)
_val(vr) = JuMP.is_fixed(vr) ? JuMP.fix_value(vr) : JuMP.start_value(vr)

# Snapshot everything we want to report, at each of the three points.
function snap()
    s = Dict{String,Any}()
    for nm in ("MainMacro", "NatMacro", "xtot", "xlab_id", "realwage_id",
               "xcap", "delVGDPEXP", "delGDPINC", "xexp_i")
        haskey(vars, nm) || continue
        v = vars[nm]
        s[nm] = v isa AbstractArray ? map(_val, collect(v)) : _val(v)
    end
    return s
end

S0 = snap()

dU  = vars["delUnity"]
lam = let v = vars["blabnat"]; v isa AbstractArray ? first(v) : v end

println("\n══ STAGE A: delUnity 0 → 1 (year forward) ══")
rA = continuation_solve!(m, vars, dU, 1.0; h0=0.25, tol=TOL, maxit=30, verbose=false)
@assert rA.reached_target "stage A failed"
println("  reached, $(rA.nsteps) step(s), ‖F‖∞=$(rA.residual)")
SA = snap()

println("\n══ STAGE B: blabnat 1 → 0.97 (3% productivity gain) ══")
rB = arclength_solve!(m, vars, lam, 0.97; ds0=0.005, dsmin=1e-7, dsmax=0.2,
                      tol=TOL, maxit=20, maxsteps=200, verbose=false)
@assert rB.reached_target "stage B failed — λ=$(rB.lam_reached)"
println("  reached, $(rB.nsteps) step(s), ‖F‖∞=$(rB.residual)")
SB = snap()

# ── macro table ─────────────────────────────────────────────────────────────
# MainMacro[macro, region], region nr+1 = National. These are indices at 1.0, so
# a percentage change is (x-1)*100 — but report against the BASE snapshot rather
# than assuming 1.0, so a non-unit benchmark cannot masquerade as a result.
pc(now, base) = abs(base) > 1e-12 ? (now / base - 1) * 100 : NaN

println("\n" * "="^78)
println("NATIONAL MACRO — % change from benchmark")
println("="^78)
@printf("%-16s %10s %10s %10s   %s\n", "", "A: year fwd", "B: full", "B−A: prod", "")
println("-"^78)

MM0, MMA, MMB = S0["MainMacro"], SA["MainMacro"], SB["MainMacro"]
natcol = nr + 1
for nm in ("RealGDP", "RealHou", "RealInv", "RealGov", "ExpVol", "ImpVolUsed",
           "AggEmploy", "realwage_io", "AggCapStock", "GDPPI", "CPI", "NomGDPexp")
    k = findfirst(==(nm), MAINMACROS); k === nothing && continue
    a = pc(MMA[k, natcol], MM0[k, natcol])
    b = pc(MMB[k, natcol], MM0[k, natcol])
    @printf("%-16s %9.3f%% %9.3f%% %9.3f%%\n", nm, a, b, b - a)
end

# ── the falsifiable one ─────────────────────────────────────────────────────
println("\n" * "="^78)
println("PREDICTION 6 (falsifiable): capital must be IDENTICAL in A and B")
println("="^78)
if haskey(S0, "xcap")
    XA, XB = SA["xcap"], SB["xcap"]
    d = maximum(abs.(XB .- XA) ./ max.(abs.(XA), 1e-12))
    @printf("max |xcap_B - xcap_A| / |xcap_A| over %d cells = %.3e\n", length(XA), d)
    println(d < 1e-6 ?
        "  ✅ capital is predetermined by the accumulation equation, exactly as\n" *
        "     `delUnity=1` implies. blabnat cannot move it, and it did not." :
        "  ❌ capital MOVED with blabnat. Under delUnity=1 the accumulation equation\n" *
        "     should pin log(xcap/CAP) = CAPADD/CAPSTOK_OLDP, independent of blabnat.\n" *
        "     This points at a wiring defect in build_dynamics! Excerpt 50.")
end

# ── labour market detail ────────────────────────────────────────────────────
println("\n" * "="^78)
println("LABOUR MARKET by occupation — % change, B−A (pure productivity effect)")
println("="^78)
if haskey(S0, "xlab_id") && haskey(S0, "realwage_id")
    LA, LB = SA["xlab_id"], SB["xlab_id"]
    WA, WB = SA["realwage_id"], SB["realwage_id"]
    @printf("%-6s %14s %14s\n", "occ", "employment", "real wage")
    for o in eachindex(LA)
        @printf("%-6d %13.3f%% %13.3f%%\n", o, pc(LB[o], LA[o]), pc(WB[o], WA[o]))
    end
end

# ── regional ────────────────────────────────────────────────────────────────
println("\n" * "="^78)
println("REGIONAL nominal GDP (expenditure side) — % change from benchmark")
println("="^78)
if haskey(S0, "delVGDPEXP")
    GES = parent(params["GDPEXPSUM"])
    base = [sum(GES[d, :]) for d in 1:nr]
    ta = [sum(SA["delVGDPEXP"][d, :]) for d in 1:nr]
    tb = [sum(SB["delVGDPEXP"][d, :]) for d in 1:nr]
    @printf("%-14s %12s %10s %10s %10s\n", "region", "bmk Rp bn", "A", "B", "B−A")
    for d in 1:nr
        a = base[d] != 0 ? ta[d] / base[d] * 100 : NaN
        b = base[d] != 0 ? tb[d] / base[d] * 100 : NaN
        @printf("%-14s %12.0f %9.3f%% %9.3f%% %9.3f%%\n", REG6[d], base[d], a, b, b - a)
    end
    println("\nNOTE: the benchmark database's regional GDP identity is off by up to 36.1%")
    println("      (worst: MalukuPapua). Regional splits inherit that; the national")
    println("      total is the reliable figure.")
end

# ── sectoral winners/losers ─────────────────────────────────────────────────
println("\n" * "="^78)
println("SECTORAL output xtot — B−A (pure productivity effect), national sum")
println("="^78)
if haskey(S0, "xtot")
    XA, XB = SA["xtot"], SB["xtot"]
    ch = [(AGGCOM[i], pc(sum(XB[i, :]), sum(XA[i, :]))) for i in 1:na]
    srt = sort(filter(r -> isfinite(r[2]), ch); by = r -> -r[2])
    println("top 6 expanding:")
    for r in first(srt, 6);  @printf("   %-14s %8.3f%%\n", r[1], r[2]); end
    println("bottom 6:")
    for r in last(srt, 6);   @printf("   %-14s %8.3f%%\n", r[1], r[2]); end
end

println("\nDone.")
