"""
Is the `pcap` homogeneity overshoot solver noise, or a real degree-of-homogeneity
error?  (Task #38)

`verify_homogeneity.jl` reports, at λ = 1.10, that two `pcap` cells scale to
1.10041 instead of 1.10 — a 4.1e-4 relative overshoot — and that the failures in
`xinvi` / `xinvitot` / `gret` / `ggro` are all downstream of it
(`gret = log(pcap) − log(pinvitot)`).  The two cells are (5,5) and (5,6), which
carry the SMALLEST capital stocks in the model: CAP = 0.0403 and 0.0328 against a
database whose flows run to 1e6.  Small denominators are exactly where a fixed
absolute residual shows up as a large relative error, so "tolerance artifact" and
"wrong price degree" both predict the observation.

THE DISCRIMINATOR IS GENERIC AND NEEDS NO PER-VARIABLE SCALE.  Solve the same λ
point to a ladder of tolerances and watch the overshoot:

  * if the overshoot FOLLOWS the residual down, it was never economics — the
    homogeneity property holds and the number was measuring the solver;
  * if it PLATEAUS at ~4e-4 while the residual falls a decade or more, the model
    genuinely violates price homogeneity in these cells and it is a translation
    defect that must be found.

A second, independent axis: repeat at a different λ.  A real degree error scales
with (λ − 1); solver noise does not care what λ is.  Both axes are reported, and
they can disagree — if they do, say so rather than picking the convenient one.

Note this script does NOT re-derive the λ point by continuation.  Prices at λ are
known analytically (×λ, quantities unchanged), so it warm-starts there, exactly
as the rewritten `verify_homogeneity.jl` does.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_homogeneity_tol.jl
Env:  LAMBDAS (default "1.05,1.10,1.20"), TOLS (default "1e-6,1e-7,1e-8,1e-9")
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, Printf

const LAMBDAS = parse.(Float64, split(get(ENV, "LAMBDAS", "1.05,1.10,1.20"), ","))
const TOLS    = parse.(Float64, split(get(ENV, "TOLS", "1e-6,1e-7,1e-8,1e-9"), ","))

# The cells verify_homogeneity.jl named, plus one healthy control. Without the
# control a ladder that improves everything equally would look like a targeted
# result; the control is what makes a difference at (5,5)/(5,6) mean something.
const WATCH = [(5, 5), (5, 6), (1, 2)]     # (industry, region); (1,2) = control

# BOTH shape tests below are blind to magnitude, and that is not a detail: the
# control cell PLATEAUS and SCALES WITH λ exactly like the failing ones, because
# accumulated round-off in a price is itself proportional to the price. On the
# first run this script therefore printed "REAL" for the control — a false
# positive that would have sent the investigation after a defect that is not
# there. Shape alone cannot separate a 4e-4 economic error from a 1e-9 numerical
# one; only size can. A converged solve here drives equilibrated rows to ~1e-9
# absolute and pcap is O(1), so relative errors up to ~1e-7 are what a correct
# model looks like at this tolerance and mean nothing economically. Cells at or
# below that are exempted from the verdict regardless of how they scale.
const HEALTHY_FLOOR = 1e-7

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

r0 = solve_newton!(m, vars; maxit=12, tol=1e-9, verbose=false)
@printf("benchmark: status=%s  ‖F‖∞=%.4g\n", r0.status, r0.residual)
@assert r0.residual <= 1e-8 "benchmark regressed — nothing below is readable"

_val(vr) = JuMP.is_fixed(vr) ? JuMP.fix_value(vr) : JuMP.start_value(vr)

num = vars["NatMacro"][findfirst(==("GDPPI"), MAINMACROS)]
freevars = [vr for (_, v) in vars for vr in (v isa AbstractArray ? vec(v) : [v])
            if !JuMP.is_fixed(vr)]
base = Dict(vr => _val(vr) for vr in freevars)

# Benchmark values of the watched cells, read once. Everything below is a ratio
# against these, so a cell that is zero at the benchmark is unusable and is
# reported as such rather than dividing by it.
pcap = vars["pcap"]
base_pcap = Dict(c => _val(pcap[c...]) for c in WATCH)
for (c, b) in base_pcap
    @printf("base pcap%s = %.10g%s\n", c, b,
            abs(b) < 1e-12 ? "   ⚠️ ZERO — ratio at this cell is meaningless" : "")
end

function warm_start!(lam::Float64)
    for vr in freevars
        b = get(base, vr, NaN); isfinite(b) && JuMP.set_start_value(vr, b)
    end
    for (nm, v) in vars
        startswith(nm, "p") || continue
        for vr in (v isa AbstractArray ? v : (v,))
            JuMP.is_fixed(vr) && continue
            b = get(base, vr, NaN); isfinite(b) && JuMP.set_start_value(vr, b * lam)
        end
    end
end

# rows[(cell, λ)] = vector of (tol, residual, overshoot) across the ladder
rows = Dict{Tuple{Tuple{Int,Int},Float64},Vector{NTuple{3,Float64}}}()

for lam in LAMBDAS
    println("\n" * "="^74)
    @printf("λ = %.4g\n", lam)
    println("="^74)
    @printf("%-10s %-12s %-6s   %s\n", "tol", "‖F‖∞", "iters",
            join([@sprintf("pcap%s overshoot", c) for c in WATCH], "   "))

    for tol in TOLS
        JuMP.fix(num, lam; force=true)
        # Reach every (λ, tol) point from the analytic prediction, NOT from the
        # previous ladder rung. Warm-starting off the previous rung would make
        # each measurement depend on the one before it, and a plateau could then
        # be nothing but the solver declining to move from where it already was.
        warm_start!(lam)
        r = solve_newton!(m, vars; maxit=40, tol=tol, verbose=false)

        overs = Float64[]
        for c in WATCH
            b = base_pcap[c]
            o = abs(b) < 1e-12 ? NaN : _val(pcap[c...]) / b / lam - 1.0
            push!(overs, o)
            push!(get!(rows, (c, lam), NTuple{3,Float64}[]), (tol, r.residual, o))
        end
        @printf("%-10.0e %-12.4g %-6d   %s\n", tol, r.residual, r.iters,
                join([@sprintf("%+.3e", o) for o in overs], "   "))
    end
end

JuMP.fix(num, 1.0; force=true)

println("\n" * "="^74)
println("READING")
println("="^74)

for c in WATCH
    lab = c == (1, 2) ? " (control)" : ""
    println("\npcap$c$lab")

    # Magnitude gate first — see HEALTHY_FLOOR. If the largest overshoot this cell
    # ever showed is at the numerical floor, the shape tests below are describing
    # round-off and must not be allowed to call it a defect.
    allov = [abs(x[3]) for lam in LAMBDAS for x in get(rows, (c, lam), NTuple{3,Float64}[])]
    peak = isempty(allov) ? NaN : maximum(allov)
    if isfinite(peak) && peak <= HEALTHY_FLOOR
        @printf("  HOMOGENEOUS. Largest overshoot over the whole sweep is %.3g, at or below\n", peak)
        @printf("  the %.0e numerical floor — this cell scales correctly. (It still plateaus\n", HEALTHY_FLOOR)
        println("  and still tracks λ, because round-off is proportional to price; that is why")
        println("  the shape tests are not consulted here.)")
        continue
    end
    isfinite(peak) && @printf("  peak overshoot %.3g — above the %.0e floor, so the shape tests below decide.\n",
                              peak, HEALTHY_FLOOR)

    for lam in LAMBDAS
        v = get(rows, (c, lam), NTuple{3,Float64}[])
        isempty(v) && continue
        res = [x[2] for x in v]; ov = [abs(x[3]) for x in v]
        # Compare the span of the residual against the span of the overshoot over
        # the same ladder. Tracking means the error was the solver's; a flat
        # overshoot against a falling residual means it was not.
        dres = maximum(res) / max(minimum(res), eps())
        dov  = maximum(ov)  / max(minimum(ov),  eps())
        verdict = !all(isfinite, ov)      ? "cell unusable (zero benchmark)" :
                  dres < 3                ? "ladder did not move the residual — inconclusive" :
                  dov  > 0.3 * dres       ? "TRACKS the residual ⇒ solver artifact" :
                  dov  < 3                ? "PLATEAU while residual fell $(round(Int,dres))× ⇒ REAL" :
                                            "partial tracking — inconclusive, tighten further"
        @printf("  λ=%-6.4g residual fell %8.3g×   |overshoot| fell %8.3g×   → %s\n",
                lam, dres, dov, verdict)
    end
    # The λ axis: a genuine degree-of-homogeneity error scales with (λ−1).
    tight = Dict(lam => last(get(rows, (c, lam), [(0.0, 0.0, NaN)]))[3] for lam in LAMBDAS)
    if length(LAMBDAS) >= 2 && all(isfinite, values(tight))
        l1, l2 = LAMBDAS[1], LAMBDAS[end]
        expected = (l2 - 1) / (l1 - 1)
        got = abs(tight[l1]) < eps() ? NaN : tight[l2] / tight[l1]
        @printf("  λ-scaling at the tightest tol: got %.4g, a real degree error predicts %.4g",
                got, expected)
        println(isfinite(got) && abs(got / expected - 1) < 0.3 ?
                "  → scales with (λ−1) ⇒ REAL" :
                "  → does NOT scale with (λ−1) ⇒ not a degree error")
    end
end

println("""

Both axes must agree before this closes task #38. If the tolerance ladder says
"artifact" and the λ axis says "real", neither answer has been earned — report
the disagreement and look for a third discriminator rather than picking one.""")
println("\nDone.")
