"""
Does ALL SIX converge to an ECONOMICALLY MEANINGFUL point, or merely a feasible one?

`test/diag_allsix_deep.jl` showed ALL SIX converges in 2 Newton iterations at a
1e-6 shock — so the closure is NOT singular. But the trace carries a number that
convergence alone does not excuse:

    it 1  max|d| = 36.61      (shock = 1e-6)
    it 2  max|d| =  2.70

A Newton step of 36.6 in response to a shock of 1e-6 is an amplification of
~4e7, and most variables in this port are log-gaps or log-levels, where a value
near 39 means e^39. "Converged" says the residual is small; it says nothing about
whether the point is economics. A solver can land exactly on a mathematically
valid equilibrium that no economy could occupy.

This names what moved. Every free variable is compared against its benchmark
value and ranked by ABSOLUTE move — deliberately not by |d/x|, which is the trap
that produced the phantom "574x runaway" earlier in this investigation: additive
variables legitimately sit at 0 (down to 1e-27), and dividing by that manufactures
enormous numbers out of nothing.

Expected if the near-null direction is the `houslack` one: the movers are
`houslack` and the `fhou2[d]` block, moving together and in opposite sign, since
that is precisely the direction the Jacobian barely resolves. If instead the
movers are ordinary economic quantities, the amplification is a different story
and the houslack explanation does not carry over from the singleton.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_allsix_movers.jl
Env:  SHOCK (default 0.999999), TOPN (default 25)
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, Printf

SHOCK = parse(Float64, get(ENV, "SHOCK", "0.999999"))
TOPN  = parse(Int,     get(ENV, "TOPN", "25"))

agg6, params = cached_pipeline(6)
bmk = benchmark_levels(params)

(m, vars) = build_model_full!(agg6, params)
initialize_model!(m, vars; bmk_levels=bmk)
apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("\nbenchmark: status=$(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= 1e-8 "benchmark regressed"

# Snapshot AFTER the benchmark solve, so the comparison is shock-vs-benchmark and
# not shock-vs-seed.
base = Dict{JuMP.VariableRef,Float64}()
for v in JuMP.all_variables(m)
    val = JuMP.is_fixed(v) ? JuMP.fix_value(v) : JuMP.start_value(v)
    base[v] = val === nothing ? NaN : val
end

sv = vars["blabnat"]
svr = sv isa AbstractArray ? first(sv) : sv
JuMP.fix(svr, SHOCK; force=true)
r = solve_newton!(m, vars; maxit=30, verbose=false)
println("shock blabnat=$SHOCK: status=$(r.status)  iters=$(r.iters)  ‖F‖∞=$(r.residual)")

movers = Tuple{String,Float64,Float64,Float64}[]   # name, base, now, |move|
for (nm, v) in vars
    for vr in (v isa AbstractArray ? vec(collect(v)) : (v,))
        b = get(base, vr, NaN)
        isfinite(b) || continue
        JuMP.is_fixed(vr) && continue
        now = JuMP.start_value(vr)
        now === nothing && continue
        push!(movers, (JuMP.name(vr), b, now, abs(now - b)))
    end
end
sort!(movers; by = x -> -x[4])

@printf("\n── top %d movers by ABSOLUTE change (not |d/x| — that metric lies on additive vars) ──\n", TOPN)
@printf("  %-30s %-16s %-16s %s\n", "variable", "benchmark", "shocked", "|move|")
for (nm, b, now, d) in first(movers, TOPN)
    @printf("  %-30s %-16.8g %-16.8g %.6g\n", nm, b, now, d)
end

big = count(x -> x[4] > 1.0, movers)
@printf("\n  %d of %d free variables moved by more than 1.0 in absolute terms\n",
        big, length(movers))

for probe in ("houslack", "fhou2[1]", "fhou2[2]", "natfhou")
    k = findfirst(x -> x[1] == probe, movers)
    k === nothing ? println("  $probe: fixed or absent") :
        @printf("  %-12s benchmark=%-14.8g shocked=%-14.8g move=%.6g\n",
                movers[k][1], movers[k][2], movers[k][3], movers[k][4])
end

println("\n", big == 0 ?
    "✅ No free variable moved by more than 1.0 — the large Newton STEP was transient\n" *
    "  (an overshoot the next iteration undid), not a large displacement. The closure\n" *
    "  lands somewhere sane." :
    "⚠️  $big variable(s) moved by more than 1.0 for a 1e-6 shock. Check the names above\n" *
    "  before treating this closure as usable: on log-gap variables a move of that size\n" *
    "  is not a small perturbation, and 'converged' does not make it economics.")
println("\nDone.")
