"""
Is the shock solve merely HARD, or is something actually WRONG?

diag_full_newton.jl + diag_cell25.jl localised 97% of the residual to ONE cell,
sector 2 / region 5 (Livestock / BaliNusa), which is the joint extreme of two
structural ratios: rank 0/150 on lowest CAP_v/INVEST_C (0.142 vs median 1.476)
and rank 6/150 on lowest fixed-factor share (0.172 vs median 0.457). With
xcap/xlnd exogenous, a small fixed-factor share forces a LARGE pcap swing to
absorb any adjustment — measured 8.6% off a 0.03% labour shock.

That is a hard nonlinear problem, not necessarily a wrong one. This script
separates the two readings:

  * If each graduated shock converges to `tol` and the last one reaches the SAME
    place as a direct 0.9997 solve, the equations are right and the remedy is
    continuation (planned as C6b) — not more solver tuning.
  * If even a 1e-6 shock stalls at ~1e-5, the residual floor is NOT a step-size
    problem and there is a genuine defect in the [2,5] block to find.

Each step starts from the previous solution, which IS the continuation method,
so a pass here doubles as the C6b prototype.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_continuation.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, LinearAlgebra

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

# SWAPS=1 runs TERM.CMF's six active swaps instead of the base static closure.
# This is not a cosmetic option: in the BASE closure `xcap` is exogenous, so the
# capital-accumulation block is inert and CAPSTOK/RNORMAL/GRETEXP never enter the
# residual at all. That is why the N.8 aggregation fix left this ladder
# bit-identical — same iteration counts, same pcap[2,5] to eight digits as the
# pre-fix run recorded in PLAN.md §I. The `xcap -> endogenous, faccum -> exogenous`
# swap is what makes those coefficients live, so the swapped ladder is the only one
# that can show whether the fix moved the fold. Run BOTH before drawing a
# conclusion about "the" fold; PLAN.md records two different ones (base 0.999752,
# swapped 0.99981).
const USE_SWAPS = get(ENV, "SWAPS", "0") != "0"
if USE_SWAPS
    apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))
    println("closure: TERM.CMF all six swaps")
else
    println("closure: base static")
end

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("\n══ benchmark ══\nstatus=$(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= 1e-8 "benchmark regressed"

sv = vars["blabnat"]
svr = sv isa AbstractArray ? first(sv) : sv

# Graduated: each step starts from the previous converged point.
# Default ladder brackets the gap where the coarse run broke: 0.9999 converged in
# 3 iterations, 0.9997 floored at 1.4e-5. If the branch has a FOLD (limit point)
# between them, the rungs converge with pcap[2,5] accelerating and then stop
# abruptly at a specific shock size; if it was merely too big a step, the fine
# ladder walks all the way to 0.9997.
LADDER = let s = get(ENV, "LADDER", "")
    isempty(s) ? (0.999999, 0.99999, 0.9999, 0.99985, 0.9998, 0.99975, 0.9997) :
        Tuple(parse.(Float64, split(s, ',')))
end
println("ladder: ", LADDER)
results = Tuple{Float64,Symbol,Int,Float64,Float64}[]
for s in LADDER
    JuMP.fix(svr, s; force=true)
    println("\n══ blabnat = $s ══")
    el = @elapsed r = solve_newton!(m, vars; maxit=40, verbose=false)
    push!(results, (s, r.status, r.iters, r.residual, el))
    println("status=$(r.status)  iters=$(r.iters)  ‖F‖∞=$(r.residual)  ($(round(el,digits=1))s)")
    println(r.residual <= 1e-8 ? "  ✅ converged" : "  ❌ did not converge")
    # Report the cell that carried the residual, so a partial pass is still legible.
    for nm in ("pcap", "gret", "xinvitot")
        haskey(vars, nm) || continue
        v = vars[nm]
        val = JuMP.is_fixed(v[2, 5]) ? JuMP.fix_value(v[2, 5]) : JuMP.start_value(v[2, 5])
        println("    $(rpad(nm,9))[2,5] = $(round(val; sigdigits=8))")
    end
    r.residual <= 1e-8 || break   # no point climbing further from a bad point
end

println("\n── ladder summary ──")
println("shock          status        iters   ‖F‖∞          secs")
for (s, st, it, res, el) in results
    println(rpad(s, 14), rpad(st, 13), rpad(it, 7), rpad(round(res; sigdigits=6), 13),
            round(el, digits=1))
end

allpass = length(results) == length(LADDER) && all(r -> r[4] <= 1e-8, results)
println("\n", allpass ?
    "✅ CONTINUATION WORKS — equations are right; the direct solve was a step-size problem." :
    """❌ ladder broke at rung $(length(results)+1) of $(length(LADDER)).
   This does NOT establish a fold. A FIXED ladder cannot distinguish "the branch is
   obstructed here" from "that rung was simply too long a step from the last one" —
   only an ADAPTIVE walk that drives h to hmin at one parameter value can. Read the
   amplification column first: if it is still linear at the last good rung, the step
   was too long (see PLAN.md §N.14, where 0.995 -> 0.99 failed with amplification
   still constant at 3.62). Re-measure with run_model!'s homotopy before concluding
   anything about the [2,5] block.""")
println("\nDone.")
