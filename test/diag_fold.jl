"""
Is the obstruction at blabnat ~ 0.99975 a FOLD, and is it an artifact of the closure?

Evidence so far (test/verify_continuation.jl, fine ladder):

    blabnat   status      pcap[2,5]     slope d(pcap dev)/d(shock)
    0.9999    converged   0.98267596
    0.99985   converged   0.97186606      216   (over 1e-4 -> 1.5e-4)
    0.9998    converged   0.95785452      280   (over 1.5e-4 -> 2e-4)
    0.99975   FAILED      0.92623464
    0.9997    FAILED      0.9263198

Near a fold the response goes like A - B*sqrt(d* - d), so the slope goes like
1/sqrt(d* - d). Two slopes determine both unknowns:
    280/216 = sqrt((d*-1.25e-4)/(d*-1.75e-4))  =>  d* = 2.485e-4, i.e. blabnat* = 0.999752
and the same fit gives pcap[2,5]* = 0.9249. The ladder broke at 0.99975 (a location NOT
used in the fit) and both failed runs park at pcap ~ 0.9262 REGARDLESS of the target
shock. Stalling at the same point no matter where you aim is what a turning point does.

Two things this script measures, because the fit alone is not proof:

(1) WHERE the branch actually ends. `continuation_solve!` halves its step on failure,
    so `s_reached` brackets the obstruction automatically and far more finely than a
    hand-written ladder. A fold predicts s_reached ~ 0.99975 no matter how small h gets.

(2) WHETHER IT IS ONE-SIDED. This is the economically decisive test. A fold at a 0.025%
    labour shock is not a property of the Indonesian economy — it is a property of a
    closure in which xcap AND xlnd are exogenous in all 25 sectors, so fixed-factor
    PRICES absorb every adjustment, with ~200x leverage in the cell whose fixed-factor
    share is smallest (0.172). If the equal-and-opposite shock (blabnat = 1.00025, MORE
    labour) walks straight through, the branch is genuinely one-sided and the model is
    usable in that direction; if it folds symmetrically at the same distance, the
    obstruction is a generic property of the closure, not of the contraction.

    Either answer redirects the next step, which is why this runs BEFORE any
    pseudo-arclength implementation. Arclength can trace a branch around a fold, but it
    cannot manufacture an equilibrium that does not exist beyond it.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_fold.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP

agg6, params = cached_pipeline(6)
bmk = benchmark_levels(params)

# maxit is deliberately LOW. A converging step near the fold took 3-4 iterations, so
# 12 is generous for a success, while a failure is what dominates wall-clock (the 40-it
# failures above cost 740-880s each). Adaptive halving needs several failures to
# bracket the wall; each must be cheap.
const MAXIT   = 12
const RESULTS = Tuple{Float64,Bool,Float64,Int,Int}[]   # target, ok, s_reached, steps, rejects

for target in (0.9997, 1.0003)
    println("\n", "="^72)
    println("══ direction: blabnat → $target  ($(target < 1 ? "LESS" : "MORE") labour) ══")
    println("="^72)

    # A FRESH MODEL per direction, not just a re-`initialize_model!`. `solve_newton!`
    # mutates the model — it pins orphan variables and DELETES dead rows — so once a
    # solve has run, re-initialising only resets start values and closure fixes; the
    # pinned orphans stay pinned at whatever the previous (shocked) solve left them at.
    # Reusing the model that way gave `no_progress` with ‖F‖∞ = NaN on the second
    # direction. The build costs ~156s; correctness is worth it.
    t = @elapsed ((m, vars) = build_model_full!(agg6, params))
    println("build_model_full!: $(round(t, digits=1))s")
    initialize_model!(m, vars; bmk_levels=bmk)
    r0 = solve_newton!(m, vars; maxit=5, verbose=false)
    println("benchmark: status=$(r0.status)  ‖F‖∞=$(r0.residual)")
    @assert r0.residual <= 1e-8 "benchmark regressed"

    sv = vars["blabnat"]
    svr = sv isa AbstractArray ? first(sv) : sv

    el = @elapsed cr = continuation_solve!(m, vars, svr, target;
                                           h0=0.15, hmax=0.5, hmin=1e-2,
                                           grow_iters=3, maxit=MAXIT,
                                           report=(("pcap", 2, 5),))
    push!(RESULTS, (target, cr.reached_target, cr.s_reached, cr.nsteps, cr.nrejects))

    println("\n  reached_target=$(cr.reached_target)  s_reached=$(cr.s_reached)")
    println("  steps=$(cr.nsteps)  rejects=$(cr.nrejects)  ($(round(el, digits=1))s)")
    println("\n  ── accepted path (shock, iters, ‖F‖∞) ──")
    for (s, it, res) in cr.path
        println("    ", rpad(round(s; sigdigits=10), 16), rpad(it, 6), round(res; sigdigits=6))
    end
end

println("\n", "="^72)
println("── verdict ──")
for (tg, ok, sr, ns, nr) in RESULTS
    dist = abs(sr - 1.0)
    println("  target $(rpad(tg,10)) reached=$(rpad(ok,6)) s_reached=$(rpad(round(sr;sigdigits=10),16)) " *
            "|Δ|=$(round(dist; sigdigits=4))  steps=$ns rejects=$nr")
end

if length(RESULTS) == 2
    down, up = RESULTS[1], RESULTS[2]
    if !down[2] && up[2]
        println("\n➡ ONE-SIDED fold: the contraction branch ends at $(down[3]), the expansion " *
                "direction runs clean to $(up[1]).\n  The model is usable for labour-INCREASE " *
                "scenarios now. Contraction beyond the fold has NO equilibrium in this closure — " *
                "that is a closure question (xcap/xlnd both exogenous), not a solver question.")
    elseif !down[2] && !up[2]
        println("\n➡ SYMMETRIC obstruction (|Δ| $(round(abs(down[3]-1);sigdigits=4)) vs " *
                "$(round(abs(up[3]-1);sigdigits=4))). A limit point in BOTH directions at " *
                "comparable distance is a generic property of the closure, not of the shock sign. " *
                "Re-examine the short-run closure before investing in arclength.")
    elseif down[2] && up[2]
        println("\n➡ BOTH directions reached target — the earlier ladder break was step-size " *
                "after all, and the adaptive driver handles it. No fold.")
    else
        println("\n➡ Contraction succeeded where expansion failed — unexpected; inspect the paths.")
    end
end
println("\nDone.")
