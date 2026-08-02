"""
Does TERM.CMF's real closure remove the fold — and can it run TERM.CMF's real scenario?

PLAN.md §J: the port applies only ONE of TERM.CMF's seven active swaps, so every
shock so far has been run against TABmate's bare automatic closure — the state
before any modelling choice. With `xcap` AND `xlnd` frozen in all 25 sectors,
fixed-factor prices carry the entire adjustment (~200x amplification in the cell
whose fixed-factor share is smallest), and the solution branch folds at
`blabnat` = 0.99975, i.e. a 0.025% shock.

`TERM.CMF:112` shocks `blabnat = -3`. That is 120x beyond the fold. It gets there
because `TERM.CMF:88` swaps `xcap` endogenous — capital accumulation on.

Three things are checked here, in order, because each is meaningless if the
previous one failed:

  1. BENCHMARK PRESERVED. `apply_swaps!` fixes each newly-exogenous variable at
     its benchmark seed, so the benchmark must remain an exact solution after the
     swap. If ‖F‖∞ jumps here, the swap is wrong and nothing below counts. This
     is the regression gate every closure change has to pass first.

  2. FOLD REMOVED. Walk to 0.9997 — just past where the base closure folds. Under
     a closure that lets capital move, this should be uneventful.

  3. THE REAL SCENARIO. Walk to TARGET (default 0.97 = TERM.CMF's -3%). Reaching
     it means the port can run the source's own simulation.

A failure at (2) or (3) is still a measurement, not a dead end: `s_reached` is the
furthest genuinely SOLVED point, so it says where the swapped branch ends.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_swapped_closure.jl
Env:  TARGET (default 0.97)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP

TARGET = parse(Float64, get(ENV, "TARGET", "0.97"))

# `test/diag_swap_bisect.jl` showed that exactly ONE of the six swaps costs the
# conditioning — `houslack = shrBoTnom2` alone reproduces the ALL-SIX failure, and
# the other five each clear a 1e-5 probe in 2-3 Newton iterations. EXCLUDE drops
# swaps by either side's name so the five good ones (crucially `xcap = faccum`)
# can be measured while that one is investigated separately.
EXCLUDE = Set(filter(!isempty, split(get(ENV, "EXCLUDE", ""), ',')))
SWAPS = [p for p in TERM_CMF_SWAPS
         if !(String(p[1]) in EXCLUDE || String(p[2]) in EXCLUDE)]
isempty(EXCLUDE) || println("EXCLUDE=$(join(EXCLUDE, ',')) → $(length(SWAPS)) of " *
                            "$(length(TERM_CMF_SWAPS)) swaps applied")

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")

bmk = benchmark_levels(params)
initialize_model!(m, vars; bmk_levels=bmk)

# Swaps are applied to a PRISTINE model, before any solve. `solve_newton!` mutates
# the model — it pins orphan variables and deletes dead rows — so swapping after a
# solve could hand `apply_swaps!` a side that is already fixed for reasons that have
# nothing to do with the closure. Applying first also means `apply_swaps!` fixes each
# newly-exogenous variable at its benchmark SEED, which is exactly the value that
# keeps the benchmark an exact solution. (The base-closure benchmark is not re-checked
# here; it is established at 1.4260876923799515e-9 by every other gate in test/.)
println("\n══ 1. applying TERM.CMF's active swaps to the benchmark closure ══")
apply_swaps!(m, vars, SWAPS)

r1 = solve_newton!(m, vars; maxit=10, verbose=false)
println("\n══ benchmark, SWAPPED closure ══  status=$(r1.status)  ‖F‖∞=$(r1.residual)")
if r1.residual > 1e-8
    println("\n❌ THE SWAP BROKE THE BENCHMARK (‖F‖∞ = $(r1.residual)).")
    println("   Nothing below would be interpretable — a closure that cannot reproduce")
    println("   the base year cannot be used to measure a shock. Fix apply_swaps! or the")
    println("   swap list before reading any continuation result.")
    exit(1)
end
println("   ✅ benchmark preserved across the swap — the closure is consistent.")

sv = vars["blabnat"]
svr = sv isa AbstractArray ? first(sv) : sv

# ── 2. past the base-closure fold ───────────────────────────────────────────
println("\n══ 2. walk to 0.9997 (base closure folds at 0.99975) ══")
el2 = @elapsed cr2 = continuation_solve!(m, vars, svr, 0.9997;
                                         h0=0.5, hmax=1.0, hmin=1e-2,
                                         grow_iters=3, maxit=15,
                                         report=(("pcap", 2, 5),))
println("  reached=$(cr2.reached_target)  s_reached=$(cr2.s_reached)  " *
        "steps=$(cr2.nsteps) rejects=$(cr2.nrejects)  ($(round(el2,digits=1))s)")
println(cr2.reached_target ?
    "  ✅ FOLD REMOVED — the base-closure fold was a closure artifact." :
    cr2.nsteps == 0 ?
    "  ⚠️  NOT A FOLD — zero steps were accepted, so the walk never left the benchmark.\n" *
    "     A fold cannot sit at s=1.0: the benchmark solves to 1e-9 there. Failing at an\n" *
    "     arbitrarily small shock means the Jacobian is near-singular AT the base year\n" *
    "     under this closure. Localise it with test/diag_swap_bisect.jl; do NOT read\n" *
    "     this as evidence about the fold in either direction." :
    "  ❌ still obstructed at $(cr2.s_reached) — the fold survives a faithful closure.")

# ── 3. TERM.CMF's own scenario ──────────────────────────────────────────────
if cr2.reached_target
    println("\n══ 3. walk to TERM.CMF's scenario: blabnat = $TARGET ══")
    el3 = @elapsed cr3 = continuation_solve!(m, vars, svr, TARGET;
                                             h0=0.05, hmax=0.25, hmin=1e-3,
                                             grow_iters=3, maxit=20,
                                             report=(("pcap", 2, 5),))
    println("\n  reached=$(cr3.reached_target)  s_reached=$(cr3.s_reached)  " *
            "steps=$(cr3.nsteps) rejects=$(cr3.nrejects)  ($(round(el3,digits=1))s)")
    println("\n  ── accepted path ──")
    for (s, it, res) in cr3.path
        println("    ", rpad(round(s; sigdigits=8), 14), rpad(it, 6), round(res; sigdigits=6))
    end
    println("\n", cr3.reached_target ?
        "✅ THE PORT RUNS TERM.CMF's OWN SCENARIO (blabnat = $TARGET)." :
        "⚠️  reached $(cr3.s_reached) of $TARGET — a genuine solution, but short of the " *
        "scenario. Report the distance honestly; do not present a partial walk as the result.")
end

println("\nDone.")
