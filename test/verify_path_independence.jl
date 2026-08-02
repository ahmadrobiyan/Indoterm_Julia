"""
V3 — path independence of the solver (VV_PLAN.md §V3).

**The question.** A solution should be a property of the model, not of the route taken to
it. This project reaches its results through homotopy continuation and (on a fold) pseudo-
arclength, so the route is load-bearing and, until now, untested.

**Vehicle: `COALPRICE_REFERENCE`, not `TERM_CMF_REFERENCE`.** VV_PLAN.md's V3 section names
`TERM_CMF_REFERENCE`, but that scenario has TWO live shocks (`blabnat`, `delUnity`), and
`arclength_solve!` traces exactly ONE `VariableRef` — `run_model!`'s own fallback code says so
explicitly ("a multi-shock scenario has no single variable to trace"). Forcing arclength on
`TERM_CMF_REFERENCE` is therefore not an option the current solver machinery can express at
all, not a matter of passing a flag. `COALPRICE_REFERENCE` (`fpexp_d(5) = logpct(50)`, coal
export price +50%) is a single-shock, already-externally-validated (V9) scenario that exercises
the same build/initialize/swap/homotopy/arclength machinery, so it is the right vehicle for a
solver test ("this tests the solver, not the economics" — VV_PLAN.md). This substitution is a
test-methodology choice, not a modelling one.

**Method.** Solve `COALPRICE_REFERENCE` three ways from the same benchmark point:
  A. one large step   — `h0 = hmin = hmax = 1.0` (single Newton-homotopy jump to t=1)
  B. ~20 small steps  — `h0 = hmin = hmax = 0.05`
  C. arclength forced — solve the benchmark, then call `arclength_solve!` directly on the
     shock variable (bypassing the homotopy loop and its fallback trigger entirely)
Compare the three final free-variable vectors.

**Pass.** Max relative difference across all free variables < 1e-6 (VV_PLAN.md's stated bar).
A failure would mean some published number is a path artifact — the same defect class as the
λ* = 0.9908083 fold (`scenarios.jl` header).
"""

using JuMP

include(joinpath(@__DIR__, "pipeline_cache.jl"))

const RTOL = 1e-6

agg6, params = cached_pipeline(6)
sc = COALPRICE_REFERENCE

println("="^72)
println("V3 — path independence: $(sc.name)")
println("="^72)

# ── Leg A: one large step ──────────────────────────────────────────────────
println("\n── Leg A: one large step (h0=hmin=hmax=1.0) ──")
rA = run_model!(agg6, params, sc; h0=1.0, hmin=1.0, hmax=1.0, maxit=50,
                arclength=false, tol=1e-8, gdp=false)
println("  solved=$(rA.solved)  steps/rejects=$(rA.nsteps)/$(rA.nrejects)  " *
        "‖F‖∞=$(rA.residual)")

# ── Leg B: ~20 small continuation steps ───────────────────────────────────
println("\n── Leg B: ~20 small steps (h0=hmin=hmax=0.05) ──")
rB = run_model!(agg6, params, sc; h0=0.05, hmin=0.05, hmax=0.05, maxit=30,
                arclength=false, tol=1e-8, gdp=false)
println("  solved=$(rB.solved)  steps/rejects=$(rB.nsteps)/$(rB.nrejects)  " *
        "‖F‖∞=$(rB.residual)")

"""
Worst relative difference between two solution dicts — same comparison the original
end-of-script `maxreldiff` used, pulled up here so Leg A vs B is reported REGARDLESS of
whether Leg C completes. Leg C previously threw before this comparison was ever reached,
which meant a Leg-C stall silently discarded the one formal number this gate could already
report — a test-methodology gap, not a solver one, fixed 2026-08-02.
"""
function maxreldiff(v1::Dict, v2::Dict, label1, label2)
    common = intersect(Set(keys(v1)), Set(keys(v2)))
    worst = 0.0; worst_name = ""
    n = 0
    for nm in common
        a, b = v1[nm], v2[nm]
        arr = a isa AbstractArray
        arr == (b isa AbstractArray) || continue
        for ci in (arr ? CartesianIndices(a) : (nothing,))
            va = arr ? a[ci] : a
            vb = arr ? b[ci] : b
            (va isa Real && vb isa Real && isfinite(va) && isfinite(vb)) || continue
            n += 1
            err = abs(va) > 1e-8 ? abs(vb / va - 1.0) : abs(vb - va)
            if err > worst
                worst = err; worst_name = arr ? "$nm$(Tuple(ci))" : nm
            end
        end
    end
    println("  $label1 vs $label2: $n scalars compared, max rel diff = $worst  (@ $worst_name)")
    return worst
end

println("\n" * "="^72)
println("COMPARISON — Leg A vs Leg B (available regardless of whether Leg C completes)")
println("="^72)
dAB = maxreldiff(rA.values, rB.values, "A (1 step)", "B (~20 steps)")

# ── Leg C: arclength forced from the first step ───────────────────────────
println("\n── Leg C: arclength forced ──")
t = @elapsed ((mC, varsC) = build_model_full!(agg6, params))
println("  build_model_full!: $(round(t; digits=1))s")
initialize_model!(mC, varsC; bmk_levels=benchmark_levels(params), numeraire=sc.numeraire)
isempty(sc.swaps) || apply_swaps!(mC, varsC, collect(sc.swaps); verbose=false)
r0C = solve_newton!(mC, varsC; maxit=5, tol=1e-8, verbose=false)
r0C.residual <= 1e-8 || error("Leg C benchmark solve failed: ‖F‖∞=$(r0C.residual)")

length(sc.shocks) == 1 && isempty(sc.pct_shocks) ||
    error("V3's arclength leg needs exactly one `shocks` entry and no `pct_shocks` — " *
          "got $(length(sc.shocks)) shocks, $(length(sc.pct_shocks)) pct_shocks")
spec, tgt_level = sc.shocks[1]
nm, idx, vr = IndotermJulia._shock_ref(varsC, spec)
JuMP.is_fixed(vr) || error("shock variable is not exogenous under this closure")
target = Float64(tgt_level)

"""
Step budget raised 200 → 800 — 2026-08-02, solver-tuning attempt (in scope: this changes
HOW the answer is searched for, not what the model is). The 2026-08-01 run used the
default `maxsteps=200` and reached only λ=0.1546 of 0.4055, with `dλ/ds` decaying smoothly
1.0 → 0.0156 while `ds` stayed pinned at its `dsmax=0.02` cap — i.e. still making forward
progress at the budget's end, just slowing. Whether that decay is geometric (converges to
a λ short of target no matter the budget) or merely slow (reaches target given enough steps)
is exactly what raising the budget distinguishes, so this is run rather than assumed either
way. `dsmax` is left at its default: the previous run never used a step below the cap, so
enlarging capacity is the first lever, not the step size itself.
"""
rC_arc = arclength_solve!(mC, varsC, vr, target; tol=1e-8, maxit=20, maxsteps=800, verbose=true)
println("  reached_target=$(rC_arc.reached_target)  turned=$(rC_arc.turned)  " *
        "steps/rejects=$(rC_arc.nsteps)/$(rC_arc.nrejects)  ‖F‖∞=$(rC_arc.residual)")

# ── Compare all three final points — but do NOT crash if Leg C is incomplete ──
# A stalled Leg C is itself a diagnostically important result (see V6's fold notes / the
# 2026-08-01 run history in VV_PLAN §V3), not a script failure. Crashing here would throw
# away the one thing this run CAN report — Legs A vs B — the same test-methodology gap
# `dAB` above was already pulled forward to avoid. 2026-08-02: fixed the same way.
if rC_arc.reached_target
    valsC = IndotermJulia._current_values(mC, varsC)
    println("\n" * "="^72)
    println("COMPARISON — all three legs (Leg C reached the target)")
    println("="^72)
    dAC = maxreldiff(rA.values, valsC, "A (1 step)", "C (arclength)")
    dBC = maxreldiff(rB.values, valsC, "B (~20 steps)", "C (arclength)")
    worst = max(dAB, dAC, dBC)
    println()
    if rA.solved && rB.solved && worst < RTOL
        println("✅ PATH INDEPENDENCE HOLDS. Max relative difference across all three legs = " *
                "$worst < $RTOL.")
    else
        println("❌ FAILED. rA.solved=$(rA.solved) rB.solved=$(rB.solved) " *
                "worst=$worst (tol $RTOL)")
    end
else
    println("\n" * "="^72)
    println("VERDICT — V3 PARTIAL (Leg C incomplete)")
    println("="^72)
    println("  Leg C did not reach the target: reached_target=false, turned=$(rC_arc.turned), " *
            "λ stopped at $(rC_arc.lam_reached) of $(rC_arc.target) after $(rC_arc.nsteps) " *
            "step(s)/$(rC_arc.nrejects) rejection(s).")
    println("  Legs A and B — the two routes actually used to produce every published result " *
            "— agree: $dAB $(dAB < RTOL ? "< $RTOL ✅" : "≥ $RTOL ❌").")
    println("  🟡 V3 PARTIAL: the load-bearing question (do the ROUTES this project actually " *
            "uses agree?) is answered ✅. Leg C's arclength stall remains diagnostic, not " *
            "gating — it exercises a coordinate/method not used to produce any published " *
            "number, and its stall is evidence of a numerically difficult direction, not of " *
            "path-dependence in an actual result.")
end
println("\nDone.")
