# Stage 2b.4 + 2b.5: 34-region COALPRICE probe (scratch, logs only)
#
# What changed vs the earlier small-shock probes: those ran `swaps = ()` — the
# bare automatic closure, which src/closures.jl and src/run_model!.jl both
# document folds on blabnat (fixed xcap+xlnd, ~200x amplification on the
# smallest fixed-factor cell). So their stalls were at least partly a closure
# artifact, NOT a 34-region fact.
#
# SCOPE (per user, 2026-09-17): the real-world scenario is a +12% coal export
# price — a 50% jump is an artificial test case (it only exists to match the
# authors' published V9 reference at 6 regions) and will not happen in practice.
# The +12% shock directly checks, and is directly checked by, the existing
# 6-region V4 round-6 reference (Kalimantan GDP ~+5.2%, all other islands ~0).
#
# The 34-region solve runs COALPRICE_SWAPS (the externally-validated closure)
# with numeraire = :exrate. CONVERGENCE TOLERANCE is 1e-4, NOT 1e-8: at 34
# regions (1,196,221 equations, kappa(J) ~ 1e10) the honest double-precision
# floor is ~1e-5 of the equilibrated residual — measured: the +12% first step
# plateaus at 1.9e-5, which a 1e-8 tolerance can never accept. The answer is
# judged economically (signs/magnitude vs the 6-region reference), not on pure
# arithmetic.
#
# Each accepted Newton step at 34 regions is ~10-30 min (LU ~260s), so this is
# a multi-hour background job. run_model! flushes per line, but the FIRST solve
# can sit silent for a long stretch; poll the process CPU alongside the log.
# Run in background and poll:
#   julia --project=. test/scratch/_solve_34_coalprice.jl
# Log: logs/solve_34_coalprice_<date>.log
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using IndotermJulia
using Printf

# Stage selection: "A" only +12%, "B" only +50% (kept for the record / future
# V9-style work), "NONE" = parse/load-only smoke test. Default is "A".
# Prefer running A and B in SEPARATE processes: each 34-region solve peaks ~16 GB
# and two live JuMP models in one process risk OOM on this 25.5 GB machine.
#
# Tolerance policy (34-region): 1e-4 default (see header). Pass a second arg to
# override, e.g. `A 1e-6`.
STAGE = length(ARGS) >= 1 ? uppercase(ARGS[1]) : "A"
TOL   = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1e-4

say(msg) = (println(msg); flush(stdout))
say("--- Stage 2b.4/2b.5: 34-region COALPRICE probe (stage=$STAGE tol=$TOL) ---")

t_pipe = @elapsed agg34, params34 = cached_pipeline(34)
say("pipeline: $(round(t_pipe; digits=1))s (cache HIT expected)")
if STAGE == "NONE"
    say("stage=NONE — parse/load smoke test only; exiting.")
    exit(0)
end

# ── National macro table for any successfully solved leg ───────────────────
function national_macro_tbl(r, tag)
    r.solved || (say("\n[$(tag)] not solved (t=$(r.t_reached)) — no macro table."); return)
    natmacro = r.values["NatMacro"]
    say("\n=== NATIONAL MACRO (% change vs benchmark) — $tag ===")
    @printf("%-14s %10s\n", "macro", "%chg")
    for (label, macname) in (
            ("Real HousCon", "RealHou"), ("Real Invest", "RealInv"),
            ("Export vol", "ExpVol"), ("Import vol", "ImpVolUsed"),
            ("Real GNE", "RealGNE"), ("Real GDP", "RealGDP"),
            ("Employment", "AggEmploy"), ("CPI", "CPI"))
        i = findfirst(==(macname), MAINMACROS)
        i === nothing && continue
        v = natmacro[i]
        @printf("%-14s %+10.3f\n", label, 100 * (v - 1.0))
    end
    # regional real GDP % change from the GDPReport
    if r.gdp !== nothing
        say("\n=== REGIONAL nominal income-side GDP levels & change — $tag ===")
        for d in 1:length(REG)
            tot = r.gdp.inc_total[d]
            base = sum(parent(params34["GDPINCSUM"])[d, :])
            pct = base > 0 ? 100 * (r.gdp.inc_change[d, :] |> sum) / base : NaN
            @printf("%-12s %12.1f  %+10.3f%%\n", REG[d], tot, pct)
        end
        # aggregate to 6 island groups for comparison vs Table 2 / V9
        say("\n=== ISLAND-GROUP (34->6) income-side GDP % change — $tag ===")
        g6 = zeros(6); g6base = zeros(6)
        for d in 1:length(REG)
            g = REG_MAP_34_to_6[d]
            g6[g] += sum(r.gdp.inc_change[d, :])
            g6base[g] += sum(parent(params34["GDPINCSUM"])[d, :])
        end
        for g in 1:6
            pct = g6base[g] > 0 ? 100 * g6[g] / g6base[g] : NaN
            @printf("%-12s %+10.3f%%\n", REG6[g], pct)
        end
    end
end

# ── (a) +12% reduced coal shock: 2b.4 shock-convergence proof at 34 regions ──
say("\n=== STAGE A: coal export price +12% (logpct(12)) at 34 regions ===")
sc_a = Scenario(
    name   = "34reg coalprice +12% (2b.4)",
    source = "constructed; COALPRICE_SWAPS closure, fpexp_d(coal) +12%",
    swaps  = COALPRICE_SWAPS,
    shocks = [("fpexp_d", 5) => logpct(12)],   # coal = sector 5
    numeraire = :exrate,
    notes = "34-region shock-convergence probe (Stage 2b.4). Reduced shock chosen so the point is reachable.",
)
if STAGE in ("A", "AB")
    describe(sc_a)
    t_a = @elapsed r_a = run_model!(agg34, params34, sc_a;
        h0 = 0.05, hmin = 1e-4, hmax = 0.25,
        tol = TOL, maxit = 30,
        linsolve = :schur, verbose = true, gdp = true)
    say("RESULT A: solved=$(r_a.solved) t_reached=$(r_a.t_reached) residual=$(r_a.residual) time=$(round(t_a; digits=1))s")
    national_macro_tbl(r_a, "A +12%")
end

# ── (b) full +50% reference: 2b.5 ──────────────────────────────────────────
say("\n=== STAGE B: COALPRICE_REFERENCE +50% at 34 regions ===")
if STAGE in ("B", "AB")
    describe(COALPRICE_REFERENCE)
    t_b = @elapsed r_b = run_model!(agg34, params34, COALPRICE_REFERENCE;
        h0 = 0.05, hmin = 1e-4, hmax = 0.25,
        tol = TOL, maxit = 30,
        linsolve = :schur, verbose = true, gdp = true)
    say("RESULT B: solved=$(r_b.solved) t_reached=$(r_b.t_reached) residual=$(r_b.residual) time=$(round(t_b; digits=1))s")
    national_macro_tbl(r_b, "B +50%")
end

say("\n--- 34-REGION COALPRICE PROBE DONE ---")
