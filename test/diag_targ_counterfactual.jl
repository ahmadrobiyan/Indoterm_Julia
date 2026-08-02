"""
Does TARG[2] cause the Livestock fold mode?  (The decisive counterfactual.)

`test/diag_pcap_rank.jl` localised the fold to a sector-2 mode — all six regions of
Livestock responding 168-1187x to a 1e-6 shock while every other cell stays under
86x — and `PLAN.md` §N.4 refuted the four obvious causes ([2,5] degeneracy, dust
capital, near-Leontief SIGMAPRIM, intermediate self-loop). §N.5 found the one
genuine sector-2 outlier:

    TARG[2] = 0.030321   the MINIMUM of all 25 sectors
                         5.1x below median (0.15506), 2.1x below next lowest
    DPRC    = 0.05       uniform

so sector 2 is the ONLY sector whose normal rate of return sits below its
depreciation rate. `TERM.TAB:2719` reads TARG as RNORMAL.

But correlation is not cause, and single-parameter reasoning has ALREADY failed
once here: sectors 4 and 7 carry strictly lower SIGMAPRIM than sector 2 and respond
20x less. So this does not argue from the anomaly — it removes it and re-measures.

Counterfactual: set RNORMAL[2] to the grid median and re-derive. RNORMAL reaches
the model by two channels, and a test that moved only one would be worthless:

  1. DIRECTLY, as a multiplier in GRETEXP (`prepare_parameters.jl:843`);
  2. INDIRECTLY, through the capital stock — `build_premod!.jl:131-139` sets
     CAPSTOK = CAP/RNORMAL precisely so the benchmark GROSSRET equals RNORMAL.
     Raising RNORMAL therefore SHRINKS the stock by oldT/newT, which is why
     Livestock currently carries a stock 46.6x its rental against 5-8x elsewhere.

Rather than hand-patch the ~12 derived arrays (GROSSRET, GROSSGRO, GROMAX, GRETEXP,
MCOEFF, CAPADD, CAPSTOK_OLDP, CAPSTOK_D, ...) and risk missing one, this patches the
two SOURCE headers in `agg` and re-runs the real `prepare_parameters!`, so every
derived coefficient is recomputed by the production code path.

Control and counterfactual each start from an independent `deepcopy` of the cached
database and go through the identical path, so the two rankings are comparable.

GUARD: the counterfactual changes the capital calibration, so the benchmark must be
re-checked. If it no longer solves to ~1e-9 the perturbed database is not a valid
base year and the comparison is meaningless — the script says so and stops rather
than reporting a ranking computed off a broken benchmark.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_targ_counterfactual.jl
Env:  SHOCK (default 0.999999), TOPN (default 12)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, Printf

SHOCK = parse(Float64, get(ENV, "SHOCK", "0.999999"))
TOPN  = parse(Int,     get(ENV, "TOPN", "12"))

agg0, _ = cached_pipeline(6)

"""Build, swap, benchmark and shock `agg`; return (ok, ranking, benchmark residual)."""
function measure(agg, label)
    println("\n" * "="^72)
    println("══ $label ══")
    println("="^72)
    params = prepare_parameters!(agg)
    bmk = benchmark_levels(params)
    (m, vars) = build_model_full!(agg, params)
    initialize_model!(m, vars; bmk_levels=bmk)
    apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))

    r0 = solve_newton!(m, vars; maxit=8, verbose=false)
    @printf("  benchmark: status=%s  ‖F‖∞=%.6g\n", r0.status, r0.residual)
    if r0.residual > 1e-8
        println("  ⛔ benchmark does NOT hold — this database is not a valid base year,")
        println("     so any ranking from it would be meaningless. Stopping this case.")
        return (false, Tuple{Int,Int,Float64}[], r0.residual)
    end

    pcap = vars["pcap"]
    na, nr = size(parent(agg["1CAP"]))
    base = [JuMP.start_value(pcap[i,d]) for i in 1:na, d in 1:nr]

    sv = vars["blabnat"]
    svr = sv isa AbstractArray ? first(sv) : sv
    JuMP.fix(svr, SHOCK; force=true)
    r = solve_newton!(m, vars; maxit=30, verbose=false)
    @printf("  shock %.7g: status=%s iters=%d ‖F‖∞=%.6g\n", SHOCK, r.status, r.iters, r.residual)
    r.status == :converged || println("  ⚠️  shock did not converge — ranking is unreliable.")

    δ = 1.0 - SHOCK
    rows = Tuple{Int,Int,Float64}[]
    for i in 1:na, d in 1:nr
        b = base[i,d]
        (b === nothing || !isfinite(b) || abs(b) < 1e-12) && continue
        now = JuMP.start_value(pcap[i,d])
        now === nothing && continue
        push!(rows, (i, d, abs(now - b) / abs(b) / δ))     # amplification per unit shock
    end
    sort!(rows; by = x -> -x[3])
    return (true, rows, r0.residual)
end

# ── control ───────────────────────────────────────────────────────────────────
(ok0, rank0, _) = measure(deepcopy(agg0), "CONTROL — TARG as calibrated")

# ── counterfactual ────────────────────────────────────────────────────────────
aggX = deepcopy(agg0)
T = parent(aggX["TARG"]); S = parent(aggX["STOC"])
na = size(T, 1)
oldT = T[2,1]
newT = sort([T[i,1] for i in 1:na])[cld(na, 2)]
@printf("\nTARG[2]: %.6g → %.6g (grid median);  CAPSTOK[2,:] scaled by %.6g\n",
        oldT, newT, oldT / newT)
for d in 1:size(T, 2)
    S[2,d] *= oldT / newT          # CAPSTOK = CAP/RNORMAL, so it scales as 1/RNORMAL
    T[2,d]  = newT
end
(okX, rankX, _) = measure(aggX, "COUNTERFACTUAL — TARG[2] set to the grid median")

# ── verdict ───────────────────────────────────────────────────────────────────
if !(ok0 && okX)
    println("\n⛔ one of the two cases failed its benchmark gate — no comparison possible.")
else
    amp(rank, i) = (k = findfirst(x -> x[1] == i, rank); k === nothing ? NaN : rank[k][3])
    println("\n" * "="^72)
    @printf("── top %d amplifications (relative pcap move per unit shock) ──\n", TOPN)
    @printf("  %-14s %-24s %s\n", "", "CONTROL", "COUNTERFACTUAL")
    for k in 1:TOPN
        a = k <= length(rank0) ? @sprintf("[%d,%d] %10.1f", rank0[k][1], rank0[k][2], rank0[k][3]) : ""
        b = k <= length(rankX) ? @sprintf("[%d,%d] %10.1f", rankX[k][1], rankX[k][2], rankX[k][3]) : ""
        @printf("  %-14d %-24s %s\n", k, a, b)
    end

    s2_0 = maximum(x -> x[1] == 2 ? x[3] : -Inf, rank0)
    s2_X = maximum(x -> x[1] == 2 ? x[3] : -Inf, rankX)
    oth0 = maximum(x -> x[1] != 2 ? x[3] : -Inf, rank0)
    othX = maximum(x -> x[1] != 2 ? x[3] : -Inf, rankX)
    @printf("\n  worst sector-2 amplification : %10.1f → %10.1f  (×%.3f)\n", s2_0, s2_X, s2_X / s2_0)
    @printf("  worst non-sector-2           : %10.1f → %10.1f  (×%.3f)\n", oth0, othX, othX / oth0)
    @printf("  sector-2 dominance ratio     : %10.1f → %10.1f\n", s2_0 / oth0, s2_X / othX)

    collapsed = s2_X < 0.5 * s2_0 && (s2_X / othX) < 0.5 * (s2_0 / oth0)
    println("\n", collapsed ?
        "➡ TARG[2] CAUSES the Livestock mode. Removing the anomaly collapses both the\n" *
        "  sector-2 amplification and its dominance over every other sector. The early\n" *
        "  fold is then a DATA DEFECT, not a property of the economy: the fix is to\n" *
        "  correct TARG in the source database (verify the INFILE header against\n" *
        "  TERM.CMF first), and pseudo-arclength would be the wrong build — it would\n" *
        "  faithfully trace a branch that should not exist." :
        "➡ TARG[2] does NOT explain the Livestock mode — it survives the counterfactual.\n" *
        "  That makes it the fifth refuted single-cause hypothesis, and it means the\n" *
        "  stiffness is emergent from the calibrated structure rather than any one\n" *
        "  number. The limit point is then a genuine property of the model and\n" *
        "  pseudo-arclength IS the right next build.")
end
println("\nDone.")
