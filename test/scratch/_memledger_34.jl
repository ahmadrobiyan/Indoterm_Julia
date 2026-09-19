# 34-region memory ledger: WHERE does the first Newton iteration die?
#
# Replaces the broken _mem_probe_34.jl (which reported "System size: 208" and
# then threw UndefVarError). Differences that matter:
#   * keeps the JuMP model + MOI evaluator ALIVE while factorizing, like the real
#     solver does — the Stage-1 probe dropped the model (`m = nothing`) and so
#     under-reported peak RSS;
#   * every stage is FLUSHED through solve_newton!'s `memlog` ledger, so a hard
#     OOM still leaves the last stage reached in the log;
#   * it drives the actual `solve_newton!` path, not a re-implementation.
#
# Run:  julia --project=. test/scratch/_memledger_34.jl
# Log:  logs/memledger_34_<date>.log   (stdout unbuffered at each mem() marker)
# Scratch only — no src/ changes beyond the opt-in memlog ledger.

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, Dates

say(msg) = (println(msg); flush(stdout))
gib(x) = round(x / 2^30; digits = 2)
peak() = gib(Sys.maxrss())

say("="^78)
say("34-region memory ledger   $(Dates.now())")
say("="^78)

t = @elapsed ((agg34, params) = cached_pipeline(34))
say("cache loaded: $(round(t; digits = 1))s   peakRSS=$(peak()) GiB")

t = @elapsed ((m, vars) = build_model_full!(agg34, params))
say("build_model_full!: $(round(t; digits = 1))s   peakRSS=$(peak()) GiB")

initialize_model!(m, vars; bmk_levels = benchmark_levels(params))
say("initialize_model!   peakRSS=$(peak()) GiB")

# ── benchmark: start point already satisfies F(x)=0 by construction, so this
#    evaluates the Jacobian once but takes no step and does no factorization.
r0 = solve_newton!(m, vars; maxit = 2, tol = 1e-8, verbose = true,
                   linsolve = :schur, memlog = true)
say("benchmark: $(r0.status)  ‖F‖∞=$(r0.residual)   peakRSS=$(peak()) GiB")

# ── the real test: one shock, then the FIRST Newton iteration on the :schur path.
JuMP.fix(vars["blabnat"], 0.99; force = true)
say("shock applied: blabnat = 0.99")

r1 = solve_newton!(m, vars; maxit = 1, tol = 1e-8, verbose = true,
                   linsolve = :schur, memlog = true)
say("shock solve: status=$(r1.status)  iters=$(r1.iters)  " *
    "‖F‖∞=$(r1.residual)   peakRSS=$(peak()) GiB")

say("="^78)
say("final peakRSS=$(peak()) GiB   (machine has ~23.8 GiB)")
say("="^78)
