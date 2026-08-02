"""
Gate for `continuation_solve!` (src/continuation.jl).

`test/verify_continuation.jl` established the fact this driver is built on: the
Newton solver converges a shock in 2–3 iterations when it starts near the answer
(blabnat 0.999999 / 0.99999 / 0.9999 all reached ‖F‖∞ ~ 1e-9), and fails only
when asked to leap. The driver automates the walk with an adaptive step, so no
ladder has to be hand-tuned per scenario.

Pass = both of:
  1. the benchmark still solves to ~1.4e-9 (no regression), and
  2. `continuation_solve!` reaches blabnat = TARGET.

A failure is still informative: `s_reached` is the furthest point actually
SOLVED (residual ≤ tol, not a stalled iterate). If the step size collapses at a
specific s while every earlier step took 2–3 iterations, the branch is
obstructed there — a fold — and no step-size policy will pass it.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_continuation_driver.jl
Env:  TARGET (default 0.9997), H0 (default 0.25)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP

TARGET = parse(Float64, get(ENV, "TARGET", "0.9997"))
H0     = parse(Float64, get(ENV, "H0", "0.25"))

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("\n══ benchmark ══  status=$(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= 1e-8 "benchmark regressed — fix that before reading anything below"

sv = vars["blabnat"]
svr = sv isa AbstractArray ? first(sv) : sv

println("\n══ continuation to blabnat = $TARGET ══")
el = @elapsed cr = continuation_solve!(m, vars, svr, TARGET;
                                       h0=H0, maxit=40,
                                       report=(("pcap", 2, 5), ("xinvitot", 2, 5)))

println("\n── path ──")
println("shock             iters   ‖F‖∞")
for (s, it, res) in cr.path
    println(rpad(round(s; sigdigits=10), 18), rpad(it, 8), round(res; sigdigits=6))
end

println("\nreached_target = $(cr.reached_target)")
println("s_reached      = $(cr.s_reached)   (target $(cr.target))")
println("steps          = $(cr.nsteps)  rejections = $(cr.nrejects)")
println("h_final        = $(cr.h_final)")
println("elapsed        = $(round(el, digits=1))s")

println("\n", cr.reached_target ?
    "✅ CONTINUATION DRIVER WORKS — shocks are now solvable; proceed to the real scenario." :
    "❌ obstructed at s=$(cr.s_reached) — the last SOLVED point. If every step before it " *
    "converged in 2-3 iterations, this is a fold and needs pseudo-arclength, not a smaller h.")
println("\nDone.")
