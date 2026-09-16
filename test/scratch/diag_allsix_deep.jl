"""
Why does ALL SIX fail, given the houslack null direction does NOT survive there?

`test/diag_houslack_null.jl` settled the singleton exactly: with `houslack`,
`fhou` and `fhou2` all free, `‖J·v‖∞ = 0.0` — identically zero, so
`swap houslack = shrBoTnom2` alone is a structurally singular closure. But under
ALL SIX, `xhouhtot = fhou` pins `fhou` and that candidate fails with residual
exactly 1.0 in the six `E_fhou` rows (the uncompensated `houslack` coefficient).
So ALL SIX stalls for some other reason, and the honest options are:

  (a) it is merely HARD, not singular — `maxit=20` in the bisect was too small,
      and the residual was still falling (1.03e-6 at iteration 20);
  (b) it is singular along a different, less obvious direction;
  (c) it is genuinely near-singular — no exact dependency, but a direction so
      weak that Newton cannot make progress in finite iterations.

This script separates (a) from (b)/(c) the cheap way, by simply letting it run:
60 iterations at a 1e-6 shock, verbose, so the per-iteration residual trajectory
is visible. A steadily falling residual that reaches tol says (a) and the whole
question dissolves. A trajectory that plateaus says (b) or (c) and earns a proper
null-space computation.

It also settles a bookkeeping point raised by the previous run, which reported
`83515 rows x 83516 free columns` under the swaps. That count is taken BEFORE
`solve_newton!` pins orphan variables, whereas the familiar 83,515 x 83,515
figure is post-squaring — so the extra column may be entirely normal. The control
case prints the same count with no swaps applied, which is the only way to tell.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_allsix_deep.jl
Env:  SHOCK (default 0.999999), MAXIT (default 60)
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP

SHOCK = parse(Float64, get(ENV, "SHOCK", "0.999999"))
MAXIT = parse(Int,     get(ENV, "MAXIT", "60"))

agg6, params = cached_pipeline(6)
bmk = benchmark_levels(params)

nfree(m, vars) = begin
    seen = Set{Int64}(); n = 0
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) && continue
        iv = vr.index.value
        iv in seen && continue
        push!(seen, iv); n += 1
    end
    n
end

println("\n══ column-count control (no swaps, before any solve) ══")
(mc, vc) = build_model_full!(agg6, params)
initialize_model!(mc, vc; bmk_levels=bmk)
nc = nfree(mc, vc)
println("  free columns with NO swaps: $nc")
println("  → the swapped runs reported 83516. Same count means the extra column is")
println("    an orphan that solve_newton! pins, and is NOT caused by the swaps.")

println("\n══ ALL SIX: $MAXIT iterations at blabnat = $SHOCK ══")
(m, vars) = build_model_full!(agg6, params)
initialize_model!(m, vars; bmk_levels=bmk)
apply_swaps!(m, vars, collect(TERM_CMF_SWAPS))
println("  free columns with ALL SIX: $(nfree(m, vars))")

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("benchmark: status=$(r0.status)  ‖F‖∞=$(r0.residual)")
@assert r0.residual <= 1e-8 "benchmark regressed"

sv = vars["blabnat"]
svr = sv isa AbstractArray ? first(sv) : sv
JuMP.fix(svr, SHOCK; force=true)

el = @elapsed r = solve_newton!(m, vars; maxit=MAXIT, verbose=true)
println("\nALL SIX @ $SHOCK: status=$(r.status)  iters=$(r.iters)  " *
        "‖F‖∞=$(r.residual)  ($(round(el, digits=1))s)")
println(r.status == :converged ?
    "\n➡ NOT (b) — convergence rules out a second exact null direction. The bisect's\n" *
    "  maxit=20 was too small, and the only structurally broken configuration is the\n" *
    "  houslack swap WITHOUT `xhouhtot = fhou`, which is not a closure anyone runs.\n" *
    "  This does NOT yet separate (a) from (c): convergence certifies a small residual,\n" *
    "  not a well-conditioned one. Read `max|d|` in the trace above against the shock\n" *
    "  size — but do NOT convert that into an 'amplification' figure here, because the\n" *
    "  step is in variable units (mostly Rp billion) while the shock is dimensionless.\n" *
    "  `test/diag_allsix_movers.jl` names the movers and their benchmark scales, and\n" *
    "  `test/diag_linearity.jl` decides (a) vs (c) by halving the shock: a genuine\n" *
    "  derivative halves with it, drift along a weak direction does not." :
    "\n➡ (b) or (c) — the residual did not reach tol in $MAXIT iterations. Read the\n" *
    "  trajectory above: a plateau means a weak or null direction that needs a real\n" *
    "  null-space computation, not more iterations.")
println("\nDone.")
