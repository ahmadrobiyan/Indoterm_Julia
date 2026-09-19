# Tight repro: Coal_Drop_20 stall (H1-H6 loop).
#
# Full test/research_coal_trap.jl costs hours per iteration — too loose for
# solver forensics. This script rebuilds the exact stall context cheaply:
# last accepted point t=0.046875 of fpexp_d[5]=logpct(-20) under
# COALPRICE_SWAPS/:exrate, then fires ONE solve_newton! at the next t_try
# with verbose=ON.
#
# RED (bug present): it 1 max|d|>2000 + 20 dead backtracks + CGNR lin_rel>0.3,
#   or a crawling corrector that trips STALL-ABORT.
# GREEN (fixed): converges in a handful of its, or fails fast with a stall
#   signal instead of burning maxit iterations.
#
# Run: julia --project=. test/scratch/_repro_coaldrop_stall.jl
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP

agg6, params = cached_pipeline(6)

COAL = findfirst(==("Coal"), AGGCOM)
@assert COAL == 5 "AGGCOM order changed"

# Same scenario object as test/research_coal_trap.jl Coal_Drop_20
sc = Scenario(
    name   = "repro — world coal export price -20%",
    source = "constructed; stall repro for H1-H6",
    swaps  = COALPRICE_SWAPS,
    shocks = [("fpexp_d", COAL) => logpct(-20.0)],
    numeraire = :exrate,
    notes = "Tight loop only — no results quoted from this run.",
)

# Build + close exactly as run_model! does, then drive t manually.
(m, vars) = build_model_full!(agg6, params)
initialize_model!(m, vars; bmk_levels=benchmark_levels(params), numeraire=:exrate)
apply_swaps!(m, vars, collect(COALPRICE_SWAPS); verbose=true)

r0 = solve_newton!(m, vars; maxit=5, tol=1e-8, verbose=true)
@assert r0.residual <= 1e-8 "benchmark must reproduce under this closure"

# Walk to the last accepted point t=0.046875 in small verified steps, then
# attempt the step that stalled in the full run.
nm, idx, vr = IndotermJulia._shock_ref(vars, ("fpexp_d", COAL))
b = JuMP.fix_value(vr); tg = logpct(-20.0)
println("shock base=$b target=$tg")
flush(stdout)

function try_t!(t_; diag=false)
    JuMP.fix(vr, b + t_ * (tg - b); force=true)
    r = solve_newton!(m, vars; maxit=30, tol=1e-8, verbose=true,
                      diagnose_worst_row=diag)
    println("t=$t_ -> $(r.status) iters=$(r.iters) residual=$(r.residual)")
    flush(stdout)
    return r
end

# Staged walk: 0.046875 in halves (each step easy, mirrors homotopy growth).
for t_ in (0.015625, 0.03125, 0.046875)
    r = try_t!(t_)
    r.residual <= 1e-8 || error("staging failed at t=$t_ — repro invalid")
end
println("staging complete — now the stall step:")
flush(stdout)
r = try_t!(0.0625; diag=true)
println("STALL-STEP verdict: $(r.status) iters=$(r.iters) residual=$(r.residual)")
println("RED = max|d|>2000 + dead backtracks + lin_rel>0.3 above, or STALL-ABORT.")
println("GREEN = converged, or fast :no_progress with stall signal.")
println("Done.")
