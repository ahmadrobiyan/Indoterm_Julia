# Stage 2a seam equivalence: :lu vs :schur must agree (scratch, logs only).
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
say(msg) = (println(msg); flush(stdout))

agg6, params = cached_pipeline(6)

function fresh()
    (m, vars) = build_model_full!(agg6, params)
    initialize_model!(m, vars; bmk_levels=benchmark_levels(params),
                      numeraire=:exrate)
    apply_swaps!(m, vars, collect(COALPRICE_SWAPS); verbose=false)
    return m, vars
end

for method in (:lu, :schur)
    m, vars = fresh()
    t = @elapsed r = solve_newton!(m, vars; maxit=5, tol=1e-8, verbose=false,
                                     linsolve=method)
    say("benchmark $method: status=$(r.status) res=$(r.residual) $(round(t;digits=1))s")
end

results = Dict()
for method in (:lu, :schur)
    t = @elapsed r = run_model!(agg6, params, COALPRICE_REFERENCE;
                                h0=0.25, tol=1e-8, gdp=false, verbose=false,
                                linsolve=method)
    nm = r.values["NatMacro"]
    hc(m) = 100 * (nm[findfirst(==(m), MAINMACROS)] - 1.0)
    results[method] = Dict(k => hc(k) for k in
        ["RealHou","RealInv","ExpVol","ImpVolUsed","RealGNE","RealGDP","AggEmploy","CPI"])
    say("coalprice $method: solved=$(r.solved) res=$(r.residual) $(round(t/60;digits=1))min")
end
mx = maximum(abs(results[:lu][k] - results[:schur][k]) for k in keys(results[:lu]))
for k in sort(collect(keys(results[:lu])))
    say("  $k  lu=$(round(results[:lu][k];digits=5))  schur=$(round(results[:schur][k];digits=5))")
end
say("max|Δ| headline metrics lu vs schur = $mx")
say(mx < 1e-6 ? "SEAM EQUIVALENT ✅" : "SEAM DIVERGED ❌")
