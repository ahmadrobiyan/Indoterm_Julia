# Seam path check: A/B step sizes under :schur, A under :lu (scratch, logs only).
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
say(msg) = (println(msg); flush(stdout))

agg6, params = cached_pipeline(6)

function runleg(h, method)
    r = run_model!(agg6, params, COALPRICE_REFERENCE;
                   h0=h, hmin=h, hmax=h, tol=1e-8, gdp=false,
                   verbose=false, linsolve=method)
    return r
end

rA = runleg(1.0, :schur)
say("A schur h=1: solved=$(rA.solved) res=$(rA.residual)")
rB = runleg(0.05, :schur)
say("B schur h=0.05: solved=$(rB.solved) res=$(rB.residual)")
rC = runleg(1.0, :lu)
say("C lu h=1: solved=$(rC.solved) res=$(rC.residual)")

function maxdiff(a, b)
    m = 0.0
    for (k, va) in a.values
        vb = b.values[k]
        if va isa AbstractArray
            m = max(m, maximum(abs.(va .- vb)))
        else
            m = max(m, abs(va - vb))
        end
    end
    return m
end
say("max|Δ| A-B (path independence under schur) = $(maxdiff(rA, rB))")
say("max|Δ| A-C (schur vs lu same path) = $(maxdiff(rA, rC))")
