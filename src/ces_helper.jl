export ces, ces_calibrate

# ── CES/CET closed-form helper ─────────────────────────────────────────────
# Ported from WayangJulia's dependency `ComputableGeneralEquilibriumHelpers`
# (~/.julia/packages/GlobalTradeAnalysisProjectModelV7/GTAPJulia/WayangJulia/),
# the same closed-form CES demand system used to convert TERM.TAB's %-change
# equations into genuine nonlinear levels equations (Step 5c, see PLAN.md's
# "Course correction" section). A negative `sigma` gives the CET (constant
# elasticity of transformation) dual of the same functional form.
#
# ces(y, p, α, σ, γ) returns the vector of cost-minimizing/revenue-maximizing
# component demands for aggregate `y` at component prices `p`, given share
# parameters `α`, elasticity `σ`, and scale `γ`.
function ces(y, p, α, σ, γ)
    if σ == 1
        y / (γ * prod((α ./ p) .^ α)) * (α ./ p)
    elseif σ == 0
        y .* α ./ γ
    elseif σ > 0
        α = Vector(α)
        c = 1 / γ * sum((α .^ σ) .* (p .^ (1.0 - σ)))^(1.0 / (1.0 - σ))
        toRet = Vector((y / γ) .* ((α .* γ .* c) ./ (p)) .^ σ)
        toRet[α.==0] .= 0
        return toRet
    else
        σ_adj = fill(σ, length(p))
        σ_adj[Vector(α).==0] .= 1
        c = 1 / γ * sum((α .^ σ_adj) .* (p .^ (1 - σ)))^(1 / (1 - σ))
        toret = (y / γ) .* ((α .* γ .* c) ./ (p)) .^ σ_adj
        return toret
    end
end

# ces_calibrate(quantities, sigma, output) finds the closed-form share (α) and
# scale (γ) parameters such that, at the benchmark (all prices = 1),
# ces(output, ones(length(quantities)), α, sigma, γ) == quantities exactly —
# i.e. the CES/CET nest reproduces the benchmark data by construction.
function ces_calibrate(quantities::AbstractVector{<:Real}, sigma::Real, output::Real)
    q = collect(Float64.(quantities))
    if output == 0
        return zeros(length(q)), 1.0
    end
    if sigma == 1
        α = q ./ sum(q)
        γ = output / prod(q .^ α)
        return α, γ
    elseif sigma == 0
        α = q ./ output
        return α, 1.0
    else
        # A zero-quantity component must get a zero share. For negative sigma
        # (the CET dual, used by the MAKE nest) `0 ^ (1/sigma)` is `0 ^ (negative)
        # = Inf`, which poisons `α = θ/sum(θ)` (Inf/Inf = NaN at the zero entries,
        # finite/Inf = 0 at the nonzero ones), silently zeroing the whole share
        # vector. That made `E_xmake!`'s `any(αv .> 0)` guard fail and take the
        # `xmake == 0` branch, breaking benchmark replication on the diagonal of
        # every diagonal-heavy MAKE column. Guarding qi==0 → θ=0 (matching the
        # `terms` line below and `ces`'s own `α.==0` handling) fixes it; for
        # positive sigma `0 ^ (positive)` is already 0, so this is a no-op there.
        θ = [qi == 0 ? 0.0 : qi ^ (1 / sigma) for qi in q]
        α = θ ./ sum(θ)
        ρ = (sigma - 1) / sigma
        terms = [qi == 0 ? 0.0 : αi * qi^ρ for (αi, qi) in zip(α, q)]
        γ = output * sum(terms)^(-1 / ρ)
        return α, γ
    end
end
