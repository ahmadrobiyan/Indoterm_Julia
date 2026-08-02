"""
Continuation that can pass a fold — arclength in a *response-weighted* metric,
with a *local coordinate* closing the system.

`continuation.jl` walks the shock parameter λ directly and runs a Newton solve at
each value. That is *natural-parameter* continuation, and PLAN.md §N.16 measured
why it stops: at `blabnat ≈ 0.9908` the branch has a genuine **fold** (limit
point). There the Jacobian `J = ∂F/∂x` is singular, `‖dx/dλ‖∞` blows up like
`1/(λ*−λ)` (measured: 1.95e6 → 3.51e7 while σ_min collapsed 3.85e-9 → 1.16e-10),
and λ *stops being a valid coordinate for the branch* — past the turn there is no
solution at that λ at all. No step size and no solver tuning can fix that; it is a
property of the equations.

The cure is to stop privileging λ: make it an unknown, and supply the missing
equation from the geometry of the curve. Two choices matter, and the first version
of this file got both wrong in ways that made it useless in practice. Both are
recorded below, because the reasons are specific to this model.

## 1. The metric — why not plain `‖Δz‖₂`

Keller's pseudo-arclength closes the system with `τᵀ(z − z_prev) = Δs`, `τ` the
unit tangent. "Unit" needs a norm, and the obvious `‖·‖₂` over the scaled state is
**catastrophic here**. Measured on the first run (PLAN.md §N.19):

    4 accepted steps, Δs = 0.75 total  ⟹  Δλ = −1.155e-7,  τ_λ = 1.54e-7

`τ_λ = 1/‖v‖₂` with `v = dz/dλ`, so `‖v‖₂ = 6.49e6` — while the largest genuine
*relative* response is `whouhtot[4]` at only **2.69** (−1.39e6 on a base of
515,617). The 6.49e6 is entirely the `del*` variables: additive adjustment terms
that legitimately sit at zero, so their column scale is `max(|x|,1) = 1`, while
their derivatives are Rp-billions. Three reporting aggregates
(`delVGDPEXP`, `delXGDPEXP`, `delGDPINC`) were consuming the whole arclength
budget and leaving λ 1.5e-7 of it. Reaching `blabnat = 0.97` would have taken
~195,000 steps ≈ two months of wall-clock. **This is the `del*`-at-zero units
trap (PLAN.md §F) reproduced inside the fix for it — the seventh instance.**

The fix is that arclength may be measured in *any* fixed positive-definite
metric — the choice affects efficiency, not validity. So the metric is derived
from the branch's own sensitivity at the start point:

    v   = dz/dλ                      (the bootstrap tangent, before normalising)
    wt  = 1 / max(|vᵢ|, εfloor)      (εfloor = 1e-6·‖v‖∞, so a near-null
                                      component cannot get infinite weight)

and the tangent is normalised so `‖wt ⊙ τ‖∞ = 1`. Every component that actually
moves along the branch then contributes O(1), no component can dominate on the
strength of its unit of measurement, and `Δs` recovers its intended meaning:
roughly "how far λ would travel if the branch were straight". Note `wt` is a
metric, kept deliberately separate from `cscale`, which still governs state
representation, row scaling and the bound logic — conflating those two is what
produced the defect.

The ∞-norm rather than the 2-norm is chosen so the measure does not carry a
`√N` penalty that would shrink with model size, and because it pairs exactly with
the local coordinate below.

## 2. The closing row — why not the dense tangent

With `τᵀ` as the last row, the augmented matrix has one **fully dense** row of
83,516 entries. The first run died on its fifth factorisation with
`OpenBLAS: malloc failed in gemm_driver` — memory exhausted in the dense kernels
UMFPACK uses on its frontal matrices, exactly the blow-up this file's earlier
notes predicted for that row.

So the system is closed instead by Rheinboldt's **local parameterisation**: pick
the index `k` of the largest weighted tangent component and advance *that*
coordinate by a fixed amount, leaving λ free.

    F(x, λ)                 = 0                (NC rows, the model)
    z[k] − (z_prev[k] + δ)  = 0                (1 row, δ = Δs·τ[k])

The border row is a single nonzero — no fill-in, no dense kernels, and a
factorisation that stays genuinely sparse. The augmented Jacobian

    A = [  J     ∂F/∂λ ]
        [ e_kᵀ     0   ]

is **nonsingular at the fold** even though `J` is singular there, provided the
null direction of `J` has a nonzero `k` component — which is guaranteed by `k`
being the argmax of the tangent, since at a fold the tangent *is* the null
direction. That is the whole point of the method, and it is why this file does not
reuse `solve_newton!` (which needs a square system in x alone with λ fixed).

Near the turn `k` migrates to whichever variable is driving the fold, and λ's
tangent component passes through zero and changes sign. `k` is printed at every
step, so the fold direction is visible live rather than inferred afterwards.

The **bordering/block-elimination** shortcut (factor `J`, solve twice, combine) is
deliberately NOT used: it factors exactly the matrix that is singular at the fold,
which would throw away the method's only advantage.

## What this can and cannot deliver

It can traverse the turn and trace the branch's far side, and it reports where the
turn is and which way the curve goes after it.

It **cannot make an unreachable target reachable.** If λ genuinely folds back at
`λ* = 0.9908`, then `blabnat = 0.97` does not exist on this branch — the honest
result is "turned at λ*, target unattainable in this closure", not a solution.
Read `turned` and `lam_fold` before reading anything else.

## Other implementation notes

* **Row scaling is computed once and frozen**, for the reason recorded at
  `solve_newton!.jl:273` — a yardstick recomputed between readings changes length
  between them, and both the convergence test and the step-acceptance test compare
  readings taken at different points.
* **Direct sparse LU, never CGNR.** §N.10: at this conditioning CGNR returned a
  step of magnitude 0.025 against a true 1221.
* The line search uses the **2-norm merit** `½‖F‖₂² + ½n²` with Armijo, not
  strict descent on `‖F‖∞`. PLAN.md "Step 2": over 83,515 rows a single
  high-curvature row (`gret[2,5]`, ~1.1e-3 of pure log-curvature at a full step)
  can veto a step that improves everything else. `‖F‖∞` remains the *convergence*
  yardstick.

The model is left at the last ACCEPTED point — a genuine solution, with `λ` fixed
back at its value there — so a caller can inspect or restart from it.
"""

"""
Result of an [`arclength_solve!`](@ref) run.

`reached_target` is the only success flag. `turned` says the branch passed a limit
point in λ (the tangent's λ component changed sign) and `lam_fold` locates it — if
`turned` is true and `reached_target` is false, the target is **not on this
branch**, and no solver will find it. `path` records `(λ, iters, ‖F‖∞_scaled)` for
every accepted point. `local_coord` names the local parameterisation coordinate at
the last step, which near a fold is the variable driving it.
"""
struct ArclengthResult
    reached_target::Bool
    lam_reached::Float64
    target::Float64
    nsteps::Int
    nrejects::Int
    residual::Float64
    turned::Bool
    lam_fold::Float64
    ds_final::Float64
    path::Vector{Tuple{Float64,Int,Float64}}
    local_coord::String
end

"""
    arclength_solve!(m, vars, shock_vr, target; kwargs...)

Trace the solution branch from the current point toward `target` by continuation
in a response-weighted arclength with a local coordinate (see this file's header),
and return an [`ArclengthResult`](@ref).

The model must already be at a converged solution with `shock_vr` **fixed** — call
`solve_newton!` first. That solve is also what squares the system (it pins orphan
variables and deletes dead rows), which this routine relies on.

Keyword arguments:
- `ds0`/`dsmin`/`dsmax`  step in the weighted metric (defaults `0.005`, `1e-7`,
  `0.02`). Because the metric is built from `dz/dλ` at the start point, `Δs` is
  approximately a step in λ *while the branch is straight*; near the fold almost
  all of it is spent moving `x` instead, which is the intended behaviour. Scale
  these to the shock, not to 1.
- `tol`     residual tolerance, on the row-scaled `‖F‖∞`, matching
  `solve_newton!`'s yardstick (default `1e-8`).
- `maxit`   corrector iterations per step (default `20`).
- `maxsteps` give up after this many accepted steps (default `200`).
- `grow_iters` grow `ds` after a step that converged in at most this many
  corrector iterations (default `4`).
- `report`  variable names (or `(name, i, j)` tuples) to print each step.
"""
function arclength_solve!(m::JuMP.Model, vars::Dict{String,Any},
                          shock_vr::JuMP.VariableRef, target::Float64;
                          ds0::Float64=0.005, dsmin::Float64=1e-7, dsmax::Float64=0.02,
                          tol::Float64=1e-8, maxit::Int=20, maxsteps::Int=200,
                          grow_iters::Int=4, report=(), verbose::Bool=true)
    # flush(stdout) is NOT optional here. Julia block-buffers stdout when it is
    # redirected to a file, so without it a run of this length is completely
    # opaque: an empty log looks identical to a hung job, and the only way to find
    # out whether it was progressing is to let it finish. That cost 2.4 hours of
    # blind waiting once already.
    say(msg) = verbose && (println(msg); flush(stdout))

    JuMP.is_fixed(shock_vr) ||
        error("shock variable must be FIXED at its current value before continuation")

    # ── index variables ─────────────────────────────────────────────────────
    allv = JuMP.all_variables(m)
    allidx = [JuMP.index(v) for v in allv]
    pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))

    nowfree = VariableRef[]
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) || push!(nowfree, vr)
    end
    N = length(nowfree)
    freecols = [pos[JuMP.index(v).value] for v in nowfree]
    col2free = zeros(Int, length(allv))
    for (i, c) in enumerate(freecols); col2free[c] = i; end
    lamcol = pos[JuMP.index(shock_vr).value]
    col2free[lamcol] == 0 ||
        error("the continuation parameter is ENDOGENOUS — it cannot also be a free unknown")

    # names for the local coordinate, so the fold direction is visible live
    freename = [JuMP.name(v) for v in nowfree]
    coordname(k) = k == N + 1 ? JuMP.name(shock_vr) :
                   (isempty(freename[k]) ? "free#$k" : freename[k])

    # ── evaluator over the (already squared) model ──────────────────────────
    nlmodel = MOI.Nonlinear.Model()
    rhs = Float64[]
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            co = JuMP.constraint_object(con)
            r = _rhs_of(co.set)
            MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r))
            push!(rhs, r)
        end
    end
    NC = length(rhs)
    NC == N || error("""
        continuation needs NC == N with λ fixed, got $NC equations vs $N free variables.
        Run solve_newton! first — it is what pins orphan variables and deletes dead
        rows. See PLAN.md "Gate #5c RESOLVED" if the difference is genuine.""")

    evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
    MOI.initialize(evaluator, [:Grad, :Jac])

    st = MOI.jacobian_structure(evaluator)
    jr_all = getindex.(st, 1); jc_all = getindex.(st, 2)
    Jbuf = zeros(length(jr_all))
    keepm = [col2free[c] != 0 for c in jc_all]
    jr = jr_all[keepm]; jc = [col2free[c] for c in jc_all[keepm]]
    lamk = findall(==(lamcol), jc_all)
    lamrows = jr_all[lamk]
    say("  arclength: $NC equations, $N free + λ = $(N+1) unknowns, " *
        "∂F/∂λ has $(length(lamk)) structural nonzero(s)")

    # ── state, in SCALED coordinates ────────────────────────────────────────
    xfull = Vector{Float64}(undef, length(allv))
    for (k, v) in enumerate(allv)
        xfull[k] = JuMP.is_fixed(v) ? JuMP.fix_value(v) :
                   (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v))
    end
    # max(|x|, 1): relative for value flows, absolute for the del* variables that
    # sit legitimately at zero and whose *ratios* are meaningless (PLAN.md §F).
    # This governs state representation, row scaling and the bound logic ONLY —
    # the arclength metric is `wt`, built separately below, because using cscale
    # as a metric is precisely the defect described in this file's header.
    cscale = [max(abs(xfull[freecols[i]]), 1.0) for i in 1:N]
    lscale = max(abs(xfull[lamcol]), 1.0)
    lb = [JuMP.has_lower_bound(v) ? JuMP.lower_bound(v) : -Inf for v in nowfree]

    z = Vector{Float64}(undef, N + 1)                 # z = (y, μ) scaled
    for i in 1:N; z[i] = xfull[freecols[i]] / cscale[i]; end
    z[N+1] = xfull[lamcol] / lscale

    function push_state!(zv)
        @inbounds for i in 1:N; xfull[freecols[i]] = zv[i] * cscale[i]; end
        xfull[lamcol] = zv[N+1] * lscale
        return nothing
    end
    lam_of(zv) = zv[N+1] * lscale

    gbuf = zeros(NC)
    function residual!(zv)
        push_state!(zv)
        MOI.eval_constraint(evaluator, gbuf, xfull)
        @inbounds for i in 1:NC; gbuf[i] -= rhs[i]; end
        return gbuf
    end

    # ── row scaling: computed ONCE at the start point, then FROZEN ───────────
    # (solve_newton!.jl:273 — a ruler that changes length between readings makes
    # rescaling look like progress.)
    push_state!(z)
    MOI.eval_constraint_jacobian(evaluator, Jbuf, xfull)
    rowscale = ones(NC)
    let mx = zeros(NC)
        @inbounds for k in eachindex(jr_all)
            c = jc_all[k]
            s = c == lamcol ? lscale : (col2free[c] != 0 ? cscale[col2free[c]] : 0.0)
            s == 0.0 && continue
            a = abs(Jbuf[k] * s)
            a > mx[jr_all[k]] && (mx[jr_all[k]] = a)
        end
        @inbounds for i in 1:NC; rowscale[i] = mx[i] > 1e-12 ? mx[i] : 1.0; end
    end

    """
    Augmented matrix `[J ∂F/∂λ; rowᵀ]`, scaled, at the current `xfull`.

    `row` is passed as its sparsity pattern — `(ridx, rval)` — never as a dense
    vector. With a dense last row the factorisation blew memory
    (`OpenBLAS: malloc failed in gemm_driver`, see the header); the local
    parameterisation used here makes it a single nonzero.
    """
    function augmented(ridx::Vector{Int}, rval::Vector{Float64})
        MOI.eval_constraint_jacobian(evaluator, Jbuf, xfull)
        nk = length(jr); nb = length(ridx)
        tot = nk + length(lamk) + nb
        I = Vector{Int}(undef, tot); J = similar(I)
        V = Vector{Float64}(undef, tot)
        vk = Jbuf[keepm]
        @inbounds for k in 1:nk
            I[k] = jr[k]; J[k] = jc[k]
            V[k] = vk[k] * cscale[jc[k]] / rowscale[jr[k]]
        end
        p = nk
        @inbounds for (q, k) in enumerate(lamk)
            r = lamrows[q]
            p += 1; I[p] = r; J[p] = N + 1; V[p] = Jbuf[k] * lscale / rowscale[r]
        end
        @inbounds for t in 1:nb
            p += 1; I[p] = NC + 1; J[p] = ridx[t]; V[p] = rval[t]
        end
        # `+` combines duplicate (i,j) entries — the Jacobian structure can list a
        # coordinate more than once, and dropping the later one would be wrong.
        return sparse(I, J, V, NC + 1, N + 1, +)
    end

    # Both norms of the scaled residual, from ONE evaluation.
    #
    # `linf` is the convergence yardstick, matching solve_newton!'s. `l2sq` is the
    # MERIT function for the line search, and the difference is not cosmetic:
    # PLAN.md "Step 2" records that a strict-descent test on ‖F‖∞ over 83,515 rows
    # lets a SINGLE high-curvature row veto a step that improves everything else —
    # measured on `gret[2,5]`, which picks up ~1.1e-3 of pure log-curvature at a
    # full step, more than the entire starting residual. Newton is meant to absorb
    # that at the next iteration; a monotone ∞-norm test never lets it.
    function snorms(zv)
        residual!(zv)
        linf = 0.0; l2sq = 0.0
        @inbounds for i in 1:NC
            a = gbuf[i] / rowscale[i]
            aa = abs(a); aa > linf && (linf = aa)
            l2sq += a * a
        end
        return linf, l2sq
    end
    sres(zv) = first(snorms(zv))

    # ── bootstrap tangent, and the metric derived from it ────────────────────
    # Closing with `μ moves at unit rate` is exactly natural-parameter
    # continuation, valid HERE because the start point is far from the fold. The
    # raw solution IS `v = dz/dλ`, which is what the metric needs.
    e = zeros(NC + 1); e[NC+1] = 1.0
    push_state!(z)
    tlu0 = @elapsed v = lu(augmented([N + 1], [1.0])) \ e
    all(isfinite, v) ||
        error("bootstrap sensitivity dz/dλ is not finite — the start point is already singular")

    vmax = maximum(abs, v)
    floorv = 1e-6 * vmax
    wt = [1.0 / max(abs(v[i]), floorv) for i in 1:(N+1)]
    mnorm(t) = maximum(abs(wt[i] * t[i]) for i in 1:(N+1))

    tau = v ./ mnorm(v)
    (tau[N+1] * (target - lam_of(z)) < 0) && (tau .*= -1)   # point it at the target
    dir0 = sign(target - lam_of(z))                         # +1 λ rises, −1 λ falls

    # The metric and the local coordinate answer DIFFERENT questions, and must be
    # computed differently — conflating them is a trap, because `wt` normalises
    # every significant component to exactly 1.0, so an argmax over the *weighted*
    # tangent is a mass tie broken by index order (it would pick the first variable
    # above the floor, which means nothing).
    #
    #  * `wt` sets the STEP SIZE, so that no variable dominates on the strength of
    #    its unit of measurement and λ actually progresses.
    #  * the local coordinate must keep `[J ∂F/∂λ; e_kᵀ]` nonsingular, which holds
    #    iff `τ[k] ≠ 0`; the safest margin is therefore the argmax of the tangent
    #    in the UNWEIGHTED scaled coordinates.
    #
    # BUT "unweighted scaled coordinates" still includes the del*-at-zero trap this
    # file's header (§1) says was fixed: `wt` only protects `ds` sizing (§1's
    # ‖wt⊙τ‖∞=1 is a property of the AGGREGATE, not of any one component), while
    # `mnorm(v)` on the bootstrap tangent is 1 by construction (wt is DERIVED from
    # v, so ‖wt⊙v‖∞=1 trivially) — so `tau = v ./ mnorm(v)` leaves `v` numerically
    # untouched. A del* reporting variable's raw dz/dλ is Rp-billions, so a plain
    # argmax over ALL components picks it every time, and the closing displacement
    # `ds*tau[kloc]` then forces that dummy accounting variable to jump by
    # `ds*(that huge derivative)` in one step — hundreds of units, in a variable
    # with no economic content — dragging the corrector to an unrelated root.
    # Measured 2026-07-31 (`test/diag_arclength_direction.jl`): `kloc` locked onto
    # `delPGDPEXP[2,1]` (bootstrap |v|=1.204e5) for every step of a 12-step run, the
    # very first accepted step moved λ the WRONG way (0.991→0.99663, target 0.97),
    # and the reported "turning point" at λ≈0.9966 was nowhere near the documented
    # fold at λ*≈0.9908083 (PLAN.md §N.19c) — because the branch being traced was
    # never the real one.
    #
    # `cscale[i] > 1.0` is NOT a fix for this: cscale is a max(|current value|, 1)
    # STATE-SCALING factor, re-evaluated at whatever point continuation is
    # currently at — a del* variable that has drifted away from its benchmark
    # zero mid-shock has |value| > 1 same as a genuine flow, so a cscale-based
    # mask lets it right back in (confirmed: `driver()`, which used exactly this
    # mask, was ALSO reporting delPGDPEXP[2,1] in the pre-fix run — the mask was
    # never excluding it). The only reliable signal is the model's own naming
    # convention: every additive reporting/adjustment variable in this codebase
    # is named with a `del` prefix (delPGDPEXP, delGDPINC, delempratio, ... —
    # see build_model!.jl/build_equations.jl/build_dynamics!.jl), independent of
    # its current magnitude. `flowmask` is therefore name-based, using `freename`
    # (already computed above for `coordname`), and restricts `locidx` candidates
    # to non-del* variables plus λ itself (index N+1, never masked — natural-
    # parameter continuation is legitimate whenever λ is the right local
    # coordinate).
    flowmask = BitVector(!startswith(freename[i], "del") for i in 1:N)
    kmask = vcat(flowmask, true)
    function locidx(t)
        best = 0; bv = -1.0
        @inbounds for i in 1:(N+1)
            kmask[i] || continue
            a = abs(t[i]); a > bv && (bv = a; best = i)
        end
        best == 0 && error("locidx: no candidate coordinate (flowmask all false and λ excluded — cannot happen)")
        return best
    end
    function driver(t)
        best = 0; bv = -1.0
        @inbounds for i in 1:N
            flowmask[i] || continue
            a = abs(t[i]); a > bv && (bv = a; best = i)
        end
        return best == 0 ? "—" : "$(coordname(best)) ($(round(bv; sigdigits=4)))"
    end

    # Report the metric, because getting this wrong is what killed the first run:
    # if `dλ/ds` is ~1e-7 again, the metric is still being eaten by something and
    # the run should be stopped rather than left to grind for weeks.
    say("  metric: bootstrap LU $(round(tlu0; digits=1))s, ‖dz/dλ‖∞=$(round(vmax; sigdigits=4)) " *
        "(2-norm $(round(norm(v); sigdigits=4)) — an unweighted metric would give " *
        "dλ/ds≈$(round(1/norm(v); sigdigits=3)))")
    say("  metric: weighted  dλ/ds=$(round(tau[N+1]*lscale; sigdigits=4))  ⟹  " *
        "≈$(ceil(Int, abs(target - lam_of(z)) / max(abs(tau[N+1]*lscale), 1e-30) / dsmax)) " *
        "step(s) at dsmax=$dsmax if the branch were straight")

    say("arclength: λ $(round(lam_of(z); sigdigits=10)) → $target  (ds0=$ds0, tol=$tol)")

    ds = clamp(ds0, dsmin, dsmax)
    nsteps = 0; nrejects = 0; res = sres(z)
    turned = false; lam_fold = NaN
    path = Tuple{Float64,Int,Float64}[]
    reached = false
    kloc = locidx(tau)

    # ── main loop ───────────────────────────────────────────────────────────
    while nsteps < maxsteps
        z_prev = copy(z); tau_prev = copy(tau); lam_prev = lam_of(z)

        # If this step would carry λ past the target, finish with an exact
        # natural-parameter landing instead: close with μ = μ*.
        lam_pred = lam_prev + ds * tau[N+1] * lscale
        land = (target - lam_prev) * (target - lam_pred) <= 0

        # Local coordinate: the largest tangent component. At a fold the tangent IS
        # the null direction of J, so its argmax is exactly the index that keeps
        # the augmented matrix nonsingular there.
        kloc = locidx(tau)

        zc = copy(z)
        if land
            zc[N+1] = target / lscale
            ridx = [N + 1]; rval = [1.0]
            nfun = zv -> zv[N+1] - target / lscale
        else
            zc .= z .+ ds .* tau                       # Euler predictor
            δ = ds * tau[kloc]
            ridx = [kloc]; rval = [1.0]
            nfun = zv -> zv[kloc] - (z_prev[kloc] + δ)
        end

        # ── corrector ───────────────────────────────────────────────────────
        # Tolerance on the closing row. In the `land` case it pins λ to the target
        # and must be tight — the answer's λ IS the deliverable. Otherwise it only
        # says how far along the curve we stopped, and being off by a whisker of a
        # step is physically meaningless, so a tight tolerance there buys nothing
        # and risks rejecting a perfectly good solution of F = 0.
        ntol = land ? 1e-12 : 1e-8 * max(abs(z_prev[kloc]), 1.0)

        ok = false; it = 0; rr = NaN
        # Kept so the tangent solve at the end of an accepted step can REUSE this
        # factorisation instead of building a fresh one: the tangent equation's
        # matrix is `[J ∂F/∂λ; e_klocᵀ]`, exactly the matrix already built here.
        # It is evaluated one corrector iterate early, but the tangent is a
        # predictor that gets renormalised, so that error is immaterial and it
        # saves a full LU per step (roughly a third of the run time).
        Flast = nothing
        for outer it in 1:maxit
            push_state!(zc)
            rr, l2 = snorms(zc)
            nres = abs(nfun(zc))
            if rr <= tol && nres <= ntol
                ok = true; break
            end
            tlu = @elapsed begin
                A = augmented(ridx, rval)
                F = try
                    lu(A)
                catch err
                    say("     corrector it=$it: factorisation failed ($(typeof(err))) — rejecting step")
                    nothing
                end
            end
            F === nothing && break
            Flast = F
            b = Vector{Float64}(undef, NC + 1)
            @inbounds for i in 1:NC; b[i] = -gbuf[i] / rowscale[i]; end
            b[NC+1] = -nfun(zc)
            dz = F \ b
            verbose && say("     it=$it  ‖F‖∞=$(round(rr; sigdigits=4))  " *
                           "n=$(round(nres; sigdigits=3))  LU $(round(tlu; digits=1))s")
            all(isfinite, dz) || break

            # Fraction-to-boundary, so a positivity-bounded variable is never
            # stepped through its bound (a domain exit shows up as NaN later, far
            # from the cause).
            #
            # The `gap <= tiny` branch is NOT optional, and mirrors
            # solve_newton!.jl:566-572 deliberately. A variable sitting exactly ON
            # its bound gives `gap = 0`, so the ratio test returns α = 0 — and with
            # 83,515 free variables at least one always does. The first version of
            # this file lacked the branch and therefore had α = 0 on every
            # iteration: ONE pinned variable vetoed the step for all the others.
            # Zeroing that single component instead keeps the rest of the Newton
            # direction usable, exactly as the working solver does.
            α = 1.0; nzeroed = 0; ibind = 0
            @inbounds for i in 1:N
                isfinite(lb[i]) || continue
                d = dz[i] * cscale[i]
                d < -1e-10 || continue
                gap = zc[i] * cscale[i] - lb[i]
                if gap > 1e-8
                    cap = 0.99 * gap / (-d)
                    if cap > 0 && cap < α; α = cap; ibind = i; end
                else
                    dz[i] = 0.0; nzeroed += 1
                end
            end

            # A collapsed α is a diagnosis, not a step size. Report the binding
            # variable rather than backtracking from a hopeless starting point.
            if α < 1e-6
                say("     it=$it: fraction-to-boundary collapsed α to " *
                    "$(round(α; sigdigits=3))" *
                    (ibind == 0 ? "" : " on $(coordname(ibind)) " *
                     "(gap $(round(zc[ibind]*cscale[ibind] - lb[ibind]; sigdigits=3)))") *
                    " — rejecting step")
                break
            end

            # Backtrack on the 2-norm merit ½‖F‖₂² + ½n², with an Armijo-style
            # sufficient-decrease factor rather than bare descent.
            #
            # `α >= αmin` is load-bearing: with α = 0 the Armijo test degenerates to
            # `base <= base`, which PASSES, so a zero-length step would be recorded
            # as a successful corrector iteration. That is precisely what happened
            # in run 2 — the residual stayed bit-identical for all 20 iterations at
            # ~2s of LU each while the corrector did nothing at all. A no-op must
            # never be able to masquerade as progress.
            αmin = 1e-10
            base = l2 + nfun(zc)^2
            zt = similar(zc); accepted = false
            for _ in 1:14
                α >= αmin || break
                zt .= zc .+ α .* dz
                rt, l2t = snorms(zt)
                mt = l2t + nfun(zt)^2
                if isfinite(mt) && mt <= (1 - 1e-4 * α) * base
                    zc .= zt; accepted = true; break
                end
                α /= 2
            end
            if !accepted
                say("     it=$it: no merit decrease down to α=$(round(α; sigdigits=3)) " *
                    "(zeroed $nzeroed at-bound component(s)) — rejecting step")
                break
            end
        end

        if ok
            nsteps += 1; z .= zc; res = rr
            lam = lam_of(z)
            push!(path, (lam, it, rr))
            # write back so the JuMP model tracks the branch
            push_state!(z)
            for (i, v2) in enumerate(nowfree); JuMP.set_start_value(v2, xfull[freecols[i]]); end
            JuMP.fix(shock_vr, lam; force=true)

            say("  ✅ λ=$(round(lam; sigdigits=10))  it=$it  " *
                "‖F‖∞=$(round(rr; sigdigits=4))  ds=$(round(ds; sigdigits=3))\n" *
                "       coord=$(coordname(kloc))   driver=$(driver(tau))   " *
                "dλ/ds=$(round(tau[N+1]*lscale; sigdigits=4))")
            for spec in report
                nm, idx = spec isa Tuple ? (spec[1], spec[2:end]) : (spec, ())
                haskey(vars, nm) || continue
                vv = vars[nm]
                vr = isempty(idx) ? (vv isa AbstractArray ? first(vv) : vv) : vv[idx...]
                val = JuMP.is_fixed(vr) ? JuMP.fix_value(vr) : JuMP.start_value(vr)
                say("       $(rpad(string(nm, idx), 18)) = $(round(val; sigdigits=8))")
            end

            if land
                reached = true
                break
            end

            # New tangent from the AUGMENTED system — nonsingular through the turn.
            # With the closing row `e_klocᵀ`, solving `A τ = e_{NC+1}` gives the
            # tangent normalised to τ[kloc] = 1; the metric normalisation follows.
            tnew = try
                if Flast === nothing
                    ttan = @elapsed Ft = lu(augmented(ridx, rval))
                    verbose && say("     tangent LU $(round(ttan; digits=1))s (fresh)")
                    Ft \ e
                else
                    Flast \ e
                end
            catch
                say("  ⛔ tangent factorisation failed at λ=$(round(lam; sigdigits=10))")
                break
            end
            mn = all(isfinite, tnew) ? mnorm(tnew) : NaN
            (isfinite(mn) && mn > 0) || begin
                say("  ⛔ tangent solve returned a degenerate vector at " *
                    "λ=$(round(lam; sigdigits=10)) — the augmented system is singular here, " *
                    "which a fold alone cannot cause (suspect a bifurcation)")
                break
            end
            tnew ./= mn
            dot(tnew, tau) < 0 && (tnew .*= -1)      # keep travelling the same way
            if !turned && tnew[N+1] * tau[N+1] < 0
                turned = true
                # The sign change is DETECTED one step past the extremum, so `lam` here
                # is already back on the far side and overstates how far the branch got.
                # λ* is the extremal λ actually attained in the direction of travel.
                lams = [p[1] for p in path]
                lam_fold = dir0 < 0 ? minimum(lams) : maximum(lams)
                say("  ↩ TURNING POINT: λ* ≈ $(round(lam_fold; sigdigits=10)) — the tangent's " *
                    "λ component changed sign. λ now moves BACK; the branch folds here.\n" *
                    "     fold direction: $(driver(tnew))   (local coordinate: $(coordname(kloc)))")
            end
            tau = tnew
            it <= grow_iters && (ds = min(dsmax, 2ds))
        else
            nrejects += 1
            z .= z_prev; tau .= tau_prev
            push_state!(z)
            JuMP.fix(shock_vr, lam_prev; force=true)
            ds /= 2
            say("  ✗ step from λ=$(round(lam_prev; sigdigits=10)) failed " *
                "(‖F‖∞=$(round(rr; sigdigits=4))) — ds → $(round(ds; sigdigits=3))")
            if ds < dsmin
                say("  ⛔ step collapsed below dsmin at λ=$(round(lam_prev; sigdigits=10)). " *
                    "A fold cannot defeat this method, so this is NOT one — suspect a " *
                    "bifurcation, a domain boundary, or a genuinely singular augmented system.")
                break
            end
        end
    end

    # restore start values to the last accepted point
    push_state!(z)
    for (i, v2) in enumerate(nowfree); JuMP.set_start_value(v2, xfull[freecols[i]]); end
    JuMP.fix(shock_vr, lam_of(z); force=true)

    if reached
        say("arclength: reached target in $nsteps step(s), $nrejects rejection(s)" *
            (turned ? " — after turning at λ ≈ $(round(lam_fold; sigdigits=10))" : ""))
    elseif turned
        say("arclength: TURNED at λ ≈ $(round(lam_fold; sigdigits=10)) without reaching $target.\n" *
            "  λ folds back there, so `$(JuMP.name(shock_vr)) = $target` does not exist on this\n" *
            "  branch. That is a property of the model and closure, not of the solver.")
    else
        say("arclength: stopped at λ=$(round(lam_of(z); sigdigits=10)) after $nsteps step(s) " *
            "without reaching $target and without turning.")
    end
    return ArclengthResult(reached, lam_of(z), target, nsteps, nrejects, res,
                           turned, lam_fold, ds, path, coordname(kloc))
end
