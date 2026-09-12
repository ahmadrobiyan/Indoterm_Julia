"""
Natural-parameter continuation for the levels model.

The Newton solver converges a shock in 2–3 iterations *provided it starts near
the answer* (`test/verify_continuation.jl`: 1e-6, 1e-5 and 1e-4 shocks in
`blabnat` all reach ‖F‖∞ ~ 1e-9 in two or three iterations). It does not
converge the same shock in one leap from the benchmark. Continuation is the
standard remedy: walk the shock parameter from its benchmark value to the target
in steps small enough that each Newton solve starts inside its own basin.

`solve_newton!` writes its solution back as each variable's start value, so the
next step warm-starts from the previous solution automatically — the driver only
has to choose the step size and roll back when a step fails.

The step size is adaptive because the required step is NOT uniform along the
path: the response accelerates sharply where a cell's fixed-factor price starts
to move (see the memory note on cell [2,5]). Steps that converge easily grow;
steps that fail are halved and retried from the last good point.

**Reading a failure.** If `h` collapses to `hmin` at a specific parameter value
while every earlier step converged in 2–3 iterations, the obstruction is at that
point on the branch, not in the stepping — that is the signature of a fold
(limit point), which natural-parameter continuation provably cannot pass. The
returned `s_reached` is then the useful number: it locates the fold. Getting
past a genuine fold needs pseudo-arclength continuation, which parameterises by
arclength instead of by the shock and so stays nonsingular through the turn.
"""

"""
Result of a [`continuation_solve!`](@ref) run.

`reached_target` is the only success flag. `s_reached` is the furthest shock
value actually solved; on failure it locates where the branch became hard.
`path` records `(s, iters, residual)` for every ACCEPTED step, so the response
can be plotted against the parameter.
"""
struct ContinuationResult
    reached_target::Bool
    s_reached::Float64
    target::Float64
    nsteps::Int
    nrejects::Int
    residual::Float64
    h_final::Float64
    path::Vector{Tuple{Float64,Int,Float64}}
end

# Snapshot / restore every variable's start value, so a failed step can be
# rolled back. solve_newton! overwrites start values in place, so a rejected
# step would otherwise leave the model sitting at a non-solution.
_snapshot_starts(m::JuMP.Model) =
    Float64[(sv = JuMP.start_value(v); sv === nothing ? 0.0 : sv)
            for v in JuMP.all_variables(m)]

function _restore_starts!(m::JuMP.Model, snap::Vector{Float64})
    for (k, v) in enumerate(JuMP.all_variables(m))
        JuMP.is_fixed(v) || JuMP.set_start_value(v, snap[k])
    end
    return nothing
end

"""
    continuation_solve!(m, vars, shock_vr, target; kwargs...)

Walk `shock_vr` from its current fixed value to `target`, running a full Newton
solve at every step, and return a [`ContinuationResult`](@ref).

Keyword arguments:
- `h0`        initial step as a fraction of the total span (default `0.25`).
- `hmin`      give up below this fraction (default `1e-4`).
- `hmax`      never step further than this fraction (default `1.0`).
- `tol`       residual tolerance passed to `solve_newton!` (default `1e-8`).
- `maxit`     Newton iterations per step (default `30`).
- `grow_iters` grow `h` when a step converged in at most this many iterations
  (default `4`). Growth is what keeps an easy stretch of the path cheap.
- `report`    optional `Vector{String}` of variable names to print each step;
  entry `(name, i, j)` tuples are also accepted for indexed variables.

The model is left at the last ACCEPTED point — on failure that is a genuine
solution of the system at `s_reached`, not a half-converged iterate, so a caller
can inspect or restart from it.
"""
function continuation_solve!(m::JuMP.Model, vars::Dict{String,Any},
                             shock_vr::JuMP.VariableRef, target::Float64;
                             h0::Float64=0.25, hmin::Float64=1e-4,
                             hmax::Float64=1.0, tol::Float64=1e-8,
                             maxit::Int=30, grow_iters::Int=4,
                             report=(), verbose::Bool=true,
                             linsolve::Symbol=:lu)
    log(msg) = verbose && println(msg)

    s = JuMP.fix_value(shock_vr)
    span = target - s
    path = Tuple{Float64,Int,Float64}[]
    if abs(span) < eps()
        return ContinuationResult(true, s, target, 0, 0, NaN, 0.0, path)
    end

    h = clamp(h0, hmin, hmax)
    nsteps = 0; nrejects = 0; res = NaN
    log("continuation: $s → $target  (h0=$h0, tol=$tol)")

    while true
        # Never overshoot: the last step lands exactly on the target.
        s_try = s + clamp(h, 0.0, 1.0) * span
        if (span > 0 && s_try > target) || (span < 0 && s_try < target)
            s_try = target
        end

        snap = _snapshot_starts(m)
        JuMP.fix(shock_vr, s_try; force=true)
        el = @elapsed r = solve_newton!(m, vars; maxit=maxit, tol=tol, verbose=false,
                                          linsolve=linsolve)

        if r.residual <= tol
            nsteps += 1; s = s_try; res = r.residual
            push!(path, (s, r.iters, r.residual))
            log("  ✅ s=$(round(s; sigdigits=10))  it=$(r.iters)  " *
                "‖F‖∞=$(round(r.residual; sigdigits=4))  h=$(round(h; sigdigits=3))  " *
                "$(round(el; digits=1))s")
            for spec in report
                nm, idx = spec isa Tuple ? (spec[1], spec[2:end]) : (spec, ())
                haskey(vars, nm) || continue
                v = vars[nm]
                vr = isempty(idx) ? (v isa AbstractArray ? first(v) : v) : v[idx...]
                val = JuMP.is_fixed(vr) ? JuMP.fix_value(vr) : JuMP.start_value(vr)
                log("       $(rpad(string(nm, idx), 18)) = $(round(val; sigdigits=8))")
            end
            s == target && break
            # Only grow after a genuinely easy step. Growing on a step that
            # merely scraped in re-triggers the failure it just escaped.
            r.iters <= grow_iters && (h = min(hmax, 2h))
        else
            nrejects += 1
            _restore_starts!(m, snap)          # back to the last real solution
            JuMP.fix(shock_vr, s; force=true)
            h /= 2
            log("  ✗ s=$(round(s_try; sigdigits=10)) failed ($(r.status), " *
                "‖F‖∞=$(round(r.residual; sigdigits=4))) — h → $(round(h; sigdigits=3))")
            if h < hmin
                log("  ⛔ step size collapsed below hmin at s=$(round(s; sigdigits=10)); " *
                    "the branch is obstructed there, not the stepping")
                return ContinuationResult(false, s, target, nsteps, nrejects, res, h, path)
            end
        end
    end

    log("continuation: reached target in $nsteps step(s), $nrejects rejection(s)")
    return ContinuationResult(true, s, target, nsteps, nrejects, res, h, path)
end
