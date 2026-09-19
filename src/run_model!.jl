"""
Gate 7 — the scenario driver.

`continuation_solve!` walks ONE `VariableRef` from its benchmark value to a
target. That is the right primitive for a single-instrument experiment (and it
is what located the fold in PLAN.md §J), but a policy scenario generally moves
several exogenous variables at once — a tariff change is one shift variable per
commodity, a regional productivity scenario is one per (sector, region). Walking
them one after another is NOT the same experiment: it forces the model through
intermediate states that no one asked about, and those states can be harder than
the target itself.

`run_model!` therefore walks every shock **together** along a scalar homotopy
`t ∈ [0, 1]`, with each shocked variable held at `base + t·(target − base)`. At
`t = 0` the model is the benchmark, so step one is always a solve from an exact
solution; at `t = 1` every shock is fully applied. Step size is adaptive with
rollback, exactly as in `continuation_solve!`.

**Closure first, benchmark second, shock third.** The driver applies the closure
swaps, then *re-solves the benchmark* before touching any shock, and refuses to
continue if that solve does not reach `tol`. A swapped closure that cannot
reproduce its own base year is a closure bug, and every number produced past
that point would be measuring it rather than the scenario. This is the same
regression gate `closures.jl` describes; making it unskippable is the point.

**Which closure.** The default is `TERM_CMF_SWAPS`, the closure TERM.CMF
actually specifies, because that is the faithful one *and* the better-behaved
one: after the flow-ratio aggregation fix it walks `blabnat` to 0.999 in one or
two Newton iterations on a linear branch, while the bare automatic closure
(`swaps = ()`) folds at `blabnat ≈ 0.999752` with `xcap`/`xlnd` frozen. Pass
`swaps = ()` only when you deliberately want the static no-accumulation closure,
and say so when quoting any number from it — the two closures read the same cell
43× apart.

Usage:

    agg, params = ...                       # Steps 0–4
    r = run_model!(agg, params;
                   shocks = ["blabnat" => pct(-3)],   # TERM.CMF's own scenario
                   report = [("pcap", 2, 5)])
    print_gdp_report(r.gdp; regions = REG6)
"""

"""
    pct(p) -> Float64

Convert a GEMPACK percentage shock to the levels-model target for a variable
whose benchmark value is 1. `Shock blabnat = -3` in TERM.CMF is `pct(-3)` = 0.97
here. Use it to keep scenario files readable in the units the source states.

For a variable whose benchmark value is NOT 1, pass `:pct` targets via
[`run_model!`](@ref)'s `pct_shocks` instead, which scales each variable's own
base value.
"""
pct(p::Real) = 1.0 + p / 100

"""
    logpct(p) -> Float64

Convert a GEMPACK percentage shock to the levels-model target for an **additive,
log-space** shifter — benchmark 0, no lower bound, entering its equation as a
bare summand alongside other log terms. `logpct(50)` = log(1.5) = 0.405465.

`pct` and `logpct` are NOT interchangeable and the model will not tell you which
one you wanted: both produce a finite, plausible-looking solve. Which applies is
decided by how the port declared the variable, and the declaration is the test:

- `@variable(m, blabnat >= 1e-6)` — lower-bounded, so `_is_ratio_type`, so
  benchmark 1, and `E_alab_o` reads `log(alab_o) = log(blabnat) + …`. The TABLO
  sum of percentage changes became a product of levels ratios ⇒ **`pct`**.
- `@variable(m, fpexp_d[1:na])` — unbounded, so benchmark 0, and `E_xexpd` keeps
  it as a summand next to `pfexp`, which the port defines outright as
  `log(ppur) − log(phi)`. The whole term is in logs ⇒ **`logpct`**.

Passing `pct(50)` = 1.5 to a log-space shifter would silently apply
`e^1.5 − 1` ≈ **+348%** instead of the +50% the source asked for.
"""
logpct(p::Real) = log1p(p / 100)

# A shock target is `name => value` or `(name, i, j, …) => value`.
_shock_ref(vars::Dict{String,Any}, spec) = begin
    nm, idx = spec isa Tuple ? (String(spec[1]), spec[2:end]) : (String(spec), ())
    haskey(vars, nm) || error("shock \"$nm\" is not a declared variable")
    v = vars[nm]
    vr = if isempty(idx)
        v isa AbstractArray || return (nm, (), v)
        length(v) == 1 || error(
            "shock \"$nm\" has $(length(v)) elements — give an index, e.g. (\"$nm\", 1, 2)")
        first(v)
    else
        v[idx...]
    end
    (nm, idx, vr)
end

_shock_label(nm, idx) = isempty(idx) ? nm : "$nm$(idx)"

# Read the model's CURRENT point into the same `Dict{String,Any}` shape that
# `solve_newton!` returns, so a caller cannot tell from `r.values` which path
# produced the result. Both `solve_newton!` and `arclength_solve!` write their
# accepted solution back into start values, so this is that solution exactly —
# but only ever call it right after an ACCEPTED point, never mid-iteration, or it
# reports a half-converged iterate as if it were a solution.
function _current_values(m::JuMP.Model, vars::Dict{String,Any})
    out = Dict{String,Any}()
    for (nm, v) in vars
        if v isa AbstractArray
            a = Array{Float64}(undef, size(v))
            for idx in eachindex(v)
                vr = v[idx]
                a[idx] = JuMP.is_fixed(vr) ? JuMP.fix_value(vr) : JuMP.start_value(vr)
            end
            out[nm] = a
        else
            out[nm] = JuMP.is_fixed(v) ? JuMP.fix_value(v) : JuMP.start_value(v)
        end
    end
    out
end

"""
Result of a [`run_model!`](@ref) run.

`solved` is the only success flag. `t_reached` is the fraction of the scenario
actually applied — on failure it says how far the branch went before it became
obstructed, which is the useful number for diagnosing a fold. `values` and `gdp`
are from the last ACCEPTED point, so they are a genuine solution of the system
at `t_reached`, never a half-converged iterate.
"""
struct ScenarioResult
    name::String
    solved::Bool
    t_reached::Float64
    nsteps::Int
    nrejects::Int
    residual::Float64
    path::Vector{Tuple{Float64,Int,Float64}}
    values::Dict{String,Any}
    gdp::Union{Nothing,GDPReport}
    bench_residual::Float64
end

"""
    run_model!(agg, params; kwargs...) -> ScenarioResult

Build the model from an aggregated database, close it, verify the benchmark,
then walk a scenario in with homotopy continuation and report GDP.

Keyword arguments:
- `shocks`     `Vector` of `spec => target`, where `spec` is a variable name or
  a `(name, indices...)` tuple and `target` is the LEVELS value to reach.
  Percentage shocks on `bmk = 1` variables read naturally as `pct(-3)`.
- `pct_shocks` `Vector` of `spec => percent`, applied relative to each
  variable's OWN base value (`base * (1 + p/100)`). Use for variables whose
  benchmark level is not 1.
- `swaps`      closure swaps (default `TERM_CMF_SWAPS`; pass `()` for the bare
  automatic closure).
- `name`       label for the log and the result.
- `h0`/`hmin`/`hmax`/`grow_iters` homotopy step control, as in
  `continuation_solve!` but in units of `t`.
- `tol`/`maxit` passed to `solve_newton!`.
- `report`     variable specs to print at every accepted step.
- `gdp`        compute a `GDPReport` from the final point (default `true`).
- `arclength`  on a step collapse, retry with `arclength_solve!` (default `true`).
  Applies only when the scenario has exactly ONE shock — arclength traces a single
  `VariableRef` and `t` is not a model variable. With several shocks the driver
  says so and stops rather than pretending; run those in stages (see the fallback
  block for the recipe §N.19e used). Set `false` to measure the unaided homotopy,
  which is what a fold diagnosis needs.

Returns a [`ScenarioResult`](@ref). The model and `vars` are left at the last
accepted point, so a caller can inspect or restart from it.
"""
function run_model!(agg, params::Dict{String,Any};
                    shocks = (), pct_shocks = (),
                    swaps = TERM_CMF_SWAPS,
                    name::AbstractString = "scenario",
                    h0::Float64 = 0.25, hmin::Float64 = 1e-4, hmax::Float64 = 1.0,
                    tol::Float64 = 1e-8, maxit::Int = 30, grow_iters::Int = 4,
                    report = (), gdp::Bool = true, verbose::Bool = true,
                    arclength::Bool = true, numeraire::Symbol = :gdppi,
                    linsolve::Symbol = :lu)
    # Flush per line. Julia block-buffers stdout when it is redirected to a file, so an
    # unflushed log leaves a multi-hour homotopy walk indistinguishable from a hang --
    # the buffer only reaches disk when the process exits, which is precisely when the
    # progress trace has stopped being useful. Costs one syscall per logged line.
    log(msg) = verbose && (println(msg); flush(stdout))

    log("="^72)
    log("run_model!  —  $name")
    log("="^72)

    t = @elapsed ((m, vars) = build_model_full!(agg, params))
    log("build_model_full!: $(round(t; digits=1))s")
    log("numeraire: " * (numeraire === :gdppi ?
        ":gdppi — phi free, NatMacro(\"GDPPI\") pinned at 1.0" :
        numeraire === :cpi ?
        ":cpi — phi free, NatMacro(\"CPI\") pinned at 1.0" :
        ":exrate — phi exogenous, GDPPI free"))
    initialize_model!(m, vars; bmk_levels=benchmark_levels(params),
                      numeraire=numeraire)

    if !isempty(swaps)
        log("closure swaps:")
        apply_swaps!(m, vars, collect(swaps); verbose=verbose)
    else
        log("closure: bare automatic (no swaps) — xcap/xlnd frozen; " *
            "capital accumulation is INERT in this closure")
    end

    # ── Gate: the closed model must still reproduce its own base year ────────
    r0 = solve_newton!(m, vars; maxit=5, tol=tol, verbose=verbose, linsolve=linsolve)
    log("benchmark: status=$(r0.status)  ‖F‖∞=$(r0.residual)")
    r0.residual <= tol || error(
        "benchmark solve failed under this closure (‖F‖∞ = $(r0.residual) > $tol). " *
        "A closure that cannot reproduce the base year makes every scenario number " *
        "a measurement of the closure bug, not of the scenario.")

    # ── Resolve the shock list into (ref, base, target) ──────────────────────
    targets = Tuple{String,JuMP.VariableRef,Float64,Float64}[]
    for (spec, tgt) in shocks
        nm, idx, vr = _shock_ref(vars, spec)
        JuMP.is_fixed(vr) || error(
            "shock $(_shock_label(nm, idx)) is ENDOGENOUS under this closure — " *
            "a shock can only be applied to an exogenous variable. Check the swaps.")
        push!(targets, (_shock_label(nm, idx), vr, JuMP.fix_value(vr), Float64(tgt)))
    end
    for (spec, p) in pct_shocks
        nm, idx, vr = _shock_ref(vars, spec)
        JuMP.is_fixed(vr) || error(
            "shock $(_shock_label(nm, idx)) is ENDOGENOUS under this closure — " *
            "a shock can only be applied to an exogenous variable. Check the swaps.")
        b = JuMP.fix_value(vr)
        push!(targets, (_shock_label(nm, idx), vr, b, b * (1 + p / 100)))
    end

    if isempty(targets)
        log("no shocks given — returning the benchmark solve")
        rep = gdp ? calculate_gdp(r0.values, params) : nothing
        return ScenarioResult(String(name), true, 0.0, 0, 0, r0.residual,
                              Tuple{Float64,Int,Float64}[], r0.values, rep, r0.residual)
    end

    log("shocks:")
    for (lbl, _, b, tg) in targets
        log("  $(rpad(lbl, 20)) $(round(b; sigdigits=8)) → $(round(tg; sigdigits=8))")
    end

    function _set_t!(t_)
        for (_, vr, b, tg) in targets
            JuMP.fix(vr, b + t_ * (tg - b); force=true)
        end
    end

    # ── Homotopy walk in t ∈ [0,1] ──────────────────────────────────────────
    s = 0.0
    h = clamp(h0, hmin, hmax)
    nsteps = 0; nrejects = 0
    res = r0.residual
    vals = r0.values
    path = Tuple{Float64,Int,Float64}[]
    # P0 cumulative stall cap (HIGH review item): consecutive STALL-SIGNATURE
    # failures with no NET t-progress. A bad predictor (H1) clears in a few
    # halvings; a genuinely near-singular t (H5) never clears. Cap at K — on
    # breach, dump diagnostics and STOP, never crawl silently.
    # NET-progress gated (review symmetry fix): the reset fires ONLY when s
    # advances by more than STALL_T_EPS since the last reset. Alternating
    # trigger types or trivial micro-steps do NOT reset — without this, an
    # H5-true system crawls indefinitely with a lower duty cycle.
    n_stall_abort = 0
    STALL_ABORT_CAP = 8
    s_at_reset = 0.0
    STALL_T_EPS = 1e-9
    log("homotopy: t = 0 → 1  (h0=$h0, tol=$tol)")

    while true
        s_try = min(1.0, s + h)
        snap = _snapshot_starts(m)
        _set_t!(s_try)
        el = @elapsed r = solve_newton!(m, vars; maxit=maxit, tol=tol, verbose=verbose,
                                          linsolve=linsolve)

        if r.residual <= tol
            nsteps += 1; s = s_try; res = r.residual; vals = r.values
            push!(path, (s, r.iters, r.residual))
            log("  ✅ t=$(round(s; sigdigits=6))  it=$(r.iters)  " *
                "‖F‖∞=$(round(r.residual; sigdigits=4))  h=$(round(h; sigdigits=3))  " *
                "$(round(el; digits=1))s")
            for spec in report
                nm, idx, vr = _shock_ref(vars, spec)
                v = JuMP.is_fixed(vr) ? JuMP.fix_value(vr) : JuMP.start_value(vr)
                log("       $(rpad(_shock_label(nm, idx), 18)) = $(round(v; sigdigits=8))")
            end
            s >= 1.0 && break
            # Grow only after a genuinely easy step — growing on a step that
            # merely scraped in re-triggers the failure it just escaped.
            # P1 (H1 fix): cap growth at 1.5x, not 2x — one cheap success must
            # not double h into a predictor outside the corrector's basin.
            r.iters <= grow_iters && (h = min(hmax, 1.5h))
            # NET-progress-gated reset (review symmetry fix): only an s-advance
            # beyond eps since the last reset clears the stall counter.
            if s - s_at_reset > STALL_T_EPS
                n_stall_abort = 0; s_at_reset = s
            end
        else
            nrejects += 1
            _restore_starts!(m, snap)      # back to the last real solution
            _set_t!(s)
            # P0 stall-SIGNATURE accounting (review CRITICAL fix): the counter
            # tracks the STALL pattern, not the :no_progress label. Signature =
            # corrector burned maxit iters ending within 10x of tol but above
            # it (plateau), OR tripped the in-loop STALL-ABORT. A tight exact
            # LU solve whose step is then rejected (lin_rel~1e-5, H1/fold
            # territory — NOT H5) must NOT increment: H5 predicts a degraded
            # LINEAR solve, and lin_rel~1e-5 contradicts its own definition.
            is_stall_sig = r.status == :no_progress ||
                (r.status == :maxit && r.residual <= 10tol && r.residual > tol)
            if is_stall_sig
                n_stall_abort += 1
            elseif r.status != :no_progress
                # A different failure regime — but still net-gated: only a real
                # s-advance clears accumulated stall evidence.
                if s - s_at_reset > STALL_T_EPS
                    n_stall_abort = 0; s_at_reset = s
                end
            end
            h /= 2
            log("  ✗ t=$(round(s_try; sigdigits=6)) failed ($(r.status), " *
                "‖F‖∞=$(round(r.residual; sigdigits=4))) — h → $(round(h; sigdigits=3))" *
                (is_stall_sig ? "  [stall-sig $n_stall_abort/$STALL_ABORT_CAP]" : ""))
            if n_stall_abort >= STALL_ABORT_CAP
                log("  ⛔ STALL CAP BREACHED at t=$(round(s; sigdigits=6)): " *
                    "$STALL_ABORT_CAP consecutive STALL-SIGNATURE failures with no net t-progress. " *
                    "This is NOT a bad predictor (H1 clears in a few halvings) — " *
                    "suspect genuine near-singular Jacobian (H5) or floor " *
                    "non-smoothness (H6). Escalating to diagnostic dump; STOPPING, " *
                    "not retrying silently. Run diagnose_jacobian(m, vars) at this " *
                    "point, partitioned by region and margin block.")
                rep = gdp ? calculate_gdp(vals, params) : nothing
                return ScenarioResult(String(name), false, s, nsteps, nrejects, res,
                                      path, vals, rep, r0.residual)
            end
            if h < hmin
                log("  ⛔ step collapsed below hmin at t=$(round(s; sigdigits=6)). " *
                    "Every earlier step converged, so the obstruction is at that point " *
                    "on the branch, not in the stepping — the signature of a fold, " *
                    "which natural-parameter continuation provably cannot pass.")

                # ── Fallback: pseudo-arclength, when the scenario admits it ──────
                # A fold is a fold in the PARAMETERISATION; arclength's augmented
                # Jacobian stays nonsingular through the turn, so it can steer where
                # stepping in t provably cannot (§N.19, §N.19c).
                #
                # LIMITATION, stated rather than hidden: arclength_solve! traces ONE
                # VariableRef. A multi-shock scenario has no single variable to trace
                # — `t` is a bookkeeping scalar here, not a model variable — so the
                # fallback only applies when exactly one shock is live. This is the
                # same one-scalar-λ shape that made a two-shock scenario invisible in
                # the first place (§N.19e), so it is reported explicitly instead of
                # silently skipped.
                #
                # For a multi-shock scenario, the working recipe is the STAGED run
                # that §N.19e actually used: apply the discrete/switch shocks first
                # (delUnity 0→1 converged in 3 steps), then arclength the continuous
                # one from there (blabnat → 0.97, 22 steps, 0 rejections). Staging is
                # not an approximation — each stage ends on a genuine solution.
                if arclength && length(targets) == 1
                    lbl, vr, b, tg = targets[1]
                    log("  ↪ falling back to pseudo-arclength on $lbl " *
                        "(current $(round(b + s*(tg-b); sigdigits=8)) → target $(round(tg; sigdigits=8)))")
                    ra = arclength_solve!(m, vars, vr, Float64(tg);
                                          tol=tol, maxit=maxit, verbose=verbose,
                                          linsolve=linsolve)
                    # Re-express λ on the branch as the scenario fraction t.
                    t_arc = tg == b ? 1.0 : (ra.lam_reached - b) / (tg - b)
                    apath = [( (tg == b ? 1.0 : (l - b) / (tg - b)), it, rr)
                             for (l, it, rr) in ra.path]
                    append!(path, apath)
                    if ra.reached_target
                        log("  ✅ arclength reached the full shock — the obstruction was " *
                            "a fold in the parameterisation, not a limit of the model.")
                        vals = _current_values(m, vars)
                        rep = gdp ? calculate_gdp(vals, params) : nothing
                        return ScenarioResult(String(name), true, 1.0,
                                              nsteps + ra.nsteps, nrejects + ra.nrejects,
                                              ra.residual, path, vals, rep, r0.residual)
                    end
                    if ra.turned
                        log("  ↩ arclength TURNED at λ ≈ $(round(ra.lam_fold; sigdigits=10)) " *
                            "without reaching the target. Past a limit point there is no " *
                            "solution at that λ on this branch, so the target does not exist " *
                            "in this closure.\n" *
                            "     ⚠️ Before reporting that as a property of the model, check " *
                            "the SCENARIO is complete — a fold at λ* = 0.9908083 was reported " *
                            "exactly that way and was an artifact of a MISSING shock " *
                            "(TERM.CMF:113 delUnity=1). See test/verify_scenario_complete.jl.")
                    else
                        log("  ❌ arclength stopped without reaching the target AND without " *
                            "turning. A fold cannot defeat arclength, so this points at a " *
                            "bifurcation, a domain boundary, or a defect in src/arclength.jl — " *
                            "not at 'the fold is impassable'.")
                    end
                    vals = _current_values(m, vars)
                    rep = gdp ? calculate_gdp(vals, params) : nothing
                    return ScenarioResult(String(name), false, max(s, t_arc),
                                          nsteps + ra.nsteps, nrejects + ra.nrejects,
                                          ra.residual, path, vals, rep, r0.residual)
                elseif arclength
                    log("  (no arclength fallback: $(length(targets)) shocks are live and " *
                        "arclength_solve! traces a single variable. Run the scenario in " *
                        "STAGES — switch-type shocks first, then arclength the continuous " *
                        "one — as §N.19e did for delUnity + blabnat.)")
                end

                rep = gdp ? calculate_gdp(vals, params) : nothing
                return ScenarioResult(String(name), false, s, nsteps, nrejects, res,
                                      path, vals, rep, r0.residual)
            end
        end
    end

    log("reached t = 1 in $nsteps step(s), $nrejects rejection(s), ‖F‖∞ = $res")
    rep = gdp ? calculate_gdp(vals, params) : nothing
    if rep !== nothing
        chg_ok, chg_d, chg_worst = gdp_change_consistency(vals)
        log("GDP model-side check: worst |wgdpdiff| = $(round(chg_worst; sigdigits=6))" *
            (chg_ok ? " ✓" : "  ⚠ region $chg_d — income and expenditure sides disagree"))
    end
    return ScenarioResult(String(name), true, 1.0, nsteps, nrejects, res,
                          path, vals, rep, r0.residual)
end
