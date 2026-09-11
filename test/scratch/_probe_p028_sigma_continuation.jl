"""
Is `P028 ×0.5` unreachable, or only unreachable along the shock path?

Every earlier attempt on the missing V8 point walked the SHOCK: fix σ_prim at half
its calibrated value, then continue `fpexp_d[coal]` from 0 to log(1.5). All three
low-σ cases (×0.5, ×0.6, ×0.75) folded on that path, at locations that are not
monotone in σ (t = 0.64, 0.14, 0.10) — which says the obstruction is a property of
the path through shock-space at low σ, not obviously of the endpoint.

This probe walks the ELASTICITY instead. Start at the ×1.0 shocked equilibrium
(solves in 3 steps, externally validated by V9), then step σ_prim down from 1.0·σ₀
towards 0.5·σ₀ with the full shock held fixed. Each step rebuilds the model with
`prepare_parameters!` re-calibrating ALPHA/GAMMA_FAC to the new σ (the same
perturb-then-recalibrate order V8 uses), re-verifies benchmark replication under
the new parameters, then Newton-solves the shocked system warm-started from the
previous σ's solution. Step in σ is adaptive with rollback, like `run_model!`'s
walk in t.

Both outcomes are informative:

  REACHED 0.5  — the 13th V8 point exists. Its headline row is printed in the
                 same 10 metrics as `analysis/output/v8_sensitivity.csv`, and
                 compared against the recorded min/max so the range can be
                 corrected from "lower bound" to complete.
  FOLDED at f  — two independent continuation paths fail to reach σ = 0.5·σ₀.
                 That upgrades "not attainable under this continuation" to
                 "no solution found along either path", and f_reached is the
                 lowest σ multiple at which the shocked equilibrium exists.

Cost: ~3 s build + one warm Newton solve per σ step. Budgeted well under an hour.

Run:  julia --project=. test/scratch/_probe_p028_sigma_continuation.jl
Log:  logs/p028_sigma_continuation_<date>.log
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, Printf, Dates

const KEY     = "P028"
const F_START = 1.0
const F_TARGET = 0.5
const DF0     = 0.10      # first step in the multiplier
const DFMIN   = 0.005     # below this the walk is declared obstructed
const DFMAX   = 0.10
const TOL     = 1e-8
const MAXIT   = 30
const COAL    = 5

# ── headline metrics, same construction as test/sensitivity.jl ────────────────────────────
const METRICS = [("Real HousCon","RealHou"), ("Real Invest","RealInv"), ("Export vol","ExpVol"),
                 ("Import vol","ImpVolUsed"), ("Real GNE","RealGNE"), ("Real GDP","RealGDP"),
                 ("Employment","AggEmploy"), ("CPI","CPI")]
const LABELS = vcat([l for (l,_) in METRICS], ["Coal output", "Coal employ"])

function headline(values, VTOT0, BMK0)
    out = Dict{String,Float64}()
    natmacro = values["NatMacro"]
    for (lbl, mac) in METRICS
        out[lbl] = 100 * (natmacro[findfirst(==(mac), MAINMACROS)] - 1.0)
    end
    xtot = values["xtot"]
    tw = sum(VTOT0[COAL, d] for d in 1:size(xtot, 2))
    out["Coal output"] = 100 * sum(VTOT0[COAL, d] * (xtot[COAL, d] - 1) for d in 1:size(xtot, 2)) / tw
    xlab = values["xlab_o"]; bl = BMK0["xlab_o"]
    out["Coal employ"] = 100 * (sum(xlab[COAL, d] for d in 1:size(xlab, 2)) /
                                sum(bl[COAL, d] for d in 1:size(bl, 2)) - 1)
    out
end

# Copy a solution Dict (as `solve_newton!`/`run_model!` return it) onto a freshly built
# model's start values. Fixed (exogenous) variables are left as the closure set them.
function warm_start!(vars, values)
    n = 0
    for (nm, v) in vars
        haskey(values, nm) || continue
        src = values[nm]
        for (k, vr) in enumerate(v isa AbstractArray ? vec(v) : [v])
            JuMP.is_fixed(vr) && continue
            JuMP.set_start_value(vr, src isa AbstractArray ? vec(src)[k] : src)
            n += 1
        end
    end
    n
end

function read_v8_ranges()
    path = joinpath(@__DIR__, "..", "..", "analysis", "output", "v8_sensitivity.csv")
    isfile(path) || return nothing
    rng = Dict{String,Tuple{Float64,Float64,Float64}}()
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        f = split(line, ',')
        rng[f[1]] = (parse(Float64, f[2]), parse(Float64, f[3]), parse(Float64, f[4]))
    end
    rng
end

say(msg) = (println(msg); flush(stdout))

function main()
    agg6, params = cached_pipeline(6)
    VTOT0 = parent(params["VTOT"])
    BMK0  = benchmark_levels(params)
    sig0  = vec(parent(agg6[KEY]))
    sc    = COALPRICE_REFERENCE
    (spec, target) = only(sc.shocks)

    say("\n" * "="^78)
    say("P028 elasticity continuation  σ_prim: $(F_START)·σ₀ → $(F_TARGET)·σ₀ at full coal shock")
    say("="^78)
    say("σ₀ (25 sectors): min $(round(minimum(sig0);digits=3))  max $(round(maximum(sig0);digits=3))  coal $(round(sig0[COAL];digits=3))")
    say("shock held fixed: $(spec) = $(round(target; sigdigits=8))   closure: $(length(sc.swaps)) swaps, numeraire :$(sc.numeraire)")
    say("step in f: df0=$DF0  dfmin=$DFMIN  tol=$TOL  maxit=$MAXIT")

    # ── anchor: the ×1.0 shocked equilibrium via the standard driver ───────────────────────
    say("\n── anchor: f = 1.0 via run_model! (the V9-validated point)")
    t0 = time()
    r1 = run_model!(agg6, params, sc; h0 = 0.25, tol = TOL, gdp = false, verbose = false)
    r1.solved || error("anchor failed: the ×1.0 point did not solve (t_reached = $(r1.t_reached)); nothing to continue from")
    say(@sprintf("   solved in %.0fs, %d steps, ‖F‖∞=%.2e", time()-t0, r1.nsteps, r1.residual))
    h1 = headline(r1.values, VTOT0, BMK0)
    say("   " * join([@sprintf("%s %+.3f", l, h1[l]) for l in LABELS], " | "))

    vals = r1.values
    f = F_START; df = DF0
    nacc = 0; nrej = 0
    path = Tuple{Float64,Int,Float64,Float64}[]   # (f, iters, residual, seconds)
    rows = Dict{Float64,Dict{String,Float64}}()
    rows[f] = h1

    while f > F_TARGET + 1e-12
        f_try = max(F_TARGET, f - df)
        say(@sprintf("\n── try f = %.4f  (df = %.4f)", f_try, df))
        t0 = time()

        agg2 = copy(agg6)
        agg2[KEY] = sig0 .* f_try
        p2 = prepare_parameters!(agg2)

        m, vars = build_model_full!(agg6, p2)
        initialize_model!(m, vars; bmk_levels = benchmark_levels(p2), numeraire = sc.numeraire)
        apply_swaps!(m, vars, collect(sc.swaps); verbose = false)

        # Benchmark must still replicate under the re-calibrated σ — this is what V8's
        # perturb-after-calibration bug broke, so it is checked at every σ, not assumed.
        r0 = solve_newton!(m, vars; maxit = 5, tol = TOL, verbose = false)
        if r0.residual > TOL
            say(@sprintf("   ✗ benchmark does NOT replicate at f=%.4f (‖F‖∞=%.2e) — calibration defect, aborting", f_try, r0.residual))
            break
        end

        _, _, vr = IndotermJulia._shock_ref(vars, spec)
        JuMP.fix(vr, Float64(target); force = true)
        nws = warm_start!(vars, vals)
        r = solve_newton!(m, vars; maxit = MAXIT, tol = TOL, verbose = false)
        dt = time() - t0

        if r.residual <= TOL
            nacc += 1; f = f_try; vals = r.values
            push!(path, (f, r.iters, r.residual, dt))
            h = headline(vals, VTOT0, BMK0); rows[f] = h
            say(@sprintf("   ✅ f=%.4f  it=%d  ‖F‖∞=%.2e  bench=%.1e  warm=%d  %.0fs", f, r.iters, r.residual, r0.residual, nws, dt))
            say("      " * join([@sprintf("%s %+.3f", l, h[l]) for l in LABELS], " | "))
            r.iters <= 4 && (df = min(DFMAX, 1.5df))
        else
            nrej += 1
            say(@sprintf("   ✗ f=%.4f failed (%s, ‖F‖∞=%.2e, it=%d) %.0fs — df → %.4f", f_try, r.status, r.residual, r.iters, dt, df/2))
            df /= 2
            if df < DFMIN
                say(@sprintf("   ⛔ step in f collapsed below %.3f at f = %.4f. The shocked equilibrium was tracked continuously down to here and no further.", DFMIN, f))
                break
            end
        end
    end

    say("\n" * "─"^78)
    say(@sprintf("path: %d accepted, %d rejected", nacc, nrej))
    say("   f        it   ‖F‖∞      s")
    for (ff, it, rr, dt) in path
        say(@sprintf("   %.4f   %2d   %.1e   %3.0f", ff, it, rr, dt))
    end

    reached = f <= F_TARGET + 1e-12
    say("\n── RESULT: " * (reached ?
        "REACHED f = 0.5 — the 13th V8 point exists; the shock-path fold was a path property" :
        @sprintf("OBSTRUCTED at f = %.4f — two independent continuation paths fail below this σ", f)))

    if reached
        h = rows[F_TARGET]; rng = read_v8_ranges()
        say("\n   P028 ×0.5 headline row, and what it does to the recorded V8 range:")
        say(@sprintf("   %-14s %9s %9s %9s %9s  %s", "metric", "x0.5", "old min", "old max", "centre", "range"))
        for l in LABELS
            v = h[l]
            if rng !== nothing && haskey(rng, l)
                c, lo, hi = rng[l]
                widen = v < lo ? "WIDENS below" : v > hi ? "WIDENS above" : "inside"
                sgn = sign(v) == sign(c) ? "" : "  ❗ SIGN FLIP vs centre"
                say(@sprintf("   %-14s %+9.4f %+9.4f %+9.4f %+9.4f  %s%s", l, v, lo, hi, c, widen, sgn))
            else
                say(@sprintf("   %-14s %+9.4f", l, v))
            end
        end
    else
        say("\n   headline at the lowest σ reached (f = $(round(f;digits=4))):")
        h = rows[f]
        say("   " * join([@sprintf("%s %+.3f", l, h[l]) for l in LABELS], " | "))
    end
    say("\nDone  $(Dates.now())")
end

main()
