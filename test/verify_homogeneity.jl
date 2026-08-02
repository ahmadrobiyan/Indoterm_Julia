"""
Price homogeneity — an independent check on the translation (PLAN task #11).

A correctly specified CGE model is homogeneous of degree zero in prices: scale
every nominal price by λ and no real quantity moves. In the levels formulation
that is directly testable, because the numeraire is an explicit variable —
`initialize_model!` applies TERM.CMF's swap `phi = NatMacro("GDPPI")`, freeing
the exchange rate and pinning the national GDP price index at 1.0. Re-pin it at
λ instead and re-solve.

Expected:
  * every quantity variable (`x…`) unchanged — ratio 1.0;
  * every price variable (`p…`) scaled by exactly λ;
  * benchmark-zero variables — see the next block, this is where run 1 went wrong.

──────────────────────────────────────────────────────────────────────────────
WHAT RUN 1 GOT WRONG, AND HOW THIS VERSION AVOIDS REPEATING IT
──────────────────────────────────────────────────────────────────────────────
Run 1 classified every variable with a zero benchmark as "additive ⇒ must not
move", and duly reported 6,858 violations. That classification is wrong. Two
quite different kinds of variable have a zero benchmark:

  * genuine ADDITIVE SHIFTERS — `delTAXint`, `delPTX`, `finv1`, `atradmar`,
    technology and tax-rate shocks. These must not move.
  * NOMINAL CHANGE VARIABLES — `delVGDPEXP`, `delPGDPEXP`, `delGDPINC`. These
    measure *the change in a rupiah value flow*. When every price rises 10%,
    nominal GDP rises 10%, so these variables MUST move, by exactly
    (λ−1)·(base value flow). Reporting them as failures is reporting the model
    for being right.

Rather than hand-classify by name — which is guessing, and guessing is what
manufactured the false verdict — this version discriminates by MEASUREMENT.
Solve at two scale factors, λ₁ and λ₂, and use the fact that

    a nominal change variable satisfies  v(λ) = K·(λ − 1)   for some constant K,
    so                                   v(λ₂)/v(λ₁) = (λ₂−1)/(λ₁−1)  exactly,

while a genuine additive shifter is zero at both. The ratio is parameter-free:
it needs no knowledge of which base flow the variable belongs to, and it is a
STRONGER test than "does not move", because a variable that scales with the
wrong degree fails it. Anything that is neither zero nor exactly proportional to
(λ−1) is a real violation.

──────────────────────────────────────────────────────────────────────────────
STRUCTURALLY INDETERMINATE PRICES — measured, not assumed
──────────────────────────────────────────────────────────────────────────────
`E_pdelivrd!` (`build_equations.jl:803`) reads

    log(pdelivrd[c,s,r,d]) = BASSHR[c,s,r,d]·log(pbasic) + Σ_m MARSHR[…]·log(psuppmar_p)

Homogeneity of this equation rests on the weights summing to one. Where the
underlying trade flow is zero, every weight is zero, the equation degenerates to
`log(pdelivrd) = 0`, and the variable is pinned at 1.0 whatever λ is. That is not
a homogeneity failure — it is a price attached to a zero quantity, which has no
economic content and multiplies nothing anywhere in the model.

The equivalent applies on the quantity side through `E_gret!` →`E_ggro!` →
`E_xinvitot!` → `E_xinvi!`: `gret = log(pcap) − log(pinvitot)` inherits any
failure of `pcap` to scale, and `xinvitot`/`xinvi` inherit it from there.

This version does not take that on faith. For every failing cell it prints the
weight sum and the base flow, and excludes a cell from the verdict ONLY when the
weight sum is measurably zero. A failing cell with a healthy weight sum stays a
failure. The exclusion count is reported, so a rule quietly excusing thousands
of cells would be visible rather than silent.

──────────────────────────────────────────────────────────────────────────────
HOW THE λ POINTS ARE REACHED — and why not by continuation
──────────────────────────────────────────────────────────────────────────────
Homogeneity is an exact symmetry, so the solution at numeraire λ is known in
closed form before any solving: prices ×λ, quantities unchanged. This version
jumps straight there and lets one Newton solve confirm it, rather than walking a
homotopy path to an answer it already knows.

That is a ~40× speedup and a stronger test. It replaced `continuation_solve!`
after a measured failure: the 1.00→1.05 leg took one step and 131s, while the
identical-size 1.05→1.10 leg ran 115 minutes without landing — because
continuation accepts a step only when `residual <= tol`, and the achievable
residual floor rises with λ. Once the floor exceeds `tol` every step is rejected
at every size and the run reports UNTESTED after hours.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_homogeneity.jl
Env:  LAMBDA (default 1.10), LAMBDA1 (default 1.05), NTOL (default 1e-8)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, Printf

const LAMBDA  = parse(Float64, get(ENV, "LAMBDA",  "1.10"))
const LAMBDA1 = parse(Float64, get(ENV, "LAMBDA1", "1.05"))
@assert 1.0 < LAMBDA1 < LAMBDA "need 1 < LAMBDA1 < LAMBDA for the two-λ discriminator"

# Newton tolerance. This is a knob the verdict is DELIBERATELY read against: run 2
# left two residual leads both sitting at ~4e-4 relative and both on the smallest
# cells in the model (pcap where CAP≈0.03, delXGDPEXP against GDP flows of ~1e6).
# That is the signature of solver tolerance, not of a wrong price degree. The
# discriminator is generic and needs no per-variable scale: tighten the tolerance
# and see whether the error follows it down. If it does, it was never economics.
const NTOL   = parse(Float64, get(ENV, "NTOL", "1e-8"))

const QTOL   = 1e-6    # relative tolerance on an expected ratio
const ZTOL   = 1e-6    # |v| below this counts as "did not move" for a zero-benchmark var
const PROPTOL = 1e-4   # tolerance on the (λ₂−1)/(λ₁−1) proportionality ratio
const WTOL   = 1e-10   # a weight sum below this is structurally zero

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

r0 = solve_newton!(m, vars; maxit=12, tol=NTOL, verbose=false)
println("\n══ benchmark ══  status=$(r0.status)  ‖F‖∞=$(r0.residual)  (NTOL=$NTOL)")
@assert r0.residual <= max(NTOL, 1e-8) "benchmark regressed — fix that before reading anything below"
# An NTOL below the solver's achievable floor is not a tighter test, it is a
# broken one. Measured floor at 25×6 is ~6e-10, so NTOL=1e-11 fails this way.
# Say so here rather than 40 minutes later.
#
# This guard is NECESSARY BUT NOT SUFFICIENT, and the gap cost 2.5 hours once: it
# measures the floor at the BENCHMARK, and the floor rises with λ. The λ loop
# below probes the floor again at each stop, which is where it actually bites.
if r0.status != :converged
    println("\n⚠️  NTOL=$NTOL is at or below the achievable residual floor " *
            "(benchmark stalled at $(r0.residual) with status=$(r0.status)).")
    println("    continuation_solve! will reject every step. Raise NTOL and rerun.")
    exit(2)
end

_val(vr) = JuMP.is_fixed(vr) ? JuMP.fix_value(vr) : JuMP.start_value(vr)

# Snapshot keyed by variable so each entry is compared in its own right rather
# than through an aggregate norm.
function snapshot()
    d = Dict{JuMP.VariableRef,Float64}()
    for v in JuMP.all_variables(m)
        x = _val(v)
        d[v] = x === nothing ? NaN : x
    end
    return d
end

base = snapshot()

k_gdppi = findfirst(==("GDPPI"), MAINMACROS)
@assert haskey(vars, "NatMacro") && k_gdppi !== nothing "no NatMacro(GDPPI) — numeraire swap not applied"
num = vars["NatMacro"][k_gdppi]

# ── reaching the λ stops: ANALYTIC WARM START, not continuation ─────────────
# Homogeneity is an exact symmetry, so if it holds the solution at numeraire λ is
# known in closed form — prices ×λ, quantities unchanged. There is nothing to
# search for. Re-pin the numeraire at λ, jump the price variables straight to
# their predicted values, and let one Newton solve confirm or refute the
# prediction from a start point that should already BE the answer.
#
# Why this replaced `continuation_solve!`. Measured: the 1.00→1.05 leg took one
# step and 131s; the identical-size 1.05→1.10 leg ran 115 minutes without
# landing. Not stiffness — this path is as smooth as paths get. Continuation
# accepts a step only when `residual <= tol` (`continuation.jl:113`), and the
# achievable floor RISES with λ. Once the floor at λ exceeds `tol`, every step is
# rejected AT EVERY SIZE, `h` collapses to `hmin` through ~14 halvings of up to
# `maxit` Newton iterations each, and the run reports "UNTESTED" after hours. The
# NTOL guard at the top of this file checks the floor only at the BENCHMARK —
# one λ short of where it bites. The floor probe below closes that gap.
#
# Why the warm start is SAFE even though it rests on a name-based price guess:
# it is only a START POINT. Newton converges to the same solution whatever it
# starts from, so a misclassified variable costs iterations and never
# correctness. This is the opposite direction from run 1's error, which used a
# name-based guess to EXCUSE observed behaviour. Guessing a start point is fine;
# guessing a verdict is not.
#
# The iteration count is itself evidence. A near-exact prediction converges in
# 1-2 iterations. Many iterations means the predicted point was NOT the solution,
# i.e. the model is not scaling the way homogeneity requires — visible before any
# of the classification machinery below runs.

freevars = [vr for (_, v) in vars for vr in (v isa AbstractArray ? vec(v) : [v])
            if !JuMP.is_fixed(vr)]

function warm_start!(lam::Float64)
    # Reset every free variable to its benchmark value, then scale the prices.
    # Resetting matters: after the λ₁ solve the model sits at the λ₁ point, and
    # each λ must be reached independently or the two stops are not independent
    # observations and the proportionality discriminator is measuring its own
    # previous answer.
    for vr in freevars
        b = get(base, vr, NaN)
        isfinite(b) && JuMP.set_start_value(vr, b)
    end
    for (nm, v) in vars
        startswith(nm, "p") || continue
        for vr in (v isa AbstractArray ? v : (v,))
            JuMP.is_fixed(vr) && continue
            b = get(base, vr, NaN)
            isfinite(b) && JuMP.set_start_value(vr, b * lam)
        end
    end
end

snaps = Dict{Float64,Dict{JuMP.VariableRef,Float64}}()
for lam in (LAMBDA1, LAMBDA)
    println("\n══ scaling the numeraire NatMacro(GDPPI) → $lam ══")
    JuMP.fix(num, lam; force=true)
    warm_start!(lam)

    # Floor probe. Solve once at a deliberately LOOSE tolerance to discover what
    # this λ can actually achieve, before demanding NTOL of it. Without this the
    # failure mode is silent and expensive: a tolerance below the floor does not
    # produce a tighter test, it produces a rejected one, and the run says
    # "UNTESTED" — which reads like a model problem and is not one.
    el = @elapsed pr = solve_newton!(m, vars; maxit=30, tol=1e-6, verbose=false)
    @printf("floor probe: iters=%d  ‖F‖∞=%.4g  (%.1fs)\n", pr.iters, pr.residual, el)
    if pr.iters <= 2
        println("  ↳ converged in $(pr.iters) iteration(s) from the predicted point — " *
                "the analytic prediction was already the solution, which is the " *
                "homogeneity result stated directly.")
    end

    if pr.residual > NTOL
        el2 = @elapsed nr = solve_newton!(m, vars; maxit=30, tol=NTOL, verbose=false)
        @printf("refine to NTOL=%g: iters=%d  ‖F‖∞=%.4g  (%.1fs)\n",
                NTOL, nr.iters, nr.residual, el2)
        if nr.residual > NTOL
            println("\n⚠️  λ=$lam cannot reach NTOL=$NTOL — floor here is ≈$(nr.residual), " *
                    "above the benchmark floor of $(r0.residual).")
            println("    This is an INSTRUMENT limit, not a homogeneity failure. The verdict " *
                    "below is still readable, but read it against $(nr.residual), not NTOL:")
            println("    any per-cell error at or under that size is solver noise.")
        end
    end

    snaps[lam] = snapshot()
end
JuMP.fix(num, 1.0; force=true)   # leave the model on its benchmark numeraire
S1, S2 = snaps[LAMBDA1], snaps[LAMBDA]

# ── the structural-indeterminacy oracle ─────────────────────────────────────
# Returns the total weight on the RHS of the equation determining this cell, or
# NaN when we have no oracle for that family. NaN means "no excuse available",
# so an unexplained failure cannot be silently dropped.
BASSHR = parent(params["BASSHR"]); MARSHR = parent(params["MARSHR"])
DELIVRD = parent(params["DELIVRD"]); CAP = parent(params["CAP"])
nm_marg = size(MARSHR, 3)

function weight_sum(nm::String, idx)
    if nm == "pdelivrd"
        c, s, r, d = idx
        return BASSHR[c,s,r,d] + sum(MARSHR[c,s,mh,r,d] for mh in 1:nm_marg)
    elseif nm in ("xinvitot", "gret", "ggro", "pcap", "pinvitot")
        i, d = idx
        return CAP[i,d]                      # no capital ⇒ no rate-of-return signal
    elseif nm == "xinvi"
        _, i, d = idx
        return CAP[i,d]
    end
    return NaN
end

# ── classify and score ──────────────────────────────────────────────────────
struct Rec
    name::String
    nm::String
    idx::Tuple
    ratio::Float64      # v(λ₂)/base, or v(λ₂) itself for zero-benchmark vars
    err::Float64
    fixed::Bool
    w::Float64          # weight sum from the oracle (NaN if none)
    r1::Float64         # v(λ₁)/base, for reading how the failure develops
end

recs = Dict{String,Vector{Rec}}()
const PROP_EXPECT = (LAMBDA - 1) / (LAMBDA1 - 1)

for (nm, v) in vars
    arr = v isa AbstractArray
    for ci in (arr ? CartesianIndices(v) : (nothing,))
        vr = arr ? v[ci] : v
        b = get(base, vr, NaN)
        isfinite(b) || continue
        v1 = get(S1, vr, NaN); v2 = get(S2, vr, NaN)
        (isfinite(v1) && isfinite(v2)) || continue
        idx = arr ? Tuple(ci) : ()

        cls, ratio, err = if abs(b) < 1e-8
            # Zero benchmark: additive shifter or nominal change variable? Let the
            # two-λ proportionality answer, instead of assuming.
            #
            # ORDER MATTERS, and run 2 got it backwards. Testing invariance first
            # gates the better test behind an ABSOLUTE magnitude cutoff, which is
            # meaningless for a model whose variables span rupiah-billions down to
            # tax-rate dust. It split one class in half: delTAXint[4,2,4,5] came in
            # at λ₁=9.927e-07, λ₂=1.985e-06 — a ratio of exactly 2.000, a perfect
            # pass — and was reported as a failure because λ₁ landed 0.7% under
            # ZTOL. Proportionality is scale-free, so ask it FIRST and let the
            # magnitude cutoff catch only what proportionality cannot score.
            p = abs(v1) > 0 ? v2 / v1 : NaN
            if isfinite(p) && abs(p / PROP_EXPECT - 1) <= PROPTOL
                ("nominal-change", p, abs(p / PROP_EXPECT - 1))
            elseif max(abs(v1), abs(v2)) < ZTOL
                ("additive-invariant", v2, abs(v2))
            elseif !isfinite(p)
                ("zero-base-unexplained", v2, abs(v2))
            else
                ("nominal-change", p, abs(p / PROP_EXPECT - 1))
            end
        elseif startswith(nm, "p")
            ("price", v2 / b, abs(v2 / b - LAMBDA))
        elseif startswith(nm, "x")
            ("quantity", v2 / b, abs(v2 / b - 1.0))
        else
            ("other", v2 / b, NaN)
        end

        tol = cls == "nominal-change" ? PROPTOL : QTOL
        err > tol || continue           # only keep the offenders; passing cells just count
        push!(get!(recs, cls, Rec[]),
              Rec(JuMP.name(vr), nm, idx, ratio, err, JuMP.is_fixed(vr),
                  weight_sum(nm, idx), isfinite(b) && abs(b) > 1e-8 ? v1 / b : v1))
    end
end

# Totals per class, counted the same way the offenders were.
totals = Dict{String,Int}()
for (nm, v) in vars
    arr = v isa AbstractArray
    for ci in (arr ? CartesianIndices(v) : (nothing,))
        vr = arr ? v[ci] : v
        b = get(base, vr, NaN); isfinite(b) || continue
        v1 = get(S1, vr, NaN); v2 = get(S2, vr, NaN)
        (isfinite(v1) && isfinite(v2)) || continue
        cls = if abs(b) < 1e-8
            p = abs(v1) > 0 ? v2 / v1 : NaN
            if isfinite(p) && abs(p / PROP_EXPECT - 1) <= PROPTOL; "nominal-change"
            elseif max(abs(v1), abs(v2)) < ZTOL; "additive-invariant"
            elseif !isfinite(p); "zero-base-unexplained"
            else; "nominal-change" end
        elseif startswith(nm, "p"); "price"
        elseif startswith(nm, "x"); "quantity"
        else; "other" end
        totals[cls] = get(totals, cls, 0) + 1
    end
end

println("\n" * "="^78)
@printf("HOMOGENEITY BY CLASS   λ₁=%.4g  λ₂=%.4g   (proportionality expects %.6g)\n",
        LAMBDA1, LAMBDA, PROP_EXPECT)
println("="^78)

failed = String[]
for cls in ("quantity", "price", "nominal-change", "additive-invariant",
            "zero-base-unexplained", "other")
    n = get(totals, cls, 0); n == 0 && continue
    bad = get(recs, cls, Rec[])
    if cls == "other"
        @printf("\n%-22s %6d variables   (no expected ratio — inspected below)\n", uppercase(cls), n)
        continue
    end
    # The verdict is scored on FREE variables only. A closure-fixed price is
    # SUPPOSED not to scale — `pfimp`/`fpexp` are foreign-currency prices, and
    # holding them while the numeraire moves is what makes this a test rather
    # than a tautology. And a structurally indeterminate cell (weight sum 0)
    # has no equation that could scale it.
    badfree = filter(r -> !r.fixed, bad)
    excused = filter(r -> isfinite(r.w) && r.w <= WTOL, badfree)
    real    = filter(r -> !(isfinite(r.w) && r.w <= WTOL), badfree)
    @printf("\n%-22s %6d variables   %d off  →  %d closure-fixed, %d zero-weight, %d REAL\n",
            uppercase(cls), n, length(bad),
            length(bad) - length(badfree), length(excused), length(real))
    isempty(real) ? println("  ✅ every free, determinate cell within tolerance") :
                    push!(failed, cls)
    for r in first(sort(bad; by = x -> -x.err), 8)
        tag = r.fixed ? "FIXED by closure" :
              (isfinite(r.w) && r.w <= WTOL ? "zero-weight ⇒ indeterminate" : "REAL")
        @printf("    %-30s λ₂=%-13.9g λ₁=%-13.9g err=%-10.4g w=%-10.4g %s\n",
                r.name, r.ratio, r.r1, r.err, r.w, tag)
    end
end

# ── evidence for every exclusion, by family ─────────────────────────────────
println("\n" * "="^78)
println("EVIDENCE FOR THE EXCLUSIONS — the numbers behind 'zero-weight'")
println("="^78)
allbad = reduce(vcat, values(recs); init=Rec[])
fams = unique(r.nm for r in allbad if !r.fixed)
for fam in sort(collect(fams))
    rs = [r for r in allbad if r.nm == fam && !r.fixed]
    ws = [r.w for r in rs if isfinite(r.w)]
    @printf("\n%-14s %4d failing free cell(s)", fam, length(rs))
    if isempty(ws)
        println("   — no structural oracle for this family, so NOT excused")
        # For the GDP change variables there IS a natural denominator: the base
        # value flow they are a change IN. Absolute rupiah tell you nothing about
        # whether a miss matters; relative to the flow, they usually tell you it
        # is noise. Print it so the reader is not left comparing a Rp-billion
        # residual against a dimensionless tolerance (measurement trap #2).
        if fam in ("delXGDPEXP", "delVGDPEXP", "delPGDPEXP")
            GES = parent(params["GDPEXPSUM"])
            rel = Float64[]
            for r in rs
                d, g = r.idx
                (d <= size(GES,1) && g <= size(GES,2) && abs(GES[d,g]) > 1e-8) || continue
                push!(rel, abs(get(S2, vars[fam][r.idx...], NaN)) / abs(GES[d,g]))
            end
            if !isempty(rel)
                @printf("     relative to base GDPEXPSUM flow: max=%.3e  median=%.3e  (n=%d)\n",
                        maximum(rel), sort(rel)[length(rel) ÷ 2 + 1], length(rel))
            end
        end
        continue
    end
    @printf("\n   weight sum over failing cells: min=%.4g  max=%.4g  (%d of %d ≤ %g)\n",
            minimum(ws), maximum(ws), count(<=(WTOL), ws), length(ws), WTOL)
    if fam == "pdelivrd"
        fl = [DELIVRD[r.idx...] for r in rs]
        @printf("   base DELIVRD flow (Rp bn):     min=%.4g  max=%.4g\n",
                minimum(fl), maximum(fl))
        println("   reading: weight sum 0 ⇒ E_pdelivrd degenerates to log(p)=0, pinning p at 1.")
        println("            The flow it prices is zero, so nothing in the model multiplies it.")
    elseif fam in ("xinvi", "xinvitot")
        idd = unique(fam == "xinvi" ? (r.idx[2], r.idx[3]) : (r.idx[1], r.idx[2]) for r in rs)
        println("   affected (industry, region) cells: $idd")
        for (i, d) in idd
            @printf("     CAP[%d,%d]=%.6g   pcap ratio λ₂=%.9g   pinvitot ratio λ₂=%.9g\n",
                    i, d, CAP[i,d],
                    S2[vars["pcap"][i,d]] / base[vars["pcap"][i,d]],
                    S2[vars["pinvitot"][i,d]] / base[vars["pinvitot"][i,d]])
        end
        println("   reading: xinvi = INVEST·xinvitot, and xinvitot is driven by")
        println("            gret = log(pcap) − log(pinvitot). Any failure of pcap to scale")
        println("            propagates here; it does not originate here.")
    end
end

# ── the "other" bucket, where an unclassified convention hides ──────────────
if get(totals, "other", 0) > 0
    rs = Float64[]
    for (nm, v) in vars
        arr = v isa AbstractArray
        for ci in (arr ? CartesianIndices(v) : (nothing,))
            vr = arr ? v[ci] : v
            b = get(base, vr, NaN)
            (isfinite(b) && abs(b) >= 1e-8) || continue
            (startswith(nm, "p") || startswith(nm, "x")) && continue
            x = get(S2, vr, NaN); isfinite(x) && push!(rs, x / b)
        end
    end
    if !isempty(rs)
        sort!(rs)
        @printf("\n  other ratios: min=%.9g  median=%.9g  max=%.9g\n",
                rs[1], rs[length(rs) ÷ 2 + 1], rs[end])
        near1 = count(r -> abs(r - 1.0) <= QTOL, rs)
        nearL = count(r -> abs(r - LAMBDA) <= QTOL, rs)
        @printf("  of %d: %d ≈ 1.0, %d ≈ λ, %d neither\n",
                length(rs), near1, nearL, length(rs) - near1 - nearL)
    end
end

println("\n" * "="^78)
println(isempty(failed) ?
    "✅ PRICE HOMOGENEITY HOLDS.\n" *
    "   Quantities invariant; prices scale by λ; nominal change variables scale\n" *
    "   exactly with (λ−1). Every remaining offender is either fixed by the\n" *
    "   closure or attached to a zero weight, with the numbers shown above." :
    "❌ HOMOGENEITY VIOLATED in: $(join(failed, ", ")) — see the cells marked REAL.")
println("\nDone.")
