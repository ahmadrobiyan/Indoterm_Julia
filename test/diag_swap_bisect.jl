"""
Which of TERM.CMF's six swaps makes the base year near-singular?

`test/verify_swapped_closure.jl` established two facts that have to be held
together:

  * the swap is CORRECT — after all six, the benchmark still solves to
    ‖F‖∞ = 1.4260876923799515e-9, bit-identical to the base closure;
  * yet a continuation walk accepted ZERO steps. It failed at h=0.5 down to
    h=0.0078 (a shock of ~1e-6), with the leftover residual falling only like
    h^1.5 and NON-monotonically in h (6.56e-5 at h=0.25, worse than 3.15e-5 at
    h=0.5).

A fold cannot explain that. A fold is a property of a point on the branch, and
s=1.0 is the benchmark, which solves exactly. Failing at an arbitrarily small
shock means the JACOBIAN is near-singular at the base year under this closure:
the solution exists, but Newton cannot find the direction to move in.

So the question is no longer "does the fold survive a faithful closure" but
"which swap destroys the conditioning". Each swap is applied ALONE here, not
incrementally, because an incremental ladder confounds the culprit with its
position in the ordering — with singletons, the k-th result depends on the k-th
swap and nothing else.

The probe is a 1e-5 shock. Under the base closure that converges in 2 Newton
iterations; anything that needs more than a handful is the suspect. `maxit` is
deliberately modest: a case that needs 20 iterations for 1e-5 has already
answered the question, and failures dominate wall-clock.

Prior suspicion (to be confirmed or killed, not assumed): `xcap = faccum` and
`finv1 = finv4` free 150 variables each against equations that live in the
DYNAMIC block (build_dynamics!.jl:64-65). In a single-period static solve those
accumulation relations may pin `xcap` only weakly, which is exactly what a
near-zero pivot looks like.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_swap_bisect.jl
Env:  SHOCK (default 0.99999), MAXIT (default 20)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP

SHOCK = parse(Float64, get(ENV, "SHOCK", "0.99999"))
MAXIT = parse(Int,     get(ENV, "MAXIT", "20"))

agg6, params = cached_pipeline(6)
bmk = benchmark_levels(params)

CASES = vcat(
    [("none (control)", Any[])],
    [("$(a) = $(b)", Any[(a, b)]) for (a, b) in TERM_CMF_SWAPS],
    [("ALL SIX", collect(TERM_CMF_SWAPS))],
)

RESULTS = Tuple{String,Symbol,Int,Float64,Float64}[]   # label, status, iters, residual, seconds

for (label, subset) in CASES
    println("\n", "="^72)
    println("══ $label ══")
    println("="^72)

    # Fresh model per case. `solve_newton!` pins orphan variables and DELETES
    # dead rows, so a model that has solved once cannot be re-closed by
    # re-initialising — the pins survive and silently change the system.
    (m, vars) = build_model_full!(agg6, params)
    initialize_model!(m, vars; bmk_levels=bmk)
    isempty(subset) || apply_swaps!(m, vars, subset)

    r0 = solve_newton!(m, vars; maxit=5, verbose=false)
    println("benchmark: status=$(r0.status)  ‖F‖∞=$(r0.residual)")
    if r0.residual > 1e-8
        println("  ⛔ this swap BREAKS THE BENCHMARK — the shock probe below would be")
        println("     meaningless, so it is skipped.")
        push!(RESULTS, (label, :benchmark_broken, 0, r0.residual, 0.0))
        continue
    end

    sv = vars["blabnat"]
    svr = sv isa AbstractArray ? first(sv) : sv
    JuMP.fix(svr, SHOCK; force=true)

    el = @elapsed r = solve_newton!(m, vars; maxit=MAXIT, verbose=false)
    println("shock blabnat=$SHOCK: status=$(r.status)  iters=$(r.iters)  " *
            "‖F‖∞=$(r.residual)  ($(round(el, digits=1))s)")
    push!(RESULTS, (label, r.status, r.iters, r.residual, el))
end

println("\n", "="^72)
println("── which swap costs the conditioning (probe: blabnat = $SHOCK) ──")
println("="^72)
for (label, st, it, res, el) in RESULTS
    mark = st == :converged ? "✅" : "❌"
    println("  $mark $(rpad(label, 26)) $(rpad(st, 18)) iters=$(rpad(it, 4)) " *
            "‖F‖∞=$(rpad(round(res; sigdigits=4), 12)) $(round(el, digits=1))s")
end

ctrl = findfirst(r -> r[1] == "none (control)", RESULTS)
if ctrl !== nothing && RESULTS[ctrl][2] != :converged
    println("\n⚠️  THE CONTROL FAILED. The base closure cannot take a $SHOCK shock either,")
    println("    so nothing above is attributable to a swap. Fix the probe (smaller SHOCK")
    println("    or larger MAXIT) before reading the singletons.")
else
    bad = [r[1] for r in RESULTS if r[2] != :converged && r[1] != "ALL SIX"]
    println("\n", isempty(bad) ?
        "➡ Every swap is individually benign. The damage is an INTERACTION between\n" *
        "  swaps — re-run incrementally to find the pair, since singletons cannot\n" *
        "  show it." :
        "➡ Individually damaging: $(join(bad, ", ")).\n" *
        "  That closure choice is where the near-singularity enters. Next: is it a\n" *
        "  genuine economic under-determination (the freed variable has no strong\n" *
        "  equation in a single-period solve) or a scaling problem in the freed block?")
end
println("\nDone.")
