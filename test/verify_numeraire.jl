"""
Numeraire invariance (VV_PLAN.md V2).

Real results must not depend on which price is chosen as the numeraire. This is
a *different* bug class than `verify_homogeneity.jl`: homogeneity scales the
CHOSEN numeraire and checks every price follows it exactly (an algebraic
identity that must hold by construction). This test changes WHICH price is
pinned — GDP price index vs. CPI — and checks that the model's real economy
does not care. That catches a variable accidentally deflated by one specific
price index instead of by the model's general price level: a bug homogeneity
cannot see, because it never varies which price is numeraire.

**Method** (VV_PLAN.md §V2). Run `TERM_CMF_REFERENCE` twice: once under the
shipped `:gdppi` swap (`phi` free, `NatMacro("GDPPI")` pinned at 1.0), once
under a `:cpi` swap (`phi` free, `NatMacro("CPI")` pinned at 1.0) added to
`initialize_model!`/`scenarios.jl` for this test. Compare every variable.

**What "CPI" means here.** `NatMacro("CPI")` is a live, solved variable —
`build_macros!.jl`'s `src_of["CPI"] = d -> pfin[1,d]` constrains
`MainMacro[k_cpi, d] == pfin[1,d]` (the regional household price index) for
every region, which is then weighted into `NatMacro[k_cpi]` via WMAIN/WNAT.
That is distinct from `build_macros!.jl`'s *other* `"CPI"` entry (in `wsrc`,
using `PUR_CS[u_hou,d]`) — that dict only builds the fixed WMAIN weight, it is
not a constraint on a variable, and pinning against it is not an option.

**Pass.**
  * every quantity variable (`x…`, following the same prefix convention as
    `verify_homogeneity.jl`) agrees between the two runs to solver tolerance —
    the real economy must be identical;
  * every free price/nominal variable (`p…`, `w…`) differs between the two
    runs by *one* common scalar (the ratio of the two runs' `phi`), not by
    several different scalars — i.e. `v_cpi / v_gdppi` is the same constant
    across all of them. A variable whose ratio deviates from that constant is
    a numeraire-dependent bug: something that should be measured relative to
    the general price level is instead pinned to one specific index.

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_numeraire.jl
Env:  NTOL (default 1e-8)
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, Printf

const NTOL = parse(Float64, get(ENV, "NTOL", "1e-8"))
const QTOL = 1e-6     # relative tolerance for quantity-variable agreement
const PROPTOL = 1e-4  # tolerance on the common-scalar ratio for price/nominal vars

agg6, params = cached_pipeline(6)

sc = TERM_CMF_REFERENCE
sc_cpi = Scenario(name = sc.name * "  [:cpi numeraire]", source = sc.source,
                   swaps = sc.swaps, shocks = sc.shocks, pct_shocks = sc.pct_shocks,
                   numeraire = :cpi, notes = sc.notes)

println("="^72)
println("V2 — numeraire invariance: TERM_CMF_REFERENCE under :gdppi vs :cpi")
println("="^72)

println("\n── run 1: numeraire = :gdppi (the shipped swap) ──")
el1 = @elapsed r1 = run_model!(agg6, params, sc; tol = NTOL, gdp = false)
@printf("solved=%s  t_reached=%.6g  steps/rejects=%d/%d  ‖F‖∞=%.3e  (%.1f min)\n",
        r1.solved, r1.t_reached, r1.nsteps, r1.nrejects, r1.residual, el1/60)
@assert r1.solved "run 1 (:gdppi) did not reach the target — nothing to compare"
@assert r1.residual <= max(NTOL, 1e-8) "run 1 (:gdppi) residual $(r1.residual) is not converged"

println("\n── run 2: numeraire = :cpi ──")
el2 = @elapsed r2 = run_model!(agg6, params, sc_cpi; tol = NTOL, gdp = false)
@printf("solved=%s  t_reached=%.6g  steps/rejects=%d/%d  ‖F‖∞=%.3e  (%.1f min)\n",
        r2.solved, r2.t_reached, r2.nsteps, r2.nrejects, r2.residual, el2/60)
@assert r2.solved "run 2 (:cpi) did not reach the target — nothing to compare"
@assert r2.residual <= max(NTOL, 1e-8) "run 2 (:cpi) residual $(r2.residual) is not converged"

V1, V2 = r1.values, r2.values
common = intersect(Set(keys(V1)), Set(keys(V2)))
@assert !isempty(common) "the two runs' r.values share no variable names — nothing to compare"

struct Rec
    name::String
    nm::String
    idx::Tuple
    v1::Float64
    v2::Float64
    ratio::Float64
    err::Float64
end

qbad = Rec[]; nbad = Rec[]; nom_ratios = Float64[]
nq = 0; nn = 0

for nm in common
    a, b = V1[nm], V2[nm]
    arr = a isa AbstractArray
    (arr == (b isa AbstractArray)) || continue
    for ci in (arr ? CartesianIndices(a) : (nothing,))
        v1 = arr ? a[ci] : a
        v2 = arr ? b[ci] : b
        (v1 isa Real && v2 isa Real && isfinite(v1) && isfinite(v2)) || continue
        idx = arr ? Tuple(ci) : ()
        if startswith(nm, "x")
            global nq += 1
            err = abs(v1) > 1e-8 ? abs(v2/v1 - 1.0) : abs(v2 - v1)
            err <= QTOL || push!(qbad, Rec(string(nm, idx), nm, idx, v1, v2, abs(v1)>1e-8 ? v2/v1 : NaN, err))
        elseif startswith(nm, "p") || startswith(nm, "w")
            abs(v1) > 1e-8 || continue    # zero-benchmark price/nominal: no ratio to read
            global nn += 1
            push!(nom_ratios, v2/v1)
        end
    end
end

println("\n" * "="^72)
@printf("QUANTITY VARIABLES (x…): %d compared, %d off by more than %.1e\n", nq, length(qbad), QTOL)
if isempty(qbad)
    println("  ✅ every quantity variable agrees to solver tolerance")
else
    for r in first(sort(qbad; by = x -> -x.err), 10)
        @printf("    %-30s gdppi=%-13.6g cpi=%-13.6g ratio=%-10.6g err=%.4g\n",
                r.name, r.v1, r.v2, r.ratio, r.err)
    end
end

med_ratio = isempty(nom_ratios) ? NaN : sort(nom_ratios)[length(nom_ratios) ÷ 2 + 1]
println("\n" * "="^72)
@printf("PRICE/NOMINAL VARIABLES (p…, w…): %d compared against median ratio %.9g\n", nn, med_ratio)

# Discriminate by MEASUREMENT, not by name, mirroring verify_homogeneity.jl's own
# design principle. A `p…`-named variable can legitimately fall in either of two
# classes here, and guessing which from the name is exactly the mistake that
# test warns against:
#   * NOMINAL price/value  — scales by the run's common ratio (v2/v1 ≈ med_ratio).
#   * REAL relative price  — some equations build log(domestic price) − log(phi)
#     (`E_pfexp!`, `build_equations.jl:711`: `pfexp = log(ppur) − log(phi)`).
#     If ppur and phi both scale by the same common factor between the two
#     numeraire runs (which is exactly what "nominal" means here), that
#     DIFFERENCE is invariant — ratio ≈ 1 — even though the variable is named
#     `p...`. That is not a numeraire-dependence bug, it is the model correctly
#     expressing a foreign-currency-relative price. Only a ratio that matches
#     NEITHER expectation is a real offender.
nombad = Rec[]; nom_real = 0; nom_nominal = 0
by_family = Dict{String,Int}()
if !isnothing(med_ratio) && isfinite(med_ratio)
    for nm in common
        a, b = V1[nm], V2[nm]
        (startswith(nm, "p") || startswith(nm, "w")) || continue
        arr = a isa AbstractArray
        (arr == (b isa AbstractArray)) || continue
        for ci in (arr ? CartesianIndices(a) : (nothing,))
            v1 = arr ? a[ci] : a
            v2 = arr ? b[ci] : b
            (v1 isa Real && v2 isa Real && isfinite(v1) && isfinite(v2) && abs(v1) > 1e-8) || continue
            idx = arr ? Tuple(ci) : ()
            ratio = v2 / v1
            err_nom  = abs(ratio / med_ratio - 1.0)
            err_real = abs(ratio - 1.0)
            if err_nom <= PROPTOL
                global nom_nominal += 1
            elseif err_real <= PROPTOL
                global nom_real += 1
            else
                err = min(err_nom, err_real)
                by_family[nm] = get(by_family, nm, 0) + 1
                push!(nombad, Rec(string(nm, idx), nm, idx, v1, v2, ratio, err))
            end
        end
    end
end
@printf("  %d nominal (ratio≈%.6g), %d real-relative (ratio≈1), %d unexplained\n",
        nom_nominal, med_ratio, nom_real, length(nombad))
if isempty(nombad)
    println("  ✅ every free price/nominal variable is explained: either scales by the " *
            "common nominal ratio or is invariant as a real relative price")
else
    println("  by family:")
    for (fam, n) in sort(collect(by_family); by = x -> -x[2])
        @printf("    %-16s %d cell(s)\n", fam, n)
    end
    @printf("  worst %d unexplained cell(s):\n", min(10, length(nombad)))
    for r in first(sort(nombad; by = x -> -x.err), 10)
        @printf("    %-30s gdppi=%-13.6g cpi=%-13.6g ratio=%-10.6g err=%.4g\n",
                r.name, r.v1, r.v2, r.ratio, r.err)
    end
end

println("\n" * "="^72)
ok = isempty(qbad) && isempty(nombad)
println(ok ?
    "✅ NUMERAIRE INVARIANCE HOLDS.\n" *
    "   Real quantities are identical under :gdppi and :cpi; every free price/nominal\n" *
    "   variable differs by exactly one common scalar." :
    "❌ NUMERAIRE INVARIANCE VIOLATED — see the REAL offenders listed above.")
println("\nDone.")
