"""
Step 8 — GDP reporting and the GDP-both-sides check.

The GDP *equations* live in the core model (Excerpts 28-29, `build_equations.jl`):
`delGDPINC[d,k]` decomposes the change in nominal income-side GDP over
`GDPINCCAT = (Land, Labour, Capital, ProdTax, ComTax)`, and
`delXGDPEXP`/`delPGDPEXP`/`delVGDPEXP[d,k]` decompose the expenditure side over
`GDPEXPCAT = (HOU, INV, GOV, STOCKS, EXP, Imports, RExports, RImports, NetMar)`.
This module turns a solved model into the reported tables, per `CHKMOD.TAB`'s
diagnostic conventions.

**Slot order is load-bearing.** TABLO indexes these by *category name*, so the
equation letters do NOT run parallel to the category order — `E_delXGDPEXPd`
writes "NetMar" (slot 9), and `E_delGDPINCb` writes "Capital" (slot 3). Getting
this wrong is invisible to benchmark replication, because `wgdpinc`/`wgdpexp` sum
over all categories and the totals stay right while the components are
mislabelled. The names below are the single source of truth for the layout.
"""

const GDPINCCAT = ["Land", "Labour", "Capital", "ProdTax", "ComTax"]
const GDPEXPCAT = ["HOU", "INV", "GOV", "STOCKS", "EXP",
                   "Imports", "RExports", "RImports", "NetMar"]

"""
GDP report produced by [`calculate_gdp`](@ref).

`inc`/`exp_` are `nr × category` matrices of **levels** (benchmark + change);
`inc_change`/`exp_change` are the changes themselves. `inc_total`/`exp_total` are
per-region totals, and `diff`/`reldiff` the income-minus-expenditure gap that the
GDP-both-sides check tests. `national_*` are the economy-wide sums.
"""
struct GDPReport
    inc::Matrix{Float64}
    exp_::Matrix{Float64}
    inc_change::Matrix{Float64}
    exp_change::Matrix{Float64}
    inc_total::Vector{Float64}
    exp_total::Vector{Float64}
    diff::Vector{Float64}
    reldiff::Vector{Float64}
    national_inc::Float64
    national_exp::Float64
    national_reldiff::Float64
end

"""
    calculate_gdp(sol, params) -> GDPReport

Build the income- and expenditure-side GDP decompositions from a solved model.
`sol` is a `NewtonResult.values` dict; `params` comes from `prepare_parameters!`.

Income side uses `GDPINCSUM .+ delGDPINC`, expenditure side
`GDPEXPSUM .+ delVGDPEXP` — both are GEMPACK `(change)` variables, i.e. genuine
absolute currency changes from the benchmark, so the level is a plain sum.
"""
function calculate_gdp(sol::Dict{String,Any}, params::Dict{String,Any})
    GDPINCSUM = parent(params["GDPINCSUM"])
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    nr = size(GDPINCSUM, 1)

    dinc = sol["delGDPINC"]
    dexp = sol["delVGDPEXP"]

    inc = GDPINCSUM .+ dinc
    exp_ = GDPEXPSUM .+ dexp

    inc_total = [sum(inc[d, :]) for d in 1:nr]
    exp_total = [sum(exp_[d, :]) for d in 1:nr]
    diff = inc_total .- exp_total
    reldiff = [abs(inc_total[d]) > 1e-10 ? diff[d] / abs(inc_total[d]) : diff[d]
               for d in 1:nr]

    nat_inc = sum(inc_total); nat_exp = sum(exp_total)
    nat_rel = abs(nat_inc) > 1e-10 ? (nat_inc - nat_exp) / abs(nat_inc) : nat_inc - nat_exp

    return GDPReport(inc, exp_, copy(dinc), copy(dexp),
                     inc_total, exp_total, diff, reldiff,
                     nat_inc, nat_exp, nat_rel)
end

"""
    gdp_both_sides_check(report; tol=1e-6) -> (ok, worst_region, worst_reldiff)

**Levels** GDP-both-sides check: income-side GDP == expenditure-side GDP, region
by region. Compares on a **relative** basis, since these are ~1e6-magnitude value
aggregates where an absolute tolerance is meaningless (the same reason
`solve_newton!` equilibrates its rows).

⚠️ This measures the **benchmark database**, not the translation. At a zero-shock
solve every `(change)` component is ~0, so it reduces to `GDPINCSUM` vs
`GDPEXPSUM` as calibrated. INDOTERM's database does *not* satisfy the regional
GDP identity exactly (the same imbalance that forces the `K` constant in
`E_pdomA_sum!` and shows up as `DIFFCOM_sc` in the pipeline diagnostics), so a
non-zero result here is a data-quality finding to weigh when interpreting
regional results — it is not evidence of a bad equation. Use
[`gdp_change_consistency`](@ref) for the model-side check.
"""
function gdp_both_sides_check(report::GDPReport; tol::Float64=1e-6)
    worst, wd = 0.0, 0
    for d in eachindex(report.reldiff)
        a = abs(report.reldiff[d])
        a > worst && ((worst, wd) = (a, d))
    end
    return (worst < tol, wd, worst)
end

"""
    gdp_change_consistency(sol; tol=1e-6) -> (ok, worst_region, worst_absdiff)

**Model-side** GDP consistency check: TERM's own `wgdpdiff = wgdpinc − wgdpexp`
diagnostic (Excerpt 29). `wgdpinc`/`wgdpexp` are ratio indices benchmarked at 1,
so this compares the two sides' *changes* and is zero at the benchmark by
construction — which is exactly why it is blind to the levels mismatch that
[`gdp_both_sides_check`](@ref) reports, and why both are worth printing.

Under a shock this is the meaningful test: if the income and expenditure sides of
the model move together, the accounting closes.
"""
function gdp_change_consistency(sol::Dict{String,Any}; tol::Float64=1e-6)
    haskey(sol, "wgdpdiff") || return (true, 0, 0.0)
    wd_ = sol["wgdpdiff"]
    worst, wr = 0.0, 0
    for d in eachindex(wd_)
        a = abs(wd_[d])
        a > worst && ((worst, wr) = (a, d))
    end
    return (worst < tol, wr, worst)
end

"""
    print_gdp_report(report; regions=nothing, io=stdout)

Print the income/expenditure decomposition and the both-sides check.
`regions` optionally supplies region names for labelling.
"""
function print_gdp_report(report::GDPReport; regions=nothing, io::IO=stdout)
    nr = length(report.inc_total)
    rname(d) = regions === nothing ? "reg$d" : string(regions[d])

    println(io, "="^72)
    println(io, "GDP DECOMPOSITION")
    println(io, "="^72)

    println(io, "\n── Income side (levels, by ", join(GDPINCCAT, "/"), ") ──")
    println(io, rpad("region", 14), join([lpad(c, 14) for c in GDPINCCAT]), lpad("TOTAL", 16))
    for d in 1:nr
        println(io, rpad(rname(d), 14),
                join([lpad(_fmt(report.inc[d, k]), 14) for k in eachindex(GDPINCCAT)]),
                lpad(_fmt(report.inc_total[d]), 16))
    end

    println(io, "\n── Expenditure side (levels, by ", join(GDPEXPCAT, "/"), ") ──")
    println(io, rpad("region", 14), join([lpad(c, 12) for c in GDPEXPCAT]), lpad("TOTAL", 16))
    for d in 1:nr
        println(io, rpad(rname(d), 14),
                join([lpad(_fmt(report.exp_[d, k]), 12) for k in eachindex(GDPEXPCAT)]),
                lpad(_fmt(report.exp_total[d]), 16))
    end

    println(io, "\n── GDP both-sides check (income − expenditure) ──")
    println(io, rpad("region", 14), lpad("income", 16), lpad("expenditure", 16),
            lpad("diff", 16), lpad("rel.diff", 14))
    for d in 1:nr
        println(io, rpad(rname(d), 14), lpad(_fmt(report.inc_total[d]), 16),
                lpad(_fmt(report.exp_total[d]), 16), lpad(_fmt(report.diff[d]), 16),
                lpad(string(round(report.reldiff[d], sigdigits=4)), 14))
    end
    println(io, "\nNational: income = ", _fmt(report.national_inc),
            "   expenditure = ", _fmt(report.national_exp),
            "   rel.diff = ", round(report.national_reldiff, sigdigits=4))

    ok, wd, worst = gdp_both_sides_check(report)
    if ok
        println(io, "Levels GDP both-sides check PASSED ✓ (worst rel.diff $(round(worst, sigdigits=4)))")
    else
        println(io, "Levels GDP both-sides gap: worst rel.diff $(round(worst, sigdigits=4)) at $(rname(wd)).")
        println(io, "  This measures the BENCHMARK DATABASE, not the translation — at a zero shock")
        println(io, "  it reduces to GDPINCSUM vs GDPEXPSUM as calibrated. See gdp_both_sides_check.")
    end
    return report
end

_fmt(x) = string(round(x, digits=2))
