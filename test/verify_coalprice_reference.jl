"""
V9 — the port against a published result (task #49).

**Why this test is different from every other one in `test/`.** Everything else
here checks the model against itself: the benchmark reproduces its own base year,
prices scale when the numeraire scales, Walras closes. Those are *verification* —
"did we build the model right". They cannot catch an error that was faithfully
translated from a misreading of the source, because the model would then be
consistently wrong and every internal check would pass.

This one is *validation*: `origin/draftreport.pdf` reports a simulation the
authors ran in GEMPACK, and `origin/coalprice.CMF` ships the exact closure and
shock that produced it. Reproducing their numbers is the only evidence available
that the translated model is the same model.

**What agreement is reasonable to demand, and why not more.** The report solves
with `method = euler; steps = 8` — an 8-step Euler linearisation. This port
solves the exact levels system by Newton. For a 50% shock those are genuinely
different numerical objects, and Euler at 8 steps carries visible truncation
error on responses this large. Demanding 3 significant figures would be a
category error; the honest tests are

  1. SIGNS, on every column. A sign flip is a translation defect, never a
     truncation artifact.
  2. The DIAGNOSTIC PAIR — Real GNE up while Real GDP is slightly down. That
     opposite-sign pair IS the Dutch-disease story the report tells, and no
     amount of Euler error manufactures or hides it.
  3. MAGNITUDES to within a tolerance that Euler truncation can plausibly
     explain, reported as a table rather than asserted, so a reader can see
     which columns agree well and which merely agree in direction.

A column that fails (3) but passes (1) is reported as a discrepancy to explain,
NOT as a pass and NOT as a failure — that distinction is the whole point of
running this.

**Two things NOT compared here, deliberately.** The regional rows of Table 2
(aggregating the report's 18 regions to our 6 needs weights we only have
approximately — see the reference file), and Table 3's Imports/Exports columns
(never recovered from the PDF; reconstructing them from prose was explicitly
ruled out).

Reference: `test/reference/draftreport_coalprice.md`
Run: julia --project=IndotermJulia IndotermJulia/test/verify_coalprice_reference.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using Printf

# ── Table 2, national row (draftreport.pdf, row 19) ────────────────────────
# (label, MAINMACROS name, reported % change)
const TABLE2_NATIONAL = [
    ("Real HousCon", "RealHou",    1.10),
    ("Real Invest",  "RealInv",    0.90),
    ("Export vol",   "ExpVol",    -3.41),
    ("Import vol",   "ImpVolUsed", 2.00),
    ("Real GNE",     "RealGNE",    0.90),
    ("Real GDP",     "RealGDP",   -0.09),
    ("Employment",   "AggEmploy", -0.27),
    ("CPI",          "CPI",        1.54),
]

# Table 3, the Coal row (sector 5). Output/Employment/Price only — the
# Imports/Exports columns of that table were never recovered.
const TABLE3_COAL = (output = 2.21, employment = 14.99, price = 51.90)

const COAL = 5          # aggregation_data.jl:42
const MAGTOL = 0.5      # percentage POINTS; see the header on what this can and
                        # cannot mean. Not a pass/fail gate on its own.

agg6, params = cached_pipeline(6)

res = run_model!(agg6, params, COALPRICE_REFERENCE; h0=0.25, tol=1e-8, gdp=false)

println("\n" * "="^78)
println("V9 — coalprice.CMF vs draftreport.pdf Table 2 (national)")
println("="^78)

if !res.solved
    println("""
    THE RUN DID NOT REACH t = 1 (got to t = $(round(res.t_reached; sigdigits=6)),
    ‖F‖∞ = $(round(res.residual; sigdigits=4))). No comparison below is meaningful and
    none is printed — a partial homotopy is not a scaled-down version of the shock,
    it is a different experiment. Fix the solve first.""")
    exit(1)
end

@printf("solved: %d steps, %d rejects, ‖F‖∞ = %.3g (benchmark %.3g)\n\n",
        res.nsteps, res.nrejects, res.residual, res.bench_residual)

natmacro = res.values["NatMacro"]
_pctchg(name) = 100 * (natmacro[findfirst(==(name), MAINMACROS)] - 1.0)

@printf("%-14s %10s %10s %10s   %s\n", "column", "report", "port", "diff", "")
println("-"^62)

signs_ok = String[]; signs_bad = String[]; mag_off = String[]
port = Dict{String,Float64}()

for (label, macname, reported) in TABLE2_NATIONAL
    got = _pctchg(macname)
    port[label] = got
    d = got - reported
    # A reported value near zero has no meaningful sign to match — Real GDP at
    # -0.09 is the obvious case. Treating |x| < 0.05 as signless keeps a rounding
    # difference on a near-zero column from being called a sign error.
    signless = abs(reported) < 0.05
    same_sign = signless || sign(got) == sign(reported)
    same_sign ? push!(signs_ok, label) : push!(signs_bad, label)
    abs(d) > MAGTOL && push!(mag_off, label)
    flag = !same_sign ? "❌ SIGN" : abs(d) > MAGTOL ? "⚠️  magnitude" : "✅"
    @printf("%-14s %10.2f %10.2f %+10.2f   %s\n", label, reported, got, d, flag)
end

println("\n" * "="^78)
println("THE DIAGNOSTIC PAIR")
println("="^78)
gne, gdp = port["Real GNE"], port["Real GDP"]
println("""
The report's central claim is Dutch disease: the terms of trade improve, so
Indonesia CONSUMES more (Real GNE +0.90) while PRODUCING slightly less
(Real GDP -0.09). Those opposite signs are the single most diagnostic pair in
Table 2 — a model that got the mechanism wrong would move them together.""")
@printf("\n  port: Real GNE %+.2f%%,  Real GDP %+.2f%%\n", gne, gdp)
pair_ok = gne > 0 && gdp < gne
println(pair_ok ?
    "  ✅ GNE rises and GDP grows by less — the Dutch-disease wedge is reproduced." :
    "  ❌ the wedge is NOT reproduced. This is a mechanism failure, not a tolerance issue.")

println("\n" * "="^78)
println("Table 3 — the Coal row")
println("="^78)
# The report's sectoral rows are NATIONAL, so each is an aggregate over our 6
# regions — and the right aggregator differs by variable type. Quantities sum:
# national coal output is the sum of regional coal output, so its % change is the
# ratio of the summed level to the summed benchmark. A PRICE does not sum; it
# must be weighted, and weighted by its OWN denominator (output), because that is
# the aggregation error already made once in this project and diagnosed — see
# memory note "INDOTERM shock fold": averaging a rate by the wrong weight
# collapsed a 1186× amplification to 11×.
#
# `xtot` and `ptot` are bmk=1 INDICES, not flows — `build_model!.jl:141` calls
# xtot "a pure real-activity index", and E_ptot! carries the benchmark value
# VTOT as a coefficient rather than as ptot's level. So neither appears in
# `benchmark_levels`, and the first version of this script reported both rows as
# "unavailable" because it looked for them there. An index aggregates as a
# VTOT-weighted mean of (level − 1); only a genuine flow like `xlab_o` sums.
bmk = benchmark_levels(params)
VTOT = parent(params["VTOT"])          # benchmark output value, na × nr

"""VTOT-weighted national % change of a bmk=1 regional index."""
function _agg_index_pct(key)
    haskey(res.values, key) || return nothing
    v = res.values[key]
    tw = sum(VTOT[COAL, d] for d in 1:size(v, 2))
    tw <= 0 && return nothing
    return 100 * sum(VTOT[COAL, d] * (v[COAL, d] - 1) for d in 1:size(v, 2)) / tw
end

function _agg_qty_pct(key)
    haskey(res.values, key) && haskey(bmk, key) || return nothing
    v, b = res.values[key], bmk[key]
    nowv = sum(v[COAL, d] for d in 1:size(v, 2))
    base = sum(b[COAL, d] for d in 1:size(b, 2))
    base <= 0 ? nothing : 100 * (nowv / base - 1)
end

t3 = [("Output",     TABLE3_COAL.output,     _agg_index_pct("xtot")),
      ("Employment", TABLE3_COAL.employment, _agg_qty_pct("xlab_o")),
      ("Price",      TABLE3_COAL.price,      _agg_index_pct("ptot"))]

@printf("%-12s %10s %10s %10s\n", "", "report", "port", "diff")
println("-"^46)
for (lbl, reported, got) in t3
    if got === nothing
        @printf("%-12s %10.2f %10s   (variable or benchmark level unavailable)\n",
                lbl, reported, "—")
    else
        @printf("%-12s %10.2f %10.2f %+10.2f\n", lbl, reported, got, got - reported)
    end
end
println("""
The Price row is the anchor: the shock is +50% on the coal export price and the
report prints +51.90. That is the one Table 3 number whose magnitude is pinned by
the shock itself rather than by model behaviour, so it is the cheapest available
check that the shock was applied in the right units — recall that `pct(50)` in
place of `logpct(50)` would have applied ≈ +348% and still solved.""")

println("\n" * "="^78)
println("VERDICT")
println("="^78)
if !isempty(signs_bad)
    println("""
    ❌ FAIL — sign mismatch on: $(join(signs_bad, ", ")).
    A sign flip cannot be blamed on the report's 8-step Euler approximation. This
    is a translation defect and must be found before any magnitude is discussed.""")
elseif !pair_ok
    println("""
    ❌ FAIL — every column has the right sign, but the GNE/GDP wedge is not
    reproduced. The mechanism, not the arithmetic, is wrong.""")
elseif isempty(mag_off)
    println("""
    ✅ PASS. Every column matches the report in sign and to within $(MAGTOL) percentage
    points, and the Dutch-disease wedge is reproduced. Given that the reference is an
    8-step Euler solve and this is an exact levels solve, this is as close as the two
    methods should be expected to come.""")
else
    println("""
    ✅ SIGNS AND MECHANISM PASS; magnitudes differ by more than $(MAGTOL) pp on:
    $(join(mag_off, ", ")).

    This is the outcome that needs a judgement rather than a label. Euler truncation
    at 8 steps over a 50% shock is a real and sufficient explanation for differences
    of this size, and the columns that differ most should be the ones with the largest
    responses if that is what is happening. Check that before concluding anything: if
    a SMALL-response column is the one that is off, truncation is not the explanation
    and something in the translation is.""")
end
println("\nDone.")
