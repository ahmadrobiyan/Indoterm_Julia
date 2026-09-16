"""
Industrial / GVC policy analysis — hilirisasi (downstreaming) of the ore-to-metal
chain, run on the 25×6 INDOTERM levels model.

Design: four runs under one closure (`COALPRICE_SWAPS`, numeraire `:exrate`),
so the legs decompose.

  S1  ban          ore export demand −30%                    (restriction only)
  S2  smelter-Java BasMetals capital +20% in Java            (build, wrong site)
  S3  smelter-East BasMetals capital +20% in Sulawesi+Maluku (build, at the ore)
  S4  full         S1 + S3                                   (the actual policy)

The interaction S4 − S1 − S3 is the question hilirisasi's case rests on: is the
export restriction a COMPLEMENT to the downstream build, or does it merely
destroy upstream rent that the build would have captured anyway?

Scenario definitions live in `src/scenarios.jl` with their policy-mapping notes;
this file only runs them and tabulates. Read those notes before quoting anything
— in particular, the OthMining export base is 72% MalukuPapua / 18% BaliNusa, so
this is the COPPER concentrate chain, not the 2020 nickel ore ban.

Run:  julia --project=IndotermJulia IndotermJulia/analysis/gvc_policy.jl
"""

include(joinpath(@__DIR__, "..", "test", "pipeline_cache.jl"))
using JuMP, Printf, Serialization

agg6, params = cached_pipeline(6)

const NA = length(AGGCOM)
const NR = length(REG6)
const OTHMINING, BASMETALS, METALMACH, TRANSPEQUIP = 7, 14, 15, 16

CAP0    = parent(params["CAP"])
XEXPD0  = parent(params["XEXPD0"])
MAKE_I0 = parent(params["MAKE_I"])

# capital actually injected by each build, in benchmark units — the denominator
# for any siting comparison (see HILIRISASI_SMELTER_* notes)
const KINJ = Dict(
    "smelter-Java" => 0.20 * CAP0[BASMETALS, 2],
    "smelter-East" => 0.20 * (CAP0[BASMETALS, 4] + CAP0[BASMETALS, 6]),
)

SCENARIOS = [
    ("ban",          HILIRISASI_BAN),
    ("smelter-Java", HILIRISASI_SMELTER_JAVA),
    ("smelter-East", HILIRISASI_SMELTER_EAST),
    ("full",         HILIRISASI_FULL),
]

results = Dict{String,Any}()
for (tag, sc) in SCENARIOS
    println("\n" * "#"^72); println("# $tag"); println("#"^72)
    r = run_model!(agg6, params, sc; h0 = 0.05, gdp = true)
    r.solved || @warn "$tag stopped at t = $(r.t_reached) — NOT a full result; do not quote it"
    results[tag] = r
    @printf("%s: solved=%s  t=%.4f  steps=%d/%d  ‖F‖∞=%.3e\n",
            tag, r.solved, r.t_reached, r.nsteps, r.nrejects, r.residual)
end

serialize(joinpath(@__DIR__, "gvc_results.jls"),
          Dict(k => (solved = v.solved, t = v.t_reached, values = v.values) for (k, v) in results))

# ── reporting helpers ───────────────────────────────────────────────────────
macpct(r, name, reg) = begin
    nm = r.values["MainMacro"]; k = findfirst(==(name), MAINMACROS)
    100 * (nm[k, reg] - 1)
end
natpct(r, name) = macpct(r, name, NR + 1)
idxpct(r, nm, i...) = 100 * (r.values[nm][i...] - 1)

TAGS = first.(SCENARIOS)
solved_tags = [t for t in TAGS if results[t].solved]

println("\n\n" * "="^78)
println("NATIONAL MACRO RESULTS  (% change from benchmark; numeraire = exchange rate)")
println("="^78)
@printf("%-16s", "macro"); for t in TAGS; @printf("%14s", t); end; println()
for mname in ["RealGDP", "RealHou", "RealInv", "ExpVol", "ImpVolUsed", "AggEmploy",
              "realwage_io", "CPI", "GDPPI", "ExportPI"]
    @printf("%-16s", mname)
    for t in TAGS; @printf("%14.4f", natpct(results[t], mname)); end
    println()
end

println("\n" * "="^78)
println("REGIONAL REAL GDP  (% change)")
println("="^78)
@printf("%-14s", "region"); for t in TAGS; @printf("%14s", t); end; println()
for d in 1:NR
    @printf("%-14s", REG6[d])
    for t in TAGS; @printf("%14.4f", macpct(results[t], "RealGDP", d)); end
    println()
end

println("\n" * "="^78)
println("SECTORAL OUTPUT, NATIONAL  (xtot is a bmk=1 index; output-weighted over regions)")
println("="^78)
@printf("%-14s", "sector"); for t in TAGS; @printf("%14s", t); end; println()
for k in 1:NA
    w = MAKE_I0[k, :]; sw = sum(w)
    sw <= 0 && continue
    row = [100 * (sum(results[t].values["xtot"][k, d] * w[d] for d in 1:NR) / sw - 1) for t in TAGS]
    maximum(abs, row) < 0.02 && !(k in (OTHMINING, BASMETALS, METALMACH, TRANSPEQUIP)) && continue
    @printf("%-14s", AGGCOM[k]); for v in row; @printf("%14.4f", v); end; println()
end

println("\n" * "="^78)
println("THE CHAIN, BY REGION  (% change in xtot)")
println("="^78)
for k in (OTHMINING, BASMETALS, METALMACH)
    println("\n$(AGGCOM[k]):")
    @printf("  %-14s", "region"); for t in TAGS; @printf("%14s", t); end; println()
    for d in 1:NR
        MAKE_I0[k, d] <= 0 && continue
        @printf("  %-14s", REG6[d])
        for t in TAGS; @printf("%14.4f", idxpct(results[t], "xtot", k, d)); end
        println()
    end
end

println("\n" * "="^78)
println("EXPORT VOLUMES, xexpd  (% change; base share of national exports in ())")
println("="^78)
for k in (OTHMINING, BASMETALS, METALMACH, TRANSPEQUIP)
    b = sum(XEXPD0[k, :])
    @printf("%-14s (%4.1f%%)", AGGCOM[k], 100 * b / sum(XEXPD0))
    for t in TAGS
        x = results[t].values["xexpd"]
        @printf("%14.4f", 100 * (sum(x[k, d] for d in 1:NR) / b - 1))
    end
    println()
end

# ── the decomposition ───────────────────────────────────────────────────────
if all(t -> results[t].solved, ["ban", "smelter-East", "full"])
    println("\n" * "="^78)
    println("INTERACTION:  full − ban − smelter-East")
    println("="^78)
    println("A positive interaction means the restriction and the build are")
    println("COMPLEMENTS — the ban makes the smelter more valuable than it is alone.")
    println("A negative one means the ban destroys more upstream than it feeds downstream.\n")
    @printf("%-16s %12s %12s %12s %12s\n", "macro", "ban", "smelt-East", "full", "interaction")
    for mname in ["RealGDP", "RealHou", "ExpVol", "AggEmploy", "realwage_io"]
        a = natpct(results["ban"], mname); b = natpct(results["smelter-East"], mname)
        c = natpct(results["full"], mname)
        @printf("%-16s %12.4f %12.4f %12.4f %12.4f\n", mname, a, b, c, c - a - b)
    end
    println("\nBy region, RealGDP:")
    @printf("%-14s %12s %12s %12s %12s\n", "region", "ban", "smelt-East", "full", "interaction")
    for d in 1:NR
        a = macpct(results["ban"], "RealGDP", d); b = macpct(results["smelter-East"], "RealGDP", d)
        c = macpct(results["full"], "RealGDP", d)
        @printf("%-14s %12.4f %12.4f %12.4f %12.4f\n", REG6[d], a, b, c, c - a - b)
    end
end

# ── siting: same policy, different place, normalised per unit of capital ────
if all(t -> results[t].solved, ["smelter-Java", "smelter-East"])
    println("\n" * "="^78)
    println("SITING:  return per unit of smelter capital injected")
    println("="^78)
    @printf("capital injected — Java %.0f units, East %.0f units (ratio %.2f)\n\n",
            KINJ["smelter-Java"], KINJ["smelter-East"],
            KINJ["smelter-Java"] / KINJ["smelter-East"])
    gdp0 = sum(parent(params["GDPEXPSUM"]))
    @printf("%-22s %14s %14s\n", "", "smelter-Java", "smelter-East")
    for mname in ["RealGDP", "RealHou", "ExpVol", "AggEmploy"]
        @printf("%-22s %14.4f %14.4f\n", mname * " (%)",
                natpct(results["smelter-Java"], mname), natpct(results["smelter-East"], mname))
    end
    for mname in ["RealGDP", "RealHou"]
        j = natpct(results["smelter-Java"], mname) / 100 * gdp0 / KINJ["smelter-Java"]
        e = natpct(results["smelter-East"], mname) / 100 * gdp0 / KINJ["smelter-East"]
        @printf("%-22s %14.4f %14.4f\n", "Δ" * mname * " per unit K", j, e)
    end
end

println("\n" * "="^78)
println("done — scenario values serialised to analysis/gvc_results.jls")
println("="^78)
