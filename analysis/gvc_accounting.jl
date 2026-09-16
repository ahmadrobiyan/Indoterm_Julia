"""
GVC accounting on the INDOTERM 25-sector national core, as a companion to
`analysis/linkages.jl` (which it re-uses wholesale — same use table, same
industry-technology conversion, same aggregation).

Three measures, all computable from a single-country table:

- **VS (vertical specialisation, Hummels–Ishii–Yi)** — imported input content of
  a sector's exports, direct plus indirect: `u' A_m (I − A_d)⁻¹ e / X`. This is
  backward GVC participation. A high VS says the sector's exports are largely
  assembled from foreign inputs.
- **DVA** — the domestic complement, 1 − VS, split into the sector's OWN value
  added and value added sourced from other domestic sectors. The second piece is
  what an industrial policy aimed at "deepening the domestic chain" is trying to
  raise.
- **Upstreamness / downstreamness (Antràs–Chor)** — output-weighted distance to
  final demand, and distance from primary factors. Note upstreamness is
  *algebraically* the Ghosh row sum, i.e. the same number `linkages.jl` reports
  as the forward linkage; it is printed here under its GVC name so the two
  literatures are visibly the same object and nobody double-counts them as
  independent evidence.

**Single-country limit.** Forward GVC participation — domestic value added
embodied in OTHER countries' exports — is not identified without a global input
output table (ADB MRIO / OECD ICIO). Nothing here should be reported as forward
participation. Where a number below looks like it is measuring re-export by
partners, it is not.

Run:  julia IndotermJulia/analysis/gvc_accounting.jl
"""

using LinearAlgebra, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "prepare_sets.jl"))
include(joinpath(ROOT, "src", "aggregation_data.jl"))

const NC = length(COM); const NI = NC
const CIDX = Dict(c => i for (i, c) in enumerate(COM))
const SIDX = Dict("dom" => 1, "imp" => 2)
const MIDX = Dict(m => i for (i, m) in enumerate(MAR))
const NA = length(AGGCOM); const MAP = SEC_MAP_185_to_25

function read_sections(path, wanted::Set{String})
    out = Dict{String,Vector{Tuple{Vector{String},Float64}}}(); cur = ""; keep = false
    for line in eachline(path)
        if startswith(line, "__header__:")
            cur = split(line[12:end], ',')[1]; keep = cur in wanted
            keep && (out[cur] = Tuple{Vector{String},Float64}[]); continue
        end
        keep || continue
        f = split(line, ','); push!(out[cur], (String.(f[1:end-1]), parse(Float64, f[end])))
    end
    out
end

@info "reading national_data.csv …"
S = read_sections(joinpath(ROOT, "data", "national_data.csv"),
                  Set(["1BAS", "1MAR", "MAKE", "4BAS", "4MAR", "1CAP", "1LAB", "1LND", "1PTX", "1OCT"]))

V1BAS = zeros(NC, 2, NI); for (k, v) in S["1BAS"]; V1BAS[CIDX[k[1]], SIDX[k[2]], CIDX[k[3]]] = v; end
V1MAR = zeros(NC, 2, NI, length(MAR))
for (k, v) in S["1MAR"]; V1MAR[CIDX[k[1]], SIDX[k[2]], CIDX[k[3]], MIDX[k[4]]] = v; end
MAKE = zeros(NC, NI); for (k, v) in S["MAKE"]; MAKE[CIDX[k[1]], CIDX[k[2]]] = v; end
V4BAS = zeros(NC); for (k, v) in S["4BAS"]; V4BAS[CIDX[k[1]]] = v; end
V4MAR = zeros(NC, length(MAR)); for (k, v) in S["4MAR"]; V4MAR[CIDX[k[1]], MIDX[k[2]]] = v; end
VA185 = zeros(NI)
for h in ("1CAP", "1LND", "1PTX", "1OCT"); for (k, v) in S[h]; VA185[CIDX[k[1]]] += v; end; end
for (k, v) in S["1LAB"]; VA185[CIDX[k[1]]] += v; end

Udom = zeros(NC, NI); Uimp = zeros(NC, NI)
for j in 1:NI, c in 1:NC
    Udom[c, j] += V1BAS[c, 1, j]; Uimp[c, j] += V1BAS[c, 2, j]
end
for j in 1:NI, m in 1:length(MAR)
    Udom[CIDX[MAR[m]], j] += sum(@view V1MAR[:, :, j, m])
end

agg2(X) = (Y = zeros(NA, NA); for j in axes(X, 2), i in axes(X, 1); Y[MAP[i], MAP[j]] += X[i, j]; end; Y)
agg1(x) = (y = zeros(NA); for i in eachindex(x); y[MAP[i]] += x[i]; end; y)

UD = agg2(Udom); UM = agg2(Uimp); MK = agg2(MAKE); VA = agg1(VA185)
# exports: basic value of the commodity, plus the margin services carrying it
EXPc = agg1(V4BAS)
for m in 1:length(MAR); EXPc[MAP[CIDX[MAR[m]]]] += sum(@view V4MAR[:, m]); end

gind = vec(sum(MK, dims = 1)); qcom = vec(sum(MK, dims = 2))
D = zeros(NA, NA); for c in 1:NA, i in 1:NA; qcom[c] > 0 && (D[i, c] = MK[c, i] / qcom[c]); end
Bd = zeros(NA, NA); Bm = zeros(NA, NA)
for j in 1:NA, c in 1:NA
    gind[j] > 0 && (Bd[c, j] = UD[c, j] / gind[j]; Bm[c, j] = UM[c, j] / gind[j])
end
Ad = D * Bd; Am = D * Bm
Ld = inv(I - Ad)
EXPi = D * EXPc                       # exports allocated to producing industries

# value-added coefficient per unit of industry output
vcoef = [gind[j] > 0 ? VA[j] / gind[j] : 0.0 for j in 1:NA]
mcoef = vec(sum(Am, dims = 1))        # imported input per unit of output

# ── VS: imported content of exports, direct + indirect ──────────────────────
VS = zeros(NA); DVA_own = zeros(NA); DVA_oth = zeros(NA)
for j in 1:NA
    EXPi[j] <= 0 && continue
    e = zeros(NA); e[j] = EXPi[j]
    x = Ld * e                                  # gross output pulled by these exports
    VS[j]      = sum(mcoef .* x) / EXPi[j]
    DVA_own[j] = vcoef[j] * x[j] / EXPi[j]
    DVA_oth[j] = (sum(vcoef .* x) - vcoef[j] * x[j]) / EXPi[j]
end

# ── Antràs–Chor upstreamness / downstreamness ──────────────────────────────
Zd = Ad .* reshape(gind, 1, :)
Bg = zeros(NA, NA); for i in 1:NA, j in 1:NA; gind[i] > 0 && (Bg[i, j] = Zd[i, j] / gind[i]); end
UP   = vec(sum(inv(I - Bg), dims = 2))
DOWN = vec(sum(Ld, dims = 2))

println("\n" * "="^92)
println("GVC ACCOUNTING — 25 sectors, INDOTERM national core")
println("="^92)
@printf("%-14s %9s %8s %8s %8s %8s %7s %7s %7s\n",
        "sector", "exports", "exp%out", "VS%", "DVAown%", "DVAoth%", "upstr", "downstr", "VA/out")
for k in sortperm(VS, rev = true)
    EXPi[k] <= 0 && continue
    @printf("%-14s %9.0f %8.1f %8.1f %8.1f %8.1f %7.3f %7.3f %7.3f\n",
            AGGCOM[k], EXPi[k], 100 * EXPi[k] / gind[k], 100VS[k], 100DVA_own[k],
            100DVA_oth[k], UP[k], DOWN[k], vcoef[k])
end

tot = sum(EXPi)
@printf("\nnational: exports %.0f, aggregate VS = %.2f%%, DVA = %.2f%%\n",
        tot, 100 * sum(VS .* EXPi) / tot, 100 * sum((DVA_own .+ DVA_oth) .* EXPi) / tot)

println("\nchain position of the three target sectors (higher upstr = further from final demand):")
for k in (7, 14, 15, 16)
    @printf("  %-13s upstreamness %.3f   downstreamness %.3f   VS %.1f%%   DVA-from-others %.1f%%\n",
            AGGCOM[k], UP[k], DOWN[k], 100VS[k], 100DVA_oth[k])
end
