"""
TKDN screening at the 185-commodity level — turning the three conditions stated
in the executive summary into an operational positive list.

The claim being tested is that a targeted local-content list can be derived from
the database itself rather than negotiated sector by sector. Three screens:

  S1  Is there a domestic supplier, and could it plausibly scale?
      MAKE row sum > 0; import penetration high enough that a mandate has
      something to bite on; domestic output large relative to the imports it
      would have to displace.

  S2  Does the multiplier stay home?
      Domestic value added generated per unit of the commodity supplied,
      direct + indirect, split into the making industry's OWN value added and
      value added pulled from OTHER domestic sectors. The second piece is the
      whole point of a local-content rule — if it is low, forcing the purchase
      onshore just relocates one assembly margin.

  S3  Is it a capital good?
      Share of domestic absorption going to capital formation, to Construction,
      and to BasMetals itself. A mandate on these items taxes the smelter build
      that hilirisasi depends on — the two policies cancel.

A commodity joins the positive list only if it passes all three.

Value basis: basic values, with margins re-attributed to the margin commodities
exactly as in `analysis/linkages.jl`, so the A matrix here is the same object.
Commodity tax wedges are not in the Leontief core, so vcoef + input coefficients
do not sum to exactly one; the residual is small and does not move the ranking.

Run:  julia IndotermJulia/analysis/tkdn_screening.jl
"""

using LinearAlgebra, Printf, Statistics, DelimitedFiles

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "prepare_sets.jl"))
include(joinpath(ROOT, "src", "aggregation_data.jl"))

const NC = length(COM); const NI = NC
const CIDX = Dict(c => i for (i, c) in enumerate(COM))
const SIDX = Dict("dom" => 1, "imp" => 2)
const MIDX = Dict(m => i for (i, m) in enumerate(MAR))
const MAP = SEC_MAP_185_to_25

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
                  Set(["1BAS", "1MAR", "2BAS", "3BAS", "4BAS", "5BAS", "6BAS",
                       "MAKE", "1CAP", "1LAB", "1LND", "1PTX", "1OCT"]))

V1BAS = zeros(NC, 2, NI); for (k, v) in S["1BAS"]; V1BAS[CIDX[k[1]], SIDX[k[2]], CIDX[k[3]]] = v; end
V2BAS = zeros(NC, 2, NI); for (k, v) in S["2BAS"]; V2BAS[CIDX[k[1]], SIDX[k[2]], CIDX[k[3]]] = v; end
V3BAS = zeros(NC, 2);     for (k, v) in S["3BAS"]; V3BAS[CIDX[k[1]], SIDX[k[2]]] = v; end
V5BAS = zeros(NC, 2);     for (k, v) in S["5BAS"]; V5BAS[CIDX[k[1]], SIDX[k[2]]] = v; end
V6BAS = zeros(NC, 2);     for (k, v) in S["6BAS"]; V6BAS[CIDX[k[1]], SIDX[k[2]]] = v; end
V4BAS = zeros(NC);        for (k, v) in S["4BAS"]; V4BAS[CIDX[k[1]]] = v; end
MAKE  = zeros(NC, NI);    for (k, v) in S["MAKE"]; MAKE[CIDX[k[1]], CIDX[k[2]]] = v; end
V1MAR = zeros(NC, 2, NI, length(MAR))
for (k, v) in S["1MAR"]; V1MAR[CIDX[k[1]], SIDX[k[2]], CIDX[k[3]], MIDX[k[4]]] = v; end

VA = zeros(NI)
for h in ("1CAP", "1LND", "1PTX", "1OCT"); for (k, v) in S[h]; VA[CIDX[k[1]]] += v; end; end
for (k, v) in S["1LAB"]; VA[CIDX[k[1]]] += v; end

# ── Leontief core, identical construction to linkages.jl ────────────────────
Udom = V1BAS[:, 1, :] |> copy
for j in 1:NI, m in 1:length(MAR)
    Udom[CIDX[MAR[m]], j] += sum(@view V1MAR[:, :, j, m])
end
gind = vec(sum(MAKE, dims = 1)); qcom = vec(sum(MAKE, dims = 2))
D = zeros(NI, NC); for c in 1:NC, i in 1:NI; qcom[c] > 0 && (D[i, c] = MAKE[c, i] / qcom[c]); end
Bd = zeros(NC, NI)
for j in 1:NI, c in 1:NC; gind[j] > 0 && (Bd[c, j] = Udom[c, j] / gind[j]); end
Ad = D * Bd
Ld = inv(I - Ad)
vcoef = [gind[j] > 0 ? VA[j] / gind[j] : 0.0 for j in 1:NI]

# ── S1: domestic supplier and room to scale ────────────────────────────────
intdom = vec(sum(V1BAS[:, 1, :], dims = 2)); intimp = vec(sum(V1BAS[:, 2, :], dims = 2))
invdom = vec(sum(V2BAS[:, 1, :], dims = 2)); invimp = vec(sum(V2BAS[:, 2, :], dims = 2))
absdom = intdom .+ invdom .+ V3BAS[:, 1] .+ V5BAS[:, 1] .+ V6BAS[:, 1]
absimp = intimp .+ invimp .+ V3BAS[:, 2] .+ V5BAS[:, 2] .+ V6BAS[:, 2]
ABS    = absdom .+ absimp                       # domestic-market absorption, excl. exports

IMPSHR   = [ABS[c] > 0 ? absimp[c] / ABS[c] : 0.0 for c in 1:NC]
CAPRATIO = [absimp[c] > 0 ? qcom[c] / absimp[c] : Inf for c in 1:NC]

# ── S2: domestic value added pulled per unit supplied ──────────────────────
DVA_tot = zeros(NC); DVA_own = zeros(NC); DVA_oth = zeros(NC)
for c in 1:NC
    qcom[c] <= 0 && continue
    x = Ld * D[:, c]                            # gross output pulled by 1 unit of commodity c
    DVA_tot[c] = sum(vcoef .* x)
    DVA_own[c] = sum(D[i, c] * vcoef[i] * Ld[i, i] for i in 1:NI)
    DVA_oth[c] = DVA_tot[c] - DVA_own[c]
end

# ── S3: capital-goods intensity ────────────────────────────────────────────
const CONSTR = [i for i in 1:NI if MAP[i] == 20]      # Construction industries
const BASMET = [i for i in 1:NI if MAP[i] == 14]      # BasMetals industries
inv_use    = invdom .+ invimp
constr_use = vec(sum(V1BAS[:, 1, CONSTR], dims = 2)) .+ vec(sum(V1BAS[:, 2, CONSTR], dims = 2))
basmet_use = vec(sum(V1BAS[:, 1, BASMET], dims = 2)) .+ vec(sum(V1BAS[:, 2, BASMET], dims = 2))
KAPINT = [ABS[c] > 0 ? (inv_use[c] + constr_use[c] + basmet_use[c]) / ABS[c] : 0.0 for c in 1:NC]

# ── thresholds ─────────────────────────────────────────────────────────────
const T_IMPSHR = 0.15     # below this a mandate has nothing to bite on
const T_CAPRAT = 1.00     # domestic output must at least match the imports displaced
const T_KAPINT = 0.25     # above this the item is mainly a capital good
elig = [qcom[c] > 0 && ABS[c] > 0 for c in 1:NC]
T_DVAOTH = median(DVA_oth[elig .& (IMPSHR .>= T_IMPSHR)])   # data-driven, not assumed

p1 = [elig[c] && IMPSHR[c] >= T_IMPSHR && CAPRATIO[c] >= T_CAPRAT for c in 1:NC]
p2 = [DVA_oth[c] >= T_DVAOTH for c in 1:NC]
p3 = [KAPINT[c] <= T_KAPINT for c in 1:NC]
pass = p1 .& p2 .& p3

println("\n" * "="^100)
println("TKDN SCREENING — 185 commodities")
println("="^100)
@printf("thresholds: import penetration >= %.0f%%, dom.output/imports >= %.2f, ", 100T_IMPSHR, T_CAPRAT)
@printf("DVA-from-others >= %.4f (median of eligible), capital intensity <= %.0f%%\n\n", T_DVAOTH, 100T_KAPINT)
@printf("commodities with a domestic supplier at all ......... %3d / %d\n", count(elig), NC)
@printf("  S1  supplier exists AND room to scale ............. %3d\n", count(p1))
@printf("  S2  multiplier stays home ......................... %3d  (of those passing S1: %d)\n",
        count(p2), count(p1 .& p2))
@printf("  S3  not a capital good ............................ %3d  (of those passing S1+S2: %d)\n",
        count(p3), count(pass))
@printf("\nPOSITIVE LIST: %d commodities = %.1f%% of the 185\n", count(pass), 100count(pass) / NC)
@printf("they carry %.1f%% of total domestic-market absorption and %.1f%% of all imports\n",
        100sum(ABS[pass]) / sum(ABS), 100sum(absimp[pass]) / sum(absimp))

println("\n" * "="^100)
println("THE POSITIVE LIST  (sorted by imports at stake)")
println("="^100)
@printf("%-42s %10s %8s %8s %8s %8s  %s\n",
        "commodity", "imports", "imp%", "dom/imp", "DVAoth", "capint", "group")
for c in sortperm(absimp, rev = true)
    pass[c] || continue
    @printf("%-42s %10.0f %7.1f%% %8.2f %8.3f %7.1f%%  %s\n",
            COM[c], absimp[c], 100IMPSHR[c], min(CAPRATIO[c], 99.99),
            DVA_oth[c], 100KAPINT[c], AGGCOM[MAP[c]])
end

println("\n" * "="^100)
println("REJECTED — and the reason. Top 25 by imports at stake.")
println("="^100)
why(c) = begin
    r = String[]
    !elig[c] && push!(r, "NO DOMESTIC SUPPLIER")
    elig[c] && IMPSHR[c] < T_IMPSHR && push!(r, @sprintf("already %.0f%% domestic", 100(1 - IMPSHR[c])))
    elig[c] && IMPSHR[c] >= T_IMPSHR && CAPRATIO[c] < T_CAPRAT &&
        push!(r, @sprintf("cannot scale (dom/imp %.2f)", CAPRATIO[c]))
    !p2[c] && push!(r, @sprintf("thin multiplier (%.3f)", DVA_oth[c]))
    !p3[c] && push!(r, @sprintf("CAPITAL GOOD (%.0f%%)", 100KAPINT[c]))
    join(r, "; ")
end
for c in first([c for c in sortperm(absimp, rev = true) if !pass[c] && absimp[c] > 0], 25)
    @printf("%-42s %10.0f  %s\n", COM[c], absimp[c], why(c))
end

println("\n" * "="^100)
println("WHERE A MANDATE CANNOT BE MET AT ALL  (no domestic production, imports > 0)")
println("="^100)
nosup = [c for c in 1:NC if qcom[c] <= 0 && absimp[c] > 0]
if isempty(nosup)
    println("none — every imported commodity has some domestic production in the 2016 base.")
else
    for c in sort(nosup, by = c -> -absimp[c])
        @printf("%-42s imports %10.0f   group %s\n", COM[c], absimp[c], AGGCOM[MAP[c]])
    end
end

println("\n" * "="^100)
println("BY 25-SECTOR GROUP  (how much of each group's import bill the positive list covers)")
println("="^100)
@printf("%-14s %5s %5s %12s %12s %7s\n", "group", "n", "pass", "imports", "on list", "cover%")
for g in 1:length(AGGCOM)
    idx = [c for c in 1:NC if MAP[c] == g]
    isempty(idx) && continue
    ti = sum(absimp[idx]); ti <= 0 && continue
    pi = sum(absimp[idx][pass[idx]])
    @printf("%-14s %5d %5d %12.0f %12.0f %6.1f%%\n",
            AGGCOM[g], length(idx), count(pass[idx]), ti, pi, 100pi / ti)
end

println("\n" * "="^100)
println("THE THREE TARGET SECTORS, COMMODITY BY COMMODITY")
println("="^100)
for g in (14, 15, 16)
    println("\n$(AGGCOM[g]):")
    @printf("  %-40s %10s %7s %8s %8s %7s %6s\n",
            "commodity", "imports", "imp%", "dom/imp", "DVAoth", "capint", "verdict")
    for c in sort([c for c in 1:NC if MAP[c] == g], by = c -> -absimp[c])
        @printf("  %-40s %10.0f %6.1f%% %8.2f %8.3f %6.1f%%  %s\n",
                COM[c], absimp[c], 100IMPSHR[c], min(CAPRATIO[c], 99.99),
                DVA_oth[c], 100KAPINT[c], pass[c] ? "PASS" : "reject")
    end
end

open(joinpath(@__DIR__, "tkdn_screening.csv"), "w") do io
    println(io, "commodity,group,absorption,imports,import_share,dom_over_imp," *
                "DVA_total,DVA_own,DVA_others,capital_intensity,S1,S2,S3,positive_list")
    for c in 1:NC
        @printf(io, "%s,%s,%.2f,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d\n",
                COM[c], AGGCOM[MAP[c]], ABS[c], absimp[c], IMPSHR[c],
                isfinite(CAPRATIO[c]) ? CAPRATIO[c] : -1.0,
                DVA_tot[c], DVA_own[c], DVA_oth[c], KAPINT[c],
                p1[c], p2[c], p3[c], pass[c])
    end
end
println("\nwritten: analysis/tkdn_screening.csv")
