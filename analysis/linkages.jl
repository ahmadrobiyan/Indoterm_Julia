#!/usr/bin/env julia
# Backward / forward linkage analysis for INDOTERM (Julia port).
# Uses the 185-sector national IO core, aggregated to the model's 25 sectors.

using LinearAlgebra, Printf, DelimitedFiles, Statistics

const ROOT = raw"C:\Users\ahmad\Downloads\INDOTERM CGE_2016_New\IndotermJulia"
include(joinpath(ROOT, "src", "prepare_sets.jl"))
include(joinpath(ROOT, "src", "aggregation_data.jl"))

const NC = length(COM); const NI = NC
const CIDX = Dict(c => i for (i, c) in enumerate(COM))
const SIDX = Dict("dom" => 1, "imp" => 2)
const MIDX = Dict(m => i for (i, m) in enumerate(MAR))
const NA = length(AGGCOM)
const MAP = SEC_MAP_185_to_25

# ── read national_data.csv sections we need ──────────────────────────────────
"""Return Dict(header => Vector{Tuple{Vector{String},Float64}}) for wanted headers."""
function read_sections(path, wanted::Set{String})
    out = Dict{String,Vector{Tuple{Vector{String},Float64}}}()
    cur = ""
    keep = false
    for line in eachline(path)
        if startswith(line, "__header__:")
            cur = split(line[12:end], ',')[1]
            keep = cur in wanted
            keep && (out[cur] = Tuple{Vector{String},Float64}[])
            continue
        end
        keep || continue
        f = split(line, ',')
        push!(out[cur], (String.(f[1:end-1]), parse(Float64, f[end])))
    end
    return out
end

wanted = Set(["1BAS", "1MAR", "1TAX", "MAKE", "1CAP", "1LAB", "1LND", "1PTX", "1OCT",
              "3BAS", "4BAS", "5BAS", "6BAS", "2BAS", "R001"])
@info "reading national_data.csv …"
S = read_sections(joinpath(ROOT, "data", "national_data.csv"), wanted)

# ── build 185-level arrays ───────────────────────────────────────────────────
V1BAS = zeros(NC, 2, NI)
for (k, v) in S["1BAS"]; V1BAS[CIDX[k[1]], SIDX[k[2]], CIDX[k[3]]] = v; end
V2BAS = zeros(NC, 2, NI)
for (k, v) in S["2BAS"]; V2BAS[CIDX[k[1]], SIDX[k[2]], CIDX[k[3]]] = v; end
V1MAR = zeros(NC, 2, NI, length(MAR))
for (k, v) in S["1MAR"]; V1MAR[CIDX[k[1]], SIDX[k[2]], CIDX[k[3]], MIDX[k[4]]] = v; end
V1TAX = zeros(NC, 2, NI)
for (k, v) in S["1TAX"]; V1TAX[CIDX[k[1]], SIDX[k[2]], CIDX[k[3]]] = v; end
MAKE = zeros(NC, NI)
for (k, v) in S["MAKE"]; MAKE[CIDX[k[1]], CIDX[k[2]]] = v; end
V1CAP = zeros(NI); for (k, v) in S["1CAP"]; V1CAP[CIDX[k[1]]] = v; end
V1LND = zeros(NI); for (k, v) in S["1LND"]; V1LND[CIDX[k[1]]] = v; end
V1PTX = zeros(NI); for (k, v) in S["1PTX"]; V1PTX[CIDX[k[1]]] = v; end
V1OCT = zeros(NI); for (k, v) in S["1OCT"]; V1OCT[CIDX[k[1]]] = v; end
V1LAB = zeros(NI); for (k, v) in S["1LAB"]; V1LAB[CIDX[k[1]]] += v; end
V3BAS = zeros(NC, 2); for (k, v) in S["3BAS"]; V3BAS[CIDX[k[1]], SIDX[k[2]]] = v; end
V4BAS = zeros(NC);    for (k, v) in S["4BAS"]; V4BAS[CIDX[k[1]]] = v; end
V5BAS = zeros(NC, 2); for (k, v) in S["5BAS"]; V5BAS[CIDX[k[1]], SIDX[k[2]]] = v; end
V6BAS = zeros(NC, 2); for (k, v) in S["6BAS"]; V6BAS[CIDX[k[1]], SIDX[k[2]]] = v; end

# ── use table at basic prices, margins re-attributed to margin commodities ───
# Udom[c,j] : domestic commodity c absorbed by industry j (incl. margin services)
# Uimp[c,j] : imported commodity c absorbed by industry j
Udom = zeros(NC, NI); Uimp = zeros(NC, NI)
for j in 1:NI, c in 1:NC
    Udom[c, j] += V1BAS[c, 1, j]
    Uimp[c, j] += V1BAS[c, 2, j]
end
# margins are domestically produced services, whatever the source of the good
for j in 1:NI, m in 1:length(MAR)
    cm = CIDX[MAR[m]]
    Udom[cm, j] += sum(@view V1MAR[:, :, j, m])
end

# ── aggregate 185 -> 25 (flows first, coefficients after) ────────────────────
agg2(X) = begin
    Y = zeros(NA, NA)
    for j in axes(X, 2), i in axes(X, 1); Y[MAP[i], MAP[j]] += X[i, j]; end
    Y
end
agg1(x) = begin
    y = zeros(NA); for i in eachindex(x); y[MAP[i]] += x[i]; end; y
end

UD = agg2(Udom); UM = agg2(Uimp); MK = agg2(MAKE)
CAP = agg1(V1CAP); LAB = agg1(V1LAB); LND = agg1(V1LND)
PTX = agg1(V1PTX); OCT = agg1(V1OCT)
TAXint = agg2(dropdims(sum(V1TAX, dims = 2), dims = 2))
HH  = agg1(V3BAS[:, 1]); EXP = agg1(V4BAS); GOV = agg1(V5BAS[:, 1])
INV = agg1(dropdims(sum(V2BAS[:, 1, :], dims = 2), dims = 2))
STK = agg1(V6BAS[:, 1])

gind = vec(sum(MK, dims = 1))   # industry output (25)
qcom = vec(sum(MK, dims = 2))   # commodity output (25)

# ── commodity x industry -> industry x industry (industry-technology) ────────
# D[i,c] = MK[c,i]/qcom[c]  (market share); A = D * B ; B[c,j] = U[c,j]/g_j
D = zeros(NA, NA)
for c in 1:NA, i in 1:NA
    qcom[c] > 0 && (D[i, c] = MK[c, i] / qcom[c])
end
Bd = zeros(NA, NA); Bm = zeros(NA, NA)
for j in 1:NA, c in 1:NA
    gind[j] > 0 && (Bd[c, j] = UD[c, j] / gind[j]; Bm[c, j] = UM[c, j] / gind[j])
end
Ad = D * Bd                 # domestic industry-by-industry technical coefficients
At = D * (Bd .+ Bm)         # total (incl. imports) — "technological" linkage
Zd = Ad .* reshape(gind, 1, :)   # domestic ind x ind transaction matrix

Ld = inv(I - Ad)
Lt = inv(I - At)

# ── Ghosh (forward) on the domestic transaction matrix ───────────────────────
Bg = zeros(NA, NA)
for i in 1:NA, j in 1:NA
    gind[i] > 0 && (Bg[i, j] = Zd[i, j] / gind[i])
end
G = inv(I - Bg)

BL   = vec(sum(Ld, dims = 1))          # backward, Leontief column sums
FL   = vec(sum(G,  dims = 2))          # forward, Ghosh row sums
BLt  = vec(sum(Lt, dims = 1))
BLi  = BL ./ (sum(BL) / NA)            # Rasmussen power of dispersion
FLi  = FL ./ (sum(FL) / NA)            # Rasmussen sensitivity of dispersion
# dispersion (coefficient of variation) — lower = more evenly spread
Vj = [ sqrt(sum((Ld[:, j] .- mean(Ld[:, j])).^2) / (NA - 1)) / mean(Ld[:, j]) for j in 1:NA ]
Vi = [ sqrt(sum((G[i, :]  .- mean(G[i, :])).^2)  / (NA - 1)) / mean(G[i, :])  for i in 1:NA ]

# ── hypothetical extraction (total linkage) ─────────────────────────────────
fd = HH .+ EXP .+ GOV .+ INV .+ STK           # domestic final demand by commodity
fdi = D * fd                                   # allocated to industries
xbase = Ld * fdi
HEM = zeros(NA)
for k in 1:NA
    Ak = copy(Ad); Ak[k, :] .= 0; Ak[:, k] .= 0
    fk = copy(fdi); fk[k] = 0
    xk = inv(I - Ak) * fk
    HEM[k] = sum(xbase) - sum(xk)
end
HEMpc = 100 .* HEM ./ sum(xbase)

# ── output ──────────────────────────────────────────────────────────────────
open(joinpath(@__DIR__, "linkage_results.csv"), "w") do io
    println(io, "sector,output,BL_dom,BLindex_dom,FL_dom,FLindex_dom,BL_total,Vj,Vi,HEM_loss_pct,share_of_output_pct")
    for k in 1:NA
        @printf(io, "%s,%.1f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                AGGCOM[k], gind[k], BL[k], BLi[k], FL[k], FLi[k], BLt[k], Vj[k], Vi[k],
                HEMpc[k], 100 * gind[k] / sum(gind))
    end
end

println("\n=== Rasmussen linkage indices (domestic, industry x industry, 25 sectors) ===")
@printf("%-14s %10s %8s %8s %8s %8s %8s\n", "sector", "output", "BL", "BLidx", "FL", "FLidx", "HEM%")
for k in sortperm(BLi, rev = true)
    @printf("%-14s %10.0f %8.3f %8.3f %8.3f %8.3f %8.2f\n",
            AGGCOM[k], gind[k], BL[k], BLi[k], FL[k], FLi[k], HEMpc[k])
end

# ── detail for the three target sectors ─────────────────────────────────────
TARGETS = [14, 15, 16]   # BasMetals, MetalMach, TranspEquip
println("\n=== Detail: top domestic suppliers (backward) and users (forward) ===")
for k in TARGETS
    println("\n--- $(AGGCOM[k])  (output $(round(gind[k], digits=0)), $(round(100*gind[k]/sum(gind), digits=2))% of gross output) ---")
    col = Zd[:, k]; ord = sortperm(col, rev = true)[1:8]
    println("  suppliers (share of $(AGGCOM[k]) total intermediate cost):")
    tot = sum(col)
    for i in ord
        col[i] <= 0 && continue
        @printf("    %-14s %8.2f%%   (A=%.4f)\n", AGGCOM[i], 100 * col[i] / tot, Ad[i, k])
    end
    @printf("    imported inputs / total intermediate = %.1f%%\n",
            100 * sum(UM[:, k]) / (sum(UM[:, k]) + sum(UD[:, k])))
    row = Zd[k, :]; ord2 = sortperm(row, rev = true)[1:8]
    println("  users (share of $(AGGCOM[k]) intermediate sales):")
    tot2 = sum(row)
    for j in ord2
        row[j] <= 0 && continue
        @printf("    %-14s %8.2f%%   (B=%.4f)\n", AGGCOM[j], 100 * row[j] / tot2, Bg[k, j])
    end
    fdshare = 100 * fdi[k] / gind[k]
    @printf("  final demand absorbs %.1f%% of output;  intermediate sales %.1f%%\n",
            fdshare, 100 * tot2 / gind[k])
end

# ── final-demand composition & import penetration (commodity basis) ────────
println("\n=== Commodity balance: where output goes, and import share of the domestic market ===")
@printf("%-13s %8s %8s %8s %8s %8s %8s %8s\n",
        "commodity", "interm%", "hhold%", "invest%", "govt%", "export%", "stock%", "impPen%")
UDc = vec(sum(UD, dims = 2))
V2c = agg1(dropdims(sum(V2BAS[:, 1, :], dims = 2), dims = 2))
V2ci = agg1(dropdims(sum(V2BAS[:, 2, :], dims = 2), dims = 2))
UMc = vec(sum(UM, dims = 2))
HHi = agg1(V3BAS[:, 2]); GOVi = agg1(V5BAS[:, 2]); STKi = agg1(V6BAS[:, 2])
for k in TARGETS
    tot = UDc[k] + HH[k] + V2c[k] + GOV[k] + EXP[k] + STK[k]
    domuse = UDc[k] + HH[k] + V2c[k] + GOV[k] + STK[k]
    impuse = UMc[k] + HHi[k] + V2ci[k] + GOVi[k] + STKi[k]
    @printf("%-13s %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n", AGGCOM[k],
            100UDc[k]/tot, 100HH[k]/tot, 100V2c[k]/tot, 100GOV[k]/tot,
            100EXP[k]/tot, 100STK[k]/tot, 100impuse/(domuse+impuse))
end

println("\n=== Domestic vs total (import-inclusive) backward linkage ===")
@printf("%-13s %10s %10s %10s\n", "sector", "BL_dom", "BL_total", "leakage")
for k in TARGETS
    @printf("%-13s %10.3f %10.3f %10.1f%%\n", AGGCOM[k], BL[k], BLt[k],
            100 * (BLt[k] - BL[k]) / BLt[k])
end

# ── the 3-sector bloc treated as one: internal vs external dependence ───────
println("\n=== Metals–machinery–transport bloc (sectors 14,15,16) ===")
blocout = sum(gind[TARGETS])
@printf("bloc gross output = %.0f  (%.2f%% of national)\n", blocout, 100 * blocout / sum(gind))
intra = sum(Zd[TARGETS, TARGETS])
@printf("intra-bloc intermediate flows = %.0f (%.2f%% of bloc intermediate purchases)\n",
        intra, 100 * intra / sum(Zd[:, TARGETS]))
# joint hypothetical extraction
Ak = copy(Ad); Ak[TARGETS, :] .= 0; Ak[:, TARGETS] .= 0
fk = copy(fdi); fk[TARGETS] .= 0
@printf("joint hypothetical extraction loss = %.2f%% of gross output (sum of individual = %.2f%%)\n",
        100 * (sum(xbase) - sum(inv(I - Ak) * fk)) / sum(xbase), sum(HEMpc[TARGETS]))

# ── regional distribution using R001 (industry share by region) ─────────────
if haskey(S, "R001")
    R001 = zeros(NI, length(REG))
    RIDX = Dict(r => i for (i, r) in enumerate(REG))
    for (k, v) in S["R001"]; R001[CIDX[k[1]], RIDX[k[2]]] = v; end
    # regional output of each agg sector = sum_i R001[i,r]*g185_i
    g185 = vec(sum(MAKE, dims = 1))
    RA = zeros(NA, length(REG))
    for r in 1:length(REG), i in 1:NI
        RA[MAP[i], r] += R001[i, r] * g185[i]
    end
    println("\n=== Regional concentration of the three sectors (top 8 provinces, % of sector output) ===")
    for k in TARGETS
        tot = sum(RA[k, :]); ord = sortperm(RA[k, :], rev = true)[1:8]
        print(rpad(AGGCOM[k], 14), ": ")
        println(join([@sprintf("%s %.1f%%", REG[r], 100 * RA[k, r] / tot) for r in ord], ", "))
    end
    writedlm(joinpath(@__DIR__, "regional_shares.csv"),
             vcat(reshape(vcat(["sector"], REG), 1, :), hcat(AGGCOM, RA)), ',')
end
