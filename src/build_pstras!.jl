"""
Translation of `pstras.tab` — check RAS results, reconcile regional IO and
trade data, recompute DISTGONE, compute FLAG for convergence.

Input: `RasResult.ras` + reg1 data (for MAKE/V1LAB/V1CAP/V1LND).
Output: `PstrasResult`.
"""

using NamedArrays

Base.@kwdef struct PstrasResult
    pstras::Dict{String,Any}
    diag::Dict{String,Any}
    converged::Bool
end

function build_pstras!(ras::Dict{String,Any}, reg1::Dict{String,Any},
                       old_distgone::Union{Dict{String,Any},Nothing}=nothing)
    T = Float64
    NC = length(COM); NI = length(IND); NS = 2; NM = length(MAR); NR = length(REG); NO = length(OCC)
    NONMAR_idx = [c for c in 1:NC if !(COM[c] in Set(MAR))]

    # ── Read inputs ─────────────────────────────────────────────────────────
    TRADE   = ras["TRAD"]  # COM×SRC×ORG×DST
    TRADMAR = ras["TMAR"]  # COM×SRC×MAR×ORG×DST
    SUPPMAR = ras["MARS"]  # MAR×ORG×DST×PRD
    BASIC_U = ras["BSCU"]
    MARGINS_U = ras["MRGU"]
    IMPORT_cr = ras["IMPS"]
    MAKE_I_cd0 = ras["COST"]

    # ── Read MAKE, V1LAB, V1CAP, V1LND from reg1 ───────────────────────────
    MAKE_r   = reg1["MAKE"]  # COM×IND×DST
    V1LAB_io = haskey(reg1, "FACT") ? nothing : nothing  # We'll recompute below
    # FAC_r from reg1: IND×FACTOR×DST
    FAC_r    = reg1["FACT"]  # IND × 4 × DST

    V1LAB_iod = zeros(T, NI, NO, NR)
    V1CAP_id  = zeros(T, NI, NR)
    V1LND_id  = zeros(T, NI, NR)
    V1PTX_id  = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        g_lab = findfirst(==("Labour"),   ["Labour","Capital","Land","ProdTax"])
        g_cap = findfirst(==("Capital"),  ["Labour","Capital","Land","ProdTax"])
        g_lnd = findfirst(==("Land"),     ["Labour","Capital","Land","ProdTax"])
        g_ptx = findfirst(==("ProdTax"),  ["Labour","Capital","Land","ProdTax"])
        for o in 1:NO; V1LAB_iod[i,o,d] = 0.0; end
        V1CAP_id[i,d]  = FAC_r[i,g_cap,d]
        V1LND_id[i,d]  = FAC_r[i,g_lnd,d]
        V1PTX_id[i,d]  = FAC_r[i,g_ptx,d]
    end

    # ── MAKE_I ──────────────────────────────────────────────────────────────
    MAKE_I_cd = zeros(T, NC, NR)
    for c in 1:NC, d in 1:NR
        MAKE_I_cd[c,d] = sum(MAKE_r[c,i,d] for i in 1:NI)
    end

    # ── TRADMAR_CS (86-89) ─────────────────────────────────────────────────
    TRADMAR_CS_mrd = zeros(T, NM, NR, NR)
    for m in 1:NM, r in 1:NR, d in 1:NR
        TRADMAR_CS_mrd[m,r,d] = sum(TRADMAR[c,s,m,r,d] for c in 1:NC, s in 1:NS)
    end

    # ── Check conditions ────────────────────────────────────────────────────
    # TRADE vs BASIC_U (93-97)
    TEMTOTcsd = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        TEMTOTcsd[c,s,d] = BASIC_U[c,s,d] - sum(TRADE[c,s,:,d])
    end

    # TRADMAR vs MARGINS_U (99-106)
    TEMTOTcsmd = zeros(T, NC, NS, NM, NR)
    for c in 1:NC, s in 1:NS, m in 1:NM, d in 1:NR
        TEMTOTcsmd[c,s,m,d] = MARGINS_U[c,s,m,d] - sum(TRADMAR[c,s,m,:,d])
    end

    # TRADE(imp) vs IMPORT (108-114)
    IMPERR = zeros(T, NC, NR)
    for c in 1:NC, r in 1:NR
        IMPERR[c,r] = IMPORT_cr[c,r] - sum(TRADE[c,2,r,d] for d in 1:NR)
    end

    # TRADE(dom) + SUPPMAR vs MAKE_I (116-124)
    DOMERR = zeros(T, NC, NR)
    for c in NONMAR_idx
        for r in 1:NR
            DOMERR[c,r] = MAKE_I_cd[c,r] - sum(TRADE[c,1,r,d] for d in 1:NR)
        end
    end
    for m in 1:NM
        for r in 1:NR
            dom_sum = sum(TRADE[m,1,r,d] for d in 1:NR)
            supp_sum = sum(SUPPMAR[m,rr,d,r] for rr in 1:NR, d in 1:NR)
            DOMERR[m,r] = MAKE_I_cd[m,r] - dom_sum - supp_sum
        end
    end

    # SUPPMAR vs TRADMAR_CS (126-136)
    TEMTOTmrd = zeros(T, NM, NR, NR)
    for m in 1:NM, r in 1:NR, d in 1:NR
        TEMTOTmrd[m,r,d] = TRADMAR_CS_mrd[m,r,d] - sum(SUPPMAR[m,r,d,p] for p in 1:NR)
    end

    # ── Adjustments (138-183) ───────────────────────────────────────────────
    # 1. Adjust IMPORT by error (139)
    IMPORT_adj = zeros(T, NC, NR)
    for c in 1:NC, r in 1:NR
        IMPORT_adj[c,r] = IMPORT_cr[c,r] - IMPERR[c,r]
    end

    # 2. Factor/MAKE adjustments (140-183)
    INDERR = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        for c in 1:NC
            if MAKE_I_cd[c,d] > 0
                INDERR[i,d] += DOMERR[c,d] * MAKE_r[c,i,d] / MAKE_I_cd[c,d]
            end
        end
    end

    FACTOR_id = zeros(T, NI, NR)
    FACTOR2_id = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        g_lab = findfirst(==("Labour"),   ["Labour","Capital","Land","ProdTax"])
        g_cap = findfirst(==("Capital"),  ["Labour","Capital","Land","ProdTax"])
        g_lnd = findfirst(==("Land"),     ["Labour","Capital","Land","ProdTax"])
        FACTOR_id[i,d] = FAC_r[i,g_cap,d] + FAC_r[i,g_lnd,d] + sum(V1LAB_iod[i,:,d])
        FACTOR2_id[i,d] = max(0.0, FACTOR_id[i,d] - INDERR[i,d])
    end

    # SCALEFAC
    MAKE_C_id = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        MAKE_C_id[i,d] = sum(MAKE_r[c,i,d] for c in 1:NC)
    end
    SCALEFAC = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        SCALEFAC[i,d] = 1.0
        if FACTOR_id[i,d] > 0.01 * MAKE_C_id[i,d]
            SCALEFAC[i,d] = FACTOR2_id[i,d] / FACTOR_id[i,d]
        end
        SCALEFAC[i,d] = clamp(SCALEFAC[i,d], 0.7, 1.5)
    end

    # Apply SCALEFAC
    for i in 1:NI, d in 1:NR
        for o in 1:NO
            V1LAB_iod[i,o,d] *= SCALEFAC[i,d]
        end
        V1CAP_id[i,d] *= SCALEFAC[i,d]
        V1LND_id[i,d] *= SCALEFAC[i,d]
    end
    delMAKE = zeros(T, NC, NI, NR)
    for c in 1:NC, i in 1:NI, d in 1:NR
        if MAKE_C_id[i,d] > 0
            delMAKE[c,i,d] = INDERR[i,d] * MAKE_r[c,i,d] / MAKE_C_id[i,d]
        end
    end
    for c in 1:NC, i in 1:NI, d in 1:NR
        if SCALEFAC[i,d] != 1.0
            MAKE_r[c,i,d] -= delMAKE[c,i,d]
        end
    end
    # Recompute MAKE_I
    for c in 1:NC, d in 1:NR
        MAKE_I_cd[c,d] = sum(MAKE_r[c,i,d] for i in 1:NI)
    end

    # ── Recompute DISTGONE (258-277) ────────────────────────────────────────
    DISTANCE = reg1["MAKE"]  # not right — need DISTANCE from reg2
    # Actually DISTANCE is in ras["DIST"]
    DISTANCE = ras["DIST"]  # ORG×DST (NR×NR)

    TRADE_R_csd = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        TRADE_R_csd[c,s,d] = sum(TRADE[c,s,:,d])
    end

    DISTGONE = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        if TRADE_R_csd[c,s,d] > 0
            num = sum(DISTANCE[r,d] * TRADE[c,s,r,d] for r in 1:NR)
            DISTGONE[c,s,d] = num / TRADE_R_csd[c,s,d]
        else
            DISTGONE[c,s,d] = sum(DISTANCE[:,d]) / NR
        end
    end

    # ── Compare with old DISTGONE ───────────────────────────────────────────
    converged = false
    if old_distgone !== nothing
        OLDDISTGONE = old_distgone["DGON"]
        DIFFDIST = zeros(T, NC, NS, NR)
        for c in 1:NC, s in 1:NS, d in 1:NR
            DIFFDIST[c,s,d] = DISTGONE[c,s,d] - OLDDISTGONE[c,s,d]
        end
        TOTTRADE = sum(TRADE_R_csd)
        TOTABSDIFF = sum(abs.(DIFFDIST))
        AVEABSDIFF = TOTTRADE > 0 ? sum(TRADE_R_csd .* abs.(DIFFDIST)) / TOTTRADE : 0.0
        converged = AVEABSDIFF < 0.001
    end

    # ── Output ──────────────────────────────────────────────────────────────
    pstras_out = Dict{String,Any}(
        "IMPS" => NamedArray(IMPORT_adj, Tuple([COM,REG]), (:COM,:ORG)),
        "MAKE" => NamedArray(MAKE_r, Tuple([COM,IND,REG]), (:COM,:IND,:DST)),
        "1LAB" => NamedArray(V1LAB_iod, Tuple([IND,OCC,REG]), (:IND,:OCC,:DST)),
        "1CAP" => NamedArray(V1CAP_id, Tuple([IND,REG]), (:IND,:DST)),
        "1LND" => NamedArray(V1LND_id, Tuple([IND,REG]), (:IND,:DST)),
        "1PTX" => NamedArray(V1PTX_id, Tuple([IND,REG]), (:IND,:DST)),
        "TRAD" => ras["TRAD"],
        "MARS" => ras["MARS"],
        "TMAR" => ras["TMAR"],
        "DIST" => DISTANCE,
        "DGON" => NamedArray(DISTGONE, Tuple([COM,SRC,REG]), (:COM,:SRC,:DST)),
        # Pass through USE and TAXES for premod
        "BSMR" => haskey(ras, "BSMR") ? ras["BSMR"] : nothing,
        "UTAX" => haskey(reg1, "UTAX") ? reg1["UTAX"] : nothing,
        "2PUR" => haskey(ras, "2PUR") ? ras["2PUR"] : nothing,
    )

    PstrasResult(pstras=pstras_out, diag=Dict{String,Any}(), converged=converged)
end
