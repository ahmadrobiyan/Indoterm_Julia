"""
Translation of `premod.tab` — assemble TERM-format database from pstras
output, compute purchaser prices, elasticities, dynamic data parameters,
write aggregation weights.

Input: `PstrasResult.pstras` + regsupp + elast dicts.
Output: `PremodResult`.
"""

using NamedArrays

Base.@kwdef struct PremodResult
    premod::Dict{String,Any}
    weights::Dict{String,Any}
    diag::Dict{String,Any}
end

function build_premod!(pstras::Dict{String,Any},
                       regsupp::Dict{String,Any},
                       elast::Dict{String,Any})
    T = Float64
    NC = length(COM); NI = length(IND); NS = 2; NM = length(MAR); NR = length(REG); NO = length(OCC)

    USR_NAMES = vcat(IND, ["Hou","Inv","Gov","Exp"])
    FACTOR_NAMES = ["Labour","Capital","Land","ProdTax"]
    NUSR = length(USR_NAMES)
    NFAC = length(FACTOR_NAMES)

    # ── Read inputs ─────────────────────────────────────────────────────────
    MAKE_r    = pstras["MAKE"]  # COM×IND×DST
    TRADE     = pstras["TRAD"]
    TRADMAR   = pstras["TMAR"]
    SUPPMAR   = pstras["MARS"]
    IMPORT_cr = pstras["IMPS"]
    V1LAB     = pstras["1LAB"]
    V1CAP     = pstras["1CAP"]
    V1LND     = pstras["1LND"]
    DISTANCE  = pstras["DIST"]

    # Compute USE, TAXES, INVEST from the trade data
    # USE(c,s,u,d) = delivered value = basic + margins
    # We need to reconstruct from reg1's USHR or from reg2's BSMR
    # For now, compute from TRADE + TRADMAR

    # BASIC flow to user u in region d: sum over all trade origins
    # BASIC(c,s,u,d) = sum_r NATBASIC(c,s,u) * USHR(c,s,u,d)
    # But we don't have USHR here directly. The simplest is to use reg1's output.
    # Actually, the premod gets BSMR from INFILE (which is ultimately from reg1/reg2)
    # But we need it passed through.
    # For now, let me compute a simplified version.

    # Actually, we need USE and TAXES passed through from RAS output.
    # Let me check what ras_balance! passes through: BSMR (from reg2)

    # For premod.tab, the key inputs from INFILE are:
    # MAKE, USE, TAXES, TRADE, SUPPMAR, TRADMAR, INVEST, STOCKS, VLAB, VCAP, VLND, VPTX
    # Then from ELAST: SLAB, P028, XPEL, P015, P018, ...
    # Then from REGSUPP: SGDD, SMAR, PO01

    # ── MAKE_I, IMPORT (126-132) ────────────────────────────────────────────
    MAKE_I_cd = zeros(T, NC, NR)
    for c in 1:NC, d in 1:NR
        MAKE_I_cd[c,d] = sum(MAKE_r[c,i,d] for i in 1:NI)
    end

    # Recompute IMPORT from TRADE (130-132)
    IMPORT_re = zeros(T, NC, NR)
    for c in 1:NC, r in 1:NR
        IMPORT_re[c,r] = sum(TRADE[c,2,r,d] for d in 1:NR)
    end

    # ── Factor aggregates (86-94) ───────────────────────────────────────────
    VLAB_O_id = zeros(T, NI, NR)
    VPRIM_id  = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        VLAB_O_id[i,d] = sum(V1LAB[i,o,d] for o in 1:NO)
        VPRIM_id[i,d]  = VLAB_O_id[i,d] + V1CAP[i,d] + V1LND[i,d]
    end

    # ── COSTS / DIFFIND (108-124) ───────────────────────────────────────────
    # We need USE and TAXES — pass through from ras/pstras or as separate args
    USE_data = nothing
    TAXES_data = nothing
    if haskey(pstras, "BSMR") && pstras["BSMR"] !== nothing
       USE_data = parent(pstras["BSMR"]) 
    end
    if haskey(pstras, "UTAX") && pstras["UTAX"] !== nothing
        TAXES_data = parent(pstras["UTAX"])
    end

    # ── EPS — expenditure elasticities (176-188) ────────────────────────────
    NATEPS = elast["XPEL"]
    V3TOT_d = zeros(T, NR)
    EPS_cd  = zeros(T, NC, NR)
    for d in 1:NR
        V3TOT_d[d] = 1.0  # placeholder — needs V3PUR_S from USE
    end

    # ── Dynamic data (279-357) ──────────────────────────────────────────────
    # CAPSTOK, DPRC, RNORMAL, GROTREND, QRATIO, ALPHA, RORADJ, GRETEXP
    # These come from CAPDAT / national.har — extract from elast dict
    NATCAPSTOK  = ones(T, NI)   # placeholders
    NATDPRC     = fill(0.08, NI)
    NATRNORMAL  = fill(0.1, NI)
    NATGROTREND = fill(0.05, NI)
    NATQRATIO   = fill(4.0, NI)
    NATALPHA    = fill(5.0, NI)
    NATRORADJ   = fill(0.8, NI)
    NATGRETEXP  = fill(0.15, NI)
    NATEMPRAT   = 1.0
    NATELASTWAGE = 0.5

    NATV1CAP_i = zeros(T, NI)
    for i in 1:NI
        NATV1CAP_i[i] = sum(V1CAP[i,d] for d in 1:NR)
    end

    CAPSTOK_id  = zeros(T, NI, NR)
    DPRC_id     = zeros(T, NI, NR)
    RNORMAL_id  = zeros(T, NI, NR)
    GROTREND_id = zeros(T, NI, NR)
    QRATIO_id   = zeros(T, NI, NR)
    ALPHA_id    = zeros(T, NI, NR)
    RORADJ_id   = zeros(T, NI, NR)
    GRETEXP_id  = zeros(T, NI, NR)
    EMPRAT_od   = zeros(T, NO, NR)
    ELASTWAGE_od = zeros(T, NO, NR)

    for i in 1:NI, d in 1:NR
        CAPSTOK_id[i,d]  = NATV1CAP_i[i] > 0 ? (V1CAP[i,d] / NATV1CAP_i[i]) * NATCAPSTOK[i] : 0.0
        DPRC_id[i,d]     = NATDPRC[i]
        RNORMAL_id[i,d]  = NATRNORMAL[i]
        GROTREND_id[i,d] = NATGROTREND[i]
        QRATIO_id[i,d]   = NATQRATIO[i]
        ALPHA_id[i,d]    = NATALPHA[i]
        RORADJ_id[i,d]   = NATRORADJ[i]
        GRETEXP_id[i,d]  = NATGRETEXP[i]
    end
    for o in 1:NO, d in 1:NR
        EMPRAT_od[o,d]   = NATEMPRAT
        ELASTWAGE_od[o,d] = NATELASTWAGE
    end

    # ── TRADE_CS for aggregation (332-334) ───────────────────────────────────
    TRADE_CS_rd = zeros(T, NR, NR)
    for r in 1:NR, d in 1:NR
        TRADE_CS_rd[r,d] = sum(TRADE[c,s,r,d] for c in 1:NC, s in 1:NS)
    end

    # ── SIGMAs from elast ───────────────────────────────────────────────────
    SIGMA1LAB  = elast["SLAB"]
    SIGMA1PRIM = elast["P028"]
    ARMSIGMA   = elast["P015"]
    SIGMAOUT   = elast["SCET"]

    # ── Read REGSUPP values ─────────────────────────────────────────────────
    SIGMADOMDOM = regsupp["SGDD"]
    SIGMAMAR    = regsupp["SMAR"]
    POP         = regsupp["PO01"]

    # ── FRISCH ────────────────────────────────────────────────────────────────
    NATFRISCH = elast["P021"] isa NamedArray ? vec(parent(elast["P021"]))[1] : Float64(elast["P021"])
    FRISCH_d = zeros(T, NR)
    for d in 1:NR
        FRISCH_d[d] = -abs(NATFRISCH)
    end

    # ── NATEPS / EPS with scaling ──────────────────────────────────────────
    EPS_cd = zeros(T, NC, NR)
    for c in 1:NC, d in 1:NR
        EPS_cd[c,d] = NATEPS[c]
    end

    # ── Weights for aggregation (248-277) ───────────────────────────────────
    NATV1PRIM_i = zeros(T, NI)
    NATV1LAB_O_i = zeros(T, NI)
    REGLAB_od = zeros(T, NO, NR)
    NATV1CAP_i_out = zeros(T, NI)
    SALE_c = zeros(T, NC)
    MARUSE_m = zeros(T, NM)
    NATEXPORT_c = zeros(T, NC)
    V3PUR_S_cd = zeros(T, NC, NR)  # placeholder
    NATIMP_c = [sum(IMPORT_re[c,r] for r in 1:NR) for c in 1:NC]

    for i in 1:NI
        NATV1PRIM_i[i] = sum(VPRIM_id[i,d] for d in 1:NR)
        NATV1LAB_O_i[i] = sum(VLAB_O_id[i,d] for d in 1:NR)
        NATV1CAP_i_out[i] = sum(V1CAP[i,d] for d in 1:NR)
    end
    for o in 1:NO, d in 1:NR
        REGLAB_od[o,d] = sum(V1LAB[i,o,d] for i in 1:NI)
    end
    for c in 1:NC
        SALE_c[c] = sum(TRADE[c,s,r,d] for s in 1:NS, r in 1:NR, d in 1:NR)
    end
    for m in 1:NM
        MARUSE_m[m] = sum(SUPPMAR[m,r,d,p] for r in 1:NR, d in 1:NR, p in 1:NR)
    end
    for c in 1:NC, d in 1:NR
        NATEXPORT_c[c] = 0.0  # placeholder — need USE data
    end

    NATCOSTS_i = zeros(T, NI)
    for i in 1:NI
        tot = 0.0
        for d in 1:NR
            tot += VPRIM_id[i,d]
            if USE_data !== nothing && TAXES_data !== nothing
                for c in 1:NC, s in 1:NS
                    tot += USE_data[c,s,i,d] + TAXES_data[c,s,i,d]
                end
            end
        end
        NATCOSTS_i[i] = tot
    end

    weights = Dict{String,Any}(
        "3PUR" => NamedArray(V3PUR_S_cd, Tuple([COM,REG]), (:COM,:DST)),
        "1PRM" => NamedArray(NATV1PRIM_i, IND, (:IND,)),
        "LABR" => NamedArray(NATV1LAB_O_i, IND, (:IND,)),
        "NIMP" => NamedArray(NATIMP_c, COM, (:COM,)),
        "NCAP" => NamedArray(NATV1CAP_i_out, IND, (:IND,)),
        "3TOT" => NamedArray(V3TOT_d, REG, (:DST,)),
        "WSGD" => NamedArray(SALE_c, COM, (:COM,)),
        "MRUS" => NamedArray(MARUSE_m, MAR, (:MAR,)),
        "1TOT" => NamedArray(NATCOSTS_i, IND, (:IND,)),
        "4TOT" => NamedArray(NATEXPORT_c, COM, (:COM,)),
        "RLAB" => NamedArray(REGLAB_od, Tuple([OCC,REG]), (:OCC,:DST)),
        "STOC" => NamedArray(CAPSTOK_id, Tuple([IND,REG]), (:IND,:DST)),
        "ATRD" => NamedArray(TRADE_CS_rd, Tuple([REG,REG]), (:ORG,:DST)),
    )

    # ── Output ──────────────────────────────────────────────────────────────
    premod_out = Dict{String,Any}(
        "MAKE" => MAKE_r,
        "TRAD" => TRADE,
        "MARS" => SUPPMAR,
        "TMAR" => TRADMAR,
        "DIST" => DISTANCE,
        "1LAB" => V1LAB,
        "1CAP" => V1CAP,
        "1LND" => V1LND,
        "SLAB" => SIGMA1LAB,
        "P028" => SIGMA1PRIM,
        "SGDD" => SIGMADOMDOM,
        "SMAR" => SIGMAMAR,
        "PO01" => POP,
        "P021" => NamedArray(FRISCH_d, REG, (:DST,)),
        "XPEL" => NamedArray(EPS_cd, Tuple([COM,REG]), (:COM,:DST)),
        "SCET" => SIGMAOUT,
        "P018" => elast["P018"],
        "DPRC" => NamedArray(DPRC_id, Tuple([IND,REG]), (:IND,:DST)),
        "TARG" => NamedArray(RNORMAL_id, Tuple([IND,REG]), (:IND,:DST)),
        "TFRO" => NamedArray(GROTREND_id, Tuple([IND,REG]), (:IND,:DST)),
        "QRAT" => NamedArray(QRATIO_id, Tuple([IND,REG]), (:IND,:DST)),
        "ALFA" => NamedArray(ALPHA_id, Tuple([IND,REG]), (:IND,:DST)),
        "RADJ" => NamedArray(RORADJ_id, Tuple([IND,REG]), (:IND,:DST)),
        "REXP" => NamedArray(GRETEXP_id, Tuple([IND,REG]), (:IND,:DST)),
        "EMPR" => NamedArray(EMPRAT_od, Tuple([OCC,REG]), (:OCC,:DST)),
        "ELWG" => NamedArray(ELASTWAGE_od, Tuple([OCC,REG]), (:OCC,:DST)),
        "STOC" => NamedArray(CAPSTOK_id, Tuple([IND,REG]), (:IND,:DST)),
        # Add LCOM for top-down (pass from nat data if available)
        "LCOM" => haskey(elast, "LCOM") ? elast["LCOM"] : nothing,
        # Pass-through for Step 4 derived parameters
        "BSMR" => haskey(pstras, "BSMR") ? pstras["BSMR"] : nothing,
        "UTAX" => haskey(pstras, "UTAX") ? pstras["UTAX"] : nothing,
        "2PUR" => haskey(pstras, "2PUR") ? pstras["2PUR"] : nothing,
        "STOK" => haskey(pstras, "STOK") ? pstras["STOK"] : nothing,
        "1PTX" => haskey(pstras, "1PTX") ? pstras["1PTX"] : nothing,
    )

    PremodResult(premod=premod_out, weights=weights, diag=Dict{String,Any}())
end
