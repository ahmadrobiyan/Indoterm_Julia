"""
Translation of `reg1.tab` — use regional shares (R001–R005) to split user
columns according to destination region.

Input: `Reg0Result.reg0` dict, `read_regsupp_data()`, `read_distgone_data()`.
Output: `Reg1Result` — struct with reg1 and diag dicts.
"""

using NamedArrays

Base.@kwdef struct Reg1Result
    reg1::Dict{String,Any}
    diag::Dict{String,Any}
end

function build_reg1!(reg0::Dict{String,Any},
                     regsupp::Dict{String,Any},
                     distgone::Dict{String,Any})

    T = Float64
    R001 = parent(reg0["R001"]); R002 = parent(reg0["R002"])
    R003 = parent(reg0["R003"]); R004 = parent(reg0["R004"]); R005 = parent(reg0["R005"])
    NATFAC   = parent(reg0["FACT"])
    NATBASIC = parent(reg0["UBAS"])
    NATMARGINS = parent(reg0["UMAR"])
    NATTAXES = parent(reg0["UTAX"])
    NATMAKE  = parent(reg0["MAKE"])
    NATSTOCKS = parent(reg0["STOK"])
    INVSHR   = parent(reg0["ISHR"])

    NC = length(COM); NI = length(IND); NS = 2; NM = length(MAR); NR = length(REG); NO = length(OCC)
    FINDEM_NAMES = ["Hou","Inv","Gov","Exp"]
    USR_NAMES = vcat(IND, FINDEM_NAMES)
    NUSR = length(USR_NAMES)
    FACTOR_NAMES = ["Labour","Capital","Land","ProdTax"]
    NFAC = length(FACTOR_NAMES)

    # ── NATIONAL aggregates (reg1.tab 55-66) ────────────────────────────────
    NATPUR = zeros(T, NC, NS, NUSR)
    for c in 1:NC, s in 1:NS, u in 1:NUSR
        NATPUR[c,s,u] = NATBASIC[c,s,u] + NATTAXES[c,s,u] + sum(NATMARGINS[c,s,u,m] for m in 1:NM)
    end

    INVPUR_I = zeros(T, NC)
    for c in 1:NC
        INVPUR_I[c] = sum(NATPUR[c,s, findfirst(==("Inv"), USR_NAMES)] for s in 1:NS)
    end
    NATINVEST = zeros(T, NC, NI)
    for c in 1:NC, i in 1:NI
        NATINVEST[c,i] = INVSHR[c,i] * INVPUR_I[c]
    end

    # ── Handle zero-share rows (86-103) ─────────────────────────────────────
    DEFREGSHR = 1.0 / NR

    for i in 1:NI
        if sum(R001[i,:]) == 0; for r in 1:NR; R001[i,r] = DEFREGSHR; end; end
        if sum(R002[i,:]) == 0; for r in 1:NR; R002[i,r] = DEFREGSHR; end; end
    end
    for c in 1:NC
        if sum(R003[c,:]) == 0; for r in 1:NR; R003[c,r] = DEFREGSHR; end; end
        if sum(R004[c,:]) == 0; for r in 1:NR; R004[c,r] = DEFREGSHR; end; end
        if sum(R005[c,:]) == 0; for r in 1:NR; R005[c,r] = DEFREGSHR; end; end
    end

    # Normalize shares (127-136)
    for i in 1:NI
        t = sum(R001[i,:]); if t > 0; R001[i,:] ./= t; end
        t = sum(R002[i,:]); if t > 0; R002[i,:] ./= t; end
    end
    for c in 1:NC
        t = sum(R003[c,:]); if t > 0; R003[c,:] ./= t; end
        t = sum(R004[c,:]); if t > 0; R004[c,:] ./= t; end
        t = sum(R005[c,:]); if t > 0; R005[c,:] ./= t; end
    end

    @assert all(i -> abs(1 - sum(R001[i,:])) < 1e-4, 1:NI) "R001 sum != 1"
    @assert all(i -> abs(1 - sum(R002[i,:])) < 1e-4, 1:NI) "R002 sum != 1"
    @assert all(c -> abs(1 - sum(R003[c,:])) < 1e-4, 1:NC) "R003 sum != 1"
    @assert all(c -> abs(1 - sum(R004[c,:])) < 1e-4, 1:NC) "R004 sum != 1"
    @assert all(c -> abs(1 - sum(R005[c,:])) < 1e-4, 1:NC) "R005 sum != 1"

    # ── USHR — user-by-region share matrix (138-167) ────────────────────────
    # USHR(c,s,u,d): COM × SRC × USR × DST
    USHR = zeros(T, NC, NS, NUSR, NR)
    for c in 1:NC, s in 1:NS, i in 1:NI, d in 1:NR
        USHR[c,s,i,d] = R001[i,d]
    end
    for c in 1:NC, s in 1:NS, d in 1:NR
        USHR[c,s, findfirst(==("Hou"), USR_NAMES), d] = R003[c,d]
        USHR[c,s, findfirst(==("Gov"), USR_NAMES), d] = R005[c,d]
        USHR[c,s, findfirst(==("Exp"), USR_NAMES), d] = R004[c,d]
    end

    # INVEST(c,i,d) — investment at purchasers prices
    INVEST = zeros(T, NC, NI, NR)
    for c in 1:NC, i in 1:NI, d in 1:NR
        INVEST[c,i,d] = R002[i,d] * NATINVEST[c,i]
    end
    INVEST_I = zeros(T, NC, NR)
    for c in 1:NC, d in 1:NR
        INVEST_I[c,d] = sum(INVEST[c,i,d] for i in 1:NI)
    end
    INVEST_ID = zeros(T, NC)
    for c in 1:NC
        INVEST_ID[c] = sum(INVEST_I[c,d] for d in 1:NR)
    end
    for c in 1:NC, s in 1:NS, d in 1:NR
        if INVEST_ID[c] > 0
            USHR[c,s, findfirst(==("Inv"), USR_NAMES), d] = INVEST_I[c,d] / INVEST_ID[c]
        else
            USHR[c,s, findfirst(==("Inv"), USR_NAMES), d] = DEFREGSHR
        end
    end

    # ── MARWGT margin weighting — destination side (170-175) ────────────────
    MARWGT_dm = zeros(T, NR, NM)
    mwgt_raw = regsupp["MWGT"]
    for d in 1:NR, m in 1:NM
        MARWGT_dm[d,m] = max(0.01, mwgt_raw[d,m])
    end

    # Adjust USHR for margin commodities by MARWGT (179-180)
    for m in 1:NM, s in 1:NS, u in 1:NUSR, d in 1:NR
        USHR[m,s,u,d] *= MARWGT_dm[d,m]
    end

    # Scale USHR to sum to 1 (182-192)
    USHRTOT = zeros(T, NC, NS, NUSR)
    for c in 1:NC, s in 1:NS, u in 1:NUSR
        USHRTOT[c,s,u] = sum(USHR[c,s,u,d] for d in 1:NR)
    end
    for c in 1:NC, s in 1:NS, u in 1:NUSR, d in 1:NR
        if USHRTOT[c,s,u] > 0
            USHR[c,s,u,d] /= USHRTOT[c,s,u]
        end
    end
    # Recompute and assert
    for c in 1:NC, s in 1:NS, u in 1:NUSR
        tot = sum(USHR[c,s,u,d] for d in 1:NR)
        @assert abs(1 - tot) < 1e-5 "USHR sum != 1: c=$c s=$s u=$u tot=$tot"
    end

    # ── Read DISTGONE from DISTGONE.HAR (196-199) ───────────────────────────
    DGON = parent(distgone["DGON"])  # COM × SRC × DST

    # ── DMAR: distance-relevant margin commodities ──────────────────────────
    DMAR_idx = collect(1:NM)

    # ── Split national flows to regions (210-258) ───────────────────────────
    BASIC_r  = zeros(T, NC, NS, NUSR, NR)
    TAXES_r  = zeros(T, NC, NS, NUSR, NR)
    MARGINS_r = zeros(T, NC, NS, NUSR, NM, NR)
    for c in 1:NC, s in 1:NS, u in 1:NUSR, d in 1:NR
        BASIC_r[c,s,u,d] = NATBASIC[c,s,u] * USHR[c,s,u,d]
        TAXES_r[c,s,u,d] = NATTAXES[c,s,u] * USHR[c,s,u,d]
        for m in 1:NM
            MARGINS_r[c,s,u,m,d] = NATMARGINS[c,s,u,m] * USHR[c,s,u,d]
        end
    end
    # Distance-related margins (DMAR) — scale by sqrt(DISTGONE) * MARWGT (217-222)
    for d_idx in DMAR_idx
        m = findfirst(==(MAR[d_idx]), MAR)
        if m === nothing; continue; end
        for c in 1:NC, s in 1:NS, u in 1:NUSR, d in 1:NR
            factor = MARWGT_dm[d,m] * sqrt(max(DGON[c,s,d], 0.0))
            MARGINS_r[c,s,u,m,d] *= factor
        end
    end
    # Rescale margins to national totals (219-222)
    MARGINS_D = zeros(T, NC, NS, NUSR, NM)
    for c in 1:NC, s in 1:NS, u in 1:NUSR, m in 1:NM
        MARGINS_D[c,s,u,m] = sum(MARGINS_r[c,s,u,m,d] for d in 1:NR)
    end
    for c in 1:NC, s in 1:NS, u in 1:NUSR, m in 1:NM, d in 1:NR
        if MARGINS_D[c,s,u,m] > 0
            MARGINS_r[c,s,u,m,d] *= NATMARGINS[c,s,u,m] / MARGINS_D[c,s,u,m]
        end
    end

    # Stocks and factor costs (223-224)
    STOCKS_r = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        STOCKS_r[i,d] = R001[i,d] * NATSTOCKS[i]
    end
    FAC_r = zeros(T, NI, NFAC, NR)
    for i in 1:NI, g in 1:NFAC, d in 1:NR
        FAC_r[i,g,d] = R001[i,d] * NATFAC[i,g]
    end

    # ── Check BADTAX (241-258) ──────────────────────────────────────────────
    for c in 1:NC, s in 1:NS, u in 1:NUSR, d in 1:NR
        if BASIC_r[c,s,u,d] <= 0
            @assert abs(TAXES_r[c,s,u,d]) < 1e-8
            TAXES_r[c,s,u,d] = 0.0
        end
    end

    # ── MAKE split (261-278) ────────────────────────────────────────────────
    COSTS_rd = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        val = sum(FAC_r[i,g,d] for g in 1:NFAC)
        for c in 1:NC, s in 1:NS
            val += BASIC_r[c,s,i,d] + TAXES_r[c,s,i,d] + sum(MARGINS_r[c,s,i,m,d] for m in 1:NM)
        end
        COSTS_rd[i,d] = val
    end
    NETPUT_rd = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        NETPUT_rd[i,d] = COSTS_rd[i,d] - STOCKS_r[i,d]
    end
    R001A = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        tot = sum(NETPUT_rd[i,:])
        R001A[i,d] = tot > 0 ? NETPUT_rd[i,d] / tot : DEFREGSHR
    end

    MAKE_r = zeros(T, NC, NI, NR)
    for c in 1:NC, i in 1:NI, d in 1:NR
        MAKE_r[c,i,d] = R001A[i,d] * NATMAKE[c,i]
    end

    # ── HEAVYEXP / LOCAL exports — override USHR for specific commodities (281-295) ──
    HEAVYEXP_COMS = ["BauxiteOre","PetrolNatGas","CopperOre","NonFerrMetal","Cocoa","NickelOre","CoalLignite"]
    for hc in HEAVYEXP_COMS
        ci = findfirst(==(hc), COM)
        if ci !== nothing
            for d in 1:NR
                USHR[ci, 1, findfirst(==("Exp"), USR_NAMES), d] = R001A[ci, d]
            end
            # rescale to 1
            tot = sum(USHR[ci, 1, findfirst(==("Exp"), USR_NAMES), :])
            if tot > 0
                USHR[ci, 1, findfirst(==("Exp"), USR_NAMES), :] ./= tot
            end
            for d in 1:NR
                BASIC_r[ci, 1, findfirst(==("Exp"), USR_NAMES), d] = NATBASIC[ci, 1, findfirst(==("Exp"), USR_NAMES)] * USHR[ci, 1, findfirst(==("Exp"), USR_NAMES), d]
            end
        end
    end

    # ── BASIC_U — sum over users (297-304) ──────────────────────────────────
    BASIC_U = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        BASIC_U[c,s,d] = sum(BASIC_r[c,s,u,d] for u in 1:NUSR)
    end

    # ── IMPORTS — MSHR from regsupp (306-332) ───────────────────────────────
    IMPUSE_cd = zeros(T, NC, NR)
    for c in 1:NC, d in 1:NR
        IMPUSE_cd[c,d] = BASIC_U[c,2,d]
    end
    NATIMP_c = [sum(IMPUSE_cd[c,d] for d in 1:NR) for c in 1:NC]

    # Default MSHR from use shares
    DEFAULTMSHR = zeros(T, NC, NR)
    for c in 1:NC, r in 1:NR
        DEFAULTMSHR[c,r] = NATIMP_c[c] > 0 ? IMPUSE_cd[c,r] / NATIMP_c[c] : 0.0
    end

    MSHR_cr = zeros(T, NC, NR)
    if haskey(regsupp, "MSHR")
        mshr_raw = regsupp["MSHR"]
        for c in 1:min(NC, size(mshr_raw,1)), r in 1:min(NR, size(mshr_raw,2))
            MSHR_cr[c,r] = mshr_raw[c,r]
        end
    end
    # If whole row zero, use default
    for c in 1:NC
        if sum(MSHR_cr[c,:]) == 0
            for r in 1:NR; MSHR_cr[c,r] = DEFAULTMSHR[c,r]; end
        end
        t = sum(MSHR_cr[c,:])
        if t > 0; MSHR_cr[c,:] ./= t; end
    end

    IMPORT_cr = zeros(T, NC, NR)
    for c in 1:NC, r in 1:NR
        IMPORT_cr[c,r] = MSHR_cr[c,r] * NATIMP_c[c]
    end

    # ── ALLDEMAND / ALLSUPPLY / DIFFCOM (334-352) ──────────────────────────
    MARDEMAND_md = zeros(T, NM, NR)
    for m in 1:NM, d in 1:NR
        MARDEMAND_md[m,d] = sum(MARGINS_r[c,s,u,m,d] for c in 1:NC, s in 1:NS, u in 1:NUSR)
    end

    ALLDEMAND = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        ALLDEMAND[c,s,d] = BASIC_U[c,s,d]
    end
    for m in 1:NM, d in 1:NR
        if m <= NC; ALLDEMAND[m,1,d] += MARDEMAND_md[m,d]; end
    end

    MAKE_I_cd = zeros(T, NC, NR)
    for c in 1:NC, d in 1:NR
        MAKE_I_cd[c,d] = sum(MAKE_r[c,i,d] for i in 1:NI)
    end

    ALLSUPPLY = zeros(T, NC, NS, NR)
    for c in 1:NC, d in 1:NR
        ALLSUPPLY[c,1,d] = MAKE_I_cd[c,d]
        ALLSUPPLY[c,2,d] = IMPORT_cr[c,d]
    end

    DIFFCOM_rsd = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        DIFFCOM_rsd[c,s,d] = ALLSUPPLY[c,s,d] - ALLDEMAND[c,s,d]
    end

    # ── INTDEMAND / FINDEMAND / DEMRATIO (395-415) ──────────────────────────
    INTDEMAND = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        INTDEMAND[c,s,d] = sum(BASIC_r[c,s,i,d] for i in 1:NI)
    end
    FINDEMAND = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        FINDEMAND[c,s,d] = sum(BASIC_r[c,s,u,d] for u in NI+1:NUSR)
    end

    DEMRATIO_csd = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        all_sup = ALLSUPPLY[c,s,d]
        all_dem = ALLDEMAND[c,s,d]
        DEMRATIO_csd[c,s,d] = (0.01 + all_dem) / (0.01 + all_sup)
    end

    # ── MARGINS_U — sum MARGINS over users (436-444) ────────────────────────
    MARGINS_U = zeros(T, NC, NS, NM, NR)
    for c in 1:NC, s in 1:NS, m in 1:NM, d in 1:NR
        MARGINS_U[c,s,m,d] = sum(MARGINS_r[c,s,u,m,d] for u in 1:NUSR)
    end

    USE_rsud = zeros(T, NC, NS, NUSR, NR)
    for c in 1:NC, s in 1:NS, u in 1:NUSR, d in 1:NR
        USE_rsud[c,s,u,d] = BASIC_r[c,s,u,d] + sum(MARGINS_r[c,s,u,m,d] for m in 1:NM)
    end

    # ── ASSEMBLE OUTPUT ─────────────────────────────────────────────────────
    reg1_out = Dict{String,Any}(
        "FACT" => NamedArray(FAC_r, Tuple([IND, FACTOR_NAMES, REG]), (:IND,:FACTOR,:DST)),
        "UTAX" => NamedArray(TAXES_r, Tuple([COM,SRC,USR_NAMES,REG]), (:COM,:SRC,:USR,:DST)),
        "MAKE" => NamedArray(MAKE_r, Tuple([COM,IND,REG]), (:COM,:IND,:DST)),
        "STOK" => NamedArray(STOCKS_r, Tuple([IND,REG]), (:IND,:DST)),
        "BSCU" => NamedArray(BASIC_U, Tuple([COM,SRC,REG]), (:COM,:SRC,:DST)),
        "MRGU" => NamedArray(MARGINS_U, Tuple([COM,SRC,MAR,REG]), (:COM,:SRC,:MAR,:DST)),
        "BSMR" => NamedArray(USE_rsud, Tuple([COM,SRC,USR_NAMES,REG]), (:COM,:SRC,:USR,:DST)),
        "IMPS" => NamedArray(IMPORT_cr, Tuple([COM,REG]), (:COM,:ORG)),
        "ADEM" => NamedArray(ALLDEMAND, Tuple([COM,SRC,REG]), (:COM,:SRC,:DST)),
        "ASUP" => NamedArray(ALLSUPPLY, Tuple([COM,SRC,REG]), (:COM,:SRC,:ORG)),
        "INVS" => NamedArray(INVEST, Tuple([COM,IND,REG]), (:COM,:IND,:DST)),
        "2PUR" => NamedArray(INVEST, Tuple([COM,IND,REG]), (:COM,:IND,:DST)),
        "OSHR" => reg0["OSHR"],
    )

    diag_out = Dict{String,Any}(
        "USHR" => USHR,
        "COST" => COSTS_rd,
        "MAKC" => [sum(MAKE_r[c,i,d] for c in 1:NC) for i in 1:NI, d in 1:NR],
        "ADEM" => ALLDEMAND,
        "ASUP" => ALLSUPPLY,
        "DCOM" => DIFFCOM_rsd,
        "DRAT" => DEMRATIO_csd,
        "MARD" => MARDEMAND_md,
    )

    Reg1Result(reg1=reg1_out, diag=diag_out)
end
