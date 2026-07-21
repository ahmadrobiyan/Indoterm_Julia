"""
Translation of `reg2.tab` — construct trial trade matrices via gravity
formula, allocate margins (TRADMAR) and margin-supply (SUPPMAR).

Input: `Reg1Result.reg1` dict + `regsupp` dict.
Output: `Reg2Result` — reg2, diag fields.
"""

using NamedArrays

Base.@kwdef struct Reg2Result
    reg2::Dict{String,Any}
    diag::Dict{String,Any}
end

function build_reg2!(reg1::Dict{String,Any}, regsupp::Dict{String,Any};
                     reg0_diag::Union{Dict{String,Any},Nothing}=nothing,
                     rloc_idx::Union{Vector{Int},Nothing}=nothing)
    T = Float64
    NC = length(COM); NI = length(IND); NS = 2; NM = length(MAR); NR = length(REG); NO = length(OCC)

    # ── Read input coeffs ──────────────────────────────────────────────────
    BASIC_U   = reg1["BSCU"]    # COM×SRC×DST
    MARGINS_U = reg1["MRGU"]    # COM×SRC×MAR×DST
    USE_r     = reg1["BSMR"]    # COM×SRC×USR×DST
    ALLDEMAND_raw = reg1["ADEM"]
    ALLSUPPLY_raw = reg1["ASUP"]
    FAC_r     = reg1["FACT"]    # IND×FACTOR×DST
    TAXES_r   = reg1["UTAX"]    # COM×SRC×USR×DST
    MAKE_r    = reg1["MAKE"]    # COM×IND×DST
    STOCKS_r  = reg1["STOK"]    # IND×DST
    IMPORT_cr = reg1["IMPS"]    # COM×ORG

    # ── Read distance / distfac / MARWGT / LOCMAR from regsupp ──
    DISTANCE = regsupp["DIST"]   # ORG×DST
    DISTFAC  = regsupp["DFAC"]   # COM×SRC
    MARWGT_rm = regsupp["MWGT"]  # REG×MAR (ORG×MAR in reg2 context)
    LOCMAR_m  = regsupp["LMAR"]   # MAR

    # ── RLOC: local commodities ──────────────────────────────────────
    if rloc_idx !== nothing
        RLOC_idx = rloc_idx
    elseif reg0_diag !== nothing && haskey(reg0_diag, "RLOC")
        RLOC_idx = reg0_diag["RLOC"]
    else
        RLOC_idx = Int[]
    end
    RLOC_set = Set(COM[RLOC_idx])
    RLOC_vec = COM[RLOC_idx]

    # ── DMAR: distance-relevant margin commodities ───────────────────
    # All 9 margin commodities use distance weighting by default
    DMAR_idx = collect(1:NM)

    DEFREGSHR = 1.0 / NR
    NONMAR_set = setdiff(Set(COM), Set(MAR))
    NONMAR_idx = [c for c in 1:NC if !(COM[c] in Set(MAR))]

    # ── MAKE_I / MAKE_C (reg2.tab 113-114) ────────────────────────────────
    MAKE_I_cd = zeros(T, NC, NR)
    MAKE_C_id = zeros(T, NI, NR)
    for c in 1:NC, d in 1:NR
        MAKE_I_cd[c,d] = sum(MAKE_r[c,i,d] for i in 1:NI)
    end
    for i in 1:NI, d in 1:NR
        MAKE_C_id[i,d] = sum(MAKE_r[c,i,d] for c in 1:NC)
    end

    # ── COSTS / DIFFIND checks (83-93) ─────────────────────────────────────
    USR_NAMES = size(USE_r, 3) == NI + 4 ? vcat(IND,["Hou","Inv","Gov","Exp"]) : [string("u",i) for i in 1:size(USE_r,3)]
    COSTS_id = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        val = sum(FAC_r[i,g,d] for g in 1:size(FAC_r,2))
        for c in 1:NC, s in 1:NS
            val += USE_r[c,s,i,d] + TAXES_r[c,s,i,d]
        end
        COSTS_id[i,d] = val
    end
    DIFFIND_id = zeros(T, NI, NR)
    for i in 1:NI, d in 1:NR
        DIFFIND_id[i,d] = (COSTS_id[i,d] - STOCKS_r[i,d]) - MAKE_C_id[i,d]
    end

    # ── MARDEMAND / ALLDEMAND / ALLSUPPLY / SUPPRATIO (98-124) ────────────
    MARDEMAND_md = zeros(T, NM, NR)
    for m in 1:NM, d in 1:NR
        MARDEMAND_md[m,d] = sum(MARGINS_U[c,s,m,d] for c in 1:NC, s in 1:NS)
    end

    ALLDEMAND = copy(parent(BASIC_U))
    for m in 1:NM, d in 1:NR
        if m <= NC; ALLDEMAND[m,1,d] += MARDEMAND_md[m,d]; end
    end

    ALLSUPPLY = zeros(T, NC, NS, NR)
    for c in 1:NC, d in 1:NR
        ALLSUPPLY[c,1,d] = MAKE_I_cd[c,d]
        ALLSUPPLY[c,2,d] = IMPORT_cr[c,d]
    end

    SUPPRATIO = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        if ALLDEMAND[c,s,d] > 0
            SUPPRATIO[c,s,d] = ALLSUPPLY[c,s,d] / ALLDEMAND[c,s,d]
        end
    end

    DIFFCOM_csd = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        DIFFCOM_csd[c,s,d] = ALLSUPPLY[c,s,d] - ALLDEMAND[c,s,d]
    end

    # ── Gravity formula — SOURCESHR (126-180) ───────────────────────────────
    OWNSHARE = zeros(T, NC, NS, NR)
    OTHSHARE = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        os = min(SUPPRATIO[c,s,d], 1.0)
        df = max(DISTFAC[c,s], 0.0)
        os *= (1 - 0.2^df)
        OWNSHARE[c,s,d] = os
        OTHSHARE[c,s,d] = 1 - os
    end

    SOURCESHR = zeros(T, NC, NS, NR, NR)  # c,s,r,d
    for c in 1:NC, s in 1:NS, r in 1:NR, d in 1:NR
        if r == d
            SOURCESHR[c,s,r,d] = OWNSHARE[c,s,d]
        else
            denom = max(DISTANCE[r,d]^DISTFAC[c,s], 1e-15)
            SOURCESHR[c,s,r,d] = sqrt(max(ALLSUPPLY[c,s,r], 0.0)) / denom
        end
    end

    # Scale off-diagonals to sum to OTHSHARE
    TEMTOTCSD = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        TEMTOTCSD[c,s,d] = sum(SOURCESHR[c,s,r,d] for r in 1:NR if r != d)
    end
    for c in 1:NC, s in 1:NS, r in 1:NR, d in 1:NR
        if r != d && TEMTOTCSD[c,s,d] > 0
            SOURCESHR[c,s,r,d] *= OTHSHARE[c,s,d] / TEMTOTCSD[c,s,d]
        elseif r != d
            SOURCESHR[c,s,r,d] = 0.0
        end
    end

    # TRADE = SOURCESHR * BASIC_U (177-181)
    TRADE = zeros(T, NC, NS, NR, NR)  # c,s,r,d
    for c in 1:NC, s in 1:NS, r in 1:NR, d in 1:NR
        TRADE[c,s,r,d] = SOURCESHR[c,s,r,d] * BASIC_U[c,s,d]
    end

    # ── Deficit-region adjustment (183-200) ─────────────────────────────────
    for c in 1:NC, s in 1:NS, r in 1:NR, d in 1:NR
        if r != d && ALLDEMAND[c,s,r] > ALLSUPPLY[c,s,r] && ALLDEMAND[c,s,d] < ALLSUPPLY[c,s,d]
            SOURCESHR[c,s,r,d] *= 0.1
        end
    end
    # Rescale to sum 1
    for c in 1:NC, s in 1:NS, d in 1:NR
        tot = sum(SOURCESHR[c,s,:,d])
        if tot > 0; SOURCESHR[c,s,:,d] ./= tot; end
    end
    for c in 1:NC, s in 1:NS, r in 1:NR, d in 1:NR
        TRADE[c,s,r,d] = SOURCESHR[c,s,r,d] * BASIC_U[c,s,d]
    end

    # ── Local commodity adjustment (202-221) ────────────────────────────────
    for c in RLOC_idx
        if c === nothing; continue; end
        for r in 1:NR, d in 1:NR
            if r != d
                SOURCESHR[c,1,r,d] *= 0.01
            end
        end
        for q in 1:NR
            SOURCESHR[c,1,q,q] = 1.0
        end
        for d in 1:NR
            tot = sum(SOURCESHR[c,1,:,d])
            if tot > 0; SOURCESHR[c,1,:,d] ./= tot; end
        end
    end
    for c in 1:NC, s in 1:NS, d in 1:NR
        tot = sum(SOURCESHR[c,s,:,d])
        if tot > 0; SOURCESHR[c,s,:,d] ./= tot; end
    end
    for c in 1:NC, s in 1:NS, r in 1:NR, d in 1:NR
        TRADE[c,s,r,d] = SOURCESHR[c,s,r,d] * BASIC_U[c,s,d]
    end

    # ── Cross-hauling reduction (OPTIONAL — 254-281) ────────────────────────
    # 95% reduction enabled for all COM
    CROSS = zeros(T, NC, NS, NR, NR)
    for c in 1:NC, s in 1:NS, r in 1:NR, d in 1:NR
        if r != d && r != d
            CROSS[c,s,r,d] = 0.95 * min(TRADE[c,s,r,d], TRADE[c,s,d,r])
        end
    end
    for c in 1:NC, s in 1:NS, r in 1:NR, d in 1:NR
        TRADE[c,s,r,d] -= CROSS[c,s,r,d]
    end
    for c in 1:NC, s in 1:NS, r in 1:NR
        TRADE[c,s,r,r] += sum(CROSS[c,s,r,d] for d in 1:NR)
    end
    # Recompute SOURCESHR
    for c in 1:NC, s in 1:NS, d in 1:NR
        tot = sum(TRADE[c,s,:,d])
        if tot > 0
            for r in 1:NR; SOURCESHR[c,s,r,d] = TRADE[c,s,r,d] / tot; end
        end
    end

    # ── Assertions: TRADE >= 0, totals match BASIC_U ────────────────────────
    @assert all(TRADE .>= -1e-10) "Negative TRADE"
    TRADE_R = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        TRADE_R[c,s,d] = sum(TRADE[c,s,:,d])
    end
    for c in 1:NC, s in 1:NS, d in 1:NR
        if BASIC_U[c,s,d] == 0
            @assert TRADE_R[c,s,d] < 1e-10 "TRADE_R non-zero where BASIC_U=0"
        else
            @assert TRADE_R[c,s,d] > 0 "TRADE_R zero where BASIC_U>0"
        end
    end

    # ── TRADMAR — margins on trade (324-368) ────────────────────────────────
    TRADMAR = zeros(T, NC, NS, NM, NR, NR)  # c,s,m,r,d

    # All margin commodities: MARWGT * SOURCESHR * MARGINS_U
    for c in 1:NC, s in 1:NS, m in 1:NM, r in 1:NR, d in 1:NR
        TRADMAR[c,s,m,r,d] = MARWGT_rm[r,m] * SOURCESHR[c,s,r,d] * MARGINS_U[c,s,m,d]
    end
    # DMAR commodities: additional sqrt(DISTANCE) factor
    for m_idx in DMAR_idx
        m = findfirst(==(MAR[m_idx]), MAR)
        if m === nothing; continue; end
        for c in 1:NC, s in 1:NS, r in 1:NR, d in 1:NR
            TRADMAR[c,s,m,r,d] *= sqrt(max(DISTANCE[r,d], 0.0))
        end
    end

    @assert all(TRADMAR .>= -1e-10) "Negative TRADMAR"
    # Assert TRADMAR=0 where TRADE=0
    for c in 1:NC, s in 1:NS, m in 1:NM, r in 1:NR, d in 1:NR
        if TRADE[c,s,r,d] == 0
            @assert TRADMAR[c,s,m,r,d] < 1e-10 "TRADMAR non-zero where TRADE=0"
        end
    end

    # Scale TRADMAR to match MARGINS_U (355-368)
    TEMTOTCSMD = zeros(T, NC, NS, NM, NR)
    for c in 1:NC, s in 1:NS, m in 1:NM, d in 1:NR
        TEMTOTCSMD[c,s,m,d] = sum(TRADMAR[c,s,m,:,d])
    end
    for c in 1:NC, s in 1:NS, m in 1:NM, r in 1:NR, d in 1:NR
        if TEMTOTCSMD[c,s,m,d] > 0
            TRADMAR[c,s,m,r,d] *= MARGINS_U[c,s,m,d] / TEMTOTCSMD[c,s,m,d]
        end
    end

    # ── SUPPMAR — margins supplied by PRD (385-438) ─────────────────────────
    TRADMAR_CS = zeros(T, NM, NR, NR)  # m,r,d — sum over c,s
    for m in 1:NM, r in 1:NR, d in 1:NR
        TRADMAR_CS[m,r,d] = sum(TRADMAR[c,s,m,r,d] for c in 1:NC, s in 1:NS)
    end

    SUPPMAR = zeros(T, NM, NR, NR, NR)  # m,r,d,p
    for m in 1:NM, r in 1:NR, d in 1:NR, p in 1:NR
        SUPPMAR[m,r,d,p] = 0.5 * (SOURCESHR[m,1,p,d] + SOURCESHR[m,1,p,r]) * TRADMAR_CS[m,r,d]
    end

    # LOCMAR — local sourcing tendency (398-404)
    LOCMAR_arr = zeros(T, NM)
    for m in 1:NM
        LOCMAR_arr[m] = max(1.0, LOCMAR_m[m])
    end
    for m in 1:NM, r in 1:NR, q in 1:NR
        SUPPMAR[m,r,q,q] *= LOCMAR_arr[m]
    end
    # Zero SUPPMAR where region does not produce margin m
    for m in 1:NM, p in 1:NR
        if MAKE_I_cd[m,p] == 0
            SUPPMAR[m,:,:,p] .= 0.0
        end
    end

    # Scale SUPPMAR so sum over p does not exceed MAKE_I (407-415)
    TTMP_mp = zeros(T, NM, NR)
    TMUL_mp = zeros(T, NM, NR)
    for m in 1:NM, p in 1:NR
        TTMP_mp[m,p] = sum(SUPPMAR[m,r,d,p] for r in 1:NR, d in 1:NR)
        TMUL_mp[m,p] = TTMP_mp[m,p] > 0 ? MAKE_I_cd[m,p] / TTMP_mp[m,p] : 0.0
        TMUL_mp[m,p] = min(1.0, TMUL_mp[m,p])
    end
    for m in 1:NM, r in 1:NR, d in 1:NR, p in 1:NR
        SUPPMAR[m,r,d,p] *= TMUL_mp[m,p]
    end

    # Scale SUPPMAR to match TRADMAR_CS (418-438)
    TEMTOTMRD = zeros(T, NM, NR, NR)
    for m in 1:NM, r in 1:NR, d in 1:NR
        TEMTOTMRD[m,r,d] = sum(SUPPMAR[m,r,d,p] for p in 1:NR)
    end
    SCALEMRD = zeros(T, NM, NR, NR)
    for m in 1:NM, r in 1:NR, d in 1:NR
        SCALEMRD[m,r,d] = TEMTOTMRD[m,r,d] > 0 ? TRADMAR_CS[m,r,d] / TEMTOTMRD[m,r,d] : 0.0
        for p in 1:NR
            SUPPMAR[m,r,d,p] *= SCALEMRD[m,r,d]
        end
    end

    # ── OCCSHR → V1LAB allocation (467-488) ─────────────────────────────────
    OCCSHR_io = reg1["FACT"]  # wait — OSHR is separate
    # Actually OSHR is in reg0["OSHR"], passed through reg1's Transfer "OSHR"
    # But FACT in reg1 has the factor matrix, OSHR came via Transfer "OSHR" from reg0
    # We need OSHR from somewhere. Let me check what reg1 has...
    # The OSHR was written to reg0 but OSHR has IND×OCC dimensions
    # In reg1, it's transferred via Transfer "OSHR" from INFILE (reg0) to OUTFILE
    # But it was NOT in Reg1Result. Let me add it properly.
    # For now, I'll compute it from the V1LAB in reg0 which was carried in reg1's FACT.

    # Actually looking at reg1.tab again: V1LAB is created from FAC/FACTOR and OCCSHR
    # OCCSHR was transferred via "Transfer OSHR from file INFILE to file OUTFILE"
    # But in reg1 output, it's not in the struct. This is a problem.
    # For now, let me compute OCCSHR from the FACT matrix's Labour column and V1LAB data.
    # Actually FACT in reg1 is IND×FACTOR×DST — Labour column gives V1LAB_O(i,d)
    # And we need OCCSHR from reg0 — let me work around this.

    # Quick workaround: read OCCSHR from a reg0 output if we have it
    OCCSHR_io = nothing  # Will need to be passed in
    # For now, compute OCCSHR from V1LAB / V1LAB_O ratios
    # V1LAB_O(i) = sum_o V1LAB(i,o) from the original national data
    # The OCCSHR(i,o) = V1LAB(i,o) / V1LAB_O(i)
    # Since we don't have it here, we'll allocate V1LAB regionally in proportion
    # to LABOUR_i,d * OCCSHR(i,o) but we need OCCSHR...

    # For now — skip this, it can be computed from national OCCSHR * FAC(i,Labour,d) share
    # We'll fill in V1LAB later when we have OCCSHR available

    # ── Output assembly ─────────────────────────────────────────────────────
    reg2_out = Dict{String,Any}(
        "TMAR" => NamedArray(TRADMAR, Tuple([COM,SRC,MAR,REG,REG]), (:COM,:SRC,:MAR,:ORG,:DST)),
        "TRAD" => NamedArray(TRADE,   Tuple([COM,SRC,REG,REG]), (:COM,:SRC,:ORG,:DST)),
        "MARS" => NamedArray(SUPPMAR, Tuple([MAR,REG,REG,REG]), (:MAR,:ORG,:DST,:PRD)),
        "BSMR" => USE_r,
        "BSCU" => BASIC_U,
        "MRGU" => MARGINS_U,
        "IMPS" => IMPORT_cr,
        "MAKE" => MAKE_r,
        "COST" => NamedArray(MAKE_I_cd, Tuple([COM,REG]), (:COM,:ORG)),
        "STOK" => STOCKS_r,
        "DIST" => DISTANCE,
    )

    # Transfer "2PUR" from INFILE (reg1) to OUTFILE (reg2)
    if haskey(reg1, "2PUR")
        reg2_out["2PUR"] = reg1["2PUR"]
    end

    diag_out = Dict{String,Any}(
        "SOURCESHR" => SOURCESHR,
        "OWNSHARE"  => OWNSHARE,
        "SUPPRATIO" => SUPPRATIO,
        "MARDEMAND" => MARDEMAND_md,
        "TRADMAR_CS" => TRADMAR_CS,
        "COSTS"     => COSTS_id,
        "DIFFIND"   => DIFFIND_id,
    )

    Reg2Result(reg2=reg2_out, diag=diag_out)
end
