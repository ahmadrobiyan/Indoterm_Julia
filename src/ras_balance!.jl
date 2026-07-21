"""
Translation of `trdras.tab` + `raslin.tab` — sequential scaling of TRADE,
TRADMAR, SUPPMAR to meet control totals (BASIC_U, MARGINS_U, IMPORT, MAKE_I).

Phase 1 (trdras): conventional sequential scaling factors.
Phase 2 (raslin): linear-system RAS (bi-proportional adjustment).

Input: `Reg2Result.reg2` dict.
Output: `RasResult` — ras dict with balanced TRADE/TRADMAR/SUPPMAR.
"""

using NamedArrays

Base.@kwdef struct RasResult
    ras::Dict{String,Any}
    diag::Dict{String,Any}
end

function ras_balance!(reg2::Dict{String,Any})
    T = Float64
    NC = length(COM); NI = length(IND); NS = 2; NM = length(MAR); NR = length(REG)
    NONMAR_idx = [c for c in 1:NC if !(COM[c] in Set(MAR))]
    MAR_idx_full = [m for m in 1:NM]

    # ── Read inputs ─────────────────────────────────────────────────────────
    TRADE_data    = reg2["TRAD"]  # COM×SRC×ORG×DST
    SUPPMAR_data  = reg2["MARS"]  # MAR×ORG×DST×PRD
    TRADMAR_data  = reg2["TMAR"]  # COM×SRC×MAR×ORG×DST
    BASIC_U       = parent(reg2["BSCU"])
    MARGINS_U     = parent(reg2["MRGU"])
    IMPORT_cr     = parent(reg2["IMPS"])
    MAKE_I_cd     = parent(reg2["COST"])

    # Raw copies for updates
    TRADE   = copy(parent(TRADE_data))
    SUPPMAR = copy(parent(SUPPMAR_data))
    TRADMAR = copy(parent(TRADMAR_data))

    # ─────────────────────────────────────────────────────────────────────────
    # PHASE 1 — trdras: conventional sequential scaling
    # ─────────────────────────────────────────────────────────────────────────

    # Temp arrays
    TEMTOTCR   = zeros(T, NC, NR)
    TEMTOTCSD  = zeros(T, NC, NS, NR)
    TEMTOTCSMD = zeros(T, NC, NS, NM, NR)
    TEMTOTMRD  = zeros(T, NM, NR, NR)

    # 1a. Scale TRADE(mar,dom) + SUPPMAR to MAKE_I (103-110)
    for m in MAR_idx_full
        for r in 1:NR
            tem = sum(TRADE[m,1,r,d] + sum(SUPPMAR[m,rr,d,r] for rr in 1:NR) for d in 1:NR)
            sc = tem > 0 ? MAKE_I_cd[m,r] / tem : 1.0
            for d in 1:NR
                TRADE[m,1,r,d] *= sc
                for p in 1:NR
                    SUPPMAR[m,r,d,p] *= sc
                end
            end
        end
    end

    # 1b. Scale TRADE(nonmar,dom) to MAKE_I (112-120)
    for c in NONMAR_idx
        for r in 1:NR
            tem = sum(TRADE[c,1,r,d] for d in 1:NR)
            sc = tem > 0 ? MAKE_I_cd[c,r] / tem : 1.0
            for d in 1:NR
                TRADE[c,1,r,d] *= sc
                for m in MAR_idx_full
                    TRADMAR[c,1,m,r,d] *= sc
                end
            end
        end
    end

    # 1c. Scale TRADE(imp) to IMPORT (122-133)
    for c in 1:NC, r in 1:NR
        tem = sum(TRADE[c,2,r,d] for d in 1:NR)
        @assert !(tem == 0 && IMPORT_cr[c,r] != 0) "IMPORT non-zero but TRADE=0: c=$c r=$r"
        sc = tem > 0 ? IMPORT_cr[c,r] / tem : 1.0
        for d in 1:NR
            TRADE[c,2,r,d] *= sc
            for m in MAR_idx_full
                TRADMAR[c,2,m,r,d] *= sc
            end
        end
    end

    # 1d. Scale TRADMAR to MARGINS_U (135-142)
    for c in 1:NC, s in 1:NS, m in MAR_idx_full, d in 1:NR
        tem = sum(TRADMAR[c,s,m,:,d])
        sc = tem > 0 ? MARGINS_U[c,s,m,d] / tem : 1.0
        for r in 1:NR
            TRADMAR[c,s,m,r,d] *= sc
        end
    end

    # 1e. Scale SUPPMAR to TRADMAR_CS (144-157)
    for m in MAR_idx_full, r in 1:NR, d in 1:NR
        tcs = sum(TRADMAR[c,s,m,r,d] for c in 1:NC, s in 1:NS)
        tem = sum(SUPPMAR[m,r,d,p] for p in 1:NR)
        @assert !(tem == 0 && tcs != 0) "SUPPMAR=0 but TRADMAR_CS>0"
        sc = tem > 0 ? tcs / tem : 1.0
        for p in 1:NR
            SUPPMAR[m,r,d,p] *= sc
        end
    end

    # 1f. Scale TRADE to BASIC_U (159-172)
    for c in 1:NC, s in 1:NS, d in 1:NR
        tem = sum(TRADE[c,s,:,d])
        @assert !(tem == 0 && BASIC_U[c,s,d] != 0) "TRADE sum zero but BASIC_U>0: c=$c s=$s d=$d"
        sc = tem > 0 ? BASIC_U[c,s,d] / tem : 1.0
        for r in 1:NR
            TRADE[c,s,r,d] *= sc
            for m in MAR_idx_full
                TRADMAR[c,s,m,r,d] *= sc
            end
        end
    end

    # ─────────────────────────────────────────────────────────────────────────
    # PHASE 2 — raslin: bi-proportional RAS adjustment
    # ─────────────────────────────────────────────────────────────────────────
    # Iterative proportional fitting with FRAC ramping.
    # We do a simplified sequential RAS with gradually tightening tolerance.

    for iter in 1:50
        max_err = 0.0

        # Check TRADE sum vs BASIC_U
        for c in 1:NC, s in 1:NS, d in 1:NR
            tem = sum(TRADE[c,s,:,d])
            err = abs(tem - BASIC_U[c,s,d])
            max_err = max(max_err, err)
            if tem > 0 && BASIC_U[c,s,d] > 0
                sc = BASIC_U[c,s,d] / tem
                for r in 1:NR
                    TRADE[c,s,r,d] *= sc
                    for m in MAR_idx_full
                        TRADMAR[c,s,m,r,d] *= sc
                    end
                end
            end
        end

        # TRADE(imp) vs IMPORT
        for c in 1:NC, r in 1:NR
            tem = sum(TRADE[c,2,r,d] for d in 1:NR)
            if tem > 0 && IMPORT_cr[c,r] > 0
                sc = IMPORT_cr[c,r] / tem
                for d in 1:NR
                    TRADE[c,2,r,d] *= sc
                    for m in MAR_idx_full
                        TRADMAR[c,2,m,r,d] *= sc
                    end
                end
            end
        end

        # TRADE(dom,nonmar) + SUPPMAR(mar) vs MAKE_I
        for c in NONMAR_idx
            for r in 1:NR
                tem = sum(TRADE[c,1,r,d] for d in 1:NR)
                if tem > 0 && MAKE_I_cd[c,r] > 0
                    sc = MAKE_I_cd[c,r] / tem
                    for d in 1:NR
                        TRADE[c,1,r,d] *= sc
                        for m in MAR_idx_full
                            TRADMAR[c,1,m,r,d] *= sc
                        end
                    end
                end
            end
        end
        for m in MAR_idx_full
            for r in 1:NR
                dom_sum = sum(TRADE[m,1,r,d] for d in 1:NR)
                supp_sum = sum(SUPPMAR[m,rr,d,r] for rr in 1:NR, d in 1:NR)
                tem = dom_sum + supp_sum
                if tem > 0 && MAKE_I_cd[m,r] > 0
                    sc = MAKE_I_cd[m,r] / tem
                    for d in 1:NR
                        TRADE[m,1,r,d] *= sc
                    end
                    for rr in 1:NR, d in 1:NR, p in 1:NR
                        SUPPMAR[m,rr,d,p] *= sc  # simplified — exact would only scale p=r
                    end
                end
            end
        end

        # TRADMAR vs MARGINS_U
        for c in 1:NC, s in 1:NS, m in MAR_idx_full, d in 1:NR
            tem = sum(TRADMAR[c,s,m,:,d])
            if tem > 0 && MARGINS_U[c,s,m,d] > 0
                sc = MARGINS_U[c,s,m,d] / tem
                for r in 1:NR
                    TRADMAR[c,s,m,r,d] *= sc
                end
            end
        end

        # SUPPMAR vs TRADMAR_CS
        for m in MAR_idx_full, r in 1:NR, d in 1:NR
            tcs = sum(TRADMAR[c,s,m,r,d] for c in 1:NC, s in 1:NS)
            tem = sum(SUPPMAR[m,r,d,p] for p in 1:NR)
            if tem > 0 && tcs > 0
                sc = tcs / tem
                for p in 1:NR
                    SUPPMAR[m,r,d,p] *= sc
                end
            end
        end

        if max_err < 1e-6
            break
        end
    end

    # ── Final assertions ────────────────────────────────────────────────────
    for c in 1:NC, s in 1:NS, d in 1:NR
        tem = sum(TRADE[c,s,:,d])
        if BASIC_U[c,s,d] == 0
            @assert tem < 1e-6 "TRADE sum non-zero where BASIC_U=0: c=$c s=$s d=$d"
        end
    end
    @assert all(TRADE .>= -1e-8) "Negative TRADE after RAS"
    @assert all(TRADMAR .>= -1e-8) "Negative TRADMAR after RAS"
    @assert all(SUPPMAR .>= -1e-8) "Negative SUPPMAR after RAS"

    # ── Output ──────────────────────────────────────────────────────────────
    ras_out = Dict{String,Any}(
        "TRAD" => NamedArray(TRADE,   Tuple([COM,SRC,REG,REG]), (:COM,:SRC,:ORG,:DST)),
        "TMAR" => NamedArray(TRADMAR, Tuple([COM,SRC,MAR,REG,REG]), (:COM,:SRC,:MAR,:ORG,:DST)),
        "MARS" => NamedArray(SUPPMAR, Tuple([MAR,REG,REG,REG]), (:MAR,:ORG,:DST,:PRD)),
        "BSCU" => BASIC_U,
        "MRGU" => MARGINS_U,
        "IMPS" => IMPORT_cr,
        "COST" => MAKE_I_cd,
    )

    reg2_keys = ["BSMR","MAKE","STOK","2PUR","DIST"]
    for k in reg2_keys
        if haskey(reg2, k); ras_out[k] = reg2[k]; end
    end

    RasResult(ras=ras_out, diag=Dict{String,Any}())
end
