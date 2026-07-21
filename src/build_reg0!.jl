"""
Translation of `reg0.tab` — reformat ORANI-G-style national database, normalize
distance matrix, make preliminary `DISTGONE` estimate.

Input: header-dicts from `read_national_data()`, `read_regsupp_data()`.
Output: `Reg0Result` — see struct fields below.
"""

using NamedArrays

Base.@kwdef struct Reg0Result
    reg0::Dict{String,Any}
    elast::Dict{String,Any}
    diag::Dict{String,Any}
end

function _na(arr, dims)
    d = copy(parent(arr))
    dimvals = Tuple(collect(keys(arr.dicts[k])) for k in 1:ndims(arr))
    if ndims(arr) == 1
        NamedArray(vec(d), dimvals, Tuple(Symbol.(dims)))
    else
        NamedArray(d, dimvals, Tuple(Symbol.(dims)))
    end
end

function build_reg0!(nat::Dict{String,Any}, regsupp::Dict{String,Any})
    T = Float64

    # ── Unpack raw data ────────────────────────────────────────────────────
    V1BAS  = _na(nat["1BAS"], ["COM","SRC","IND"])
    V2BAS  = _na(nat["2BAS"], ["COM","SRC","IND"])
    V3BAS  = _na(nat["3BAS"], ["COM","SRC"])
    V4BASA = vec(nat["4BAS"].array)
    V5BAS  = _na(nat["5BAS"], ["COM","SRC"])
    V6BAS  = _na(nat["6BAS"], ["COM","SRC"])

    V1TAX  = _na(nat["1TAX"], ["COM","SRC","IND"])
    V2TAX  = _na(nat["2TAX"], ["COM","SRC","IND"])
    V3TAX  = _na(nat["3TAX"], ["COM","SRC"])
    V4TAXA = vec(nat["4TAX"].array)
    V5TAX  = _na(nat["5TAX"], ["COM","SRC"])

    V1MAR  = _na(nat["1MAR"], ["COM","SRC","IND","MAR"])
    V2MAR  = _na(nat["2MAR"], ["COM","SRC","IND","MAR"])
    V3MAR  = _na(nat["3MAR"], ["COM","SRC","MAR"])
    V4MARA = _na(nat["4MAR"], ["COM","MAR"])
    V5MAR  = _na(nat["5MAR"], ["COM","SRC","MAR"])

    V1CAP  = vec(nat["1CAP"].array)
    V1LAB  = _na(nat["1LAB"], ["IND","OCC"])
    V1LND  = vec(nat["1LND"].array)
    V1OCT  = vec(nat["1OCT"].array)
    V1PTX  = vec(nat["1PTX"].array)
    MAKE   = _na(nat["MAKE"], ["COM","IND"])
    V0TAR  = vec(nat["0TAR"].array)

    R001   = _na(nat["R001"], ["IND","REG"])
    R002   = _na(nat["R002"], ["IND","REG"])
    R003   = _na(nat["R003"], ["COM","REG"])
    R004   = _na(nat["R004"], ["COM","REG"])
    R005   = _na(nat["R005"], ["COM","REG"])

    # ── Derived sets ───────────────────────────────────────────────────────
    FINDEM_NAMES = ["Hou","Inv","Gov","Exp"]
    USR_NAMES    = vcat(IND, FINDEM_NAMES)
    FACTOR_NAMES = ["Labour","Capital","Land","ProdTax"]
    NUSR = length(USR_NAMES)
    NFAC = length(FACTOR_NAMES)
    NC = length(COM); NI = length(IND); NS = 2; NM = length(MAR); NR = length(REG); NO = length(OCC)

    # ── Add import dimension to exports (104-110) ──────────────────────────
    V4BAS = zeros(T, NC, NS)
    V4MAR = zeros(T, NC, NS, NM)
    V4TAX = zeros(T, NC, NS)
    for c in 1:NC
        V4BAS[c,1] = V4BASA[c]
        V4TAX[c,1] = V4TAXA[c]
        for m in 1:NM
            V4MAR[c,1,m] = V4MARA[c,m]
        end
    end

    # ── MAKE_I / MAKE_C (129-133) ──────────────────────────────────────────
    MAKE_I = vec(sum(MAKE.array, dims=2))
    MAKE_C = vec(sum(MAKE.array, dims=1))

    # ── Capital floor (135-143) ────────────────────────────────────────────
    CAPFRAC = 1e-9
    CAPADJ = zeros(T, NI)
    for i in 1:NI
        if V1CAP[i] == 0
            CAPADJ[i] = CAPFRAC * MAKE_C[i]
            V1CAP[i]  = CAPADJ[i]
        end
    end

    # ── Remove stocks from MAKE (146-153) ──────────────────────────────────
    STOCKMAKE = zeros(T, NC, NI)
    for c in 1:NC
        mc = MAKE_I[c]
        if mc > 0
            for i in 1:NI
                STOCKMAKE[c,i] = V6BAS[c,1] * MAKE[c,i] / mc
            end
        end
    end
    STOCKS = [sum(STOCKMAKE[:,i]) for i in 1:NI]
    for i in 1:NI, c in 1:NC
        MAKE[c,i] -= STOCKMAKE[c,i]
    end
    for c in 1:NC; MAKE_I[c] = sum(MAKE[c,:]); end
    for i in 1:NI; MAKE_C[i] = sum(MAKE[:,i]); end

    # ── BASIC flows (155-161) ──────────────────────────────────────────────
    # BASIC(c,s,u): COM × SRC × USR
    BASIC_a = zeros(T, NC, NS, NUSR)
    for c in 1:NC, s in 1:NS, i in 1:NI
        BASIC_a[c,s,i] = V1BAS[c,s,i]
    end
    for c in 1:NC, s in 1:NS
        BASIC_a[c,s, findfirst(==("Hou"), USR_NAMES)] = V3BAS[c,s]
        BASIC_a[c,s, findfirst(==("Exp"), USR_NAMES)] = V4BAS[c,s]
        BASIC_a[c,s, findfirst(==("Gov"), USR_NAMES)] = V5BAS[c,s]
        BASIC_a[c,s, findfirst(==("Inv"), USR_NAMES)] = sum(V2BAS[c,s,i] for i in 1:NI)
    end

    # ── Factor costs (164-172) ─────────────────────────────────────────────
    V1LAB_O = vec(sum(V1LAB.array, dims=2))
    OCCSHR  = zeros(T, NI, NO)
    for i in 1:NI
        lo = V1LAB_O[i]
        for o in 1:NO
            OCCSHR[i,o] = lo > 0 ? V1LAB[i,o] / lo : 0.0
        end
    end

    FAC_a = zeros(T, NI, NFAC)
    for i in 1:NI
        FAC_a[i, findfirst(==("Labour"),  FACTOR_NAMES)] = V1LAB_O[i]
        FAC_a[i, findfirst(==("Capital"), FACTOR_NAMES)] = V1CAP[i]
        FAC_a[i, findfirst(==("Land"),    FACTOR_NAMES)] = V1LND[i]
        FAC_a[i, findfirst(==("ProdTax"), FACTOR_NAMES)] = V1OCT[i] + V1PTX[i]
    end

    # ── Commodity taxes TAX (174-180) ──────────────────────────────────────
    TAX_a = zeros(T, NC, NS, NUSR)
    for c in 1:NC, s in 1:NS, i in 1:NI
        TAX_a[c,s,i] = V1TAX[c,s,i]
    end
    for c in 1:NC, s in 1:NS
        TAX_a[c,s, findfirst(==("Inv"), USR_NAMES)] = sum(V2TAX[c,s,i] for i in 1:NI)
        TAX_a[c,s, findfirst(==("Hou"), USR_NAMES)] = V3TAX[c,s]
        TAX_a[c,s, findfirst(==("Exp"), USR_NAMES)] = V4TAX[c,s]
        TAX_a[c,s, findfirst(==("Gov"), USR_NAMES)] = V5TAX[c,s]
    end

    # ── Tariff adjustment (186-199) ────────────────────────────────────────
    IMPTOT = [sum(BASIC_a[c,2,:]) for c in 1:NC]
    TARF_a = zeros(T, NC, NUSR)
    for c in 1:NC
        if IMPTOT[c] > 0
            for u in 1:NUSR
                TARF_a[c,u] = V0TAR[c] * BASIC_a[c,2,u] / IMPTOT[c]
                TAX_a[c,2,u]  += TARF_a[c,u]
                BASIC_a[c,2,u] -= TARF_a[c,u]
            end
        end
    end

    # ── MARGINS (202-209) ──────────────────────────────────────────────────
    MARGINS_a = zeros(T, NC, NS, NUSR, NM)
    for c in 1:NC, s in 1:NS, i in 1:NI, m in 1:NM
        MARGINS_a[c,s,i,m] = V1MAR[c,s,i,m]
    end
    for c in 1:NC, s in 1:NS, m in 1:NM
        MARGINS_a[c,s, findfirst(==("Inv"), USR_NAMES), m] = sum(V2MAR[c,s,i,m] for i in 1:NI)
        MARGINS_a[c,s, findfirst(==("Hou"), USR_NAMES), m] = V3MAR[c,s,m]
        MARGINS_a[c,s, findfirst(==("Exp"), USR_NAMES), m] = V4MAR[c,s,m]
        MARGINS_a[c,s, findfirst(==("Gov"), USR_NAMES), m] = V5MAR[c,s,m]
    end

    MAR_OGL = copy(MARGINS_a)

    # ── Zero margins on zero BASIC + rescale (220-244) ─────────────────────
    MARTOT_um = zeros(T, NUSR, NM)
    for u in 1:NUSR, m in 1:NM
        MARTOT_um[u,m] = sum(MARGINS_a[c,s,u,m] for c in 1:NC, s in 1:NS)
    end

    for c in 1:NC, s in 1:NS, u in 1:NUSR, m in 1:NM
        if BASIC_a[c,s,u] == 0
            MARGINS_a[c,s,u,m] = 0.0
        end
    end

    MARTOT2_um = zeros(T, NUSR, NM)
    for u in 1:NUSR, m in 1:NM
        MARTOT2_um[u,m] = sum(MARGINS_a[c,s,u,m] for c in 1:NC, s in 1:NS)
    end
    for u in 1:NUSR, m in 1:NM
        if MARTOT2_um[u,m] > 0
            sf = MARTOT_um[u,m] / MARTOT2_um[u,m]
            for c in 1:NC, s in 1:NS
                MARGINS_a[c,s,u,m] *= sf
            end
        end
    end

    # ── Bad tax on zero basic flows (247-258) ──────────────────────────────
    TRNSTX = zeros(T, NI)
    for c in 1:NC, s in 1:NS, u in 1:NUSR
        if BASIC_a[c,s,u] <= 0
            if u <= NI
                TRNSTX[u] += TAX_a[c,s,u]
            end
            TAX_a[c,s,u] = 0.0
        end
    end
    for i in 1:NI
        FAC_a[i, findfirst(==("ProdTax"), FACTOR_NAMES)] += TRNSTX[i]
    end

    # ── Enforce same margin rate for dom/imp (260-268) ─────────────────────
    MARGINS_S_cum = zeros(T, NC, NUSR, NM)
    BASIC_S_cu    = zeros(T, NC, NUSR)
    for c in 1:NC, u in 1:NUSR, m in 1:NM
        MARGINS_S_cum[c,u,m] = sum(MARGINS_a[c,s,u,m] for s in 1:NS)
        BASIC_S_cu[c,u]      = sum(BASIC_a[c,s,u] for s in 1:NS)
    end
    for c in 1:NC, s in 1:NS, u in 1:NUSR, m in 1:NM
        if BASIC_S_cu[c,u] > 0
            MARGINS_a[c,s,u,m] = MARGINS_S_cum[c,u,m] * BASIC_a[c,s,u] / BASIC_S_cu[c,u]
        end
    end

    # ── Assert: TAX=0, MARGINS=0 where BASIC=0 ─────────────────────────────
    for c in 1:NC, s in 1:NS, u in 1:NUSR
        @assert BASIC_a[c,s,u] <= 0 ? abs(TAX_a[c,s,u]) < 1e-10 : true
        for m in 1:NM
            @assert BASIC_a[c,s,u] <= 0 ? abs(MARGINS_a[c,s,u,m]) < 1e-10 : true
        end
    end

    # ── COSTS / MAKE_C / DIFFIND / DIFFCOM (319-354) ──────────────────────
    COSTS_i = zeros(T, NI)
    for i in 1:NI
        val = sum(FAC_a[i,:])
        for c in 1:NC, s in 1:NS
            val += BASIC_a[c,s,i] + TAX_a[c,s,i] + sum(MARGINS_a[c,s,i,m] for m in 1:NM)
        end
        COSTS_i[i] = val
    end
    DIFFIND_i = [COSTS_i[i] - MAKE_C[i] - STOCKS[i] for i in 1:NI]

    SALES_c = zeros(T, NC)
    for c in 1:NC
        SALES_c[c] = sum(BASIC_a[c,1,u] for u in 1:NUSR) + V6BAS[c,1]
    end
    MARSALES_m = zeros(T, NM)
    for m in 1:NM
        MARSALES_m[m] = sum(MARGINS_a[c,s,u,m] for c in 1:NC, s in 1:NS, u in 1:NUSR)
    end
    for mi in 1:NM
        m_com = findfirst(==(MAR[mi]), COM)
        if m_com !== nothing
            SALES_c[m_com] += MARSALES_m[mi]
        end
    end
    DIFFCOM_c = [SALES_c[c] - MAKE_I[c] for c in 1:NC]

    # Scaled diagnostics
    DIFFIND_sc = [COSTS_i[i] > 0 ? DIFFIND_i[i] / COSTS_i[i] : 0.0 for i in 1:NI]
    DIFFCOM_sc = [SALES_c[c] > 0 ? DIFFCOM_c[c] / SALES_c[c] : 0.0 for c in 1:NC]
    maxDIFFIND = maximum(abs.(DIFFIND_sc))
    maxDIFFCOM = maximum(abs.(DIFFCOM_sc))
    println("  max DIFFIND_sc = $maxDIFFIND, max DIFFCOM_sc = $maxDIFFCOM")
    worst = sortperm(abs.(DIFFCOM_sc), rev=true)[1:min(5,NC)]
    for ci in worst
        println("    COM[$ci] ($(COM[ci])): MAKE_I=$(MAKE_I[ci]), SALES=$(SALES_c[ci]), DIFFCOM_sc=$(DIFFCOM_sc[ci])")
    end
    if maxDIFFIND >= 0.05
        @warn "DIFFIND large: $maxDIFFIND"
    end
    if maxDIFFCOM >= 0.50
        @warn "DIFFCOM large: $maxDIFFCOM"
    end

    # ── Investment shares (292-306) ────────────────────────────────────────
    V2PUR_csi = zeros(T, NC, NS, NI)
    for c in 1:NC, s in 1:NS, i in 1:NI
        V2PUR_csi[c,s,i] = V2BAS[c,s,i] + V2TAX[c,s,i] + sum(V2MAR[c,s,i,m] for m in 1:NM)
    end
    V2PUR_S_ci = [sum(V2PUR_csi[c,s,i] for s in 1:NS) for c in 1:NC, i in 1:NI]
    V2PUR_SI_c = [sum(V2PUR_S_ci[c,i] for i in 1:NI) for c in 1:NC]
    INVSHR_ci  = zeros(T, NC, NI)
    for c in 1:NC, i in 1:NI
        if V2PUR_SI_c[c] > 0
            INVSHR_ci[c,i] = V2PUR_S_ci[c,i] / V2PUR_SI_c[c]
        end
    end

    # ── Elasticities (308-317) ─────────────────────────────────────────────
    elast_dict = Dict{String,Any}()
    elast_dict["P015"] = nat["1ARM"]
    for h in ["P028","P018","P021","XPEL","SLAB","SCET"]
        if haskey(nat,h); elast_dict[h] = nat[h]; end
    end

    # ── Regional shares: floor & normalize (369-430) ───────────────────────
    DEFREGSHR_r = zeros(T, NR)
    for r in 1:NR
        DEFREGSHR_r[r] = sum(R001[i,r] * sum(FAC_a[i,:]) for i in 1:NI)
    end
    t = sum(DEFREGSHR_r); if t > 0; DEFREGSHR_r ./= t; end
    MINSHR_r = [DEFREGSHR_r[r] / 1e5 for r in 1:NR]

    for i in 1:NI, r in 1:NR
        if R001[i,r] < MINSHR_r[r]; R001[i,r] = MINSHR_r[r]; end
        if R002[i,r] < MINSHR_r[r]; R002[i,r] = MINSHR_r[r]; end
    end
    for c in 1:NC, r in 1:NR
        if R003[c,r] < MINSHR_r[r]; R003[c,r] = MINSHR_r[r]; end
        if R004[c,r] < MINSHR_r[r]; R004[c,r] = MINSHR_r[r]; end
        if R005[c,r] < MINSHR_r[r]; R005[c,r] = MINSHR_r[r]; end
    end

    for i in 1:NI
        t = sum(R001[i,:]); if t > 0; R001[i,:] ./= t; end
        t = sum(R002[i,:]); if t > 0; R002[i,:] ./= t; end
    end
    for c in 1:NC
        t = sum(R003[c,:]); if t > 0; R003[c,:] ./= t; end
        t = sum(R004[c,:]); if t > 0; R004[c,:] ./= t; end
        t = sum(R005[c,:]); if t > 0; R005[c,:] ./= t; end
    end

    # ── R001 adjustment for local commodities (432-467) ────────────────────
    LCOM_arr = nat["LCOM"]
    local_flags = vec(parent(LCOM_arr))
    RLOC_idx = [c for c in 1:NC if local_flags[c] > 0]
    RLOC_vec = COM[RLOC_idx]
    LOCRAT_i = zeros(T, NI)
    for i in 1:NI
        if MAKE_C[i] > 0
            LOCRAT_i[i] = sum(MAKE[c,i] for c in RLOC_idx) / MAKE_C[i]
        end
    end
    LOCIND_set = Set(i for i in 1:NI if LOCRAT_i[i] > 0.9)

    # TOTDEMREG(c,r)
    TOTDEMREG_cr = zeros(T, NC, NR)
    for c in 1:NC, r in 1:NR
        val = sum(R001[i,r] * V1BAS[c,1,i] + R002[i,r] * V2BAS[c,1,i] for i in 1:NI)
        val += R003[c,r] * V3BAS[c,1] + R004[c,r] * V4BAS[c,1] + R005[c,r] * V5BAS[c,1]
        TOTDEMREG_cr[c,r] = val
    end
    # margins demand for margin commodities
    for cm in 1:NM
        for r in 1:NR
            val = sum(
                R003[cc,r] * V3MAR[cc,s,cm] + R004[cc,r] * V4MAR[cc,s,cm] + R005[cc,r] * V5MAR[cc,s,cm]
                + sum(R001[i,r] * V1MAR[cc,s,i,cm] + R002[i,r] * V2MAR[cc,s,i,cm] for i in 1:NI)
                for cc in 1:NC, s in 1:NS
            )
            TOTDEMREG_cr[cm,r] += val
        end
    end

    # TOTDEMINDREG & adjust R001 for local industries
    INDTOT_i = zeros(T, NI)
    for i in 1:NI
        if i in LOCIND_set
            tot = 0.0
            for c in 1:NC, r in 1:NR
                if MAKE_I[c] > 0
                    tot += TOTDEMREG_cr[c,r] * MAKE[c,i] / MAKE_I[c]
                end
            end
            INDTOT_i[i] = tot
        end
    end
    for i in 1:NI, r in 1:NR
        val = sum(TOTDEMREG_cr[c,r] * MAKE[c,i] / MAKE_I[c] for c in 1:NC if MAKE_I[c] > 0)
        if i in LOCIND_set && INDTOT_i[i] > 0
            R001[i,r] = val / INDTOT_i[i]
        end
    end

    # ── Distance matrix (476-624) ──────────────────────────────────────────
    DIST_raw = regsupp["DIST"]
    D = copy(DIST_raw.array)
    for r in 1:NR, d in 1:NR
        if D[r,d] < 0; D[r,d] = 0; end
    end
    for r in 1:NR, d in 1:NR
        if D[r,d] == 0 && D[d,r] > 0; D[r,d] = D[d,r]; end
    end
    D = max.(D, D')  # symmetrise

    NCELL     = sum(D .> 0)
    DISTTOT   = sum(D[D .> 0])
    AVEDIST   = DISTTOT / max(NCELL, 1)
    MAXERR    = AVEDIST / 1000

    for r in 1:NR, d in 1:NR
        if r != d && D[r,d] <= 0; D[r,d] = DISTTOT; end
    end

    CLOSESTDIST = zeros(T, NR)
    for r in 1:NR
        m = minimum(D[r, setdiff(1:NR, r)])
        CLOSESTDIST[r] = m
    end
    for r in 1:NR
        if D[r,r] == 0; D[r,r] = 0.5 * CLOSESTDIST[r]; end
    end

    # Floyd (12 iters)
    ND = copy(D)
    for _ in 1:12
        TMP = [minimum(ND[r,:] + ND[:,d]) for r in 1:NR, d in 1:NR]
        ND  = min.(ND, TMP)
    end

    FLOYD_SAV = D - ND
    @assert all(abs.(FLOYD_SAV) .< MAXERR)
    @assert all(abs.(FLOYD_SAV - FLOYD_SAV') .< MAXERR)

    # Scale to avg 10
    scale = 10 / (sum(ND) / (NR*NR))
    ND .*= scale
    DISTANCE = NamedArray(ND, (REG, REG), (:ORG,:DST))

    # ── Preliminary DISTGONE (626-639) ─────────────────────────────────────
    DISTAVE_r = zeros(T, NR)
    for r in 1:NR
        DISTAVE_r[r] = sum(ND[q,r] for q in 1:NR) / NR
        DISTAVE_r[r] = (DISTAVE_r[r] + ND[r,r]) / 2
    end
    DGON_csd = zeros(T, NC, NS, NR)
    for c in 1:NC, s in 1:NS, d in 1:NR
        DGON_csd[c,s,d] = DISTAVE_r[d]
    end

    # ── Assemble output dicts ──────────────────────────────────────────────
    reg0_out = Dict{String,Any}(
        "FACT" => NamedArray(FAC_a,         Tuple([IND, FACTOR_NAMES]),  (:IND,:FACTOR)),
        "UBAS" => NamedArray(BASIC_a,       Tuple([COM,SRC,USR_NAMES]), (:COM,:SRC,:USR)),
        "UMAR" => NamedArray(MARGINS_a,     Tuple([COM,SRC,USR_NAMES,MAR]), (:COM,:SRC,:USR,:MAR)),
        "UTAX" => NamedArray(TAX_a,         Tuple([COM,SRC,USR_NAMES]), (:COM,:SRC,:USR)),
        "MAKE" => MAKE,
        "STOK" => NamedArray(STOCKS,        IND,  (:IND,)),
        "1LAB" => V1LAB,
        "OSHR" => NamedArray(OCCSHR,        Tuple([IND,OCC]), (:IND,:OCC)),
        "R001" => R001,
        "R002" => R002,
        "R003" => R003,
        "R004" => R004,
        "R005" => R005,
        "ISHR" => NamedArray(INVSHR_ci,     Tuple([COM,IND]), (:COM,:IND)),
        "DIST" => DISTANCE,
    )

    diag_out = Dict{String,Any}(
        "CADJ"  => CAPADJ,
        "TARF"  => TARF_a,
        "TRTX"  => TRNSTX,
        "RLOC"  => RLOC_idx,
        "DIND"  => DIFFIND_i,
        "DCOM"  => DIFFCOM_c,
        "SIND"  => DIFFIND_sc,
        "SCOM"  => DIFFCOM_sc,
        "LIND"  => LOCRAT_i,
        "DERR"  => FLOYD_SAV,
        "DIST"  => ND,
    )

    Reg0Result(reg0=reg0_out, elast=elast_dict, diag=diag_out)
end
