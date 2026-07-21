function prepare_parameters!(agg::Dict{String,Any})
    T = Float64
    na = size(agg["MAKE"], 1)
    nr = size(agg["MAKE"], 3)
    ns = 2
    nm = haskey(agg, "TMAR") ? size(agg["TMAR"], 3) : length(MAR)
    no = size(agg["1LAB"], 2)

    MAKE   = parent(agg["MAKE"])
    TRADE  = parent(agg["TRAD"])
    TMAR   = haskey(agg,"TMAR") ? parent(agg["TMAR"]) : zeros(T, na, ns, nm, nr, nr)
    MARS   = haskey(agg,"MARS") ? parent(agg["MARS"]) : zeros(T, nm, nr, nr, nr)
    DIST   = haskey(agg,"DIST") ? parent(agg["DIST"]) : zeros(T, nr, nr)
    V1LAB  = parent(agg["1LAB"])
    V1CAP  = parent(agg["1CAP"])
    V1LND  = parent(agg["1LND"])

    BSMR   = haskey(agg,"BSMR") ? parent(agg["BSMR"]) : nothing
    UTAX   = haskey(agg,"UTAX") ? parent(agg["UTAX"]) : nothing
    V2PUR  = haskey(agg,"2PUR") ? parent(agg["2PUR"]) : nothing
    VSTOK  = haskey(agg,"STOK") ? parent(agg["STOK"]) : nothing
    V1PTX  = haskey(agg,"1PTX") ? parent(agg["1PTX"]) : nothing

    SLAB   = haskey(agg,"SLAB") ? vec(parent(agg["SLAB"])) : zeros(T, na)
    P028   = haskey(agg,"P028") ? vec(parent(agg["P028"])) : zeros(T, na)
    SGDD   = haskey(agg,"SGDD") ? vec(parent(agg["SGDD"])) : zeros(T, na)
    SMAR_v = haskey(agg,"SMAR") ? vec(parent(agg["SMAR"])) : zeros(T, nm)
    PO01   = haskey(agg,"PO01") ? vec(parent(agg["PO01"])) : zeros(T, nr)
    P021_v = if haskey(agg,"P021")
        raw = agg["P021"] isa NamedArray ? parent(agg["P021"]) : agg["P021"]
        v = vec(Float64.(raw))
        length(v) == nr ? v : fill(v[1], nr)
    else
        fill(-2.0, nr)
    end
    XPEL   = haskey(agg,"XPEL") ? parent(agg["XPEL"]) : zeros(T, na, nr)
    SCET   = haskey(agg,"SCET") ? vec(parent(agg["SCET"])) : zeros(T, na)
    P018   = haskey(agg,"P018") ? vec(parent(agg["P018"])) : zeros(T, na)

    p = Dict{String,Any}()

    # ── LAB_O: total labour bill ──────────────────────────────
    LAB_O = zeros(T, na, nr)
    for i in 1:na, d in 1:nr
        LAB_O[i,d] = sum(V1LAB[i,o,d] for o in 1:no)
    end
    p["LAB_O"] = LAB_O
    p["CAP"]  = V1CAP
    p["LND"]  = V1LND
    p["LAB"]  = parent(agg["1LAB"])
    if V1PTX !== nothing
        p["PRODTAX"] = V1PTX
    else
        p["PRODTAX"] = zeros(T, na, nr)
    end
    if VSTOK !== nothing
        p["STOCKS"] = VSTOK
    else
        p["STOCKS"] = zeros(T, na, nr)
    end
    if V2PUR !== nothing
        p["INVEST"] = V2PUR
    else
        p["INVEST"] = zeros(T, na, na, nr)
    end

    # ── PRIM: primary factor composite ─────────────────────────
    PRIM = LAB_O + V1CAP + V1LND
    p["PRIM"] = PRIM

    PRIMCOST = cat(V1LND, LAB_O, V1CAP; dims=ndims(V1LND)+1)  # na×nr×3
    p["PRIMCOST"] = PRIMCOST

    # ── MAKE aggregates ────────────────────────────────────────
    MAKE_C = zeros(T, na, nr)
    MAKE_I = zeros(T, na, nr)
    for c in 1:na, i in 1:na, d in 1:nr
        MAKE_C[i,d] += MAKE[c,i,d]
        MAKE_I[c,d] += MAKE[c,i,d]
    end
    p["MAKE_C"] = MAKE_C
    p["MAKE_I"] = MAKE_I

    MAKE_D = zeros(T, na, na)
    for c in 1:na, i in 1:na
        MAKE_D[c,i] = sum(MAKE[c,i,d] for d in 1:nr)
    end
    p["MAKE_D"] = MAKE_D

    MAKESHR1 = zeros(T, na, na, nr)
    MAKESHR2 = zeros(T, na, na, nr)
    for c in 1:na, i in 1:na, d in 1:nr
        MAKESHR1[c,i,d] = MAKE_C[i,d] > 0 ? MAKE[c,i,d] / MAKE_C[i,d] : 0.0
        MAKESHR2[c,i,d] = MAKE_I[c,d] > 0 ? MAKE[c,i,d] / MAKE_I[c,d] : 0.0
    end
    p["MAKESHR1"] = MAKESHR1
    p["MAKESHR2"] = MAKESHR2

    # ── Trade sums ─────────────────────────────────────────────
    TRDIAG = zeros(T, na, ns, nr)
    for c in 1:na, s in 1:ns, r in 1:nr
        TRDIAG[c,s,r] = TRADE[c,s,r,r]
    end
    p["TRDIAG"] = TRDIAG

    TRADE_D = zeros(T, na, ns, nr)
    TRADE_R = zeros(T, na, ns, nr)
    for c in 1:na, s in 1:ns, r in 1:nr, d in 1:nr
        TRADE_D[c,s,r] += TRADE[c,s,r,d]
        TRADE_R[c,s,d] += TRADE[c,s,r,d]
    end
    p["TRADE_D"] = TRADE_D
    p["TRADE_R"] = TRADE_R

    TRADE_RD = zeros(T, na, ns)
    for c in 1:na, s in 1:ns
        TRADE_RD[c,s] = sum(TRADE_D[c,s,r] for r in 1:nr)
    end
    p["TRADE_RD"] = TRADE_RD

    # ── DELIVRD, BASSHR, MARSHR (Excerpt 19) ──────────────────
    DELIVRD = zeros(T, na, ns, nr, nr)
    BASSHR  = zeros(T, na, ns, nr, nr)
    MARSHR  = zeros(T, na, ns, nm, nr, nr)
    for c in 1:na, s in 1:ns, r in 1:nr, d in 1:nr
        tm_sum = sum(TMAR[c,s,m,r,d] for m in 1:nm)
        DELIVRD[c,s,r,d] = TRADE[c,s,r,d] + tm_sum
        dv = DELIVRD[c,s,r,d]
        BASSHR[c,s,r,d] = dv > 0 ? TRADE[c,s,r,d] / dv : 0.0
        for m in 1:nm
            MARSHR[c,s,m,r,d] = dv > 0 ? TMAR[c,s,m,r,d] / dv : 0.0
        end
    end
    p["MAKE"] = MAKE
    p["DELIVRD"] = DELIVRD
    p["BASSHR"] = BASSHR
    p["MARSHR"] = MARSHR

    # ── DELIVRD_R (Excerpt 20) ─────────────────────────────────
    DELIVRD_R = zeros(T, na, ns, nr)
    for c in 1:na, s in 1:ns, d in 1:nr
        DELIVRD_R[c,s,d] = sum(DELIVRD[c,s,r,d] for r in 1:nr)
    end
    p["DELIVRD_R"] = DELIVRD_R
    p["SGDD"] = SGDD

    # ── Margin aggregates (Excerpt 21) ─────────────────────────
    TRADMAR_CS = zeros(T, nm, nr, nr)
    for m in 1:nm, r in 1:nr, d in 1:nr
        TRADMAR_CS[m,r,d] = sum(TMAR[c,s,m,r,d] for c in 1:na, s in 1:ns)
    end
    p["TRADMAR_CS"] = TRADMAR_CS

    SUPPMAR_P = zeros(T, nm, nr, nr)
    SUPPMAR_D = zeros(T, nm, nr, nr)
    SUPPMAR_R = zeros(T, nm, nr, nr)
    for m in 1:nm, r in 1:nr, d in 1:nr, p_ in 1:nr
        SUPPMAR_P[m,r,d] += MARS[m,r,d,p_]
        SUPPMAR_D[m,r,p_] += MARS[m,r,d,p_]
        SUPPMAR_R[m,d,p_] += MARS[m,r,d,p_]
    end
    p["SUPPMAR_P"] = SUPPMAR_P
    p["SUPPMAR_D"] = SUPPMAR_D
    p["SUPPMAR_R"] = SUPPMAR_R

    SUPPMAR_RD = zeros(T, nm, nr)
    for m in 1:nm, p_ in 1:nr
        SUPPMAR_RD[m,p_] = sum(SUPPMAR_D[m,r,p_] for r in 1:nr)
    end
    p["SUPPMAR_RD"] = SUPPMAR_RD

    # ── PUR/PUR_S/PUR_CS/SRCSHR (Excerpt 7) — needs BSMR+UTAX ─
    if BSMR !== nothing && UTAX !== nothing
        nu = size(BSMR, 3)  # na + 4
        USE = BSMR
        TAX = UTAX
        PUR = USE + TAX
        p["USE"] = USE
        p["TAX"] = TAX
        p["PUR"] = PUR

        PUR_S = zeros(T, na, nu, nr)
        for c in 1:na, u in 1:nu, d in 1:nr
            PUR_S[c,u,d] = PUR[c,1,u,d] + PUR[c,2,u,d]
        end
        p["PUR_S"] = PUR_S

        PUR_CS = zeros(T, nu, nr)
        for u in 1:nu, d in 1:nr
            PUR_CS[u,d] = sum(PUR_S[c,u,d] for c in 1:na)
        end
        p["PUR_CS"] = PUR_CS

        SRCSHR = zeros(T, na, ns, nu, nr)
        for c in 1:na, s in 1:ns, u in 1:nu, d in 1:nr
            ps = PUR_S[c,u,d]
            SRCSHR[c,s,u,d] = ps > 0 ? PUR[c,s,u,d] / ps : 0.5
        end
        p["SRCSHR"] = SRCSHR

        PUR_D = zeros(T, na, ns, nu)
        for c in 1:na, s in 1:ns, u in 1:nu
            PUR_D[c,s,u] = sum(PUR[c,s,u,d] for d in 1:nr)
        end
        p["PUR_D"] = PUR_D

        USE_U = zeros(T, na, ns, nr)
        USE_I = zeros(T, na, ns, nr)
        for c in 1:na, s in 1:ns, d in 1:nr
            for u in 1:nu
                USE_U[c,s,d] += USE[c,s,u,d]
                if u <= na
                    USE_I[c,s,d] += USE[c,s,u,d]
                end
            end
        end
        p["USE_U"] = USE_U
        p["USE_I"] = USE_I

        # LOCUSE (Excerpt 18)
        u_exp = na + 4
        LOCUSE = zeros(T, na, ns, nr)
        for c in 1:na, s in 1:ns, d in 1:nr
            LOCUSE[c,s,d] = USE_U[c,s,d] - USE[c,s,u_exp,d]
        end
        p["LOCUSE"] = LOCUSE

        LOCUSE_S = zeros(T, na, nr)
        LOCUSE_SD = zeros(T, na)
        IMPSHR = zeros(T, na, nr)
        for c in 1:na, d in 1:nr
            LOCUSE_S[c,d] = LOCUSE[c,1,d] + LOCUSE[c,2,d]
            IMPSHR[c,d] = LOCUSE_S[c,d] > 0 ? LOCUSE[c,2,d] / LOCUSE_S[c,d] : 0.0
        end
        for c in 1:na
            LOCUSE_SD[c] = sum(LOCUSE_S[c,d] for d in 1:nr)
        end
        p["LOCUSE_S"] = LOCUSE_S
        p["LOCUSE_SD"] = LOCUSE_SD
        p["IMPSHR"] = IMPSHR

        # TAXRATE
        TAXRATE = zeros(T, na, ns, nu, nr)
        for c in 1:na, s in 1:ns, u in 1:nu, d in 1:nr
            TAXRATE[c,s,u,d] = USE[c,s,u,d] > 0 ? TAX[c,s,u,d] / USE[c,s,u,d] : 0.0
        end
        p["TAXRATE"] = TAXRATE

        # HOUPUR / BUDGSHR / EPS / SLUX / BLUX (Excerpt 13)
        u_hou = na + 1
        HOUPUR = zeros(T, na, 1, nr)
        for c in 1:na, d in 1:nr
            HOUPUR[c,1,d] = PUR_S[c,u_hou,d]
        end
        p["HOUPUR"] = HOUPUR

        HOUPUR_C = zeros(T, 1, nr)
        for d in 1:nr
            HOUPUR_C[1,d] = sum(HOUPUR[c,1,d] for c in 1:na)
        end
        p["HOUPUR_C"] = HOUPUR_C

        BUDGSHR = zeros(T, na, 1, nr)
        for c in 1:na, d in 1:nr
            BUDGSHR[c,1,d] = HOUPUR_C[1,d] > 0 ? HOUPUR[c,1,d] / HOUPUR_C[1,d] : 0.0
        end
        p["BUDGSHR"] = BUDGSHR

        EPSH = zeros(T, na, 1, nr)
        for c in 1:na, d in 1:nr
            EPSH[c,1,d] = XPEL[c,d]
        end
        EPSAVE = zeros(T, 1, nr)
        for d in 1:nr
            esum = sum(EPSH[c,1,d] * BUDGSHR[c,1,d] for c in 1:na)
            EPSAVE[1,d] = esum
            if esum > 0
                for c in 1:na
                    EPSH[c,1,d] /= esum
                end
            end
        end
        p["EPSH"] = EPSH
        p["EPSAVE"] = EPSAVE

        FRISCH = -abs.(P021_v)
        BLUX = zeros(T, na, 1, nr)
        SLUX = zeros(T, na, 1, nr)
        for c in 1:na, d in 1:nr
            BLUX[c,1,d] = abs(EPSH[c,1,d] / max(abs(FRISCH[d]), 1e-10))
            SLUX[c,1,d] = EPSH[c,1,d] * BUDGSHR[c,1,d]
        end
        p["BLUX"] = BLUX
        p["SLUX"] = SLUX

        # HOUSHR (single household: identity)
        HOUPUR_H = zeros(T, na, nr)
        HOUSHR = zeros(T, na, 1, nr)
        for c in 1:na, d in 1:nr
            HOUPUR_H[c,d] = HOUPUR[c,1,d]
            HOUSHR[c,1,d] = 1.0
        end
        p["HOUPUR_H"] = HOUPUR_H
        p["HOUSHR"] = HOUSHR

        # VMAINUSE (Excerpt 44)
        u_inv = na + 2
        u_gov = na + 3
        VMAINUSE = zeros(T, na, ns, 4, nr)
        for c in 1:na, s in 1:ns, d in 1:nr
            int_sum = sum(USE[c,s,i,d] for i in 1:na)
            VMAINUSE[c,s,1,d] = int_sum  # INT
            VMAINUSE[c,s,2,d] = USE[c,s,u_hou,d]  # HOU
            VMAINUSE[c,s,3,d] = USE[c,s,u_inv,d]  # INV
            VMAINUSE[c,s,4,d] = USE[c,s,u_gov,d]  # GOV
        end
        p["VMAINUSE"] = VMAINUSE

        # VARCST / VCST / VTOT / COSTMAT (Excerpt 12)
        VARCST = zeros(T, na, nr)
        VCST   = zeros(T, na, nr)
        for i in 1:na, d in 1:nr
            VARCST[i,d] = LAB_O[i,d] + PUR_CS[i,d]
            VCST[i,d]   = PRIM[i,d] + PUR_CS[i,d]
        end
        p["VARCST"] = VARCST
        p["VCST"]   = VCST

        VTOT = zeros(T, na, nr)
        PTXRATE = zeros(T, na, nr)
        COSTMAT_arr = zeros(T, na, 7, nr)
        for i in 1:na, d in 1:nr
            ptx_val = V1PTX !== nothing ? V1PTX[i,d] : 0.0
            VTOT[i,d] = VCST[i,d] + ptx_val
            PTXRATE[i,d] = VCST[i,d] > 0 ? ptx_val / VCST[i,d] : 0.0
            COSTMAT_arr[i,1,d] = sum(BSMR[c,1,i,d] for c in 1:na)  # IntDom
            COSTMAT_arr[i,2,d] = sum(BSMR[c,2,i,d] for c in 1:na)  # IntImp
            COSTMAT_arr[i,3,d] = sum(UTAX[c,s,i,d] for c in 1:na, s in 1:ns)  # ComTax
            COSTMAT_arr[i,4,d] = LAB_O[i,d]   # LAB
            COSTMAT_arr[i,5,d] = V1CAP[i,d]   # CAP
            COSTMAT_arr[i,6,d] = V1LND[i,d]   # LND
            COSTMAT_arr[i,7,d] = ptx_val      # PRODTAX
        end
        p["VTOT"] = VTOT
        p["PTXRATE"] = PTXRATE
        p["COSTMAT"] = COSTMAT_arr

        # NATVTOT / LAB_OD / CAP_D (Excerpt 35)
        NATVTOT = zeros(T, na)
        LAB_OD = zeros(T, na)
        CAP_D = zeros(T, na)
        for i in 1:na
            ptx_val_i = V1PTX !== nothing ? V1PTX[i,:] : zeros(T, nr)
            NATVTOT[i] = sum(VCST[i,d] + ptx_val_i[d] for d in 1:nr)
            LAB_OD[i]  = sum(LAB_O[i,d] for d in 1:nr)
            CAP_D[i]   = sum(V1CAP[i,d] for d in 1:nr)
        end
        p["NATVTOT"] = NATVTOT
        p["LAB_OD"]  = LAB_OD
        p["CAP_D"]   = CAP_D

        # INVEST_I / INVEST_C (Excerpt 14)
        if V2PUR !== nothing
            INVEST_I = zeros(T, na, nr)
            INVEST_C = zeros(T, na, nr)
            for c in 1:na, i in 1:na, d in 1:nr
                INVEST_I[c,d] += V2PUR[c,i,d]
                INVEST_C[i,d] += V2PUR[c,i,d]
            end
            p["INVEST_I"] = INVEST_I
            p["INVEST_C"] = INVEST_C
        end
    end

    # ── Factor aggregates (Excerpt 27) ─────────────────────────
    LAB_I = zeros(T, no, nr)
    for o in 1:no, d in 1:nr
        LAB_I[o,d] = sum(V1LAB[i,o,d] for i in 1:na)
    end
    p["LAB_I"] = LAB_I

    LAB_IO = zeros(T, nr)
    LND_I = zeros(T, nr)
    CAP_I = zeros(T, nr)
    for d in 1:nr
        LAB_IO[d] = sum(LAB_I[o,d] for o in 1:no)
        LND_I[d]  = sum(V1LND[i,d] for i in 1:na)
        CAP_I[d]  = sum(V1CAP[i,d] for i in 1:na)
    end
    p["LAB_IO"] = LAB_IO
    p["LND_I"]  = LND_I
    p["CAP_I"]  = CAP_I

    SLAB_I_arr = zeros(T, na, no, nr)
    for i in 1:na, o in 1:no, d in 1:nr
        SLAB_I_arr[i,o,d] = LAB_I[o,d] > 0 ? V1LAB[i,o,d] / LAB_I[o,d] : 0.0
    end
    p["SLAB_I"] = SLAB_I_arr

    # ── PRIM_I / PRIMSHR (Excerpt 25) ──────────────────────────
    PRIM_I_v = zeros(T, nr)
    for i in 1:na, d in 1:nr
        PRIM_I_v[d] += PRIM[i,d]
    end
    PRIMSHR_m = zeros(T, na, nr)
    for i in 1:na, d in 1:nr
        PRIMSHR_m[i,d] = PRIM_I_v[d] > 0 ? PRIM[i,d] / PRIM_I_v[d] : 0.0
    end
    p["PRIM_I"]   = PRIM_I_v
    p["PRIMSHR"]  = PRIMSHR_m
    PRIM_D_v = zeros(T, na)
    for i in 1:na, d in 1:nr
        PRIM_D_v[i] += PRIM[i,d]
    end
    p["PRIM_D"]   = PRIM_D_v

    # ── Regional macro aggregates (Excerpt 31) ─────────────────
    TRADE_CR = zeros(T, ns, nr)
    for s in 1:ns, d in 1:nr
        TRADE_CR[s,d] = sum(TRADE_R[c,s,d] for c in 1:na)
    end
    p["TRADE_CR"]   = TRADE_CR
    p["IMPUSED_C"]  = [TRADE_CR[2,d] for d in 1:nr]
    p["IMPLANDED_C"] = [sum(TRADE_D[c,2,d] for c in 1:na) for d in 1:nr]

    # ── GDP income & expenditure breakdowns (Excerpts 28-29) ───
    if haskey(p, "PUR_CS")
        PUR_CS_p = p["PUR_CS"]  # nu × nr
        # GDPINCCAT: Land, Labour, Capital, PRODTAX, ComTax
        GDPINCSUM = zeros(T, nr, 5)
        for d in 1:nr
            GDPINCSUM[d,1] = p["LND_I"][d]
            GDPINCSUM[d,2] = p["LAB_IO"][d]
            GDPINCSUM[d,3] = p["CAP_I"][d]
            GDPINCSUM[d,4] = V1PTX !== nothing ? sum(V1PTX[i,d] for i in 1:na) : 0.0
            if BSMR !== nothing && UTAX !== nothing
                tax_sum = sum(UTAX[c,s,u,d] for c in 1:na, s in 1:ns, u in 1:nu)
                GDPINCSUM[d,5] = tax_sum
            end
        end
        GDPINC = [sum(GDPINCSUM[d,:]) for d in 1:nr]
        p["GDPINCSUM"] = GDPINCSUM
        p["GDPINC"]    = GDPINC

        # GDPEXPCAT: HOU, INV, GOV, STOCKS, EXP, Imports, RExports, RImports, NetMar
        u_hou_i = na + 1
        u_inv_i = na + 2
        u_gov_i = na + 3
        u_exp_i = na + 4
        GDPEXPSUM = zeros(T, nr, 9)
        for d in 1:nr
            GDPEXPSUM[d,1] = PUR_CS_p[u_hou_i,d]  # HOU
            GDPEXPSUM[d,2] = PUR_CS_p[u_inv_i,d]  # INV
            GDPEXPSUM[d,3] = PUR_CS_p[u_gov_i,d]  # GOV
            GDPEXPSUM[d,4] = VSTOK !== nothing ? sum(VSTOK[i,d] for i in 1:na) : 0.0  # STOCKS
            GDPEXPSUM[d,5] = PUR_CS_p[u_exp_i,d]  # EXP
            GDPEXPSUM[d,6] = -sum(p["TRADE_D"][c,2,d] for c in 1:na)  # Imports
            GDPEXPSUM[d,7] = sum(c -> p["TRADE_D"][c,1,d] - p["TRDIAG"][c,1,d], 1:na)  # RExports
            GDPEXPSUM[d,8] = -sum(c -> p["TRADE_R"][c,1,d] - p["TRDIAG"][c,1,d], 1:na)  # RImports
            netmar = 0.0
            for m in 1:nm
                netmar += sum(p["SUPPMAR_D"][m,r,d] - p["SUPPMAR_P"][m,r,d] for r in 1:nr)
            end
            GDPEXPSUM[d,9] = netmar  # NetMar
        end
        GDPEXP = [sum(GDPEXPSUM[d,:]) for d in 1:nr]
        p["GDPEXPSUM"] = GDPEXPSUM
        p["GDPEXP"]    = GDPEXP
    end

    # ── LOCSHR (Excerpt 44) ────────────────────────────────────
    if haskey(p, "USE_U")
        LOCSHR_m = zeros(T, na, ns, nr)
        for c in 1:na, s in 1:ns, r in 1:nr
            uu = p["USE_U"][c,s,r]
            LOCSHR_m[c,s,r] = uu > 0 ? TRDIAG[c,s,r] / uu : 0.0
        end
        p["LOCSHR"] = LOCSHR_m

        LOCSHR_R_m = zeros(T, na, ns, nr)
        for c in 1:na, s in 1:ns, r in 1:nr
            uu = p["USE_U"][c,s,r]
            LOCSHR_R_m[c,s,r] = uu > 0 ? TRADE_R[c,s,r] / uu : 0.0
        end
        p["LOCSHR_R"] = LOCSHR_R_m
    end

    # ── ROWDEM / EXPSHR (Excerpt 40) ───────────────────────────
    if BSMR !== nothing && haskey(p, "USE_U")
        u_exp_ix = na + 4
        ROWDEM = zeros(T, na, nr, nr)
        for c in 1:na, r in 1:nr, d in 1:nr
            uu = p["USE_U"][c,1,d]
            if uu > 0
                ROWDEM[c,r,d] = TRADE[c,1,r,d] * BSMR[c,1,u_exp_ix,d] / uu
            end
        end
        p["ROWDEM"] = ROWDEM

        EXPSHR_arr = zeros(T, na, nr)
        for c in 1:na, r in 1:nr
            rowdem_sum = sum(ROWDEM[c,r,d] for d in 1:nr)
            EXPSHR_arr[c,r] = rowdem_sum / (0.001 + MAKE_I[c,r])
        end
        p["EXPSHR"] = EXPSHR_arr
    end

    # ── Accounting checks (Excerpt 41) ─────────────────────────
    if haskey(p, "VTOT")
        EP = 1e-5
        CHECKA = zeros(T, na, nr)
        CKRATA = zeros(T, na, nr)
        for i in 1:na, d in 1:nr
            stok_val = VSTOK !== nothing ? VSTOK[i,d] : 0.0
            vtot_val = p["VTOT"][i,d]
            CHECKA[i,d] = (vtot_val - stok_val) - MAKE_C[i,d]
            CKRATA[i,d] = (EP + vtot_val - stok_val) / (EP + MAKE_C[i,d])
        end
        p["CHECKA"] = CHECKA
        p["CKRATA"] = CKRATA
    end

    if haskey(p, "USE_U") && haskey(p, "DELIVRD_R")
        EP = 1e-5
        CHECKB = zeros(T, na, ns, nr)
        CKRATB = zeros(T, na, ns, nr)
        for c in 1:na, s in 1:ns, d in 1:nr
            uu = p["USE_U"][c,s,d]
            dr = p["DELIVRD_R"][c,s,d]
            CHECKB[c,s,d] = uu - dr
            CKRATB[c,s,d] = (EP + uu) / (EP + dr)
        end
        p["CHECKB"] = CHECKB
        p["CKRATB"] = CKRATB
    end

    return p
end
