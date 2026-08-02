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
    P015   = haskey(agg,"P015") ? vec(parent(agg["P015"])) : fill(5.0, na)
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
        acc = 0.0
        @inbounds for o in 1:no
            acc += V1LAB[i,o,d]
        end
        LAB_O[i,d] = acc
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

    # ── Levels-form CES calibration (Excerpt 10-11): labour composition
    # and labour/capital/land factor nests. Calibrated once here at the
    # benchmark (P=1, Q=value-flow convention) so build_model!.jl's CES
    # equations reproduce the benchmark data by construction — see
    # ces_calibrate() in ces_helper.jl and PLAN.md's "Course correction".
    ALPHA_LAB = zeros(T, na, no, nr)
    GAMMA_LAB = ones(T, na, nr)
    for i in 1:na, d in 1:nr
        αv, γv = ces_calibrate(V1LAB[i, :, d], SLAB[i], LAB_O[i, d])
        ALPHA_LAB[i, :, d] = αv
        GAMMA_LAB[i, d] = γv
    end
    p["ALPHA_LAB"] = ALPHA_LAB
    p["GAMMA_LAB"] = GAMMA_LAB

    # Factor nest order matches E_pprim!'s RHS: [LAB_O, CAP, LND]
    ALPHA_FAC = zeros(T, na, 3, nr)
    GAMMA_FAC = ones(T, na, nr)
    for i in 1:na, d in 1:nr
        αv, γv = ces_calibrate([LAB_O[i, d], V1CAP[i, d], V1LND[i, d]], P028[i], PRIM[i, d])
        ALPHA_FAC[i, :, d] = αv
        GAMMA_FAC[i, d] = γv
    end
    p["ALPHA_FAC"] = ALPHA_FAC
    p["GAMMA_FAC"] = GAMMA_FAC

    # ── MAKE/CET calibration (Excerpt 22): each multi-product industry i
    # transforms its total output xtot[i,d] (benchmark=1) into commodity-
    # specific outputs xmake[c,i,d] (benchmark=MAKE[c,i,d]) via a CET nest.
    # A negative sigma passed to ces_calibrate/ces gives the CET dual (see
    # ces_helper.jl), so SCET (a positive output-mix elasticity) is negated.
    ALPHA_MAKE = zeros(T, na, na, nr)   # ALPHA_MAKE[i, c, d]
    GAMMA_MAKE = ones(T, na, nr)
    for i in 1:na, d in 1:nr
        αv, γv = ces_calibrate(MAKE[:, i, d], -SCET[i], 1.0)
        ALPHA_MAKE[i, :, d] = αv
        GAMMA_MAKE[i, d] = γv
    end
    p["ALPHA_MAKE"] = ALPHA_MAKE
    p["GAMMA_MAKE"] = GAMMA_MAKE

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
        acc = 0.0
        @inbounds for d in 1:nr
            acc += MAKE[c,i,d]
        end
        MAKE_D[c,i] = acc
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
        tm_sum = 0.0
        @inbounds for m in 1:nm
            tm_sum += TMAR[c,s,m,r,d]
        end
        DELIVRD[c,s,r,d] = TRADE[c,s,r,d] + tm_sum
        dv = DELIVRD[c,s,r,d]
        BASSHR[c,s,r,d] = dv > 0 ? TRADE[c,s,r,d] / dv : 0.0
        for m in 1:nm
            MARSHR[c,s,m,r,d] = dv > 0 ? TMAR[c,s,m,r,d] / dv : 0.0
        end
    end
    p["MAKE"] = MAKE
    p["TRADE"] = TRADE
    p["TMAR"] = TMAR
    p["MARS"] = MARS
    p["DELIVRD"] = DELIVRD
    p["BASSHR"] = BASSHR
    p["MARSHR"] = MARSHR

    # ── DELIVRD_R (Excerpt 20) ─────────────────────────────────
    DELIVRD_R = zeros(T, na, ns, nr)
    for c in 1:na, s in 1:ns, d in 1:nr
        acc = 0.0
        @inbounds for r in 1:nr
            acc += DELIVRD[c,s,r,d]
        end
        DELIVRD_R[c,s,d] = acc
    end
    p["DELIVRD_R"] = DELIVRD_R
    p["SGDD"] = SGDD

    # ── Elasticities wired from national HAR data (Excerpt 6-27 sigmas) ────
    # Raw HAR header codes, matching the "SGDD" naming precedent above:
    #   SLAB=SIGMALAB (labour CES), P028=SIGMAPRIM (primary factor CES),
    #   P015=SIGMADOMIMP/ARMSIGMA (Armington dom/imp CES),
    #   SMAR=SIGMAMAR (margin substitution), PO01=POP (population),
    #   SCET=SIGMAOUT (CET output-mix), P018=EXP_ELAST (export demand)
    p["SLAB"] = SLAB
    p["P028"] = P028
    p["P015"] = P015
    p["SMAR"] = SMAR_v
    p["PO01"] = PO01
    p["SCET"] = SCET
    p["P018"] = P018

    # ── Margin aggregates (Excerpt 21) ─────────────────────────
    TRADMAR_CS = zeros(T, nm, nr, nr)
    for m in 1:nm, r in 1:nr, d in 1:nr
        acc = 0.0
        @inbounds for c in 1:na, s in 1:ns
            acc += TMAR[c,s,m,r,d]
        end
        TRADMAR_CS[m,r,d] = acc
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
        acc = 0.0
        @inbounds for r in 1:nr
            acc += SUPPMAR_D[m,r,p_]
        end
        SUPPMAR_RD[m,p_] = acc
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
            acc = 0.0
            @inbounds for c in 1:na
                acc += PUR_S[c,u,d]
            end
            PUR_CS[u,d] = acc
        end
        p["PUR_CS"] = PUR_CS

        SRCSHR = zeros(T, na, ns, nu, nr)
        for c in 1:na, s in 1:ns, u in 1:nu, d in 1:nr
            ps = PUR_S[c,u,d]
            SRCSHR[c,s,u,d] = ps > 0 ? PUR[c,s,u,d] / ps : 0.5
        end
        p["SRCSHR"] = SRCSHR

        # ── Levels-form CES calibration (Excerpt 8): Armington dom/imp
        # nest, one per (c,u,d) triple, using benchmark purchaser values
        # PUR[c,:,u,d] as the benchmark "quantities" (P=1 convention).
        # ces_calibrate returns (zeros, 1.0) wherever PUR_S==0 (no flow),
        # which the equation-writing code below skips via `any(αv .> 0)`.
        ALPHA_ARMINT = zeros(T, na, ns, nu, nr)
        GAMMA_ARMINT = ones(T, na, nu, nr)
        for c in 1:na, u in 1:nu, d in 1:nr
            αv, γv = ces_calibrate(PUR[c, :, u, d], P015[c], PUR_S[c, u, d])
            ALPHA_ARMINT[c, :, u, d] = αv
            GAMMA_ARMINT[c, u, d] = γv
        end
        p["ALPHA_ARMINT"] = ALPHA_ARMINT
        p["GAMMA_ARMINT"] = GAMMA_ARMINT

        PUR_D = zeros(T, na, ns, nu)
        for c in 1:na, s in 1:ns, u in 1:nu
            acc = 0.0
            @inbounds for d in 1:nr
                acc += PUR[c,s,u,d]
            end
            PUR_D[c,s,u] = acc
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
            acc = 0.0
            @inbounds for d in 1:nr
                acc += LOCUSE_S[c,d]
            end
            LOCUSE_SD[c] = acc
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
            acc = 0.0
            @inbounds for c in 1:na
                acc += HOUPUR[c,1,d]
            end
            HOUPUR_C[1,d] = acc
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
            esum = 0.0
            @inbounds for c in 1:na
                esum += EPSH[c,1,d] * BUDGSHR[c,1,d]
            end
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

        # XSUB0/XLUX0/WLUX0/ALUX0 (Excerpt 13, levels benchmark decomposition):
        # ELES splits benchmark household purchases HOUPUR(c,d) into a
        # "subsistence" part (1-BLUX)*HOUPUR and a "luxury"/supernumerary part
        # BLUX*HOUPUR. WLUX0(d) is the aggregate benchmark supernumerary
        # expenditure and ALUX0(c,d) is alux's true (non-unity) benchmark
        # share XLUX0(c,d)/WLUX0(d) — see build_equations.jl's E_xsub!/E_xlux!
        # header comment for the full derivation.
        XSUB0 = zeros(T, na, 1, nr)
        XLUX0 = zeros(T, na, 1, nr)
        for c in 1:na, d in 1:nr
            XSUB0[c,1,d] = (1.0 - BLUX[c,1,d]) * HOUPUR[c,1,d]
            XLUX0[c,1,d] = BLUX[c,1,d] * HOUPUR[c,1,d]
        end
        p["XSUB0"] = XSUB0
        p["XLUX0"] = XLUX0

        WLUX0 = zeros(T, 1, nr)
        for d in 1:nr
            acc = 0.0
            @inbounds for c in 1:na
                acc += XLUX0[c,1,d]
            end
            WLUX0[1,d] = acc
        end
        p["WLUX0"] = WLUX0

        ALUX0 = zeros(T, na, 1, nr)
        for c in 1:na, d in 1:nr
            ALUX0[c,1,d] = WLUX0[1,d] > 0 ? XLUX0[c,1,d] / WLUX0[1,d] : 0.0
        end
        p["ALUX0"] = ALUX0

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
            int_sum = 0.0
            @inbounds for i in 1:na
                int_sum += USE[c,s,i,d]
            end
            VMAINUSE[c,s,1,d] = int_sum  # INT
            VMAINUSE[c,s,2,d] = USE[c,s,u_hou,d]  # HOU
            VMAINUSE[c,s,3,d] = USE[c,s,u_inv,d]  # INV
            VMAINUSE[c,s,4,d] = USE[c,s,u_gov,d]  # GOV
        end
        p["VMAINUSE"] = VMAINUSE

        # XGOV0 / XEXPD0 (Excerpt 16 levels benchmarks)
        u_exp = na + 4
        XGOV0 = zeros(T, na, ns, nr)
        for c in 1:na, s in 1:ns, d in 1:nr
            XGOV0[c,s,d] = PUR[c,s,u_gov,d]
        end
        p["XGOV0"] = XGOV0

        XEXPD0 = zeros(T, na, nr)
        for c in 1:na, d in 1:nr
            acc = 0.0
            @inbounds for s in 1:ns
                acc += PUR[c,s,u_exp,d]
            end
            XEXPD0[c,d] = acc
        end
        p["XEXPD0"] = XEXPD0

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
            acc1 = 0.0; acc2 = 0.0; acc3 = 0.0
            @inbounds for c in 1:na
                acc1 += BSMR[c,1,i,d]
                acc2 += BSMR[c,2,i,d]
            end
            @inbounds for c in 1:na, s in 1:ns
                acc3 += UTAX[c,s,i,d]
            end
            COSTMAT_arr[i,1,d] = acc1  # IntDom
            COSTMAT_arr[i,2,d] = acc2  # IntImp
            COSTMAT_arr[i,3,d] = acc3  # ComTax
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
            acc1 = 0.0; acc2 = 0.0; acc3 = 0.0
            @inbounds for d in 1:nr
                acc1 += VCST[i,d] + ptx_val_i[d]
                acc2 += LAB_O[i,d]
                acc3 += V1CAP[i,d]
            end
            NATVTOT[i] = acc1
            LAB_OD[i]  = acc2
            CAP_D[i]   = acc3
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
        acc = 0.0
        @inbounds for i in 1:na
            acc += V1LAB[i,o,d]
        end
        LAB_I[o,d] = acc
    end
    p["LAB_I"] = LAB_I

    LAB_IO = zeros(T, nr)
    LND_I = zeros(T, nr)
    CAP_I = zeros(T, nr)
    for d in 1:nr
        acc1 = 0.0; acc2 = 0.0; acc3 = 0.0
        @inbounds for o in 1:no
            acc1 += LAB_I[o,d]
        end
        @inbounds for i in 1:na
            acc2 += V1LND[i,d]
            acc3 += V1CAP[i,d]
        end
        LAB_IO[d] = acc1
        LND_I[d]  = acc2
        CAP_I[d]  = acc3
    end
    p["LAB_IO"] = LAB_IO
    p["LND_I"]  = LND_I
    p["CAP_I"]  = CAP_I

    SLAB_I_arr = zeros(T, na, no, nr)
    for i in 1:na, o in 1:no, d in 1:nr
        SLAB_I_arr[i,o,d] = LAB_I[o,d] > 0 ? V1LAB[i,o,d] / LAB_I[o,d] : 0.0
    end
    p["SLAB_I"] = SLAB_I_arr

    # Sum-over-industry-and-region labour aggregates (Excerpt 27, "_id" family)
    LAB_ID = zeros(T, no)
    for o in 1:no
        acc = 0.0
        @inbounds for d in 1:nr
            acc += LAB_I[o,d]
        end
        LAB_ID[o] = acc
    end
    p["LAB_ID"] = LAB_ID

    SLAB_ID_arr = zeros(T, no, nr)
    for o in 1:no, d in 1:nr
        SLAB_ID_arr[o,d] = LAB_ID[o] > 0 ? LAB_I[o,d] / LAB_ID[o] : 0.0
    end
    p["SLAB_ID"] = SLAB_ID_arr

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
        acc = 0.0
        @inbounds for c in 1:na
            acc += TRADE_R[c,s,d]
        end
        TRADE_CR[s,d] = acc
    end
    p["TRADE_CR"]   = TRADE_CR
    p["IMPUSED_C"]  = [TRADE_CR[2,d] for d in 1:nr]
    IMPLANDED_C = zeros(T, nr)
    for d in 1:nr
        acc = 0.0
        @inbounds for c in 1:na
            acc += TRADE_D[c,2,d]
        end
        IMPLANDED_C[d] = acc
    end
    p["IMPLANDED_C"] = IMPLANDED_C

    # ── GDP income & expenditure breakdowns (Excerpts 28-29) ───
    if haskey(p, "PUR_CS")
        PUR_CS_p = p["PUR_CS"]  # nu × nr
        # GDPINCCAT: Land, Labour, Capital, PRODTAX, ComTax
        GDPINCSUM = zeros(T, nr, 5)
        for d in 1:nr
            GDPINCSUM[d,1] = p["LND_I"][d]
            GDPINCSUM[d,2] = p["LAB_IO"][d]
            GDPINCSUM[d,3] = p["CAP_I"][d]
            acc = 0.0
            if V1PTX !== nothing
                @inbounds for i in 1:na
                    acc += V1PTX[i,d]
                end
            end
            GDPINCSUM[d,4] = acc
            if BSMR !== nothing && UTAX !== nothing
                tax_sum = 0.0
                @inbounds for c in 1:na, s in 1:ns, u in 1:nu
                    tax_sum += UTAX[c,s,u,d]
                end
                GDPINCSUM[d,5] = tax_sum
            end
        end
        GDPINC = zeros(T, nr)
        for d in 1:nr
            acc = 0.0
            @inbounds for k in 1:5
                acc += GDPINCSUM[d,k]
            end
            GDPINC[d] = acc
        end
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
            acc = 0.0
            if VSTOK !== nothing
                @inbounds for i in 1:na
                    acc += VSTOK[i,d]
                end
            end
            GDPEXPSUM[d,4] = acc  # STOCKS
            GDPEXPSUM[d,5] = PUR_CS_p[u_exp_i,d]  # EXP
            acc = 0.0
            @inbounds for c in 1:na
                acc += p["TRADE_D"][c,2,d]
            end
            GDPEXPSUM[d,6] = -acc  # Imports
            acc = 0.0
            @inbounds for c in 1:na, s in 1:ns
                acc += p["TRADE_D"][c,s,d] - p["TRDIAG"][c,s,d]
            end
            GDPEXPSUM[d,7] = acc  # RExports
            acc = 0.0
            @inbounds for c in 1:na, s in 1:ns
                acc += p["TRADE_R"][c,s,d] - p["TRDIAG"][c,s,d]
            end
            GDPEXPSUM[d,8] = -acc  # RImports
            netmar = 0.0
            for m in 1:nm
                acc = 0.0
                @inbounds for r in 1:nr
                    acc += p["SUPPMAR_D"][m,r,d] - p["SUPPMAR_P"][m,r,d]
                end
                netmar += acc
            end
            GDPEXPSUM[d,9] = netmar  # NetMar
        end
        GDPEXP = zeros(T, nr)
        for d in 1:nr
            acc = 0.0
            @inbounds for k in 1:9
                acc += GDPEXPSUM[d,k]
            end
            GDPEXP[d] = acc
        end
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
            acc = 0.0
            @inbounds for d in 1:nr
                acc += ROWDEM[c,r,d]
            end
            EXPSHR_arr[c,r] = acc / (0.001 + MAKE_I[c,r])
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

    # ── Dynamic extension coefficients (Excerpts 50-51, 54) ────────────────
    # Only derived when the dynamic headers survived aggregation; the static
    # model never reads them. Names mirror TERM.TAB's own coefficient names,
    # except ALPHA → ALPHA_DYN to avoid colliding with the CES ALPHA_* shares.
    if haskey(agg, "STOC") && haskey(p, "INVEST_C")
        CAPSTOK = parent(agg["STOC"])
        DPRC    = parent(agg["DPRC"])
        RNORMAL = parent(agg["TARG"])
        GROTREND= parent(agg["TFRO"])
        QRATIO  = parent(agg["QRAT"])
        ALPHA_D = parent(agg["ALFA"])
        RORADJ  = parent(agg["RADJ"])
        INVEST_C = p["INVEST_C"]; CAPv = p["CAP"]

        p["CAPSTOK"] = CAPSTOK; p["DPRC"] = DPRC; p["RNORMAL"] = RNORMAL
        p["GROTREND"] = GROTREND; p["QRATIO"] = QRATIO
        p["ALPHA_DYN"] = ALPHA_D; p["RORADJ"] = RORADJ

        GROSSRET = zeros(T, na, nr)   # PK/PI
        GROSSGRO = zeros(T, na, nr)   # investment/capital ratio
        GROMAX   = zeros(T, na, nr)   # max investment/capital ratio
        GRETEXP  = zeros(T, na, nr)   # expected gross rate of return
        MCOEFF   = zeros(T, na, nr)   # mratio coefficient
        CAPADD   = zeros(T, na, nr)   # addition to CAPSTOK from last year's investment
        for i in 1:na, d in 1:nr
            cs = CAPSTOK[i,d]
            GROSSRET[i,d] = cs > 1e-10 ? CAPv[i,d] / cs : 0.0
            GROSSGRO[i,d] = cs > 1e-10 ? INVEST_C[i,d] / cs : 0.0
            GROMAX[i,d]   = QRATIO[i,d] * GROTREND[i,d]
            # TERM.TAB floors DENOM at 0.001 to avoid raising a negative number
            # to a fractional power (its own comment, Excerpt 51).
            den = GROSSGRO[i,d] > 1e-10 ? (GROMAX[i,d] / GROSSGRO[i,d]) - 1.0 : 0.001
            den <= 0 && (den = 0.001)
            GRETEXP[i,d] = ALPHA_D[i,d] > 1e-10 ?
                RNORMAL[i,d] * ((QRATIO[i,d] - 1.0) / den)^(1.0 / ALPHA_D[i,d]) : 0.0
            MCOEFF[i,d] = GROMAX[i,d] > 1e-10 ?
                ALPHA_D[i,d] * (1.0 - GROSSGRO[i,d] / GROMAX[i,d]) : 0.0
            CAPADD[i,d] = INVEST_C[i,d] - DPRC[i,d] * cs
        end
        p["GROSSRET"] = GROSSRET; p["GROSSRET0"] = copy(GROSSRET)
        p["GROSSGRO"] = GROSSGRO; p["GROMAX"] = GROMAX
        p["GRETEXP"]  = GRETEXP;  p["GRETEXP0"] = copy(GRETEXP)
        p["MCOEFF"]   = MCOEFF;   p["CAPADD"]  = CAPADD
        p["CAPSTOK_OLDP"] = copy(CAPSTOK)
        CAPSTOK_D = zeros(T, na)
        INVEST_CD = zeros(T, na)
        for i in 1:na
            acc1 = 0.0; acc2 = 0.0
            @inbounds for d in 1:nr
                acc1 += CAPSTOK[i,d]
                acc2 += INVEST_C[i,d]
            end
            CAPSTOK_D[i] = acc1
            INVEST_CD[i] = acc2
        end
        p["CAPSTOK_D"] = CAPSTOK_D
        p["INVEST_CD"] = INVEST_CD
    end
    if haskey(agg, "EMPR")
        EMPRAT = parent(agg["EMPR"]); ELASTWAGE = parent(agg["ELWG"])
        # TERM.TAB indexes both by OCC only; the pipeline carries them OCC×REG,
        # so collapse to OCC with a plain mean (they are intensive rates).
        p["EMPRAT"]    = ndims(EMPRAT)    == 2 ? vec(sum(EMPRAT, dims=2))    ./ nr : EMPRAT
        p["ELASTWAGE"] = ndims(ELASTWAGE) == 2 ? vec(sum(ELASTWAGE, dims=2)) ./ nr : ELASTWAGE
        p["EMPRAT0"]   = copy(p["EMPRAT"])
        p["WAGERATE"]  = ones(T, length(p["EMPRAT"]))   # index rebased each period
    end

    return p
end
