function build_model!(agg::Dict{String,Any}, params::Dict{String,Any})
    T = Float64
    na = size(agg["MAKE"], 1)
    nr = size(agg["MAKE"], 3)
    ns = 2
    nm = haskey(agg, "TMAR") ? size(agg["TMAR"], 3) : 9
    no = size(agg["1LAB"], 2)
    nu = na + 4

    MAKE  = parent(agg["MAKE"])
    TRADE = parent(agg["TRAD"])
    TMAR  = haskey(agg,"TMAR") ? parent(agg["TMAR"]) : zeros(T, na, ns, nm, nr, nr)
    MARS  = haskey(agg,"MARS") ? parent(agg["MARS"]) : zeros(T, nm, nr, nr, nr)
    DIST  = haskey(agg,"DIST") ? parent(agg["DIST"]) : zeros(T, nr, nr)
    V1LAB = parent(agg["1LAB"])
    V1CAP = parent(agg["1CAP"])
    V1LND = parent(agg["1LND"])
    V1PTX = haskey(agg,"1PTX") ? parent(agg["1PTX"]) : zeros(T, na, nr)

    # Derived parameters
    LAB_O  = parent(params["LAB_O"])
    PRIM   = parent(params["PRIM"])
    PUR_S  = parent(params["PUR_S"])
    PUR_CS = parent(params["PUR_CS"])
    SRCSHR = parent(params["SRCSHR"])
    DELIVRD = parent(params["DELIVRD"])
    DELIVRD_R = parent(params["DELIVRD_R"])
    BASSHR = parent(params["BASSHR"])
    MARSHR = parent(params["MARSHR"])
    MAKE_C = parent(params["MAKE_C"])
    MAKE_I = parent(params["MAKE_I"])
    TRADE_D = parent(params["TRADE_D"])
    TRADE_R = parent(params["TRADE_R"])
    LAB_I   = parent(params["LAB_I"])
    LAB_IO  = parent(params["LAB_IO"])
    CAP_I   = parent(params["CAP_I"])
    LND_I   = parent(params["LND_I"])
    PRIM_I  = parent(params["PRIM_I"])
    VCST    = parent(params["VCST"])
    VARCST  = parent(params["VARCST"])
    VTOT    = parent(params["VTOT"])
    USE_U   = parent(params["USE_U"])
    USE_I   = parent(params["USE_I"])
    BUDGSHR = parent(params["BUDGSHR"])
    BLUX    = parent(params["BLUX"])
    SLUX    = parent(params["SLUX"])
    LOCSHR  = parent(params["LOCSHR"])
    HOUPUR  = parent(params["HOUPUR"])
    INVEST_C = parent(params["INVEST_C"])
    INVEST_I = parent(params["INVEST_I"])

    # Tax / margin data
    USE_d  = zeros(T, na, ns, nu, nr)
    TAX_d  = zeros(T, na, ns, nu, nr)
    if haskey(agg, "BSMR")
        bsmr = parent(agg["BSMR"])
        szu = min(size(bsmr,3), nu)
        for c in 1:na, s in 1:ns, u in 1:szu, d in 1:nr
            USE_d[c,s,u,d] = bsmr[c,s,u,d]
        end
    end
    if haskey(agg, "UTAX")
        utax = parent(agg["UTAX"])
        szu = min(size(utax,3), nu)
        for c in 1:na, s in 1:ns, u in 1:szu, d in 1:nr
            TAX_d[c,s,u,d] = utax[c,s,u,d]
        end
    end

    PUR_d  = USE_d
    SUPP_d = zeros(T, nm, nr, nr, nr)

    # ── JuMP Model ────────────────────────────────────────────────
    m = JuMP.Model(Ipopt.Optimizer)
    JuMP.set_silent(m)
    JuMP.set_attribute(m, "tol", 1.0e-6)
    JuMP.set_attribute(m, "print_level", 0)

    # ── %-change variables (lowercase ≡ 100 × d(log(X))) ──────────
    # Prices (Excerpt 4, 6)
    @variable(m, pdom[1:na, 1:nr])
    @variable(m, pimp[1:na, 1:nr])
    @variable(m, phi)
    @variable(m, pfimp[1:na])
    @variable(m, pbasic[1:na, 1:ns, 1:nr])

    # Purchaser prices (Excerpt 7)
    @variable(m, puse[1:na, 1:ns, 1:nr])
    @variable(m, tuser[1:na, 1:ns, 1:nu, 1:nr])
    @variable(m, tuser_ud[1:na, 1:ns])
    @variable(m, tuser_su[1:na, 1:nr])
    @variable(m, tuser_sud[1:na])
    @variable(m, ppur[1:na, 1:ns, 1:nu, 1:nr])
    @variable(m, ppur_s[1:na, 1:nu, 1:nr])

    # Armington (Excerpt 8)
    @variable(m, phou[1:na, 1:nr])
    @variable(m, pinvest[1:na, 1:nr])
    @variable(m, xint[1:na, 1:ns, 1:na, 1:nr])
    @variable(m, xhou[1:na, 1:ns, 1:nr])
    @variable(m, xinv[1:na, 1:ns, 1:nr])
    @variable(m, xint_s[1:na, 1:na, 1:nr])
    @variable(m, xhou_s[1:na, 1:nr])
    @variable(m, xinv_s[1:na, 1:nr])

    # Intermediate (Excerpt 9)
    @variable(m, atot[1:na, 1:nr])
    @variable(m, aint_s[1:na, 1:na, 1:nr])
    @variable(m, bint_scd[1:na])
    @variable(m, bint_s[1:na, 1:na, 1:nr])
    @variable(m, pint[1:na, 1:nr])

    # Labour (Excerpt 10)
    @variable(m, xlab[1:na, 1:no, 1:nr])
    @variable(m, plab[1:na, 1:no, 1:nr])
    @variable(m, xlab_o[1:na, 1:nr])
    @variable(m, plab_o[1:na, 1:nr])
    @variable(m, wlab_o[1:na, 1:nr])

    # Factor demands (Excerpt 11)
    @variable(m, xcap[1:na, 1:nr])
    @variable(m, pcap[1:na, 1:nr])
    @variable(m, xlnd[1:na, 1:nr])
    @variable(m, plnd[1:na, 1:nr])
    @variable(m, xprim[1:na, 1:nr])
    @variable(m, pprim[1:na, 1:nr])
    @variable(m, alab_o[1:na, 1:nr])
    @variable(m, acap[1:na, 1:nr])
    @variable(m, alnd[1:na, 1:nr])
    @variable(m, aprim[1:na, 1:nr])
    @variable(m, wprim[1:na, 1:nr])
    @variable(m, bprimnat)
    @variable(m, bprim_d[1:na])
    @variable(m, bprim[1:na, 1:nr])
    @variable(m, blabnat)
    @variable(m, blab_d[1:na])
    @variable(m, blab[1:na, 1:nr])

    # Output (Excerpt 12)
    @variable(m, xtot[1:na, 1:nr])
    @variable(m, ptot[1:na, 1:nr])
    @variable(m, pvar[1:na, 1:nr])
    @variable(m, pcst[1:na, 1:nr])
    @variable(m, delPTX[1:na, 1:nr])
    @variable(m, delPTXRATE[1:na, 1:nr])

    # Trade variables (Excerpt 4, 17, 19-23)
    @variable(m, xtrad[1:na, 1:ns, 1:nr, 1:nr])
    @variable(m, xtrad_d[1:na, 1:ns, 1:nr])
    @variable(m, xtrad_r[1:na, 1:ns, 1:nr])
    @variable(m, xuse[1:na, 1:ns, 1:nr])
    @variable(m, xint_i[1:na, 1:ns, 1:nr])
    @variable(m, pdelivrd[1:na, 1:ns, 1:nr, 1:nr])
    @variable(m, xcom[1:na, 1:nr])
    @variable(m, xmake[1:na, 1:na, 1:nr])
    @variable(m, pmake[1:na, 1:na, 1:nr])
    @variable(m, xtradmar[1:na, 1:ns, 1:nm, 1:nr, 1:nr])
    @variable(m, xsuppmar[1:nm, 1:nr, 1:nr, 1:nr])
    @variable(m, xsuppmar_p[1:nm, 1:nr, 1:nr])
    @variable(m, psuppmar_p[1:nm, 1:nr, 1:nr])
    @variable(m, xsuppmar_d[1:nm, 1:nr, 1:nr])
    @variable(m, xsuppmar_rd[1:nm, 1:nr])
    @variable(m, atrad[1:na, 1:ns, 1:nr, 1:nr])
    @variable(m, atradmar[1:na, 1:ns, 1:nm, 1:nr, 1:nr])
    @variable(m, asuppmar[1:nm, 1:nr, 1:nr, 1:nr])
    @variable(m, srctwist[1:na, 1:ns, 1:nr, 1:nr])
    @variable(m, avesrctwist[1:na, 1:ns, 1:nr])

    # Household (Excerpt 13) — single representative household (h=1)
    @variable(m, nhou[1:nr])
    @variable(m, xhoutot[1:nr])
    @variable(m, phoutot[1:nr])
    @variable(m, xhouhtot[1:nr])
    @variable(m, phouhtot[1:nr])
    @variable(m, whouhtot[1:nr])
    @variable(m, xlux[1:na, 1:nr])
    @variable(m, xsub[1:na, 1:nr])
    @variable(m, wlux[1:nr])
    @variable(m, alux[1:na, 1:nr])
    @variable(m, asub[1:na, 1:nr])
    @variable(m, ahou_s[1:na, 1:nr])

    # Investment / Government / Export (Excerpt 14-16)
    @variable(m, xinvitot[1:na, 1:nr])
    @variable(m, pinvitot[1:na, 1:nr])
    @variable(m, xgov[1:na, 1:ns, 1:nr])
    @variable(m, xgov_s[1:na, 1:nr])
    @variable(m, fgov[1:na, 1:ns, 1:nr])
    @variable(m, fgov_s[1:na, 1:nr])
    @variable(m, fgovtot[1:nr])
    @variable(m, fgovtot2[1:nr])
    @variable(m, fgovtot3[1:nr])
    @variable(m, fgovgen)
    @variable(m, xexp[1:na, 1:ns, 1:nr])
    @variable(m, xexp_s[1:na, 1:nr])
    @variable(m, xexpd[1:na, 1:nr])
    @variable(m, pfexp[1:na, 1:nr])
    @variable(m, fqexp[1:na, 1:nr])
    @variable(m, fpexp[1:na, 1:nr])
    @variable(m, fqexp_d[1:na])
    @variable(m, fpexp_d[1:na])
    @variable(m, natfpexp)
    @variable(m, natfqexp)
    @variable(m, xstocks[1:na, 1:nr])
    @variable(m, fxstocks[1:na, 1:nr])
    @variable(m, xinvi[1:na, 1:na, 1:nr])

    # Investment rule variables (Excerpt 15)
    @variable(m, gret[1:na, 1:nr])
    @variable(m, ggro[1:na, 1:nr])
    @variable(m, finv1[1:na, 1:nr])
    @variable(m, finv2[1:na, 1:nr])
    @variable(m, invslack)
    @variable(m, capslack)
    @variable(m, fgret[1:na, 1:nr])

    # Labour market / closures (Excerpt 27, 38)
    @variable(m, xlab_i[1:no, 1:nr])
    @variable(m, plab_i[1:no, 1:nr])
    @variable(m, wlab_i[1:no, 1:nr])
    @variable(m, rlab_i[1:no, 1:nr])
    @variable(m, xlab_io[1:nr])
    @variable(m, plab_io[1:nr])
    @variable(m, wlab_io[1:nr])
    @variable(m, realwage[1:na, 1:no, 1:nr])
    @variable(m, realwage_i[1:no, 1:nr])
    @variable(m, realwage_io[1:nr])
    @variable(m, flabsupA[1:no, 1:nr])
    @variable(m, labslack[1:no])
    @variable(m, flab_i[1:no, 1:nr])
    @variable(m, flab[1:na, 1:no, 1:nr])
    @variable(m, flab_io[1:nr])
    @variable(m, flab_iod)
    @variable(m, flab_id[1:no])
    @variable(m, flabsup_id[1:no])
    @variable(m, xlab_id[1:no])
    @variable(m, plab_id[1:no])
    @variable(m, wlab_id[1:no])
    @variable(m, rlab_id[1:no])
    @variable(m, realwage_id[1:no])
    @variable(m, rlab_io[1:nr])
    @variable(m, xcap_i[1:nr])
    @variable(m, wcap_i[1:nr])
    @variable(m, xlnd_i[1:nr])
    @variable(m, wlnd_i[1:nr])
    @variable(m, wprim_i[1:nr])

    # Final demand / GDP (Excerpt 24-30)
    @variable(m, pfin[1:4, 1:nr])      # HOU=1, INV=2, GOV=3, EXP=4
    @variable(m, xfin[1:4, 1:nr])
    @variable(m, wfin[1:4, 1:nr])
    @variable(m, delXGDPEXP[1:nr, 1:9]) # 9 GDPEXPCAT categories
    @variable(m, delPGDPEXP[1:nr, 1:9])
    @variable(m, delVGDPEXP[1:nr, 1:9])
    @variable(m, xgdpexp[1:nr])
    @variable(m, pgdpexp[1:nr])
    @variable(m, wgdpexp[1:nr])
    @variable(m, wgdpinc[1:nr])
    @variable(m, wgdpdiff[1:nr])
    @variable(m, xgne[1:nr])
    @variable(m, pgne[1:nr])
    @variable(m, wgne[1:nr])
    @variable(m, delGDPINC[1:nr, 1:5])  # Land, Capital, Labour, ProdTax, ComTax
    @variable(m, delINDTAX[1:nr])
    @variable(m, delBUDG1[1:nr])
    @variable(m, delBUDG2[1:nr])

    # Household closure (Excerpt 39)
    @variable(m, fhou[1:nr])
    @variable(m, fhou2[1:nr])
    @variable(m, houslack)
    @variable(m, natfhou)

    # Tax revenue change variables (Excerpt 4, 26)
    @variable(m, delTAXint[1:na, 1:ns, 1:na, 1:nr])
    @variable(m, delTAXhou[1:na, 1:ns, 1:nr])
    @variable(m, delTAXinv[1:na, 1:ns, 1:nr])
    @variable(m, delTAXgov[1:na, 1:ns, 1:nr])
    @variable(m, delTAXexp[1:na, 1:ns, 1:nr])

    # ── Equation Blocks ──────────────────────────────────────────
    # Helper: SIGMADOMIMP (default 5.0)
    sigmadomimp = fill(5.0, na)
    sigmalab = fill(0.5, na)
    sigmaprim = fill(0.5, na)
    sigmaout = fill(0.5, na)
    exp_elast = fill(2.0, na)

    # ... equation blocks will be added via functions below

    vars = Dict{String,Any}()
    for (nm, vl) in [
        ("pdom", pdom), ("pimp", pimp), ("phi", phi), ("pfimp", pfimp),
        ("pbasic", pbasic), ("puse", puse), ("ppur", ppur), ("ppur_s", ppur_s),
        ("tuser", tuser), ("tuser_ud", tuser_ud), ("tuser_su", tuser_su),
        ("tuser_sud", tuser_sud), ("phou", phou), ("pinvest", pinvest),
        ("xint", xint), ("xhou", xhou), ("xinv", xinv),
        ("xint_s", xint_s), ("xhou_s", xhou_s), ("xinv_s", xinv_s),
        ("atot", atot), ("aint_s", aint_s), ("bint_scd", bint_scd), ("bint_s", bint_s),
        ("pint", pint),
        ("xlab", xlab), ("plab", plab), ("xlab_o", xlab_o), ("plab_o", plab_o),
        ("wlab_o", wlab_o),
        ("xcap", xcap), ("pcap", pcap), ("xlnd", xlnd), ("plnd", plnd),
        ("xprim", xprim), ("pprim", pprim),
        ("alab_o", alab_o), ("acap", acap), ("alnd", alnd), ("aprim", aprim),
        ("wprim", wprim),
        ("bprimnat", bprimnat), ("bprim_d", bprim_d), ("bprim", bprim),
        ("blabnat", blabnat), ("blab_d", blab_d), ("blab", blab),
        ("xtot", xtot), ("ptot", ptot), ("pvar", pvar), ("pcst", pcst),
        ("delPTX", delPTX), ("delPTXRATE", delPTXRATE),
        ("xtrad", xtrad), ("xtrad_d", xtrad_d), ("xtrad_r", xtrad_r),
        ("xuse", xuse), ("xint_i", xint_i),
        ("pdelivrd", pdelivrd),
        ("xcom", xcom), ("xmake", xmake), ("pmake", pmake),
        ("xtradmar", xtradmar), ("xsuppmar", xsuppmar),
        ("xsuppmar_p", xsuppmar_p), ("psuppmar_p", psuppmar_p),
        ("xsuppmar_d", xsuppmar_d), ("xsuppmar_rd", xsuppmar_rd),
        ("atrad", atrad), ("atradmar", atradmar), ("asuppmar", asuppmar),
        ("srctwist", srctwist), ("avesrctwist", avesrctwist),
        ("nhou", nhou), ("xhoutot", xhoutot), ("phoutot", phoutot),
        ("xhouhtot", xhouhtot), ("phouhtot", phouhtot), ("whouhtot", whouhtot),
        ("xlux", xlux), ("xsub", xsub), ("wlux", wlux),
        ("alux", alux), ("asub", asub), ("ahou_s", ahou_s),
        ("xinvitot", xinvitot), ("pinvitot", pinvitot),
        ("xgov", xgov), ("xgov_s", xgov_s),
        ("fgov", fgov), ("fgov_s", fgov_s), ("fgovtot", fgovtot),
        ("fgovtot2", fgovtot2), ("fgovtot3", fgovtot3), ("fgovgen", fgovgen),
        ("xexp", xexp), ("xexp_s", xexp_s), ("xexpd", xexpd),
        ("pfexp", pfexp), ("fqexp", fqexp), ("fpexp", fpexp),
        ("fqexp_d", fqexp_d), ("fpexp_d", fpexp_d),
        ("natfpexp", natfpexp), ("natfqexp", natfqexp),
        ("xstocks", xstocks), ("fxstocks", fxstocks),
        ("xinvi", xinvi),
        ("xlab_i", xlab_i), ("plab_i", plab_i),
        ("wlab_i", wlab_i), ("rlab_i", rlab_i),
        ("xlab_io", xlab_io), ("plab_io", plab_io),
        ("wlab_io", wlab_io), ("rlab_io", rlab_io),
        ("realwage", realwage), ("realwage_i", realwage_i),
        ("realwage_io", realwage_io),
        ("flabsupA", flabsupA), ("labslack", labslack),
        ("flab_i", flab_i), ("flab", flab),
        ("flab_io", flab_io), ("flab_iod", flab_iod),
        ("flab_id", flab_id), ("flabsup_id", flabsup_id),
        ("xlab_id", xlab_id), ("plab_id", plab_id),
        ("wlab_id", wlab_id), ("rlab_id", rlab_id),
        ("realwage_id", realwage_id),
        ("xcap_i", xcap_i), ("wcap_i", wcap_i),
        ("xlnd_i", xlnd_i), ("wlnd_i", wlnd_i),
        ("wprim_i", wprim_i),
        ("gret", gret), ("ggro", ggro),
        ("finv1", finv1), ("finv2", finv2),
        ("invslack", invslack), ("capslack", capslack), ("fgret", fgret),
        ("pfin", pfin), ("xfin", xfin), ("wfin", wfin),
        ("delXGDPEXP", delXGDPEXP), ("delPGDPEXP", delPGDPEXP),
        ("delVGDPEXP", delVGDPEXP),
        ("xgdpexp", xgdpexp), ("pgdpexp", pgdpexp),
        ("wgdpexp", wgdpexp), ("wgdpinc", wgdpinc),
        ("wgdpdiff", wgdpdiff),
        ("xgne", xgne), ("pgne", pgne), ("wgne", wgne),
        ("delGDPINC", delGDPINC),
        ("delINDTAX", delINDTAX), ("delBUDG1", delBUDG1), ("delBUDG2", delBUDG2),
        ("fhou", fhou), ("fhou2", fhou2), ("houslack", houslack), ("natfhou", natfhou),
        ("delTAXint", delTAXint), ("delTAXhou", delTAXhou),
        ("delTAXinv", delTAXinv), ("delTAXgov", delTAXgov), ("delTAXexp", delTAXexp)]
        vars[nm] = vl
    end

    return m, vars
end

function build_model_full!(agg, params)
    m, vars = build_model!(agg, params)
    na = size(agg["MAKE"], 1); nr = size(agg["MAKE"], 3)
    ns = 2; no = size(agg["1LAB"], 2); nm = haskey(agg, "TMAR") ? size(agg["TMAR"], 3) : 9; nu = na + 4

    MAKE  = parent(agg["MAKE"]); TRADE = parent(agg["TRAD"])
    V1LAB = parent(agg["1LAB"]); V1CAP = parent(agg["1CAP"])
    V1LND = parent(agg["1LND"]); V1PTX = haskey(agg,"1PTX") ? parent(agg["1PTX"]) : zeros(Float64, na, nr)
    DIST  = haskey(agg,"DIST") ? parent(agg["DIST"]) : zeros(Float64, nr, nr)

    LAB_O  = parent(params["LAB_O"]); PRIM   = parent(params["PRIM"])
    PUR_S  = parent(params["PUR_S"]); PUR_CS = parent(params["PUR_CS"])
    SRCSHR = parent(params["SRCSHR"]); DELIVRD = parent(params["DELIVRD"])
    DELIVRD_R = parent(params["DELIVRD_R"]); BASSHR = parent(params["BASSHR"])
    MARSHR = parent(params["MARSHR"]); MAKE_C = parent(params["MAKE_C"])
    TRADE_D = parent(params["TRADE_D"]); TRADE_R = parent(params["TRADE_R"])

    sigmadomimp = fill(5.0, na); sigmalab = fill(0.5, na); sigmaprim = fill(0.5, na)
    sigmaout = fill(0.5, na); exp_elast = fill(2.0, na)

    # ── Excerpt 6: Basic prices ──────────────────────────────────
    E_pimp!(m, vars, na, nr, params) # pimp = pfimp + phi
    E_pbasic!(m, vars, na, nr, ns)   # pbasic = dom·pdom + imp·pimp

    # ── Excerpt 7: Purchaser prices ──────────────────────────────
    E_ppur!(m, vars, na, nr, ns, nu)
    E_tuser!(m, vars, na, nr, ns, nu)

    # ── Excerpt 8: Armington ─────────────────────────────────────
    E_ppur_s!(m, vars, na, nr, ns, nu, params)
    E_phou!(m, vars, na, nr)
    E_xint!(m, vars, na, nr, ns, sigmadomimp)
    E_xhou!(m, vars, na, nr, ns, sigmadomimp)
    E_xinv!(m, vars, na, nr, ns, sigmadomimp)

    # ── Excerpt 9: Intermediate ──────────────────────────────────
    E_aint_s!(m, vars, na, nr)
    E_xint_s!(m, vars, na, nr)
    E_pint!(m, vars, na, nr, params)

    # ── Excerpt 10: Labour ───────────────────────────────────────
    E_xlab!(m, vars, na, nr, no, sigmalab)
    E_plab_o!(m, vars, na, nr, no, params)
    E_wlab_o!(m, vars, na, nr, no, params)

    # ── Excerpt 11: Factor demands ───────────────────────────────
    E_xlab_o!(m, vars, na, nr, sigmaprim)
    E_pcap!(m, vars, na, nr, sigmaprim)
    E_plnd!(m, vars, na, nr, sigmaprim)
    E_pprim!(m, vars, na, nr, no, params)
    E_xprim!(m, vars, na, nr)
    E_aprim!(m, vars, na, nr)
    E_alab_o!(m, vars, na, nr)
    E_wprim!(m, vars, na, nr, no, params)

    # ── Excerpt 12: Output prices ────────────────────────────────
    E_pvar!(m, vars, na, nr, params)
    E_pcst!(m, vars, na, nr, params)
    E_delPTX!(m, vars, na, nr, params)
    E_ptot!(m, vars, na, nr, params)

    # ── Excerpt 13: Household ────────────────────────────────────
    E_xsub!(m, vars, na, nr)
    E_xlux!(m, vars, na, nr)
    E_xhouh_s_agg!(m, vars, na, nr, params)
    E_alux!(m, vars, na, nr, params)
    E_asub!(m, vars, na, nr, params)
    E_wlux!(m, vars, na, nr, params)
    E_phouhtot!(m, vars, na, nr, params)
    E_whouhtot!(m, vars, na, nr)
    E_xhoutot!(m, vars, na, nr, params)
    E_phoutot!(m, vars, na, nr, params)

    # ── Excerpt 14: Investment demands ───────────────────────────
    E_xinvi!(m, vars, na, nr, params)
    E_pinvest!(m, vars, na, nr, nu)
    E_pinvitot!(m, vars, na, nr, params)

    # ── Excerpt 15: Investment rule ──────────────────────────────
    E_gret!(m, vars, na, nr)
    E_xinvitot!(m, vars, na, nr)
    E_finv2!(m, vars, na, nr)

    # ── Excerpt 16: Government / Export / Stocks ─────────────────
    E_xgov!(m, vars, na, nr, ns)
    E_xgov_s!(m, vars, na, nr, ns, params)
    E_fgovtot2!(m, vars, na, nr)
    E_fgovtot3!(m, vars, na, nr)
    E_pfexp!(m, vars, na, nr)
    E_xexpd!(m, vars, na, nr, exp_elast)
    E_xexp!(m, vars, na, nr, ns)
    E_xexp_s!(m, vars, na, nr, ns, params)
    E_xstocks!(m, vars, na, nr)

    # ── Excerpt 17: Total regional demand ────────────────────────
    E_xint_i!(m, vars, na, nr, ns, params)
    E_xuse!(m, vars, na, nr, ns, params)

    # ── Excerpt 19: Delivering goods (margins) ────────────────────
    E_pdelivrd!(m, vars, na, nr, ns, nm, params)
    E_xtradmar_na!(m, vars, na, nr, ns, nm)

    # ── Excerpt 20: Regional sourcing ────────────────────────────
    E_puse!(m, vars, na, nr, ns, params)
    E_avesrctwist!(m, vars, na, nr, ns, params)
    E_xtrad!(m, vars, na, nr, ns, params)

    # ── Excerpt 21: Margin supply ────────────────────────────────
    E_xsuppmar_p!(m, vars, na, nr, ns, nm, params)
    E_psuppmar_p!(m, vars, na, nr, nm, params)
    E_xsuppmar!(m, vars, na, nr, nm)
    E_xsuppmar_d!(m, vars, na, nr, ns, nm, params)
    E_xsuppmar_rd!(m, vars, na, nr, ns, nm, params)

    # ── Excerpt 22: MAKE / CET ───────────────────────────────────
    E_xmake!(m, vars, na, nr, sigmaout, params)
    E_xtotA_B!(m, vars, na, nr, params)
    E_xcomA_B!(m, vars, na, nr, params)
    E_pmake!(m, vars, na, nr)

    # ── Excerpt 23: Market clearing ──────────────────────────────
    E_xtrad_d!(m, vars, na, nr, ns, params)
    E_xtrad_r!(m, vars, na, nr, ns, params)
    E_pdomA_sum!(m, vars, na, nr)

    # ── Excerpt 24: Final demand aggregates ──────────────────────
    E_pfin!(m, vars, na, nr, nu, params)
    E_xfina!(m, vars, na, nr)
    E_xfinb!(m, vars, na, nr, params)
    E_xfinc!(m, vars, na, nr, params)
    E_xfind!(m, vars, na, nr, ns, params)
    E_wfin!(m, vars, nr)

    # ── Excerpt 26: Commodity tax revenues ───────────────────────
    E_delTAXint!(m, vars, na, nr, ns, params)
    E_delTAXhou!(m, vars, na, nr, ns, params)
    E_delTAXinv!(m, vars, na, nr, ns, params)
    E_delTAXgov!(m, vars, na, nr, ns, params)
    E_delTAXexp!(m, vars, na, nr, ns, params)

    # ── Excerpt 27: Factor aggregates ────────────────────────────
    E_wlnd_i!(m, vars, na, nr, params)
    E_wcap_i!(m, vars, na, nr, params)
    E_wprim_i!(m, vars, na, nr, params)
    E_xlnd_i!(m, vars, na, nr, params)
    E_xcap_i!(m, vars, na, nr, params)
    E_plab!(m, vars, na, nr, no)
    E_plab_i!(m, vars, na, nr, no, params)
    E_realwage_i!(m, vars, na, nr, no, params)
    E_xlab_i!(m, vars, na, nr, no, params)
    E_wlab_i!(m, vars, na, nr, no)
    E_rlab_i!(m, vars, na, nr, no)
    E_plab_io!(m, vars, na, nr, no, params)
    E_xlab_io!(m, vars, na, nr, no, params)
    E_realwage_io!(m, vars, na, nr, no, params)
    E_wlab_io!(m, vars, na, nr)
    E_rlab_io!(m, vars, na, nr)

    # ── Excerpt 28: Income-side GDP ──────────────────────────────
    E_delGDPINCa!(m, vars, nr, params)
    E_delGDPINCb!(m, vars, nr, params)
    E_delGDPINCc!(m, vars, nr, params)
    E_delGDPINCd!(m, vars, na, nr)
    E_delGDPINCe!(m, vars, na, nr, ns)
    E_wgdpinc!(m, vars, nr, params)

    # ── Excerpt 29: Expenditure-side GDP ─────────────────────────
    E_delXGDPEXPa_setup!(m, vars, na, nr, params)
    E_delXGDPEXPb!(m, vars, na, nr, params)
    E_delXGDPEXPc!(m, vars, na, nr, ns, params)
    E_xgdpexp!(m, vars, nr, params)
    E_delPGDPEXPa!(m, vars, na, nr, nu, params)
    E_delPGDPEXPb!(m, vars, na, nr, params)
    E_delPGDPEXPc!(m, vars, na, nr, ns, params)
    E_pgdpexp!(m, vars, nr, params)
    E_wgdpexp!(m, vars, nr)
    E_wgdpdiff!(m, vars, nr)
    E_xgne!(m, vars, nr, params)
    E_pgne!(m, vars, nr, params)
    E_wgne!(m, vars, nr)
    E_delINDTAX!(m, vars, nr)
    E_delBUDG1!(m, vars, nr)
    E_delBUDG2!(m, vars, nr)
    E_delVGDPEXP!(m, vars, nr, 9)

    # ── Labour market closure (Excerpt 38) ────────────────────────
    E_labslack!(m, vars, no)
    E_flab_i!(m, vars, na, nr, no)
    E_realwage!(m, vars, na, nr, no)

    # ── Household closure (Excerpt 39) ───────────────────────────
    E_fhou!(m, vars, na, nr, ns, nu)
    E_natfhou!(m, vars, na, nr, ns, nu, params)

    return m, vars
end
