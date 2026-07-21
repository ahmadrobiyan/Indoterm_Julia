export E_pimp!, E_pbasic!, E_ppur!, E_tuser!
export E_ppur_s!, E_phou!, E_xint!, E_xhou!, E_xinv!
export E_aint_s!, E_xint_s!, E_pint!
export E_xlab!, E_plab_o!, E_wlab_o!
export E_xlab_o!, E_pcap!, E_plnd!, E_pprim!, E_xprim!, E_aprim!, E_alab_o!, E_wprim!
export E_pvar!, E_pcst!, E_delPTX!, E_ptot!
export E_xsub!, E_xlux!, E_xhouh_s_agg!, E_alux!, E_asub!, E_wlux!, E_phouhtot!, E_whouhtot!
export E_xhoutot!, E_phoutot!
export E_xinvi!, E_pinvest!, E_pinvitot!
export E_gret!, E_xinvitot!, E_finv2!
export E_xgov!, E_xgov_s!, E_fgovtot2!, E_fgovtot3!
export E_pfexp!, E_xexpd!, E_xexp!, E_xexp_s!, E_xstocks!
export E_xint_i!, E_xuse!
export E_pdelivrd!, E_xtradmar_na!
export E_puse!, E_xtrad!
export E_xsuppmar_p!, E_psuppmar_p!, E_xsuppmar!, E_xsuppmar_d!, E_xsuppmar_rd!
export E_xmake!, E_xtotA_B!, E_xcomA_B!, E_pmake!
export E_xtrad_d!, E_xtrad_r!, E_pdomA_sum!
export E_pfin!, E_xfina!, E_xfinb!, E_xfinc!, E_xfind!, E_wfin!
export E_delTAXint!, E_delTAXhou!, E_delTAXinv!, E_delTAXgov!, E_delTAXexp!
export E_wlnd_i!, E_wcap_i!, E_wprim_i!, E_xlnd_i!, E_xcap_i!
export E_plab!, E_plab_i!, E_realwage_i!, E_xlab_i!, E_wlab_i!, E_rlab_i!
export E_plab_io!, E_xlab_io!, E_realwage_io!, E_wlab_io!, E_rlab_io!
export E_delGDPINCa!, E_delGDPINCb!, E_delGDPINCc!, E_delGDPINCd!, E_delGDPINCe!, E_wgdpinc!
export E_delXGDPEXPa_setup!, E_delXGDPEXPb!, E_delXGDPEXPc!
export E_xgdpexp!, E_delPGDPEXPa!, E_delPGDPEXPb!, E_delPGDPEXPc!
export E_pgdpexp!, E_wgdpexp!, E_wgdpdiff!
export E_xgne!, E_pgne!, E_wgne!
export E_delINDTAX!, E_delBUDG1!, E_delBUDG2!, E_delVGDPEXP!
export E_labslack!, E_flab_i!, E_realwage!
export E_fhou!, E_natfhou!
export E_plab_o_setup!, INVEST_setup!, USE_IS_setup!, USE_usc_setup!
export TAX_PUR_setup!, TRADE_setup!, TRADMAR_setup!
export SUPPMAR_setup!, SUPPMAR_D_setup!, STOCKS_setup!

# ── Individual Equation Functions ──────────────────────────────────────────
# Each follows the TERM.TAB %-change linear form.
# Lowercase variables = 100 × d(log(X)).  All equations are linear.
# Uses JuMP array constraint syntax for efficiency (compiles once per block).

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 6 — Basic prices
# ═══════════════════════════════════════════════════════════════════════════
function E_pimp!(m, vars, na, nr, params)
    pimp = vars["pimp"]; pfimp = vars["pfimp"]; phi = vars["phi"]
    @constraint(m, [c=1:na, r=1:nr], pimp[c,r] == pfimp[c] + phi)
end

function E_pbasic!(m, vars, na, nr, ns)
    pbasic = vars["pbasic"]; pdom = vars["pdom"]; pimp = vars["pimp"]
    for s_dom in 1:ns
        p = s_dom == 1 ? pdom : pimp
        @constraint(m, [c=1:na, r=1:nr], pbasic[c,s_dom,r] == p[c,r])
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 7 — Purchaser prices
# ═══════════════════════════════════════════════════════════════════════════
function E_tuser!(m, vars, na, nr, ns, nu)
    tuser = vars["tuser"]; tuser_ud = vars["tuser_ud"]; tuser_su = vars["tuser_su"]
    @constraint(m, [c=1:na, s=1:ns, u=1:nu, d=1:nr],
        tuser[c,s,u,d] == tuser_ud[c,s] + tuser_su[c,d])
end

function E_ppur!(m, vars, na, nr, ns, nu)
    ppur = vars["ppur"]; puse = vars["puse"]
    tuser = vars["tuser"]; tuser_sud = vars["tuser_sud"]
    @constraint(m, [c=1:na, s=1:ns, u=1:nu, d=1:nr],
        ppur[c,s,u,d] == puse[c,s,d] + tuser[c,s,u,d] + tuser_sud[c])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 8 — Armington (dom/imp substitution)
# ═══════════════════════════════════════════════════════════════════════════
function E_ppur_s!(m, vars, na, nr, ns, nu, params)
    ppur_s = vars["ppur_s"]; ppur = vars["ppur"]
    SRCSHR = parent(params["SRCSHR"])
    for c in 1:na, u in 1:nu, d in 1:nr
        lhs = sum(SRCSHR[c,s,u,d] * ppur[c,s,u,d] for s in 1:ns)
        @constraint(m, ppur_s[c,u,d] == lhs)
    end
end

function E_phou!(m, vars, na, nr)
    phou = vars["phou"]; ppur_s = vars["ppur_s"]
    u_hou = na + 1
    @constraint(m, [c=1:na, d=1:nr], phou[c,d] == ppur_s[c,u_hou,d])
end

function E_xint!(m, vars, na, nr, ns, sigmadomimp)
    xint = vars["xint"]; xint_s = vars["xint_s"]
    ppur = vars["ppur"]; ppur_s = vars["ppur_s"]
    @constraint(m, [c=1:na, s=1:ns, i=1:na, d=1:nr],
        xint[c,s,i,d] == xint_s[c,i,d] - sigmadomimp[c] * (ppur[c,s,i,d] - ppur_s[c,i,d]))
end

function E_xhou!(m, vars, na, nr, ns, sigmadomimp)
    xhou = vars["xhou"]; xhou_s = vars["xhou_s"]
    ppur = vars["ppur"]; phou = vars["phou"]
    u_hou = na + 1
    @constraint(m, [c=1:na, s=1:ns, d=1:nr],
        xhou[c,s,d] == xhou_s[c,d] - sigmadomimp[c] * (ppur[c,s,u_hou,d] - phou[c,d]))
end

function E_xinv!(m, vars, na, nr, ns, sigmadomimp)
    xinv = vars["xinv"]; xinv_s = vars["xinv_s"]
    ppur = vars["ppur"]; pinvest = vars["pinvest"]
    u_inv = na + 2
    @constraint(m, [c=1:na, s=1:ns, d=1:nr],
        xinv[c,s,d] == xinv_s[c,d] - sigmadomimp[c] * (ppur[c,s,u_inv,d] - pinvest[c,d]))
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 9 — Intermediate demands
# ═══════════════════════════════════════════════════════════════════════════
function E_aint_s!(m, vars, na, nr)
    aint_s = vars["aint_s"]; bint_scd = vars["bint_scd"]; bint_s = vars["bint_s"]
    @constraint(m, [c=1:na, i=1:na, d=1:nr],
        aint_s[c,i,d] == bint_scd[i] + bint_s[c,i,d])
end

function E_xint_s!(m, vars, na, nr)
    xint_s = vars["xint_s"]; atot = vars["atot"]
    aint_s = vars["aint_s"]; xtot = vars["xtot"]
    @constraint(m, [c=1:na, i=1:na, d=1:nr],
        xint_s[c,i,d] == atot[i,d] + aint_s[c,i,d] + xtot[i,d])
end

function E_pint!(m, vars, na, nr, params)
    pint = vars["pint"]; ppur_s = vars["ppur_s"]; aint_s = vars["aint_s"]
    PUR_S = parent(params["PUR_S"]); PUR_CS = parent(params["PUR_CS"])
    for i in 1:na, d in 1:nr
        PUR_CS[i,d] > 1e-10 || continue
        rhs = sum(PUR_S[c,i,d] * (ppur_s[c,i,d] + aint_s[c,i,d]) for c in 1:na)
        @constraint(m, PUR_CS[i,d] * pint[i,d] == rhs)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 10 — Labour composition
# ═══════════════════════════════════════════════════════════════════════════
function E_xlab!(m, vars, na, nr, no, sigmalab)
    xlab = vars["xlab"]; xlab_o = vars["xlab_o"]
    plab = vars["plab"]; plab_o = vars["plab_o"]
    @constraint(m, [i=1:na, o=1:no, d=1:nr],
        xlab[i,o,d] == xlab_o[i,d] - sigmalab[i] * (plab[i,o,d] - plab_o[i,d]))
end

function E_plab_o!(m, vars, na, nr, no, params)
    plab_o = vars["plab_o"]; plab = vars["plab"]
    LAB_O = parent(params["LAB_O"])
    for i in 1:na, d in 1:nr
        LAB_O[i,d] > 1e-10 || continue
        rhs = sum(get(V1LAB_idx, (i,o,d), 0.0) * plab[i,o,d] for o in 1:no)
        @constraint(m, LAB_O[i,d] * plab_o[i,d] == rhs)
    end
end

const V1LAB_idx = Dict{Tuple{Int,Int,Int}, Float64}()
function E_plab_o_setup!(na, no, nr, V1LAB)
    empty!(V1LAB_idx)
    for i in 1:na, o in 1:no, d in 1:nr
        V1LAB_idx[(i,o,d)] = V1LAB[i,o,d]
    end
end

function E_wlab_o!(m, vars, na, nr, no, params)
    wlab_o = vars["wlab_o"]; plab = vars["plab"]; xlab = vars["xlab"]
    LAB_O = parent(params["LAB_O"])
    for i in 1:na, d in 1:nr
        LAB_O[i,d] > 1e-10 || continue
        rhs = sum(get(V1LAB_idx, (i,o,d), 0.0) * (plab[i,o,d] + xlab[i,o,d]) for o in 1:no)
        @constraint(m, LAB_O[i,d] * wlab_o[i,d] == rhs)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 11 — Factor demands (CES between labour, capital, land)
# ═══════════════════════════════════════════════════════════════════════════
function E_xlab_o!(m, vars, na, nr, sigmaprim)
    xlab_o = vars["xlab_o"]; alab_o = vars["alab_o"]
    xprim = vars["xprim"]; plab_o = vars["plab_o"]; pprim = vars["pprim"]
    @constraint(m, [i=1:na, d=1:nr],
        xlab_o[i,d] - alab_o[i,d] == xprim[i,d] - sigmaprim[i] * (plab_o[i,d] + alab_o[i,d] - pprim[i,d]))
end

function E_pcap!(m, vars, na, nr, sigmaprim)
    xcap = vars["xcap"]; acap = vars["acap"]
    xprim = vars["xprim"]; pcap = vars["pcap"]; pprim = vars["pprim"]
    @constraint(m, [i=1:na, d=1:nr],
        xcap[i,d] - acap[i,d] == xprim[i,d] - sigmaprim[i] * (pcap[i,d] + acap[i,d] - pprim[i,d]))
end

function E_plnd!(m, vars, na, nr, sigmaprim)
    xlnd = vars["xlnd"]; alnd = vars["alnd"]
    xprim = vars["xprim"]; plnd = vars["plnd"]; pprim = vars["pprim"]
    @constraint(m, [i=1:na, d=1:nr],
        xlnd[i,d] - alnd[i,d] == xprim[i,d] - sigmaprim[i] * (plnd[i,d] + alnd[i,d] - pprim[i,d]))
end

function E_pprim!(m, vars, na, nr, no, params)
    pprim = vars["pprim"]; plab_o = vars["plab_o"]; alab_o = vars["alab_o"]
    pcap = vars["pcap"]; acap = vars["acap"]
    plnd = vars["plnd"]; alnd = vars["alnd"]
    PRIM = parent(params["PRIM"]); CAP = parent(params["CAP"])
    LND = parent(params["LND"]); LAB_O = parent(params["LAB_O"])
    for i in 1:na, d in 1:nr
        PRIM[i,d] > 1e-10 || continue
        rhs = LAB_O[i,d] * (plab_o[i,d] + alab_o[i,d]) + CAP[i,d] * (pcap[i,d] + acap[i,d]) + LND[i,d] * (plnd[i,d] + alnd[i,d])
        @constraint(m, PRIM[i,d] * pprim[i,d] == rhs)
    end
end

function E_xprim!(m, vars, na, nr)
    xprim = vars["xprim"]; xtot = vars["xtot"]
    atot = vars["atot"]; aprim = vars["aprim"]
    @constraint(m, [i=1:na, d=1:nr], xprim[i,d] == xtot[i,d] + atot[i,d] + aprim[i,d])
end

function E_aprim!(m, vars, na, nr)
    aprim = vars["aprim"]; bprimnat = vars["bprimnat"]
    bprim_d = vars["bprim_d"]; bprim = vars["bprim"]
    @constraint(m, [i=1:na, d=1:nr], aprim[i,d] == bprimnat + bprim_d[i] + bprim[i,d])
end

function E_alab_o!(m, vars, na, nr)
    alab_o = vars["alab_o"]; blabnat = vars["blabnat"]
    blab_d = vars["blab_d"]; blab = vars["blab"]
    @constraint(m, [i=1:na, d=1:nr], alab_o[i,d] == blabnat + blab_d[i] + blab[i,d])
end

function E_wprim!(m, vars, na, nr, no, params)
    wprim = vars["wprim"]; pcap = vars["pcap"]; xcap = vars["xcap"]
    plnd = vars["plnd"]; xlnd = vars["xlnd"]
    plab = vars["plab"]; xlab = vars["xlab"]
    PRIM = parent(params["PRIM"]); CAP = parent(params["CAP"]); LND = parent(params["LND"])
    for i in 1:na, d in 1:nr
        PRIM[i,d] > 1e-10 || continue
        rhs = CAP[i,d] * (pcap[i,d] + xcap[i,d]) + LND[i,d] * (plnd[i,d] + xlnd[i,d])
        rhs += sum(get(V1LAB_idx, (i,o,d), 0.0) * (plab[i,o,d] + xlab[i,o,d]) for o in 1:no)
        @constraint(m, PRIM[i,d] * wprim[i,d] == rhs)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 12 — Output prices
# ═══════════════════════════════════════════════════════════════════════════
function E_pvar!(m, vars, na, nr, params)
    pvar = vars["pvar"]; atot = vars["atot"]
    plab_o = vars["plab_o"]; alab_o = vars["alab_o"]; pint = vars["pint"]
    VARCST = parent(params["VARCST"]); LAB_O = parent(params["LAB_O"])
    PUR_CS = parent(params["PUR_CS"])
    for i in 1:na, d in 1:nr
        VARCST[i,d] > 1e-10 || continue
        rhs = LAB_O[i,d] * (plab_o[i,d] + alab_o[i,d]) + PUR_CS[i,d] * pint[i,d]
        @constraint(m, VARCST[i,d] * (pvar[i,d] - atot[i,d]) == rhs)
    end
end

function E_pcst!(m, vars, na, nr, params)
    pcst = vars["pcst"]; atot = vars["atot"]
    aprim = vars["aprim"]; pprim = vars["pprim"]; pint = vars["pint"]
    VCST = parent(params["VCST"]); PRIM = parent(params["PRIM"])
    PUR_CS = parent(params["PUR_CS"])
    for i in 1:na, d in 1:nr
        VCST[i,d] > 1e-10 || continue
        rhs = PRIM[i,d] * (aprim[i,d] + pprim[i,d]) + PUR_CS[i,d] * pint[i,d]
        @constraint(m, VCST[i,d] * (pcst[i,d] - atot[i,d]) == rhs)
    end
end

function E_delPTX!(m, vars, na, nr, params)
    delPTX = vars["delPTX"]; xtot = vars["xtot"]; pcst = vars["pcst"]
    delPTXRATE = vars["delPTXRATE"]
    VCST = parent(params["VCST"])
    PRODTAX_v = parent(params["PRODTAX"])
    @constraint(m, [i=1:na, d=1:nr],
        delPTX[i,d] == 0.01 * PRODTAX_v[i,d] * (xtot[i,d] + pcst[i,d]) + VCST[i,d] * delPTXRATE[i,d])
end

function E_ptot!(m, vars, na, nr, params)
    ptot = vars["ptot"]; xtot = vars["xtot"]; pcst = vars["pcst"]; delPTX = vars["delPTX"]
    VTOT = parent(params["VTOT"]); VCST = parent(params["VCST"])
    for i in 1:na, d in 1:nr
        VTOT[i,d] > 1e-10 || continue
        @constraint(m, VTOT[i,d] * (ptot[i,d] + xtot[i,d]) == VCST[i,d] * (pcst[i,d] + xtot[i,d]) + 100.0 * delPTX[i,d])
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 13 — Household demands (LES/Stone-Geary)
# ═══════════════════════════════════════════════════════════════════════════
function E_xsub!(m, vars, na, nr)
    xsub = vars["xsub"]; nhou = vars["nhou"]; asub = vars["asub"]
    @constraint(m, [c=1:na, d=1:nr], xsub[c,d] == nhou[d] + asub[c,d])
end

function E_xlux!(m, vars, na, nr)
    xlux = vars["xlux"]; phou = vars["phou"]; wlux = vars["wlux"]; alux = vars["alux"]
    @constraint(m, [c=1:na, d=1:nr], xlux[c,d] + phou[c,d] == wlux[d] + alux[c,d])
end

function E_xhouh_s_agg!(m, vars, na, nr, params)
    xhou_s = vars["xhou_s"]; xlux = vars["xlux"]; xsub = vars["xsub"]
    BLUX_v = parent(params["BLUX"])
    for c in 1:na, d in 1:nr
        @constraint(m, xhou_s[c,d] == BLUX_v[c,1,d] * xlux[c,d] + (1.0 - BLUX_v[c,1,d]) * xsub[c,d])
    end
end

function E_alux!(m, vars, na, nr, params)
    alux = vars["alux"]; asub = vars["asub"]
    SLUX_v = parent(params["SLUX"])
    for c in 1:na, d in 1:nr
        sum_slux = sum(SLUX_v[k,1,d] * asub[k,d] for k in 1:na)
        @constraint(m, alux[c,d] == asub[c,d] - sum_slux)
    end
end

function E_asub!(m, vars, na, nr, params)
    asub = vars["asub"]; ahou_s = vars["ahou_s"]
    BUDGSHR_v = parent(params["BUDGSHR"])
    for c in 1:na, d in 1:nr
        sum_budg = sum(BUDGSHR_v[k,1,d] * ahou_s[k,d] for k in 1:na)
        @constraint(m, asub[c,d] == ahou_s[c,d] - sum_budg)
    end
end

function E_wlux!(m, vars, na, nr, params)
    xhoutot = vars["xhoutot"]; xhou_s = vars["xhou_s"]
    BUDGSHR_v = parent(params["BUDGSHR"])
    for d in 1:nr
        @constraint(m, xhoutot[d] == sum(BUDGSHR_v[c,1,d] * xhou_s[c,d] for c in 1:na))
    end
end

function E_phouhtot!(m, vars, na, nr, params)
    phouhtot = vars["phouhtot"]; phou = vars["phou"]
    BUDGSHR_v = parent(params["BUDGSHR"])
    for d in 1:nr
        @constraint(m, phouhtot[d] == sum(BUDGSHR_v[c,1,d] * phou[c,d] for c in 1:na))
    end
end

function E_whouhtot!(m, vars, na, nr)
    whouhtot = vars["whouhtot"]; phouhtot = vars["phouhtot"]; xhoutot = vars["xhoutot"]
    @constraint(m, [d=1:nr], whouhtot[d] == phouhtot[d] + xhoutot[d])
end

function E_xhoutot!(m, vars, na, nr, params)
    xhoutot = vars["xhoutot"]; xhouhtot = vars["xhouhtot"]
    HOUPUR_v = parent(params["HOUPUR"])
    for d in 1:nr
        @constraint(m, HOUPUR_v[1,1,d] * (xhoutot[d] - xhouhtot[d]) == 0.0)
    end
end

function E_phoutot!(m, vars, na, nr, params)
    phoutot = vars["phoutot"]; phouhtot = vars["phouhtot"]
    HOUPUR_v = parent(params["HOUPUR"])
    for d in 1:nr
        @constraint(m, HOUPUR_v[1,1,d] * (phoutot[d] - phouhtot[d]) == 0.0)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 14 — Investment demands
# ═══════════════════════════════════════════════════════════════════════════
function E_xinvi!(m, vars, na, nr, params)
    xinvi = vars["xinvi"]; xinvitot = vars["xinvitot"]
    @constraint(m, [c=1:na, i=1:na, d=1:nr], xinvi[c,i,d] == xinvitot[i,d])
end

function E_pinvest!(m, vars, na, nr, nu)
    pinvest = vars["pinvest"]; ppur_s = vars["ppur_s"]
    u_inv = na + 2
    @constraint(m, [c=1:na, d=1:nr], pinvest[c,d] == ppur_s[c,u_inv,d])
end

function E_pinvitot!(m, vars, na, nr, params)
    pinvitot = vars["pinvitot"]; pinvest = vars["pinvest"]
    INVEST_C = parent(params["INVEST_C"])
    for i in 1:na, d in 1:nr
        INVEST_C[i,d] > 1e-10 || continue
        rhs = sum(get(INVEST_C_idx, (c,i,d), 0.0) * pinvest[c,d] for c in 1:na)
        @constraint(m, INVEST_C[i,d] * pinvitot[i,d] == rhs)
    end
end

const INVEST_C_idx = Dict{Tuple{Int,Int,Int}, Float64}()
function INVEST_setup!(na, nr, V2PUR)
    empty!(INVEST_C_idx)
    for c in 1:na, i in 1:na, d in 1:nr
        INVEST_C_idx[(c,i,d)] = V2PUR[c,i,d]
    end
end

function E_xinv_s!(m, vars, na, nr, params)
    xinv_s = vars["xinv_s"]; xinvi = vars["xinvi"]
    INVEST_I = parent(params["INVEST_I"])
    for c in 1:na, d in 1:nr
        INVEST_I[c,d] > 1e-10 || continue
        rhs = sum(get(INVEST_C_idx, (c,i,d), 0.0) * xinvi[c,i,d] for i in 1:na)
        @constraint(m, INVEST_I[c,d] * xinv_s[c,d] == rhs)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 15 — Industry-specific investment
# ═══════════════════════════════════════════════════════════════════════════
function E_gret!(m, vars, na, nr)
    gret2 = vars["gret"]; pcap = vars["pcap"]; pinvitot = vars["pinvitot"]
    @constraint(m, [i=1:na, d=1:nr], gret2[i,d] == pcap[i,d] - pinvitot[i,d])
end

function E_xinvitot!(m, vars, na, nr)
    xinvitot = vars["xinvitot"]; xcap = vars["xcap"]; ggro = vars["ggro"]
    @constraint(m, [i=1:na, d=1:nr], xinvitot[i,d] - xcap[i,d] == ggro[i,d])
end

function E_finv2!(m, vars, na, nr)
    xinvitot = vars["xinvitot"]; finv2 = vars["finv2"]; xgdpexp = vars["xgdpexp"]
    @constraint(m, [i=1:na, d=1:nr], xinvitot[i,d] == finv2[i,d] + xgdpexp[d])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 16 — Government / Export / Inventories
# ═══════════════════════════════════════════════════════════════════════════
function E_xgov!(m, vars, na, nr, ns)
    xgov = vars["xgov"]; fgovtot = vars["fgovtot"]
    fgov = vars["fgov"]; fgov_s = vars["fgov_s"]; fgovgen = vars["fgovgen"]
    @constraint(m, [c=1:na, s=1:ns, d=1:nr],
        xgov[c,s,d] == fgovtot[d] + fgov[c,s,d] + fgov_s[c,d] + fgovgen)
end

function E_xgov_s!(m, vars, na, nr, ns, params)
    xgov_s = vars["xgov_s"]; xgov = vars["xgov"]
    SRCSHR = parent(params["SRCSHR"]); u_gov = na + 3
    for c in 1:na, d in 1:nr
        @constraint(m, xgov_s[c,d] == sum(SRCSHR[c,s,u_gov,d] * xgov[c,s,d] for s in 1:ns))
    end
end

function E_fgovtot2!(m, vars, na, nr)
    fgovtot = vars["fgovtot"]; fgovtot2 = vars["fgovtot2"]; xhoutot = vars["xhoutot"]
    @constraint(m, [d=1:nr], fgovtot[d] == fgovtot2[d] + xhoutot[d])
end

function E_fgovtot3!(m, vars, na, nr)
    fgovtot = vars["fgovtot"]; fgovtot3 = vars["fgovtot3"]; xgdpexp = vars["xgdpexp"]
    @constraint(m, [d=1:nr], fgovtot[d] == fgovtot3[d] + xgdpexp[d])
end

function E_pfexp!(m, vars, na, nr)
    pfexp = vars["pfexp"]; ppur = vars["ppur"]; phi = vars["phi"]
    u_exp = na + 4
    @constraint(m, [c=1:na, d=1:nr], pfexp[c,d] == ppur[c,1,u_exp,d] - phi)
end

function E_xexpd!(m, vars, na, nr, exp_elast)
    xexpd = vars["xexpd"]; pfexp = vars["pfexp"]
    fpexp = vars["fpexp"]; fpexp_d = vars["fpexp_d"]
    fqexp = vars["fqexp"]; fqexp_d = vars["fqexp_d"]
    natfqexp = vars["natfqexp"]; natfpexp = vars["natfpexp"]
    @constraint(m, [c=1:na, d=1:nr],
        xexpd[c,d] == natfqexp + fqexp[c,d] + fqexp_d[c] - exp_elast[c] * (pfexp[c,d] - fpexp[c,d] - fpexp_d[c] - natfpexp))
end

function E_xexp!(m, vars, na, nr, ns)
    xexp = vars["xexp"]; xexpd = vars["xexpd"]
    @constraint(m, [c=1:na, d=1:nr], xexp[c,1,d] == xexpd[c,d])
    @constraint(m, [c=1:na, d=1:nr], xexp[c,2,d] == 0.0)
end

function E_xexp_s!(m, vars, na, nr, ns, params)
    xexp_s = vars["xexp_s"]; xexp = vars["xexp"]
    SRCSHR = parent(params["SRCSHR"]); u_exp = na + 4
    for c in 1:na, d in 1:nr
        @constraint(m, xexp_s[c,d] == sum(SRCSHR[c,s,u_exp,d] * xexp[c,s,d] for s in 1:ns))
    end
end

function E_xstocks!(m, vars, na, nr)
    xstocks = vars["xstocks"]; xtot = vars["xtot"]; fxstocks = vars["fxstocks"]
    @constraint(m, [i=1:na, d=1:nr], xstocks[i,d] == xtot[i,d] + fxstocks[i,d])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 17 — Total regional demand
# ═══════════════════════════════════════════════════════════════════════════
function E_xint_i!(m, vars, na, nr, ns, params)
    xint_i = vars["xint_i"]; xint = vars["xint"]
    USE_I = parent(params["USE_I"])
    for c in 1:na, s in 1:ns, d in 1:nr
        USE_I[c,s,d] > 1e-10 || continue
        rhs = sum(get(USE_IS_idx, (c,s,i,d), 0.0) * xint[c,s,i,d] for i in 1:na)
        @constraint(m, USE_I[c,s,d] * xint_i[c,s,d] == rhs)
    end
end

const USE_IS_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
function USE_IS_setup!(na, nr, ns, USE)
    empty!(USE_IS_idx)
    for c in 1:na, s in 1:ns, i in 1:na, d in 1:nr
        USE_IS_idx[(c,s,i,d)] = USE[c,s,i,d]
    end
end

function E_xuse!(m, vars, na, nr, ns, params)
    xuse = vars["xuse"]; xint_i = vars["xint_i"]
    xhou = vars["xhou"]; xinv = vars["xinv"]
    xgov = vars["xgov"]; xexp = vars["xexp"]
    USE_U = parent(params["USE_U"]); USE_I = parent(params["USE_I"])
    u_hou = na + 1; u_inv = na + 2; u_gov = na + 3; u_exp = na + 4
    for c in 1:na, s in 1:ns, d in 1:nr
        USE_U[c,s,d] > 1e-10 || continue
        rhs = USE_I[c,s,d] * xint_i[c,s,d]
        rhs += get(USE_usc_idx, (c,s,u_hou,d), 0.0) * xhou[c,s,d]
        rhs += get(USE_usc_idx, (c,s,u_inv,d), 0.0) * xinv[c,s,d]
        rhs += get(USE_usc_idx, (c,s,u_gov,d), 0.0) * xgov[c,s,d]
        rhs += get(USE_usc_idx, (c,s,u_exp,d), 0.0) * xexp[c,s,d]
        @constraint(m, USE_U[c,s,d] * xuse[c,s,d] == rhs)
    end
end

const USE_usc_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
function USE_usc_setup!(na, nr, ns, nu, USE)
    empty!(USE_usc_idx)
    for c in 1:na, s in 1:ns, u in 1:nu, d in 1:nr
        USE_usc_idx[(c,s,u,d)] = USE[c,s,u,d]
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 19 — Delivering goods (margins, Leontief)
# ═══════════════════════════════════════════════════════════════════════════
function E_xtradmar_na!(m, vars, na, nr, ns, nm)
    xtradmar = vars["xtradmar"]; xtrad = vars["xtrad"]; atradmar = vars["atradmar"]
    @constraint(m, [c=1:na, s=1:ns, m_ix=1:nm, r=1:nr, d=1:nr],
        xtradmar[c,s,m_ix,r,d] == xtrad[c,s,r,d] + atradmar[c,s,m_ix,r,d])
end

function E_pdelivrd!(m, vars, na, nr, ns, nm, params)
    pdelivrd = vars["pdelivrd"]; pbasic = vars["pbasic"]
    psuppmar_p = vars["psuppmar_p"]; atradmar = vars["atradmar"]
    BASSHR = parent(params["BASSHR"]); MARSHR = parent(params["MARSHR"])
    for c in 1:na, s in 1:ns, r in 1:nr, d in 1:nr
        lhs = BASSHR[c,s,r,d] * pbasic[c,s,r]
        for mh in 1:nm
            lhs += MARSHR[c,s,mh,r,d] * (psuppmar_p[mh,r,d] + atradmar[c,s,mh,r,d])
        end
        @constraint(m, pdelivrd[c,s,r,d] == lhs)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 20 — Regional sourcing (CES between origins)
# ═══════════════════════════════════════════════════════════════════════════
function E_puse!(m, vars, na, nr, ns, params)
    puse = vars["puse"]; pdelivrd = vars["pdelivrd"]; atrad = vars["atrad"]
    DELIVRD = parent(params["DELIVRD"]); DELIVRD_R = parent(params["DELIVRD_R"])
    for c in 1:na, s in 1:ns, d in 1:nr
        DELIVRD_R[c,s,d] > 1e-10 || continue
        rhs = sum(DELIVRD[c,s,r,d] * (pdelivrd[c,s,r,d] + atrad[c,s,r,d]) for r in 1:nr)
        @constraint(m, DELIVRD_R[c,s,d] * puse[c,s,d] == rhs)
    end
end

function E_avesrctwist!(m, vars, na, nr, ns, params)
    avesrctwist = vars["avesrctwist"]; srctwist = vars["srctwist"]
    DELIVRD = parent(params["DELIVRD"]); DELIVRD_R = parent(params["DELIVRD_R"])
    for c in 1:na, s in 1:ns, d in 1:nr
        DELIVRD_R[c,s,d] > 1e-10 || continue
        rhs = sum(DELIVRD[c,s,r,d] * srctwist[c,s,r,d] for r in 1:nr)
        @constraint(m, DELIVRD_R[c,s,d] * avesrctwist[c,s,d] == rhs)
    end
end

function E_xtrad!(m, vars, na, nr, ns, params)
    xtrad = vars["xtrad"]; xuse = vars["xuse"]
    pdelivrd = vars["pdelivrd"]; puse = vars["puse"]; atrad = vars["atrad"]
    srctwist = vars["srctwist"]; avesrctwist = vars["avesrctwist"]
    SGDD = parent(params["SGDD"])  # SIGMADOMDOM(c) — CES elasticity of substitution across regional sources
    @constraint(m, [c=1:na, s=1:ns, r=1:nr, d=1:nr],
        xtrad[c,s,r,d] - atrad[c,s,r,d] == xuse[c,s,d] + srctwist[c,s,r,d] - avesrctwist[c,s,d]
            - SGDD[c] * (pdelivrd[c,s,r,d] + atrad[c,s,r,d] - puse[c,s,d]))
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 21 — Margin supply
# ═══════════════════════════════════════════════════════════════════════════
function E_xsuppmar_p!(m, vars, na, nr, ns, nm, params)
    xsuppmar_p = vars["xsuppmar_p"]; xtradmar = vars["xtradmar"]
    for m_ix in 1:nm, r in 1:nr, d in 1:nr
        id01_val = sum(abs, get(TRADMAR_idx, (c,s,m_ix,r,d), 0.0) for c in 1:na, s in 1:ns)
        id01_val > 1e-10 || continue
        rhs = sum(get(TRADMAR_idx, (c,s,m_ix,r,d), 0.0) * xtradmar[c,s,m_ix,r,d] for c in 1:na, s in 1:ns)
        @constraint(m, id01_val * xsuppmar_p[m_ix,r,d] == rhs)
    end
end

const TRADMAR_idx = Dict{Tuple{Int,Int,Int,Int,Int}, Float64}()
function TRADMAR_setup!(na, nr, ns, nm, TMAR)
    empty!(TRADMAR_idx)
    for c in 1:na, s in 1:ns, m in 1:nm, r in 1:nr, d in 1:nr
        TRADMAR_idx[(c,s,m,r,d)] = TMAR[c,s,m,r,d]
    end
end

function E_psuppmar_p!(m, vars, na, nr, nm, params)
    psuppmar_p = vars["psuppmar_p"]; pdom = vars["pdom"]; asuppmar = vars["asuppmar"]
    for m_ix in 1:nm, r in 1:nr, d in 1:nr
        id01_val = sum(abs, get(SUPPMAR_idx, (m_ix,r,d,p), 0.0) for p in 1:nr)
        id01_val > 1e-10 || continue
        rhs = sum(get(SUPPMAR_idx, (m_ix,r,d,p), 0.0) * (pdom[m_ix,p] + asuppmar[m_ix,r,d,p]) for p in 1:nr)
        @constraint(m, id01_val * psuppmar_p[m_ix,r,d] == rhs)
    end
end

const SUPPMAR_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
function SUPPMAR_setup!(nm, nr, MARS, DIST)
    empty!(SUPPMAR_idx)
    for m in 1:nm, r in 1:nr, d in 1:nr, p in 1:nr
        SUPPMAR_idx[(m,r,d,p)] = MARS[m,r,d,p] * DIST[r,d]
    end
end

function E_xsuppmar!(m, vars, na, nr, nm)
    xsuppmar = vars["xsuppmar"]; xsuppmar_p = vars["xsuppmar_p"]
    asuppmar = vars["asuppmar"]
    @constraint(m, [m_ix=1:nm, r=1:nr, d=1:nr, p=1:nr],
        xsuppmar[m_ix,r,d,p] == xsuppmar_p[m_ix,r,d] + asuppmar[m_ix,r,d,p])
end

function E_xsuppmar_d!(m, vars, na, nr, ns, nm, params)
    xsuppmar_d = vars["xsuppmar_d"]; xsuppmar = vars["xsuppmar"]
    for m_ix in 1:nm, r in 1:nr, p in 1:nr
        id01_val = sum(abs, get(SUPPMAR_idx, (m_ix,r,d,p), 0.0) for d in 1:nr)
        id01_val > 1e-10 || continue
        rhs = sum(get(SUPPMAR_idx, (m_ix,r,d,p), 0.0) * xsuppmar[m_ix,r,d,p] for d in 1:nr)
        @constraint(m, id01_val * xsuppmar_d[m_ix,r,p] == rhs)
    end
end

function E_xsuppmar_rd!(m, vars, na, nr, ns, nm, params)
    xsuppmar_rd = vars["xsuppmar_rd"]; xsuppmar_d = vars["xsuppmar_d"]
    for m_ix in 1:nm, p in 1:nr
        id01_val = sum(abs, get(SUPPMAR_D_idx, (m_ix,r,p), 0.0) for r in 1:nr)
        id01_val > 1e-10 || continue
        rhs = sum(get(SUPPMAR_D_idx, (m_ix,r,p), 0.0) * xsuppmar_d[m_ix,r,p] for r in 1:nr)
        @constraint(m, id01_val * xsuppmar_rd[m_ix,p] == rhs)
    end
end

const SUPPMAR_D_idx = Dict{Tuple{Int,Int,Int}, Float64}()
function SUPPMAR_D_setup!(nm, nr)
    empty!(SUPPMAR_D_idx)
    for m in 1:nm, r in 1:nr, p in 1:nr
        SUPPMAR_D_idx[(m,r,p)] = sum(get(SUPPMAR_idx, (m,r,d,p), 0.0) for d in 1:nr)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 22 — MAKE / CET (multi-product industries)
# ═══════════════════════════════════════════════════════════════════════════
function E_xmake!(m, vars, na, nr, sigmaout, params)
    xmake = vars["xmake"]; xtot = vars["xtot"]; pmake = vars["pmake"]; ptot = vars["ptot"]
    @constraint(m, [c=1:na, i=1:na, d=1:nr],
        xmake[c,i,d] == xtot[i,d] + sigmaout[i] * (pmake[c,i,d] - ptot[i,d]))
end

function E_xtotA_B!(m, vars, na, nr, params)
    ptot = vars["ptot"]; pmake = vars["pmake"]
    MAKE = parent(params["MAKE"])
    for i in 1:na, d in 1:nr
        sum_mk = sum(MAKE[c,i,d] for c in 1:na)
        sum_mk > 0 || continue
        rhs = sum((MAKE[c,i,d] / sum_mk) * pmake[c,i,d] for c in 1:na)
        @constraint(m, ptot[i,d] == rhs)
    end
end

function E_xcomA_B!(m, vars, na, nr, params)
    xcom = vars["xcom"]; xmake = vars["xmake"]
    MAKE = parent(params["MAKE"])
    for c in 1:na, d in 1:nr
        mk_i = sum(MAKE[c,i,d] for i in 1:na)
        mk_i > 0 || continue
        rhs = sum(MAKE[c,i,d] * xmake[c,i,d] for i in 1:na)
        @constraint(m, mk_i * xcom[c,d] == rhs)
    end
end

function E_pmake!(m, vars, na, nr)
    pmake = vars["pmake"]; pdom = vars["pdom"]; xmake = vars["xmake"]; xcom = vars["xcom"]
    @constraint(m, [c=1:na, i=1:na, d=1:nr],
        pmake[c,i,d] == pdom[c,d] - 0.05 * (xmake[c,i,d] - xcom[c,d]))
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 23 — Market clearing
# ═══════════════════════════════════════════════════════════════════════════
function E_xtrad_d!(m, vars, na, nr, ns, params)
    xtrad_d = vars["xtrad_d"]; xtrad = vars["xtrad"]
    TRADE_D = parent(params["TRADE_D"])
    for c in 1:na, s in 1:ns, r in 1:nr
        TRADE_D[c,s,r] > 1e-10 || continue
        rhs = sum(get(TRADE_idx, (c,s,r,d), 0.0) * xtrad[c,s,r,d] for d in 1:nr)
        @constraint(m, TRADE_D[c,s,r] * xtrad_d[c,s,r] == rhs)
    end
end

const TRADE_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
function TRADE_setup!(na, nr, ns, TRADE)
    empty!(TRADE_idx)
    for c in 1:na, s in 1:ns, r in 1:nr, d in 1:nr
        TRADE_idx[(c,s,r,d)] = TRADE[c,s,r,d]
    end
end

function E_xtrad_r!(m, vars, na, nr, ns, params)
    xtrad_r = vars["xtrad_r"]; xtrad = vars["xtrad"]
    TRADE_R = parent(params["TRADE_R"])
    for c in 1:na, s in 1:ns, d in 1:nr
        TRADE_R[c,s,d] > 1e-10 || continue
        rhs = sum(get(TRADE_idx, (c,s,r,d), 0.0) * xtrad[c,s,r,d] for r in 1:nr)
        @constraint(m, TRADE_R[c,s,d] * xtrad_r[c,s,d] == rhs)
    end
end

function E_pdomA_sum!(m, vars, na, nr)
    xcom = vars["xcom"]; xtrad_d = vars["xtrad_d"]
    @constraint(m, [c=1:na, r=1:nr], xcom[c,r] == xtrad_d[c,1,r])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 24 — Final demand aggregates
# ═══════════════════════════════════════════════════════════════════════════
function E_pfin!(m, vars, na, nr, nu, params)
    pfin = vars["pfin"]; ppur_s = vars["ppur_s"]
    PUR_S = parent(params["PUR_S"]); PUR_CS = parent(params["PUR_CS"])
    u_list = [na+1, na+2, na+3, na+4]
    for (fi, u) in enumerate(u_list), d in 1:nr
        PUR_CS[u,d] > 1e-10 || continue
        @constraint(m, PUR_CS[u,d] * pfin[fi,d] == sum(PUR_S[c,u,d] * ppur_s[c,u,d] for c in 1:na))
    end
end

function E_xfina!(m, vars, na, nr)
    xfin = vars["xfin"]; xhoutot = vars["xhoutot"]
    @constraint(m, [d=1:nr], xfin[1,d] == xhoutot[d])
end

function E_xfinb!(m, vars, na, nr, params)
    xfin = vars["xfin"]; xinv_s = vars["xinv_s"]
    PUR_CS = parent(params["PUR_CS"]); PUR_S = parent(params["PUR_S"])
    u_inv = na + 2
    for d in 1:nr
        PUR_CS[u_inv,d] > 1e-10 || continue
        @constraint(m, PUR_CS[u_inv,d] * xfin[2,d] == sum(PUR_S[c,u_inv,d] * xinv_s[c,d] for c in 1:na))
    end
end

function E_xfinc!(m, vars, na, nr, params)
    xfin = vars["xfin"]; xgov_s = vars["xgov_s"]
    PUR_CS = parent(params["PUR_CS"]); PUR_S = parent(params["PUR_S"])
    u_gov = na + 3
    for d in 1:nr
        PUR_CS[u_gov,d] > 1e-10 || continue
        @constraint(m, PUR_CS[u_gov,d] * xfin[3,d] == sum(PUR_S[c,u_gov,d] * xgov_s[c,d] for c in 1:na))
    end
end

function E_xfind!(m, vars, na, nr, ns, params)
    xfin = vars["xfin"]; xexp = vars["xexp"]
    PUR_CS = parent(params["PUR_CS"])
    u_exp = na + 4
    for d in 1:nr
        PUR_CS[u_exp,d] > 1e-10 || continue
        rhs = sum(get(PUR_src_idx, (c,s,u_exp,d), 0.0) * xexp[c,s,d] for c in 1:na, s in 1:ns)
        @constraint(m, PUR_CS[u_exp,d] * xfin[4,d] == rhs)
    end
end

const PUR_src_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
function PUR_src_setup!(na, nr, ns, nu, PUR)
    empty!(PUR_src_idx)
    for c in 1:na, s in 1:ns, u in 1:nu, d in 1:nr
        PUR_src_idx[(c,s,u,d)] = PUR[c,s,u,d]
    end
end

function E_wfin!(m, vars, nr)
    wfin = vars["wfin"]; xfin = vars["xfin"]; pfin = vars["pfin"]
    @constraint(m, [f=1:4, d=1:nr], wfin[f,d] == xfin[f,d] + pfin[f,d])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 26 — Commodity tax revenues
# ═══════════════════════════════════════════════════════════════════════════
function E_delTAXint!(m, vars, na, nr, ns, params)
    delTAXint = vars["delTAXint"]; xint = vars["xint"]; puse = vars["puse"]; tuser = vars["tuser"]
    @constraint(m, [c=1:na, s=1:ns, i=1:na, d=1:nr],
        delTAXint[c,s,i,d] == 0.01 * get(TAX_idx, (c,s,i,d), 0.0) * (xint[c,s,i,d] + puse[c,s,d]) + 0.01 * get(PUR_idx, (c,s,i,d), 0.0) * tuser[c,s,i,d])
end

const TAX_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
const PUR_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
function TAX_PUR_setup!(na, nr, ns, nu, TAX_data, PUR_data)
    empty!(TAX_idx); empty!(PUR_idx)
    for c in 1:na, s in 1:ns, u in 1:nu, d in 1:nr
        TAX_idx[(c,s,u,d)] = TAX_data[c,s,u,d]
        PUR_idx[(c,s,u,d)] = PUR_data[c,s,u,d]
    end
end

function E_delTAXhou!(m, vars, na, nr, ns, params)
    delTAXhou = vars["delTAXhou"]; xhou = vars["xhou"]; puse = vars["puse"]; tuser = vars["tuser"]
    u_hou = na + 1
    @constraint(m, [c=1:na, s=1:ns, d=1:nr],
        delTAXhou[c,s,d] == 0.01 * get(TAX_idx, (c,s,u_hou,d), 0.0) * (xhou[c,s,d] + puse[c,s,d]) + 0.01 * get(PUR_idx, (c,s,u_hou,d), 0.0) * tuser[c,s,u_hou,d])
end

function E_delTAXinv!(m, vars, na, nr, ns, params)
    delTAXinv = vars["delTAXinv"]; xinv = vars["xinv"]; puse = vars["puse"]; tuser = vars["tuser"]
    u_inv = na + 2
    @constraint(m, [c=1:na, s=1:ns, d=1:nr],
        delTAXinv[c,s,d] == 0.01 * get(TAX_idx, (c,s,u_inv,d), 0.0) * (xinv[c,s,d] + puse[c,s,d]) + 0.01 * get(PUR_idx, (c,s,u_inv,d), 0.0) * tuser[c,s,u_inv,d])
end

function E_delTAXgov!(m, vars, na, nr, ns, params)
    delTAXgov = vars["delTAXgov"]; xgov = vars["xgov"]; puse = vars["puse"]; tuser = vars["tuser"]
    u_gov = na + 3
    @constraint(m, [c=1:na, s=1:ns, d=1:nr],
        delTAXgov[c,s,d] == 0.01 * get(TAX_idx, (c,s,u_gov,d), 0.0) * (xgov[c,s,d] + puse[c,s,d]) + 0.01 * get(PUR_idx, (c,s,u_gov,d), 0.0) * tuser[c,s,u_gov,d])
end

function E_delTAXexp!(m, vars, na, nr, ns, params)
    delTAXexp = vars["delTAXexp"]; xexp = vars["xexp"]; puse = vars["puse"]; tuser = vars["tuser"]
    u_exp = na + 4
    @constraint(m, [c=1:na, s=1:ns, d=1:nr],
        delTAXexp[c,s,d] == 0.01 * get(TAX_idx, (c,s,u_exp,d), 0.0) * (xexp[c,s,d] + puse[c,s,d]) + 0.01 * get(PUR_idx, (c,s,u_exp,d), 0.0) * tuser[c,s,u_exp,d])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 27 — Primary Factor Aggregates
# ═══════════════════════════════════════════════════════════════════════════
function E_wlnd_i!(m, vars, na, nr, params)
    wlnd_i = vars["wlnd_i"]; plnd = vars["plnd"]; xlnd = vars["xlnd"]
    LND_I = parent(params["LND_I"]); LND = parent(params["LND"])
    for d in 1:nr
        LND_I[d] > 1e-10 || continue
        @constraint(m, LND_I[d] * wlnd_i[d] == sum(LND[i,d] * (plnd[i,d] + xlnd[i,d]) for i in 1:na))
    end
end

function E_wcap_i!(m, vars, na, nr, params)
    wcap_i = vars["wcap_i"]; pcap = vars["pcap"]; xcap = vars["xcap"]
    CAP_I = parent(params["CAP_I"]); CAP = parent(params["CAP"])
    for d in 1:nr
        CAP_I[d] > 1e-10 || continue
        @constraint(m, CAP_I[d] * wcap_i[d] == sum(CAP[i,d] * (pcap[i,d] + xcap[i,d]) for i in 1:na))
    end
end

function E_wprim_i!(m, vars, na, nr, params)
    wprim_i = vars["wprim_i"]; wprim = vars["wprim"]
    PRIM_I = parent(params["PRIM_I"]); PRIM = parent(params["PRIM"])
    for d in 1:nr
        PRIM_I[d] > 1e-10 || continue
        @constraint(m, PRIM_I[d] * wprim_i[d] == sum(PRIM[i,d] * wprim[i,d] for i in 1:na))
    end
end

function E_xlnd_i!(m, vars, na, nr, params)
    xlnd_i = vars["xlnd_i"]; xlnd = vars["xlnd"]
    LND_I = parent(params["LND_I"]); LND = parent(params["LND"])
    for d in 1:nr
        LND_I[d] > 1e-10 || continue
        @constraint(m, LND_I[d] * xlnd_i[d] == sum(LND[i,d] * xlnd[i,d] for i in 1:na))
    end
end

function E_xcap_i!(m, vars, na, nr, params)
    xcap_i = vars["xcap_i"]; xcap = vars["xcap"]
    CAP_I = parent(params["CAP_I"]); CAP = parent(params["CAP"])
    for d in 1:nr
        CAP_I[d] > 1e-10 || continue
        @constraint(m, CAP_I[d] * xcap_i[d] == sum(CAP[i,d] * xcap[i,d] for i in 1:na))
    end
end

function E_plab!(m, vars, na, nr, no)
    realwage = vars["realwage"]; plab = vars["plab"]; pfin = vars["pfin"]
    @constraint(m, [i=1:na, o=1:no, d=1:nr], realwage[i,o,d] == plab[i,o,d] - pfin[1,d])
end

function E_plab_i!(m, vars, na, nr, no, params)
    plab_i = vars["plab_i"]; plab = vars["plab"]
    for o in 1:no, d in 1:nr
        lab_sum = sum(get(V1LAB_idx, (i,o,d), 0.0) for i in 1:na)
        lab_sum > 0 || continue
        @constraint(m, lab_sum * plab_i[o,d] == sum(get(V1LAB_idx, (i,o,d), 0.0) * plab[i,o,d] for i in 1:na))
    end
end

function E_realwage_i!(m, vars, na, nr, no, params)
    realwage_i = vars["realwage_i"]; realwage = vars["realwage"]
    for o in 1:no, d in 1:nr
        rw_sum = sum(get(V1LAB_idx, (i,o,d), 0.0) for i in 1:na)
        rw_sum > 0 || continue
        @constraint(m, rw_sum * realwage_i[o,d] == sum(get(V1LAB_idx, (i,o,d), 0.0) * realwage[i,o,d] for i in 1:na))
    end
end

function E_xlab_i!(m, vars, na, nr, no, params)
    xlab_i = vars["xlab_i"]; xlab = vars["xlab"]
    for o in 1:no, d in 1:nr
        lab_sum = sum(get(V1LAB_idx, (i,o,d), 0.0) for i in 1:na)
        lab_sum > 0 || continue
        @constraint(m, lab_sum * xlab_i[o,d] == sum(get(V1LAB_idx, (i,o,d), 0.0) * xlab[i,o,d] for i in 1:na))
    end
end

function E_wlab_i!(m, vars, na, nr, no)
    wlab_i = vars["wlab_i"]; xlab_i = vars["xlab_i"]; plab_i = vars["plab_i"]
    @constraint(m, [o=1:no, d=1:nr], wlab_i[o,d] == xlab_i[o,d] + plab_i[o,d])
end

function E_rlab_i!(m, vars, na, nr, no)
    rlab_i = vars["rlab_i"]; xlab_i = vars["xlab_i"]; realwage_i = vars["realwage_i"]
    @constraint(m, [o=1:no, d=1:nr], rlab_i[o,d] == xlab_i[o,d] + realwage_i[o,d])
end

function E_plab_io!(m, vars, na, nr, no, params)
    plab_io = vars["plab_io"]; plab_i = vars["plab_i"]
    LAB_IO = parent(params["LAB_IO"]); LAB_I = parent(params["LAB_I"])
    for d in 1:nr
        LAB_IO[d] > 1e-10 || continue
        @constraint(m, LAB_IO[d] * plab_io[d] == sum(LAB_I[o,d] * plab_i[o,d] for o in 1:no))
    end
end

function E_xlab_io!(m, vars, na, nr, no, params)
    xlab_io = vars["xlab_io"]; xlab_i = vars["xlab_i"]
    LAB_IO = parent(params["LAB_IO"]); LAB_I = parent(params["LAB_I"])
    for d in 1:nr
        LAB_IO[d] > 1e-10 || continue
        @constraint(m, LAB_IO[d] * xlab_io[d] == sum(LAB_I[o,d] * xlab_i[o,d] for o in 1:no))
    end
end

function E_realwage_io!(m, vars, na, nr, no, params)
    realwage_io = vars["realwage_io"]; realwage_i = vars["realwage_i"]
    LAB_IO = parent(params["LAB_IO"]); LAB_I = parent(params["LAB_I"])
    for d in 1:nr
        LAB_IO[d] > 1e-10 || continue
        @constraint(m, LAB_IO[d] * realwage_io[d] == sum(LAB_I[o,d] * realwage_i[o,d] for o in 1:no))
    end
end

function E_wlab_io!(m, vars, na, nr)
    wlab_io = vars["wlab_io"]; xlab_io = vars["xlab_io"]; plab_io = vars["plab_io"]
    @constraint(m, [d=1:nr], wlab_io[d] == xlab_io[d] + plab_io[d])
end

function E_rlab_io!(m, vars, na, nr)
    rlab_io = vars["rlab_io"]; xlab_io = vars["xlab_io"]; realwage_io = vars["realwage_io"]
    @constraint(m, [d=1:nr], rlab_io[d] == xlab_io[d] + realwage_io[d])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 28 — Income-side GDP
# ═══════════════════════════════════════════════════════════════════════════
function E_delGDPINCa!(m, vars, nr, params)
    delGDPINC = vars["delGDPINC"]; wlnd_i = vars["wlnd_i"]
    LND_I = parent(params["LND_I"])
    @constraint(m, [d=1:nr], delGDPINC[d,1] == 0.01 * LND_I[d] * wlnd_i[d])
end

function E_delGDPINCb!(m, vars, nr, params)
    delGDPINC = vars["delGDPINC"]; wcap_i = vars["wcap_i"]
    CAP_I = parent(params["CAP_I"])
    @constraint(m, [d=1:nr], delGDPINC[d,2] == 0.01 * CAP_I[d] * wcap_i[d])
end

function E_delGDPINCc!(m, vars, nr, params)
    delGDPINC = vars["delGDPINC"]; wlab_io = vars["wlab_io"]
    LAB_IO = parent(params["LAB_IO"])
    @constraint(m, [d=1:nr], delGDPINC[d,3] == 0.01 * LAB_IO[d] * wlab_io[d])
end

function E_delGDPINCd!(m, vars, na, nr)
    delGDPINC = vars["delGDPINC"]; delPTX = vars["delPTX"]
    for d in 1:nr
        @constraint(m, delGDPINC[d,4] == sum(delPTX[i,d] for i in 1:na))
    end
end

function E_delGDPINCe!(m, vars, na, nr, ns)
    delGDPINC = vars["delGDPINC"]
    delTAXint = vars["delTAXint"]; delTAXhou = vars["delTAXhou"]
    delTAXinv = vars["delTAXinv"]; delTAXgov = vars["delTAXgov"]; delTAXexp = vars["delTAXexp"]
    for d in 1:nr
        rhs = sum(delTAXint[c,s,i,d] for c in 1:na, s in 1:ns, i in 1:na)
        for c in 1:na, s in 1:ns
            rhs += delTAXhou[c,s,d] + delTAXinv[c,s,d] + delTAXgov[c,s,d] + delTAXexp[c,s,d]
        end
        @constraint(m, delGDPINC[d,5] == rhs)
    end
end

function E_wgdpinc!(m, vars, nr, params)
    wgdpinc = vars["wgdpinc"]; delGDPINC = vars["delGDPINC"]
    GDPINCSUM = parent(params["GDPINCSUM"])
    for d in 1:nr
        gdp_inc = sum(GDPINCSUM[d,k] for k in 1:5)
        gdp_inc > 1e-10 || continue
        @constraint(m, gdp_inc * wgdpinc[d] == 100.0 * sum(delGDPINC[d,k] for k in 1:5))
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 29 — Expenditure-side GDP
# ═══════════════════════════════════════════════════════════════════════════
function E_delXGDPEXPa_setup!(m, vars, na, nr, params)
    delXGDPEXP = vars["delXGDPEXP"]; xfin = vars["xfin"]
    PUR_CS = parent(params["PUR_CS"])
    u_list = [na+1, na+2, na+3, na+4]
    for d in 1:nr, (fi, u) in enumerate(u_list)
        @constraint(m, delXGDPEXP[d,fi] == 0.01 * PUR_CS[u,d] * xfin[fi,d])
    end
end

function E_delXGDPEXPb!(m, vars, na, nr, params)
    delXGDPEXP = vars["delXGDPEXP"]; xstocks = vars["xstocks"]
    for d in 1:nr
        @constraint(m, delXGDPEXP[d,5] == 0.01 * sum(get(STOCKS_idx, (i,d), 0.0) * xstocks[i,d] for i in 1:na))
    end
end

const STOCKS_idx = Dict{Tuple{Int,Int}, Float64}()
function STOCKS_setup!(na, nr, STOK)
    empty!(STOCKS_idx)
    for i in 1:na, d in 1:nr
        STOCKS_idx[(i,d)] = STOK[i,d]
    end
end

function E_delXGDPEXPc!(m, vars, na, nr, ns, params)
    delXGDPEXP = vars["delXGDPEXP"]; xtrad_d = vars["xtrad_d"]
    TRADE_D = parent(params["TRADE_D"])
    for q in 1:nr
        @constraint(m, delXGDPEXP[q,6] == -0.01 * sum(TRADE_D[c,2,q] * xtrad_d[c,2,q] for c in 1:na))
    end
end

function E_xgdpexp!(m, vars, nr, params)
    xgdpexp = vars["xgdpexp"]; delXGDPEXP = vars["delXGDPEXP"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for d in 1:nr
        gdp_exp = sum(GDPEXPSUM[d,k] for k in 1:9)
        gdp_exp > 1e-10 || continue
        @constraint(m, gdp_exp * xgdpexp[d] == 100.0 * sum(delXGDPEXP[d,k] for k in 1:9))
    end
end

function E_delPGDPEXPa!(m, vars, na, nr, nu, params)
    delPGDPEXP = vars["delPGDPEXP"]; pfin = vars["pfin"]
    PUR_CS = parent(params["PUR_CS"])
    u_list = [na+1, na+2, na+3, na+4]
    for d in 1:nr, (fi, u) in enumerate(u_list)
        @constraint(m, delPGDPEXP[d,fi] == 0.01 * PUR_CS[u,d] * pfin[fi,d])
    end
end

function E_delPGDPEXPb!(m, vars, na, nr, params)
    delPGDPEXP = vars["delPGDPEXP"]; ptot = vars["ptot"]
    for d in 1:nr
        @constraint(m, delPGDPEXP[d,5] == 0.01 * sum(get(STOCKS_idx, (i,d), 0.0) * ptot[i,d] for i in 1:na))
    end
end

function E_delPGDPEXPc!(m, vars, na, nr, ns, params)
    delPGDPEXP = vars["delPGDPEXP"]; pimp = vars["pimp"]
    TRADE_D = parent(params["TRADE_D"])
    for q in 1:nr
        @constraint(m, delPGDPEXP[q,6] == -0.01 * sum(TRADE_D[c,2,q] * pimp[c,q] for c in 1:na))
    end
end

function E_pgdpexp!(m, vars, nr, params)
    pgdpexp = vars["pgdpexp"]; delPGDPEXP = vars["delPGDPEXP"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for d in 1:nr
        gdp_exp = sum(GDPEXPSUM[d,k] for k in 1:9)
        gdp_exp > 1e-10 || continue
        @constraint(m, gdp_exp * pgdpexp[d] == 100.0 * sum(delPGDPEXP[d,k] for k in 1:9))
    end
end

function E_wgdpexp!(m, vars, nr)
    wgdpexp = vars["wgdpexp"]; xgdpexp = vars["xgdpexp"]; pgdpexp = vars["pgdpexp"]
    @constraint(m, [d=1:nr], wgdpexp[d] == xgdpexp[d] + pgdpexp[d])
end

function E_wgdpdiff!(m, vars, nr)
    wgdpdiff = vars["wgdpdiff"]; wgdpinc = vars["wgdpinc"]; wgdpexp = vars["wgdpexp"]
    @constraint(m, [d=1:nr], wgdpdiff[d] == wgdpinc[d] - wgdpexp[d])
end

function E_xgne!(m, vars, nr, params)
    xgne = vars["xgne"]; delXGDPEXP = vars["delXGDPEXP"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for d in 1:nr
        gne_gdp = sum(GDPEXPSUM[d,k] for k in [1,2,3,4,5])
        gne_gdp > 1e-10 || continue
        @constraint(m, gne_gdp * xgne[d] == 100.0 * sum(delXGDPEXP[d,k] for k in [1,2,3,4,5]))
    end
end

function E_pgne!(m, vars, nr, params)
    pgne = vars["pgne"]; delPGDPEXP = vars["delPGDPEXP"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for d in 1:nr
        gne_gdp = sum(GDPEXPSUM[d,k] for k in [1,2,3,4,5])
        gne_gdp > 1e-10 || continue
        @constraint(m, gne_gdp * pgne[d] == 100.0 * sum(delPGDPEXP[d,k] for k in [1,2,3,4,5]))
    end
end

function E_wgne!(m, vars, nr)
    wgne = vars["wgne"]; xgne = vars["xgne"]; pgne = vars["pgne"]
    @constraint(m, [d=1:nr], wgne[d] == xgne[d] + pgne[d])
end

function E_delINDTAX!(m, vars, nr)
    delINDTAX = vars["delINDTAX"]; delGDPINC = vars["delGDPINC"]
    @constraint(m, [d=1:nr], delINDTAX[d] == delGDPINC[d,4] + delGDPINC[d,5])
end

function E_delBUDG1!(m, vars, nr)
    delBUDG1 = vars["delBUDG1"]; delINDTAX = vars["delINDTAX"]; delVGDPEXP = vars["delVGDPEXP"]
    @constraint(m, [d=1:nr], delBUDG1[d] == delINDTAX[d] - delVGDPEXP[d,3])
end

function E_delBUDG2!(m, vars, nr)
    delBUDG2 = vars["delBUDG2"]; delGDPINC = vars["delGDPINC"]; delVGDPEXP = vars["delVGDPEXP"]
    @constraint(m, [d=1:nr], delBUDG2[d] == delGDPINC[d,5] - delVGDPEXP[d,3])
end

function E_delVGDPEXP!(m, vars, nr, ngdpexp)
    delVGDPEXP = vars["delVGDPEXP"]; delPGDPEXP = vars["delPGDPEXP"]; delXGDPEXP = vars["delXGDPEXP"]
    @constraint(m, [d=1:nr, k=1:ngdpexp], delVGDPEXP[d,k] == delPGDPEXP[d,k] + delXGDPEXP[d,k])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 38 — Labour market closure
# ═══════════════════════════════════════════════════════════════════════════
function E_labslack!(m, vars, no)
    realwage_id = vars["realwage_id"]; xlab_id = vars["xlab_id"]; flabsup_id = vars["flabsup_id"]
    @constraint(m, [o=1:no], realwage_id[o] == 2.0 * xlab_id[o] + flabsup_id[o])
end

function E_flab_i!(m, vars, na, nr, no)
    realwage_i = vars["realwage_i"]; xlab_i = vars["xlab_i"]
    flabsupA = vars["flabsupA"]; labslack = vars["labslack"]
    @constraint(m, [o=1:no, d=1:nr],
        realwage_i[o,d] == (1.0/1.1) * xlab_i[o,d] + flabsupA[o,d] + labslack[o])
end

function E_realwage!(m, vars, na, nr, no)
    realwage = vars["realwage"]; flab = vars["flab"]
    flab_i = vars["flab_i"]; flab_io = vars["flab_io"]
    flab_id = vars["flab_id"]; flab_iod = vars["flab_iod"]
    @constraint(m, [i=1:na, o=1:no, d=1:nr],
        realwage[i,o,d] == flab[i,o,d] + flab_i[o,d] + flab_io[d] + flab_id[o] + flab_iod)
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 39 — Household closure
# ═══════════════════════════════════════════════════════════════════════════
function E_fhou!(m, vars, na, nr, ns, nu)
    whouhtot = vars["whouhtot"]; wlab_io = vars["wlab_io"]
    fhou = vars["fhou"]; houslack = vars["houslack"]
    @constraint(m, [d=1:nr], whouhtot[d] == wlab_io[d] + fhou[d] + houslack)
end

function E_natfhou!(m, vars, na, nr, ns, nu, params)
end
