export E_pimp!, E_pbasic!, E_ppur!, E_tuser!
export E_ppur_s!, E_phou!, E_xint!, E_xhou!, E_xinv!
export E_aint_s!, E_xint_s!, E_pint!
export E_xlab!, E_plab_o!, E_wlab_o!
export E_xlab_o!, E_pcap!, E_plnd!, E_pprim!, E_xprim!, E_aprim!, E_alab_o!, E_wprim!
export E_pvar!, E_pcst!, E_delPTX!, E_ptot!
export E_xsub!, E_xlux!, E_xhouh_s_agg!, E_alux!, E_asub!, E_wlux!, E_phouhtot!, E_whouhtot!
export E_xhoutot!, E_phoutot!
export E_xinvi!, E_pinvest!, E_pinvitot!
export E_gret!, E_xinvitot!, E_ggro!, E_fgret!, E_finv2!
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
export E_plab_id!, E_realwage_id!, E_xlab_id!, E_wlab_id!, E_rlab_id!
export E_plab_io!, E_xlab_io!, E_realwage_io!, E_wlab_io!, E_rlab_io!
export E_delGDPINCa!, E_delGDPINCb!, E_delGDPINCc!, E_delGDPINCd!, E_delGDPINCe!, E_wgdpinc!
export E_delXGDPEXPa_setup!, E_delXGDPEXPb!, E_delXGDPEXPc!
export E_xgdpexp!, E_delPGDPEXPa!, E_delPGDPEXPb!, E_delPGDPEXPc!
export E_pgdpexp!, E_wgdpexp!, E_wgdpdiff!
export E_xgne!, E_pgne!, E_wgne!
export E_delINDTAX!, E_delBUDG1!, E_delBUDG2!, E_delVGDPEXP!
export E_labslack!, E_flab_i!, E_realwage!
export E_fhou!, E_fhou2!, E_natfhou!
export E_plab_o_setup!, INVEST_setup!
export TAX_PUR_setup!, STOCKS_setup!

# ── Individual Equation Functions ──────────────────────────────────────────
# LEVELS form (Step 5c course correction, see PLAN.md). Each function writes
# the genuine nonlinear identity whose first-order log-differential around
# the benchmark reproduces TERM.TAB's original %-change linear equation
# (verified by hand for every equation below). Benchmark-level convention:
#   - Multiplicative "powers"/shifters/price INDEXES benchmark = 1
#     (pdom,pimp,phi,pfimp,pbasic,puse,tuser*,ppur,ppur_s,phou,pinvest,
#      atot,aint_s,bint_*,alab_o,acap,alnd,aprim,blab*,bprim*,
#      plab,plab_o,wlab_o,pcap,plnd,pprim,wprim,pint,pvar/atot-ratio parts)
#   - Genuine flow QUANTITIES benchmark = own benchmark data value
#     (xint,xhou,xinv,xint_s,xhou_s,xinv_s,xlab,xlab_o,xcap,xlnd,xprim)
#   - `xtot` is the exception: a pure real-activity INDEX (benchmark = 1)
#     shared across multiple aggregates with different weight coefficients
#     (VTOT, VCST, VARCST, PUR_S-ratios, PRIM-ratio) — so `ptot`,`pcst`,
#     `pvar` instead carry the corresponding benchmark VALUE (VTOT,VCST,
#     VARCST) rather than 1. `delPTX` is a genuine tax-revenue LEVEL
#     (benchmark = PRODTAX); `delPTXRATE` stays additive (benchmark = 0,
#     a decimal ad-valorem rate shift, unchanged from the %-change form).
# `ces()`/`ces_calibrate()` (ces_helper.jl) implement genuine CES/CET nests;
# equations combining fixed-technology (Leontief, no substitution
# elasticity) aggregates use direct value/product identities instead.

const _TINY = 1e-9

# Smallest benchmark flow for which a `log(price*flow + _TINY)` equation still
# determines its price. In such a row the derivative w.r.t. the PRICE is
# X/(p*X + _TINY) ≈ 1, while the derivative w.r.t. the QUANTITY is
# p/(p*X + _TINY) ≈ 1/X. Row equilibration divides by the latter, so the price
# column is left with relative weight ≈ X and the Newton step for that price is
# amplified by 1/X. Anything within a few orders of magnitude of _TINY is
# therefore a price the system cannot see.
#
# Measured (25x6, blabnat shock): PUR_S[5,16,6] = 9.95e-8 Rp bn — about 100
# rupiah of annual flow, against a median non-zero cell of 163.6 Rp bn —
# produced a Newton step of -769 on a price index whose benchmark is 1. That one
# cell throttled fraction-to-boundary to α ≈ 1.2e-3 (next limiter: 0.344) and
# drove its own log argument negative, so ‖F‖ was NaN for every α ≥ 0.01.
#
# 1e5*_TINY caps the amplification at ~1e4 and pins 26 of 4350 PUR_S cells, all
# carrying ≤ 1e-4 Rp bn (≤ 100,000 rupiah) — economically nil, numerically fatal.
const _MIN_PRICED_FLOW = 1e5 * _TINY

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 6 — Basic prices
# ═══════════════════════════════════════════════════════════════════════════
function E_pimp!(m, vars, na, nr, params)
    pimp = vars["pimp"]; pfimp = vars["pfimp"]; phi = vars["phi"]
    @constraint(m, [c=1:na, r=1:nr], log(pimp[c,r]) == log(pfimp[c]) + log(phi))
end

function E_pbasic!(m, vars, na, nr, ns)
    pbasic = vars["pbasic"]; pdom = vars["pdom"]; pimp = vars["pimp"]
    for s_dom in 1:ns
        p = s_dom == 1 ? pdom : pimp
        @constraint(m, [c=1:na, r=1:nr], pbasic[c,s_dom,r] == p[c,r])
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 7 — Purchaser prices (multiplicative tax "powers")
# ═══════════════════════════════════════════════════════════════════════════
function E_tuser!(m, vars, na, nr, ns, nu)
    tuser = vars["tuser"]; tuser_ud = vars["tuser_ud"]; tuser_su = vars["tuser_su"]
    @constraint(m, [c=1:na, s=1:ns, u=1:nu, d=1:nr],
        log(tuser[c,s,u,d]) == log(tuser_ud[c,s]) + log(tuser_su[c,d]))
end

function E_ppur!(m, vars, na, nr, ns, nu)
    ppur = vars["ppur"]; puse = vars["puse"]
    tuser = vars["tuser"]; tuser_sud = vars["tuser_sud"]
    @constraint(m, [c=1:na, s=1:ns, u=1:nu, d=1:nr],
        log(ppur[c,s,u,d]) == log(puse[c,s,d]) + log(tuser[c,s,u,d]) + log(tuser_sud[c]))
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 8 — Armington (genuine CES dom/imp substitution)
# ═══════════════════════════════════════════════════════════════════════════
# Principle-B value identity: aggregate price × CURRENT aggregate quantity
# == sum of component price × CURRENT component quantity (Shephard's lemma /
# Euler's theorem — exact for any CES nest, not just at the linearization
# point). One identity per user category (industries, hou, inv, gov, exp).
function E_ppur_s!(m, vars, na, nr, ns, nu, params)
    ppur_s = vars["ppur_s"]; ppur = vars["ppur"]
    xint = vars["xint"]; xint_s = vars["xint_s"]
    xhou = vars["xhou"]; xhou_s = vars["xhou_s"]
    xinv = vars["xinv"]; xinv_s = vars["xinv_s"]
    xgov = vars["xgov"]; xgov_s = vars["xgov_s"]
    xexp = vars["xexp"]; xexp_s = vars["xexp_s"]
    u_hou, u_inv, u_gov, u_exp = na + 1, na + 2, na + 3, na + 4
    PUR_S = parent(params["PUR_S"]); SRCSHR = parent(params["SRCSHR"])
    # Zero- AND dust-flow cells leave ppur_s UNDETERMINED even though the equation
    # is emitted: d/d(ppur_s) log(ppur_s*X + _TINY) = X/(ppur_s*X + _TINY), which is
    # exactly 0 when the benchmark flow X is 0 and — after row equilibration —
    # negligible whenever X is merely small (see _MIN_PRICED_FLOW). The constraint
    # still HOLDS (both sides collapse to log(_TINY)) but contributes a numerically
    # empty Jacobian column — a null direction that makes SHOCK results non-unique
    # (it does not affect benchmark replication). Pin those cells at the benchmark
    # price instead. Distinct from the `|| continue` guards elsewhere: there the
    # equation is skipped, here it exists but carries no information.
    #
    # The threshold was 1e-10, which was INERT: the smallest strictly positive
    # PUR_S in the 25x6 database is 2.7e-8, so `> 1e-10` pinned only the cells that
    # were already exactly zero and never once did what this comment describes.
    # Either branch supplies one equation for one variable, so the system stays
    # square for any threshold.
    #
    # THE GUARDED BRANCH MUST NOT PIN THE PRICE AT A CONSTANT. It used to read
    # `fix(ppur_s, 1.0)`, and that is a degree-of-homogeneity breaker: scale the
    # numeraire by λ and every other price becomes λ×itself while these stay at
    # 1.0. On the 594 cells with PUR_S == 0 that is harmless — nothing multiplies
    # them — but 26 cells sit in the live band 0 < PUR_S <= _MIN_PRICED_FLOW,
    # where the frozen price carries a small but real weight into the cost
    # aggregations above it. Measured cost: pcap(5,5) and pcap(5,6) failed price
    # homogeneity by 3.7e-4 and 4.1e-4 at λ=1.10, matching (λ−1)×(pinned share of
    # purchases) to within 15% (test/diag_pinned_price_share.jl, task #38).
    #
    # The replacement is the exact degenerate case of the identity it stands in
    # for. Holding the source split at its benchmark shares — SRCSHR sums to 1
    # over s by construction, including on dead cells where it is 0.5/0.5 — the
    # value identity ppur_s·X_s = Σ_s ppur·X collapses to a share-weighted average
    # of the component prices. That is homogeneous of degree 1 in prices, so it
    # scales with λ as it must; it is linear and O(1)-scaled, so it still gives
    # the Newton step a clean row instead of the dust-flow log row this branch
    # exists to avoid; and at the benchmark every ppur = 1, so it returns
    # ppur_s = 1 — exactly what the old `fix` set. Benchmark replication is
    # unchanged; only the response to a price rescaling differs.
    for c in 1:na, i in 1:na, d in 1:nr
        if PUR_S[c,i,d] > _MIN_PRICED_FLOW
            @constraint(m, log(ppur_s[c,i,d] * xint_s[c,i,d] + _TINY) ==
                log(sum(xint[c,s,i,d] * ppur[c,s,i,d] for s in 1:ns) + _TINY))
        else
            @constraint(m, ppur_s[c,i,d] ==
                sum(SRCSHR[c,s,i,d] * ppur[c,s,i,d] for s in 1:ns))
        end
    end
    for (u, xq, xs) in ((u_hou, xhou, xhou_s), (u_inv, xinv, xinv_s),
                        (u_gov, xgov, xgov_s), (u_exp, xexp, xexp_s))
        for c in 1:na, d in 1:nr
            if PUR_S[c,u,d] > _MIN_PRICED_FLOW
                @constraint(m, log(ppur_s[c,u,d] * xs[c,d] + _TINY) ==
                    log(sum(xq[c,s,d] * ppur[c,s,u,d] for s in 1:ns) + _TINY))
            else
                @constraint(m, ppur_s[c,u,d] ==
                    sum(SRCSHR[c,s,u,d] * ppur[c,s,u,d] for s in 1:ns))
            end
        end
    end
end

function E_phou!(m, vars, na, nr)
    phou = vars["phou"]; ppur_s = vars["ppur_s"]
    u_hou = na + 1
    @constraint(m, [c=1:na, d=1:nr], phou[c,d] == ppur_s[c,u_hou,d])
end

function E_xint!(m, vars, na, nr, ns, sigmadomimp, params)
    xint = vars["xint"]; xint_s = vars["xint_s"]; ppur = vars["ppur"]
    ALPHA = params["ALPHA_ARMINT"]; GAMMA = params["GAMMA_ARMINT"]
    for c in 1:na, i in 1:na, d in 1:nr
        αv = ALPHA[c,:,i,d]
        if any(αv .> 0)
            dem = ces(xint_s[c,i,d], [ppur[c,s,i,d] for s in 1:ns], αv, sigmadomimp[c], GAMMA[c,i,d])
            @constraint(m, [s=1:ns], xint[c,s,i,d] == dem[s])
        else
            @constraint(m, [s=1:ns], xint[c,s,i,d] == 0.0)
        end
    end
end

function E_xhou!(m, vars, na, nr, ns, sigmadomimp, params)
    xhou = vars["xhou"]; xhou_s = vars["xhou_s"]; ppur = vars["ppur"]
    u_hou = na + 1
    ALPHA = params["ALPHA_ARMINT"]; GAMMA = params["GAMMA_ARMINT"]
    for c in 1:na, d in 1:nr
        αv = ALPHA[c,:,u_hou,d]
        if any(αv .> 0)
            dem = ces(xhou_s[c,d], [ppur[c,s,u_hou,d] for s in 1:ns], αv, sigmadomimp[c], GAMMA[c,u_hou,d])
            @constraint(m, [s=1:ns], xhou[c,s,d] == dem[s])
        else
            @constraint(m, [s=1:ns], xhou[c,s,d] == 0.0)
        end
    end
end

function E_xinv!(m, vars, na, nr, ns, sigmadomimp, params)
    xinv = vars["xinv"]; xinv_s = vars["xinv_s"]; ppur = vars["ppur"]
    u_inv = na + 2
    ALPHA = params["ALPHA_ARMINT"]; GAMMA = params["GAMMA_ARMINT"]
    for c in 1:na, d in 1:nr
        αv = ALPHA[c,:,u_inv,d]
        if any(αv .> 0)
            dem = ces(xinv_s[c,d], [ppur[c,s,u_inv,d] for s in 1:ns], αv, sigmadomimp[c], GAMMA[c,u_inv,d])
            @constraint(m, [s=1:ns], xinv[c,s,d] == dem[s])
        else
            @constraint(m, [s=1:ns], xinv[c,s,d] == 0.0)
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 9 — Intermediate demands (Leontief: no substitution across c)
# ═══════════════════════════════════════════════════════════════════════════
function E_aint_s!(m, vars, na, nr)
    aint_s = vars["aint_s"]; bint_scd = vars["bint_scd"]; bint_s = vars["bint_s"]
    @constraint(m, [c=1:na, i=1:na, d=1:nr],
        log(aint_s[c,i,d]) == log(bint_scd[i]) + log(bint_s[c,i,d]))
end

# NOTE: TERM.TAB's E_xint_s carries an extra -0.15*{ppur_s+aint_s-pint}
# intermediate-substitution term ("works faster if pint is NOT substituted")
# that the existing %-change Julia port already omits; preserved here to
# match the established (already benchmark-verified) simplification.
function E_xint_s!(m, vars, na, nr, params)
    xint_s = vars["xint_s"]; atot = vars["atot"]
    aint_s = vars["aint_s"]; xtot = vars["xtot"]
    PUR_S = parent(params["PUR_S"])
    @constraint(m, [c=1:na, i=1:na, d=1:nr],
        xint_s[c,i,d] == PUR_S[c,i,d] * atot[i,d] * aint_s[c,i,d] * xtot[i,d])
end

function E_pint!(m, vars, na, nr, params)
    pint = vars["pint"]; ppur_s = vars["ppur_s"]; aint_s = vars["aint_s"]
    PUR_S = parent(params["PUR_S"]); PUR_CS = parent(params["PUR_CS"])
    for i in 1:na, d in 1:nr
        PUR_CS[i,d] > 1e-10 || continue
        rhs = sum(PUR_S[c,i,d] * ppur_s[c,i,d] * aint_s[c,i,d] for c in 1:na)
        @constraint(m, PUR_CS[i,d] * pint[i,d] == rhs)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 10 — Labour composition (genuine CES across occupations)
# ═══════════════════════════════════════════════════════════════════════════
function E_xlab!(m, vars, na, nr, no, sigmalab, params)
    xlab = vars["xlab"]; xlab_o = vars["xlab_o"]; plab = vars["plab"]
    ALPHA = params["ALPHA_LAB"]; GAMMA = params["GAMMA_LAB"]
    for i in 1:na, d in 1:nr
        αv = ALPHA[i,:,d]
        if any(αv .> 0)
            dem = ces(xlab_o[i,d], [plab[i,o,d] for o in 1:no], αv, sigmalab[i], GAMMA[i,d])
            @constraint(m, [o=1:no], xlab[i,o,d] == dem[o])
        else
            @constraint(m, [o=1:no], xlab[i,o,d] == 0.0)
        end
    end
end

# Principle-B: composite labour price/quantity value == sum of
# per-occupation values, using CURRENT (post-CES) xlab (not fixed V1LAB
# benchmark weights) since occupations substitute with elasticity sigmalab.
function E_plab_o!(m, vars, na, nr, no, params)
    plab_o = vars["plab_o"]; plab = vars["plab"]
    xlab = vars["xlab"]; xlab_o = vars["xlab_o"]
    @constraint(m, [i=1:na, d=1:nr],
        log(plab_o[i,d] * xlab_o[i,d] + _TINY) ==
        log(sum(plab[i,o,d] * xlab[i,o,d] for o in 1:no) + _TINY))
end

function E_wlab_o!(m, vars, na, nr, no, params)
    wlab_o = vars["wlab_o"]; plab_o = vars["plab_o"]; xlab_o = vars["xlab_o"]
    LAB_O = parent(params["LAB_O"])
    for i in 1:na, d in 1:nr
        LAB_O[i,d] > 1e-10 || continue
        @constraint(m, LAB_O[i,d] * wlab_o[i,d] == plab_o[i,d] * xlab_o[i,d])
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 11 — Factor demands (genuine CES between labour, capital, land)
# ═══════════════════════════════════════════════════════════════════════════
# ALPHA_FAC/GAMMA_FAC nest order [LAB_O, CAP, LND] (matches prepare_parameters.jl).
# alab_o/acap/alnd are augmenting-technical-change shifters (benchmark=1)
# dividing out the "effective" quantity fed into the CES demand system.
function E_xlab_o!(m, vars, na, nr, sigmaprim, params)
    xlab_o = vars["xlab_o"]; alab_o = vars["alab_o"]; acap = vars["acap"]; alnd = vars["alnd"]
    xprim = vars["xprim"]; plab_o = vars["plab_o"]; pcap = vars["pcap"]; plnd = vars["plnd"]
    ALPHA = params["ALPHA_FAC"]; GAMMA = params["GAMMA_FAC"]
    @constraint(m, [i=1:na, d=1:nr],
        xlab_o[i,d] == alab_o[i,d] * ces(xprim[i,d],
            [plab_o[i,d]*alab_o[i,d], pcap[i,d]*acap[i,d], plnd[i,d]*alnd[i,d]],
            ALPHA[i,:,d], sigmaprim[i], GAMMA[i,d])[1])
end

function E_pcap!(m, vars, na, nr, sigmaprim, params)
    xcap = vars["xcap"]; alab_o = vars["alab_o"]; acap = vars["acap"]; alnd = vars["alnd"]
    xprim = vars["xprim"]; plab_o = vars["plab_o"]; pcap = vars["pcap"]; plnd = vars["plnd"]
    ALPHA = params["ALPHA_FAC"]; GAMMA = params["GAMMA_FAC"]
    @constraint(m, [i=1:na, d=1:nr],
        xcap[i,d] == acap[i,d] * ces(xprim[i,d],
            [plab_o[i,d]*alab_o[i,d], pcap[i,d]*acap[i,d], plnd[i,d]*alnd[i,d]],
            ALPHA[i,:,d], sigmaprim[i], GAMMA[i,d])[2])
end

function E_plnd!(m, vars, na, nr, sigmaprim, params)
    xlnd = vars["xlnd"]; alab_o = vars["alab_o"]; acap = vars["acap"]; alnd = vars["alnd"]
    xprim = vars["xprim"]; plab_o = vars["plab_o"]; pcap = vars["pcap"]; plnd = vars["plnd"]
    ALPHA = params["ALPHA_FAC"]; GAMMA = params["GAMMA_FAC"]
    LND = parent(params["LND"])
    # Sectors with no benchmark land use have a zero land CES share, so land
    # demand is 0 whatever plnd does → zero Jacobian column, plnd undetermined
    # (and the row itself is dead, since xlnd is exogenous and fixed at 0 here).
    # Dropping the row and pinning the variable together keeps the system square.
    for i in 1:na, d in 1:nr
        if LND[i,d] > 1e-10
            @constraint(m,
                xlnd[i,d] == alnd[i,d] * ces(xprim[i,d],
                    [plab_o[i,d]*alab_o[i,d], pcap[i,d]*acap[i,d], plnd[i,d]*alnd[i,d]],
                    ALPHA[i,:,d], sigmaprim[i], GAMMA[i,d])[3])
        else
            fix(plnd[i,d], 1.0; force=true)
        end
    end
end

# Principle-B: factor composite value == sum of component values (labour
# composite, capital, land), using CURRENT quantities — the "a" augmentation
# shifters cancel exactly in this value identity (Euler's theorem).
function E_pprim!(m, vars, na, nr, no, params)
    pprim = vars["pprim"]; xprim = vars["xprim"]
    plab_o = vars["plab_o"]; xlab_o = vars["xlab_o"]
    pcap = vars["pcap"]; xcap = vars["xcap"]
    plnd = vars["plnd"]; xlnd = vars["xlnd"]
    @constraint(m, [i=1:na, d=1:nr],
        log(pprim[i,d] * xprim[i,d] + _TINY) == log(
            plab_o[i,d]*xlab_o[i,d] + pcap[i,d]*xcap[i,d] + plnd[i,d]*xlnd[i,d] + _TINY))
end

function E_xprim!(m, vars, na, nr, params)
    xprim = vars["xprim"]; xtot = vars["xtot"]
    atot = vars["atot"]; aprim = vars["aprim"]
    PRIM = parent(params["PRIM"])
    @constraint(m, [i=1:na, d=1:nr],
        xprim[i,d] == PRIM[i,d] * xtot[i,d] * atot[i,d] * aprim[i,d])
end

function E_aprim!(m, vars, na, nr)
    aprim = vars["aprim"]; bprimnat = vars["bprimnat"]
    bprim_d = vars["bprim_d"]; bprim = vars["bprim"]
    @constraint(m, [i=1:na, d=1:nr],
        log(aprim[i,d]) == log(bprimnat) + log(bprim_d[i]) + log(bprim[i,d]))
end

function E_alab_o!(m, vars, na, nr)
    alab_o = vars["alab_o"]; blabnat = vars["blabnat"]
    blab_d = vars["blab_d"]; blab = vars["blab"]
    @constraint(m, [i=1:na, d=1:nr],
        log(alab_o[i,d]) == log(blabnat) + log(blab_d[i]) + log(blab[i,d]))
end

# wprim: factor-cost value index normalized by the FIXED benchmark PRIM
# (distinct from pprim, which is normalized by the CURRENT xprim).
function E_wprim!(m, vars, na, nr, no, params)
    wprim = vars["wprim"]; pcap = vars["pcap"]; xcap = vars["xcap"]
    plnd = vars["plnd"]; xlnd = vars["xlnd"]
    plab = vars["plab"]; xlab = vars["xlab"]
    PRIM = parent(params["PRIM"])
    for i in 1:na, d in 1:nr
        PRIM[i,d] > 1e-10 || continue
        rhs = pcap[i,d]*xcap[i,d] + plnd[i,d]*xlnd[i,d]
        rhs += sum(plab[i,o,d]*xlab[i,o,d] for o in 1:no)
        @constraint(m, PRIM[i,d] * wprim[i,d] == rhs)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 12 — Output prices (Leontief value-added/intermediate combination)
# ═══════════════════════════════════════════════════════════════════════════
# pvar/pcst are unit-cost prices (scale-independent — no xtot term, matching
# TERM.TAB's own linear form), so the levels identity is a per-unit Leontief
# cost-dual: fixed-weight (VARCST/VCST) combination of component prices,
# scaled by the output-augmenting shifter atot.
function E_pvar!(m, vars, na, nr, params)
    pvar = vars["pvar"]; atot = vars["atot"]
    plab_o = vars["plab_o"]; alab_o = vars["alab_o"]; pint = vars["pint"]
    VARCST = parent(params["VARCST"]); LAB_O = parent(params["LAB_O"])
    PUR_CS = parent(params["PUR_CS"])
    for i in 1:na, d in 1:nr
        VARCST[i,d] > 1e-10 || continue
        rhs = atot[i,d] * (LAB_O[i,d]*plab_o[i,d]*alab_o[i,d] + PUR_CS[i,d]*pint[i,d])
        @constraint(m, pvar[i,d] * VARCST[i,d] == rhs)
    end
end

function E_pcst!(m, vars, na, nr, params)
    pcst = vars["pcst"]; atot = vars["atot"]
    aprim = vars["aprim"]; pprim = vars["pprim"]; pint = vars["pint"]
    VCST = parent(params["VCST"]); PRIM = parent(params["PRIM"])
    PUR_CS = parent(params["PUR_CS"])
    for i in 1:na, d in 1:nr
        VCST[i,d] > 1e-10 || continue
        rhs = atot[i,d] * (PRIM[i,d]*aprim[i,d]*pprim[i,d] + PUR_CS[i,d]*pint[i,d])
        @constraint(m, pcst[i,d] * VCST[i,d] == rhs)
    end
end

# delPTX is the ORDINARY CHANGE in production-tax revenue (TERM.TAB:269, and the
# update TERM.TAB:318 `(change) PRODTAX = delPTX`), benchmark = 0 — NOT a level.
# It equals (new tax) − (benchmark PRODTAX) where new tax = (PTXRATE+delPTXRATE)·
# (ex-tax cost value pcst·xtot·VCST). Since PTXRATE·VCST = PRODTAX, that is
# PRODTAX·(pcst·xtot − 1) + VCST·pcst·xtot·delPTXRATE — the exact nonlinear form
# of TERM.TAB:547-548's linearization, and 0 at the benchmark. (The earlier port
# wrongly "redefined delPTX as a level = PRODTAX", which forced delPTX=PRODTAX/VCST
# here while E_ptot below needed delPTX=0 — a self-inconsistent pair.)
function E_delPTX!(m, vars, na, nr, params)
    delPTX = vars["delPTX"]; xtot = vars["xtot"]; pcst = vars["pcst"]
    delPTXRATE = vars["delPTXRATE"]
    VCST = parent(params["VCST"])
    PRODTAX_v = parent(params["PRODTAX"])
    @constraint(m, [i=1:na, d=1:nr],
        delPTX[i,d] == PRODTAX_v[i,d] * (pcst[i,d]*xtot[i,d] - 1.0)
                       + VCST[i,d] * pcst[i,d] * xtot[i,d] * delPTXRATE[i,d])
end

# TERM.TAB:549-551 E_ptot (value identity VTOT = VCST + PRODTAX): the tax-inclusive
# output value equals the ex-tax cost value plus the (benchmark PRODTAX + its
# change delPTX). At the benchmark VTOT = VCST + PRODTAX ⇒ ptot = 1.
function E_ptot!(m, vars, na, nr, params)
    ptot = vars["ptot"]; xtot = vars["xtot"]; pcst = vars["pcst"]; delPTX = vars["delPTX"]
    VTOT = parent(params["VTOT"]); VCST = parent(params["VCST"])
    PRODTAX_v = parent(params["PRODTAX"])
    for i in 1:na, d in 1:nr
        VTOT[i,d] > 1e-10 || continue
        @constraint(m, ptot[i,d] * xtot[i,d] * VTOT[i,d] ==
            pcst[i,d] * xtot[i,d] * VCST[i,d] + PRODTAX_v[i,d] + delPTX[i,d])
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 13 — Household demands (LES/Stone-Geary)
#
# ELES splits benchmark household purchases HOUPUR(c,d) into a "subsistence"
# quantity (1-BLUX)*HOUPUR and a "luxury"/supernumerary quantity BLUX*HOUPUR
# (BLUX/SLUX are the standard Lluch(1973) ELES shares, prepared in
# prepare_parameters.jl). That decomposition means xsub/xlux do NOT share
# xhou_s's benchmark directly — they're genuine ELES-additive components with
# their OWN benchmarks XSUB0/XLUX0 = (1-BLUX)*HOUPUR / BLUX*HOUPUR, and
# xhou_s = xlux + xsub is a plain addition (Principle A), not the
# BLUX-weighted-average form the old %-change linearization used (that form
# is only the correct first-order expansion of this exact identity).
#
# alux is a second exception to the "shifter benchmark = 1" convention: it
# must absorb per-commodity variation that the per-region scalar wlux can't
# carry in a plain multiplicative identity, so its true benchmark is the
# nontrivial calibrated share ALUX0(c,d) = XLUX0(c,d)/WLUX0(d) (WLUX0(d) =
# aggregate benchmark supernumerary expenditure), both prepared above.
# wlux itself is a genuine nominal value with benchmark WLUX0(d) (not an
# index) — analogous to WLUX0 sharing HOUPUR's units.
#
# ahou_s is TABLO-`Omit`ted (origin/TERM.TAB:2575) — permanently fixed at its
# benchmark (1.0) and never shocked — so asub/alux are, in practice, always
# constant; the equations below still need to be genuinely nonlinear so they
# reduce correctly to that constant at ahou_s=1 and remain correct if that
# ever changes.
function E_asub!(m, vars, na, nr, params)
    asub = vars["asub"]; ahou_s = vars["ahou_s"]
    BUDGSHR_v = parent(params["BUDGSHR"])
    for d in 1:nr
        geomean = sum(BUDGSHR_v[k,1,d] * log(ahou_s[k,d] + _TINY) for k in 1:na)
        @constraint(m, [c=1:na], asub[c,d] == ahou_s[c,d] / exp(geomean))
    end
end

function E_alux!(m, vars, na, nr, params)
    alux = vars["alux"]; asub = vars["asub"]
    SLUX_v = parent(params["SLUX"])
    ALUX0_v = parent(params["ALUX0"])
    for d in 1:nr
        geomean = sum(SLUX_v[k,1,d] * log(asub[k,d] + _TINY) for k in 1:na)
        @constraint(m, [c=1:na], alux[c,d] == ALUX0_v[c,1,d] * asub[c,d] / exp(geomean))
    end
end

function E_xsub!(m, vars, na, nr, params)
    xsub = vars["xsub"]; nhou = vars["nhou"]; asub = vars["asub"]
    XSUB0_v = parent(params["XSUB0"])
    @constraint(m, [c=1:na, d=1:nr], xsub[c,d] == XSUB0_v[c,1,d] * nhou[d] * asub[c,d])
end

function E_xlux!(m, vars, na, nr)
    xlux = vars["xlux"]; phou = vars["phou"]; wlux = vars["wlux"]; alux = vars["alux"]
    @constraint(m, [c=1:na, d=1:nr], xlux[c,d] * phou[c,d] == wlux[d] * alux[c,d])
end

function E_xhouh_s_agg!(m, vars, na, nr)
    xhou_s = vars["xhou_s"]; xlux = vars["xlux"]; xsub = vars["xsub"]
    @constraint(m, [c=1:na, d=1:nr], xhou_s[c,d] == xlux[c,d] + xsub[c,d])
end

# xhouhtot(d) — real (constant-benchmark-price) volume index of total
# household consumption, benchmark=1. TABLO names this equation "E_wlux" (it
# is the equation used to back-solve the wlux/xhou_s/xhouh_s block jointly),
# but its literal content defines xhouhtot, so this Julia function is kept
# under that historical name for traceability to origin/TERM.TAB:659-660.
function E_wlux!(m, vars, na, nr, params)
    xhoutot = vars["xhoutot"]; xhou_s = vars["xhou_s"]
    HOUPUR_C_v = parent(params["HOUPUR_C"])
    for d in 1:nr
        @constraint(m, xhoutot[d] == sum(xhou_s[c,d] for c in 1:na) / HOUPUR_C_v[1,d])
    end
end

# phouhtot(d) — Laspeyres-style CPI: fixed-budget-share arithmetic mean of
# the (benchmark=1) sourcing prices phou(c,d) (Principle A).
function E_phouhtot!(m, vars, na, nr, params)
    phouhtot = vars["phouhtot"]; phou = vars["phou"]
    BUDGSHR_v = parent(params["BUDGSHR"])
    for d in 1:nr
        @constraint(m, phouhtot[d] == sum(BUDGSHR_v[c,1,d] * phou[c,d] for c in 1:na))
    end
end

# whouhtot(d) — total nominal household expenditure, a genuine VALUE
# (benchmark = HOUPUR_C(d)) built as benchmark x price-index x quantity-index.
function E_whouhtot!(m, vars, na, nr, params)
    whouhtot = vars["whouhtot"]; phouhtot = vars["phouhtot"]; xhoutot = vars["xhoutot"]
    HOUPUR_C_v = parent(params["HOUPUR_C"])
    @constraint(m, [d=1:nr], whouhtot[d] == HOUPUR_C_v[1,d] * phouhtot[d] * xhoutot[d])
end

# xhoutot/phoutot — national aggregation over household TYPES; with the
# single representative household type already collapsed (no h index),
# these coincide exactly with xhouhtot/phouhtot.
function E_xhoutot!(m, vars, na, nr)
    xhoutot = vars["xhoutot"]; xhouhtot = vars["xhouhtot"]
    @constraint(m, [d=1:nr], xhoutot[d] == xhouhtot[d])
end

function E_phoutot!(m, vars, na, nr)
    phoutot = vars["phoutot"]; phouhtot = vars["phouhtot"]
    @constraint(m, [d=1:nr], phoutot[d] == phouhtot[d])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 14 — Investment demands
#
# INVEST_C(c,i,d) is a fixed benchmark VALUE (industry i's investment
# commodity-composition), so xinvi is Leontief in levels: each commodity
# tracks industry i's total investment index xinvitot (benchmark=1, a real
# volume index — see Excerpt 15) in fixed proportion. pinvest/pinvitot stay
# benchmark=1 price indices exactly as before (Principle A/B already hold in
# levels form unchanged since all component prices share benchmark=1).
# ═══════════════════════════════════════════════════════════════════════════
function E_xinvi!(m, vars, na, nr, params)
    xinvi = vars["xinvi"]; xinvitot = vars["xinvitot"]
    INVEST_v = parent(params["INVEST"])  # raw V2PUR[c,i,d], not the c-summed INVEST_C[i,d]
    @constraint(m, [c=1:na, i=1:na, d=1:nr],
        xinvi[c,i,d] == INVEST_v[c,i,d] * xinvitot[i,d])
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

# TERM.TAB E_xinv_s (lines 738-740) is value-weighted:
#   INVEST_I(c,d)·xinv_s% = Σ_i INVEST(c,i,d)·xinvi% .
# The two sides have DIFFERENT benchmarks: xinv_s's is PUR_S[inv] (the purchaser
# source-composite it feeds in the Armington nest E_xinv and in E_xfinb, line
# 1217), while Σ_i xinvi's is INVEST_I = Σ_i INVEST (the investment-by-industry
# matrix, 2PUR). These come from different source data and do not reconcile —
# the 2PUR-vs-PUR[inv] gap. A literal xinv_s = Σ xinvi therefore violates the
# benchmark by (INVEST_I − PUR_S[inv]) (≈33570 worst). The faithful levels form
# is the ratio identity xinv_s/PUR_S[inv] = (Σ_i xinvi)/INVEST_I, i.e.
#   INVEST_I·xinv_s = PUR_S[inv]·Σ_i xinvi ,
# exact at the benchmark and linearizing back to TERM.TAB's E_xinv_s.
function E_xinv_s!(m, vars, na, nr, params)
    xinv_s = vars["xinv_s"]; xinvi = vars["xinvi"]
    INVEST_I = parent(params["INVEST_I"]); PUR_S = parent(params["PUR_S"])
    u_inv = na + 2
    for c in 1:na, d in 1:nr
        ii = INVEST_I[c,d]; ps = PUR_S[c,u_inv,d]
        if ii > 1e-10 && ps > 1e-10
            @constraint(m, xinv_s[c,d] / ps == sum(xinvi[c,i,d] for i in 1:na) / ii)
        else
            # No investment flow on one basis (or both): both benchmarks are 0,
            # so the plain sum is already exact — keep the literal form.
            @constraint(m, xinv_s[c,d] == sum(xinvi[c,i,d] for i in 1:na))
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 15 — Industry-specific investment (flexible-accelerator closure)
#
# gret/ggro/fgret/finv2 are rate/gap variables (rate-of-return gap, capital
# growth-rate gap) with no natural multiplicative-index reading — like
# delPTXRATE (Excerpt 12), they stay ADDITIVE (benchmark=0), built from
# log() of whichever levels index/quantity feeds them. xinvitot becomes a
# genuine real investment-volume INDEX (benchmark=1, log-additive with xcap's
# own %-deviation from its real benchmark XCAP0).
# ═══════════════════════════════════════════════════════════════════════════
function E_gret!(m, vars, na, nr)
    gret2 = vars["gret"]; pcap = vars["pcap"]; pinvitot = vars["pinvitot"]
    @constraint(m, [i=1:na, d=1:nr], gret2[i,d] == log(pcap[i,d] + _TINY) - log(pinvitot[i,d] + _TINY))
end

function E_xinvitot!(m, vars, na, nr, params)
    xinvitot = vars["xinvitot"]; xcap = vars["xcap"]; ggro = vars["ggro"]
    CAP_v = parent(params["CAP"])
    @constraint(m, [i=1:na, d=1:nr],
        log(xinvitot[i,d] + _TINY) - (log(xcap[i,d] + _TINY) - log(CAP_v[i,d] + _TINY)) == ggro[i,d])
end

function E_ggro!(m, vars, na, nr)
    ggro = vars["ggro"]; gret = vars["gret"]
    finv1 = vars["finv1"]; invslack = vars["invslack"]
    @constraint(m, [i=1:na, d=1:nr], ggro[i,d] == finv1[i,d] + 0.33 * (2.0 * gret[i,d] - invslack))
end

function E_fgret!(m, vars, na, nr)
    gret = vars["gret"]; fgret = vars["fgret"]; capslack = vars["capslack"]
    @constraint(m, [i=1:na, d=1:nr], gret[i,d] == fgret[i,d] + capslack)
end

function E_finv2!(m, vars, na, nr)
    xinvitot = vars["xinvitot"]; finv2 = vars["finv2"]; xgdpexp = vars["xgdpexp"]
    @constraint(m, [i=1:na, d=1:nr],
        log(xinvitot[i,d] + _TINY) == finv2[i,d] + log(xgdpexp[d] + _TINY))
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 16 — Government / Export / Inventories
#
# fgovtot/fgov/fgov_s/fgovgen are pure exogenous policy-shock instruments
# (closure-fixed at benchmark=0, cf. delPTXRATE) — xgov is a genuine
# quantity (benchmark=own data value) driven multiplicatively by their sum:
# XGOV = XGOV0 * exp(shocks). Same pattern for export demand (xexpd) and its
# constant-elasticity price response.
# ═══════════════════════════════════════════════════════════════════════════
function E_xgov!(m, vars, na, nr, ns, params)
    xgov = vars["xgov"]; fgovtot = vars["fgovtot"]
    fgov = vars["fgov"]; fgov_s = vars["fgov_s"]; fgovgen = vars["fgovgen"]
    XGOV0_v = parent(params["XGOV0"])
    @constraint(m, [c=1:na, s=1:ns, d=1:nr],
        xgov[c,s,d] == XGOV0_v[c,s,d] * exp(fgovtot[d] + fgov[c,s,d] + fgov_s[c,d] + fgovgen))
end

function E_xgov_s!(m, vars, na, nr, ns)
    xgov_s = vars["xgov_s"]; xgov = vars["xgov"]
    @constraint(m, [c=1:na, d=1:nr], xgov_s[c,d] == sum(xgov[c,s,d] for s in 1:ns))
end

function E_fgovtot2!(m, vars, na, nr)
    fgovtot = vars["fgovtot"]; fgovtot2 = vars["fgovtot2"]; xhoutot = vars["xhoutot"]
    @constraint(m, [d=1:nr], fgovtot[d] == fgovtot2[d] + log(xhoutot[d] + _TINY))
end

function E_fgovtot3!(m, vars, na, nr)
    fgovtot = vars["fgovtot"]; fgovtot3 = vars["fgovtot3"]; xgdpexp = vars["xgdpexp"]
    @constraint(m, [d=1:nr], fgovtot[d] == fgovtot3[d] + log(xgdpexp[d] + _TINY))
end

function E_pfexp!(m, vars, na, nr)
    pfexp = vars["pfexp"]; ppur = vars["ppur"]; phi = vars["phi"]
    u_exp = na + 4
    @constraint(m, [c=1:na, d=1:nr],
        pfexp[c,d] == log(ppur[c,1,u_exp,d] + _TINY) - log(phi + _TINY))
end

function E_xexpd!(m, vars, na, nr, exp_elast, params)
    xexpd = vars["xexpd"]; pfexp = vars["pfexp"]
    fpexp = vars["fpexp"]; fpexp_d = vars["fpexp_d"]
    fqexp = vars["fqexp"]; fqexp_d = vars["fqexp_d"]
    natfqexp = vars["natfqexp"]; natfpexp = vars["natfpexp"]
    XEXPD0_v = parent(params["XEXPD0"])
    @constraint(m, [c=1:na, d=1:nr],
        log(xexpd[c,d] + _TINY) - log(XEXPD0_v[c,d] + _TINY) ==
            natfqexp + fqexp[c,d] + fqexp_d[c] -
            exp_elast[c] * (pfexp[c,d] - fpexp[c,d] - fpexp_d[c] - natfpexp))
end

function E_xexp!(m, vars, na, nr, ns)
    xexp = vars["xexp"]; xexpd = vars["xexpd"]
    @constraint(m, [c=1:na, d=1:nr], xexp[c,1,d] == xexpd[c,d])
    @constraint(m, [c=1:na, d=1:nr], xexp[c,2,d] == 0.0)
end

function E_xexp_s!(m, vars, na, nr, ns)
    xexp_s = vars["xexp_s"]; xexp = vars["xexp"]
    @constraint(m, [c=1:na, d=1:nr], xexp_s[c,d] == sum(xexp[c,s,d] for s in 1:ns))
end

# xstocks is inventory investment — can be zero or negative at the benchmark,
# so (like delPTX) it is NOT wrapped in exp(); it stays the same additive
# log-index quantity the old %-change model used, just built from log(xtot).
function E_xstocks!(m, vars, na, nr)
    xstocks = vars["xstocks"]; xtot = vars["xtot"]; fxstocks = vars["fxstocks"]
    @constraint(m, [i=1:na, d=1:nr], xstocks[i,d] == log(xtot[i,d] + _TINY) + fxstocks[i,d])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 17 — Total regional demand
# ═══════════════════════════════════════════════════════════════════════════
# xint/xhou/xinv/xgov/xexp are all genuine quantities (benchmark = own data
# value), so their Divisia aggregates are now plain physical-quantity sums
# (Principle A, cf. xinv_s in Excerpt 14) — no value-weight machinery needed.
function E_xint_i!(m, vars, na, nr, ns)
    xint_i = vars["xint_i"]; xint = vars["xint"]
    @constraint(m, [c=1:na, s=1:ns, d=1:nr], xint_i[c,s,d] == sum(xint[c,s,i,d] for i in 1:na))
end

# E_xuse (Excerpt 17, TERM.TAB:855-861): USE_U·xuse% = USE_I·xint_i% +
# USE(hou)·xhou% + USE(inv)·xinv% + USE(gov)·xgov% + USE(exp)·xexp%, i.e. a
# BASIC-value-weighted aggregation of purchaser-basis final-user quantities.
# The literal levels form `xuse = Σ_user xuser` is the "literal-equality" bug
# class (cf. E_pdomA/E_xinv_s): it silently assumes every benchmark matches,
# but xuse's bmk is basic (TRADE_R) while each xuser's bmk is purchaser (PUR).
# The faithful levels form scales each var by its own benchmark and weights by
# the TERM.TAB basic coefficient:
#   (USE_U/TRADE_R)·xuse = Σ_user (USE_user/PUR_user)·xuser.
# At the benchmark LHS = USE_U and RHS = USE_I + Σ_fin USE(fin) = USE_U (the
# basic-value identity USE_U = Σ_u USE), so it is EXACT with no imbalance
# constant, AND its linearization reproduces TERM.TAB coefficient-for-coefficient.
function E_xuse!(m, vars, na, nr, ns, params)
    xuse = vars["xuse"]; xint_i = vars["xint_i"]
    xhou = vars["xhou"]; xinv = vars["xinv"]
    xgov = vars["xgov"]; xexp = vars["xexp"]
    USE     = parent(params["USE"])      # na×ns×nu×nr, basic values
    USE_U   = parent(params["USE_U"])    # na×ns×nr
    USE_I   = parent(params["USE_I"])    # na×ns×nr
    PUR     = parent(params["PUR"])      # na×ns×nu×nr, purchaser values
    TRADE_R = parent(params["TRADE_R"])  # na×ns×nr  (= xuse benchmark)
    u_hou = na + 1; u_inv = na + 2; u_gov = na + 3; u_exp = na + 4
    ϵ = 1e-10
    for c in 1:na, s in 1:ns, d in 1:nr
        tr = TRADE_R[c,s,d]
        if tr <= ϵ || USE_U[c,s,d] <= ϵ
            # No basic absorption of this good/source in d: pin the (zero-flow)
            # composite to the raw sum — harmless, keeps the variable determined.
            @constraint(m, xuse[c,s,d] ==
                xint_i[c,s,d] + xhou[c,s,d] + xinv[c,s,d] + xgov[c,s,d] + xexp[c,s,d])
            continue
        end
        α = USE_U[c,s,d] / tr
        xinti0 = sum(PUR[c,s,i,d] for i in 1:na)            # xint_i benchmark
        β_int = xinti0 > ϵ ? USE_I[c,s,d] / xinti0 : 0.0
        ph = PUR[c,s,u_hou,d]; pv = PUR[c,s,u_inv,d]
        pg = PUR[c,s,u_gov,d]; pe = PUR[c,s,u_exp,d]
        β_hou = ph > ϵ ? USE[c,s,u_hou,d] / ph : 0.0
        β_inv = pv > ϵ ? USE[c,s,u_inv,d] / pv : 0.0
        β_gov = pg > ϵ ? USE[c,s,u_gov,d] / pg : 0.0
        β_exp = pe > ϵ ? USE[c,s,u_exp,d] / pe : 0.0
        @constraint(m, α * xuse[c,s,d] ==
            β_int*xint_i[c,s,d] + β_hou*xhou[c,s,d] + β_inv*xinv[c,s,d] +
            β_gov*xgov[c,s,d] + β_exp*xexp[c,s,d])
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 19 — Delivering goods (margins, Leontief)
#
# xtradmar/xtrad are genuine quantities (benchmark TMAR/TRADE); atradmar is an
# additive benchmark=0 technology-shift shock. pdelivrd/pbasic/psuppmar_p are
# benchmark=1 price levels. Direct log-linear translation of TERM.TAB's
# E_xtradmar/E_pdelivrd (cf. E_xexpd! precedent, Excerpt 16): genuine
# quantity/price terms become log(X/X0), shocks stay additive.
# ═══════════════════════════════════════════════════════════════════════════
function E_xtradmar_na!(m, vars, na, nr, ns, nm, params)
    xtradmar = vars["xtradmar"]; xtrad = vars["xtrad"]; atradmar = vars["atradmar"]
    TMAR = parent(params["TMAR"]); TRADE = parent(params["TRADE"])
    for c in 1:na, s in 1:ns, m_ix in 1:nm, r in 1:nr, d in 1:nr
        if TMAR[c,s,m_ix,r,d] > 1e-10
            @constraint(m,
                log(xtradmar[c,s,m_ix,r,d] + _TINY) - log(TMAR[c,s,m_ix,r,d] + _TINY) ==
                log(xtrad[c,s,r,d] + _TINY) - log(TRADE[c,s,r,d] + _TINY) + atradmar[c,s,m_ix,r,d])
        else
            JuMP.fix(xtradmar[c,s,m_ix,r,d], 0.0; force=true)
        end
    end
end

function E_pdelivrd!(m, vars, na, nr, ns, nm, params)
    pdelivrd = vars["pdelivrd"]; pbasic = vars["pbasic"]
    psuppmar_p = vars["psuppmar_p"]; atradmar = vars["atradmar"]
    BASSHR = parent(params["BASSHR"]); MARSHR = parent(params["MARSHR"])
    for c in 1:na, s in 1:ns, r in 1:nr, d in 1:nr
        lhs = BASSHR[c,s,r,d] * log(pbasic[c,s,r] + _TINY)
        for mh in 1:nm
            lhs += MARSHR[c,s,mh,r,d] * (log(psuppmar_p[mh,r,d] + _TINY) + atradmar[c,s,mh,r,d])
        end
        @constraint(m, log(pdelivrd[c,s,r,d] + _TINY) == lhs)
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
        rhs = sum(DELIVRD[c,s,r,d] * (log(pdelivrd[c,s,r,d] + _TINY) + atrad[c,s,r,d]) for r in 1:nr)
        @constraint(m, DELIVRD_R[c,s,d] * log(puse[c,s,d] + _TINY) == rhs)
    end
end

# avesrctwist/srctwist are both additive benchmark=0 shocks — unchanged from
# the original %-change form (fixed-weight plain average of additive terms).
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
    TRADE = parent(params["TRADE"]); TRADE_R = parent(params["TRADE_R"])
    for c in 1:na, s in 1:ns, r in 1:nr, d in 1:nr
        if TRADE[c,s,r,d] > 1e-10
            @constraint(m,
                log(xtrad[c,s,r,d] + _TINY) - log(TRADE[c,s,r,d] + _TINY) - atrad[c,s,r,d] ==
                    log(xuse[c,s,d] + _TINY) - log(TRADE_R[c,s,d] + _TINY) + srctwist[c,s,r,d] - avesrctwist[c,s,d]
                    - SGDD[c] * (log(pdelivrd[c,s,r,d] + _TINY) + atrad[c,s,r,d] - log(puse[c,s,d] + _TINY)))
        else
            JuMP.fix(xtrad[c,s,r,d], 0.0; force=true)
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 21 — Margin supply
#
# TERM.TAB confirms (lines 1039-1068) the aggregation weights for
# xsuppmar_p/xsuppmar_d/xsuppmar_rd are exactly their components' own
# benchmarks (TRADMAR==TMAR, SUPPMAR==MARS) — Principle A plain-sum applies,
# no Dict/setup machinery needed. E_psuppmar_p/E_xsuppmar are genuine price-
# index/CES-demand duals (direct log-linear translation, matching E_xtrad!).
# NOTE: the prior %-change Julia port's E_xsuppmar! omitted TERM.TAB's
# SIGMAMAR elasticity term entirely (a pure Leontief, no CES response) — this
# conversion restores it, since it's the true TERM.TAB E_xsuppmar equation.
# ═══════════════════════════════════════════════════════════════════════════
function E_xsuppmar_p!(m, vars, na, nr, ns, nm)
    xsuppmar_p = vars["xsuppmar_p"]; xtradmar = vars["xtradmar"]
    @constraint(m, [m_ix=1:nm, r=1:nr, d=1:nr],
        xsuppmar_p[m_ix,r,d] == sum(xtradmar[c,s,m_ix,r,d] for c in 1:na, s in 1:ns))
end

function E_psuppmar_p!(m, vars, na, nr, nm, params)
    psuppmar_p = vars["psuppmar_p"]; pdom = vars["pdom"]; asuppmar = vars["asuppmar"]
    MARS = parent(params["MARS"]); SUPPMAR_P = parent(params["SUPPMAR_P"])
    for m_ix in 1:nm, r in 1:nr, d in 1:nr
        if SUPPMAR_P[m_ix,r,d] > 1e-10
            rhs = sum(MARS[m_ix,r,d,p_] * (log(pdom[m_ix,p_] + _TINY) + asuppmar[m_ix,r,d,p_]) for p_ in 1:nr)
            @constraint(m, SUPPMAR_P[m_ix,r,d] * log(psuppmar_p[m_ix,r,d] + _TINY) == rhs)
        else
            # No margin supplied on this route → the equation above is skipped and
            # psuppmar_p would be a free variable with no defining equation.
            fix(psuppmar_p[m_ix,r,d], 1.0; force=true)
        end
    end
end

function E_xsuppmar!(m, vars, na, nr, nm, params)
    xsuppmar = vars["xsuppmar"]; xsuppmar_p = vars["xsuppmar_p"]
    asuppmar = vars["asuppmar"]; pdom = vars["pdom"]; psuppmar_p = vars["psuppmar_p"]
    MARS = parent(params["MARS"]); SUPPMAR_P = parent(params["SUPPMAR_P"])
    SMAR = parent(params["SMAR"])  # SIGMAMAR(m) — CES elasticity of margin sourcing across p
    for m_ix in 1:nm, r in 1:nr, d in 1:nr, p_ in 1:nr
        if MARS[m_ix,r,d,p_] > 1e-10
            @constraint(m,
                log(xsuppmar[m_ix,r,d,p_] + _TINY) - log(MARS[m_ix,r,d,p_] + _TINY) - asuppmar[m_ix,r,d,p_] ==
                    log(xsuppmar_p[m_ix,r,d] + _TINY) - log(SUPPMAR_P[m_ix,r,d] + _TINY)
                    - SMAR[m_ix] * (log(pdom[m_ix,p_] + _TINY) + asuppmar[m_ix,r,d,p_] - log(psuppmar_p[m_ix,r,d] + _TINY)))
        else
            JuMP.fix(xsuppmar[m_ix,r,d,p_], 0.0; force=true)
        end
    end
end

function E_xsuppmar_d!(m, vars, na, nr, ns, nm)
    xsuppmar_d = vars["xsuppmar_d"]; xsuppmar = vars["xsuppmar"]
    @constraint(m, [m_ix=1:nm, r=1:nr, p_=1:nr],
        xsuppmar_d[m_ix,r,p_] == sum(xsuppmar[m_ix,r,d,p_] for d in 1:nr))
end

function E_xsuppmar_rd!(m, vars, na, nr, nm)
    xsuppmar_rd = vars["xsuppmar_rd"]; xsuppmar_d = vars["xsuppmar_d"]
    @constraint(m, [m_ix=1:nm, p_=1:nr],
        xsuppmar_rd[m_ix,p_] == sum(xsuppmar_d[m_ix,r,p_] for r in 1:nr))
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 22 — MAKE / CET (multi-product industries)
#
# xmake is a genuine CET nest (industry i transforms xtot[i,d], benchmark=1,
# into commodity-specific xmake[c,i,d], benchmark=MAKE[c,i,d]) — calibrated
# via ALPHA_MAKE/GAMMA_MAKE (ces_calibrate with -SCET, the CET dual, see
# prepare_parameters.jl). E_xtotA_B! is its Principle-B revenue-identity dual
# (envelope theorem: exact at any point, cf. E_pprim!/E_ppur_s!), NOT a
# log-geometric average. E_xcomA_B! is a plain physical-quantity sum (the
# OTHER direction of the MAKE matrix — MAKESHR2's weight is xmake's own
# benchmark exactly, Principle A). E_pmake! stays a direct log-linear
# translation of the ad hoc 0.05-coefficient reduced form (not a calibrated
# CES nest).
# ═══════════════════════════════════════════════════════════════════════════
function E_xmake!(m, vars, na, nr, params)
    xmake = vars["xmake"]; xtot = vars["xtot"]; pmake = vars["pmake"]
    ALPHA = params["ALPHA_MAKE"]; GAMMA = params["GAMMA_MAKE"]
    SCET = parent(params["SCET"]); MAKE = parent(params["MAKE"])
    for i in 1:na, d in 1:nr
        αv = ALPHA[i,:,d]
        if any(αv .> 0)
            dem = ces(xtot[i,d], [pmake[c,i,d] for c in 1:na], αv, -SCET[i], GAMMA[i,d])
            for c in 1:na
                if MAKE[c,i,d] > 1e-10
                    @constraint(m, xmake[c,i,d] == dem[c])
                else
                    JuMP.fix(xmake[c,i,d], 0.0; force=true)
                end
            end
        else
            @constraint(m, [c=1:na], xmake[c,i,d] == 0.0)
        end
    end
end

function E_xtotA_B!(m, vars, na, nr, params)
    ptot = vars["ptot"]; pmake = vars["pmake"]; xtot = vars["xtot"]; xmake = vars["xmake"]
    MAKE_C = parent(params["MAKE_C"])
    for i in 1:na, d in 1:nr
        MAKE_C[i,d] > 1e-10 || continue
        @constraint(m, ptot[i,d] * xtot[i,d] * MAKE_C[i,d] ==
            sum(pmake[c,i,d] * xmake[c,i,d] for c in 1:na))
    end
end

function E_xcomA_B!(m, vars, na, nr)
    xcom = vars["xcom"]; xmake = vars["xmake"]
    @constraint(m, [c=1:na, d=1:nr], xcom[c,d] == sum(xmake[c,i,d] for i in 1:na))
end

function E_pmake!(m, vars, na, nr, params)
    pmake = vars["pmake"]; pdom = vars["pdom"]; xmake = vars["xmake"]; xcom = vars["xcom"]
    MAKE = parent(params["MAKE"]); MAKE_I = parent(params["MAKE_I"])
    @constraint(m, [c=1:na, i=1:na, d=1:nr],
        log(pmake[c,i,d] + _TINY) == log(pdom[c,d] + _TINY)
            - 0.05 * (log(xmake[c,i,d] + _TINY) - log(MAKE[c,i,d] + _TINY)
                      - log(xcom[c,d] + _TINY) + log(MAKE_I[c,d] + _TINY)))
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 23 — Market clearing
#
# xtrad_d/xtrad_r aggregate xtrad using weight TRADE — exactly xtrad's own
# benchmark (TERM.TAB lines 1189-1194) — so both collapse to plain sums
# (Principle A). Market clearing (E_pdomA_sum! below) now implements BOTH of
# TERM.TAB's split forms: E_pdomA for non-margins (ratio form, exact at the
# MAKE_I-vs-TRADE_D imbalance) and E_pdomB for the margin commodities (Trade /
# Transport / InfoComm — xcom balances domestic demand plus margin supply,
# with the benchmark imbalance carried as an additive constant).
# ═══════════════════════════════════════════════════════════════════════════
function E_xtrad_d!(m, vars, na, nr, ns)
    xtrad_d = vars["xtrad_d"]; xtrad = vars["xtrad"]
    @constraint(m, [c=1:na, s=1:ns, r=1:nr], xtrad_d[c,s,r] == sum(xtrad[c,s,r,d] for d in 1:nr))
end

function E_xtrad_r!(m, vars, na, nr, ns)
    xtrad_r = vars["xtrad_r"]; xtrad = vars["xtrad"]
    @constraint(m, [c=1:na, s=1:ns, d=1:nr], xtrad_r[c,s,d] == sum(xtrad[c,s,r,d] for r in 1:nr))
end

function E_pdomA_sum!(m, vars, na, nr, nm, params)
    xcom = vars["xcom"]; xtrad_d = vars["xtrad_d"]; xsuppmar_rd = vars["xsuppmar_rd"]
    MAKE_I = parent(params["MAKE_I"]); TRADE_D = parent(params["TRADE_D"])
    SUPPMAR_RD = parent(params["SUPPMAR_RD"])

    # Map each of the 9 original margins (MAR) to its aggregated-commodity index
    # and group them by target. TMAR's margin axis is NOT aggregated, so nm stays
    # 9 while xcom is indexed by the 25 aggregated commodities: the margins land
    # on commodity 20 (Trade), 21 (Transport), 22 (InfoComm). Those commodities'
    # domestic output balances direct domestic demand PLUS the margin services
    # they supply (TERM.TAB E_pdomB, lines 1200-1203); every other commodity uses
    # E_pdomA (lines 1196-1198).
    ncom = length(COM)
    margins_of = Dict{Int,Vector{Int}}()
    for mi in 1:nm
        cpos = findfirst(==(MAR[mi]), COM)
        cpos === nothing && continue
        aggc = ncom == na ? cpos : SEC_MAP_185_to_25[cpos]
        push!(get!(margins_of, aggc, Int[]), mi)
    end

    for c in 1:na, r in 1:nr
        if haskey(margins_of, c)
            # E_pdomB in levels. The %-change form MAKE_I·xcom = TRADE_D·xtrad_d +
            # SUPPMAR_RD·xsuppmar_rd is the linearization of the flow balance
            # xcom = xtrad_d + Σ xsuppmar_rd. The benchmark data is NOT
            # output = domestic-demand + margin-supply balanced, so we carry the
            # imbalance K = MAKE_I − TRADE_D − Σ SUPPMAR_RD as an additive
            # constant; then the benchmark (xcom=MAKE_I, xtrad_d=TRADE_D,
            # xsuppmar_rd=SUPPMAR_RD) replicates exactly and the equation still
            # linearizes to TERM.TAB's E_pdomB.
            K = MAKE_I[c,r] - TRADE_D[c,1,r] - sum(SUPPMAR_RD[mi,r] for mi in margins_of[c])
            @constraint(m, xcom[c,r] ==
                xtrad_d[c,1,r] + sum(xsuppmar_rd[mi,r] for mi in margins_of[c]) + K,
                base_name = "E_pdomA_sum_$(c)_$(r)")
        else
            # E_pdomA in levels. TERM.TAB's xcom% = xtrad_d% is, in value-flow
            # variables (xcom bmk=MAKE_I, xtrad_d bmk=TRADE_D), the ratio identity
            # xcom/MAKE_I = xtrad_d/TRADE_D — exact at the benchmark whatever the
            # MAKE_I vs TRADE_D imbalance (e.g. Crops, whose surplus output is RoW
            # export/inventory the domestic trade matrix does not carry) and
            # linearizing back to xcom = xtrad_d. The literal equality xcom=xtrad_d
            # used before wrongly forced MAKE_I = TRADE_D, breaking benchmark
            # replication for Crops/Trade/Transport/InfoComm.
            mi_v = MAKE_I[c,r]; td_v = TRADE_D[c,1,r]
            nm_ = "E_pdomA_sum_$(c)_$(r)"
            if mi_v > 1e-10 && td_v > 1e-10
                @constraint(m, xcom[c,r] / mi_v == xtrad_d[c,1,r] / td_v, base_name = nm_)
            elseif mi_v > 1e-10
                @constraint(m, xtrad_d[c,1,r] == td_v, base_name = nm_)   # produced, no domestic trade
            elseif td_v > 1e-10
                @constraint(m, xcom[c,r] == mi_v, base_name = nm_)        # traded, not produced
            else
                @constraint(m, xcom[c,r] == xtrad_d[c,1,r], base_name = nm_)  # both zero
            end
        end
    end
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

# xfin is a bmk=1 QUANTITY INDEX (TERM.TAB:1209 "Final user quantity indices"),
# not a flow: TERM.TAB:1217-1224 divides the PUR_S-weighted component sum by the
# category total PUR_CS(u,d). xinv_s/xgov_s/xexp are genuine flow quantities
# (bmk = their PUR_S/PUR share), so Σ_c PUR_S·(xinv_s/PUR_S) = Σ_c xinv_s, and the
# equation is PUR_CS·xfin = Σ_c xinv_s ⇒ xfin bmk = PUR_CS/PUR_CS = 1. (An earlier
# port dropped the PUR_CS(u,d) LHS factor, wrongly making xfin a flow — which then
# forced compensating PUR_CS factors into E_wfin/E_natxfin/E_natwfin/E_delXGDPEXPa
# and left xfin("hou") inconsistent, since E_xfina keeps it a bmk=1 index.)
function E_xfinb!(m, vars, na, nr, params)
    xfin = vars["xfin"]; xinv_s = vars["xinv_s"]; PUR_CS = parent(params["PUR_CS"])
    u = na + 2
    for d in 1:nr
        PUR_CS[u,d] > 1e-10 || continue
        @constraint(m, PUR_CS[u,d] * xfin[2,d] == sum(xinv_s[c,d] for c in 1:na))
    end
end

function E_xfinc!(m, vars, na, nr, params)
    xfin = vars["xfin"]; xgov_s = vars["xgov_s"]; PUR_CS = parent(params["PUR_CS"])
    u = na + 3
    for d in 1:nr
        PUR_CS[u,d] > 1e-10 || continue
        @constraint(m, PUR_CS[u,d] * xfin[3,d] == sum(xgov_s[c,d] for c in 1:na))
    end
end

function E_xfind!(m, vars, na, nr, ns, params)
    xfin = vars["xfin"]; xexp = vars["xexp"]; PUR_CS = parent(params["PUR_CS"])
    u = na + 4
    for d in 1:nr
        PUR_CS[u,d] > 1e-10 || continue
        @constraint(m, PUR_CS[u,d] * xfin[4,d] == sum(xexp[c,s,d] for c in 1:na, s in 1:ns))
    end
end

const PUR_src_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
function PUR_src_setup!(na, nr, ns, nu, PUR)
    empty!(PUR_src_idx)
    for c in 1:na, s in 1:ns, u in 1:nu, d in 1:nr
        PUR_src_idx[(c,s,u,d)] = PUR[c,s,u,d]
    end
end

# wfin is a bmk=1 nominal-expenditure INDEX: TERM.TAB:1226 E_wfin is the
# %-change identity wfin = xfin + pfin, i.e. levels wfin = xfin·pfin (both
# bmk=1 indices). No PUR_CS factor — that only appeared in the earlier port
# because xfin had been mistranslated as a flow (see E_xfinb! header).
function E_wfin!(m, vars, na, nr, params)
    wfin = vars["wfin"]; xfin = vars["xfin"]; pfin = vars["pfin"]
    @constraint(m, [fi=1:4, d=1:nr], wfin[fi,d] == xfin[fi,d] * pfin[fi,d])
end

# National final-demand aggregates (natpfin/natxfin/natwfin). natxfin is a
# genuine national quantity level (Principle A: PUR_CS(u,d) exactly matches
# xfin(u,d)'s own benchmark, so the PUR_CS-weighted average collapses to a
# plain sum); natpfin stays a genuine PUR_CS-weighted average of the bmk=1
# price pfin (Category 2); natwfin pairs them multiplicatively (Principle B).
function E_natpfin!(m, vars, na, nr, params)
    natpfin = vars["natpfin"]; pfin = vars["pfin"]
    PUR_CS = parent(params["PUR_CS"])
    u_list = [na+1, na+2, na+3, na+4]
    for (fi, u) in enumerate(u_list)
        tot = sum(PUR_CS[u,d] for d in 1:nr)
        tot > 1e-10 || continue
        @constraint(m, tot * natpfin[fi] == sum(PUR_CS[u,d] * pfin[fi,d] for d in 1:nr))
    end
end

# natxfin is the PUR_CS-weighted national average of the bmk=1 index xfin
# (TERM.TAB:1234 Σ_d PUR_CS·[natxfin-xfin]=0 ⇒ natxfin bmk=1), exactly mirroring
# the already-correct E_natpfin. The earlier port's plain Σ_d xfin was the flow
# form (bmk=nr) — wrong once xfin is restored to a bmk=1 index.
function E_natxfin!(m, vars, na, nr, params)
    natxfin = vars["natxfin"]; xfin = vars["xfin"]
    PUR_CS = parent(params["PUR_CS"])
    u_list = [na+1, na+2, na+3, na+4]
    for (fi, u) in enumerate(u_list)
        tot = sum(PUR_CS[u,d] for d in 1:nr)
        tot > 1e-10 || continue
        @constraint(m, tot * natxfin[fi] == sum(PUR_CS[u,d] * xfin[fi,d] for d in 1:nr))
    end
end

function E_natwfin!(m, vars, na, nr, params)
    natwfin = vars["natwfin"]; natpfin = vars["natpfin"]; natxfin = vars["natxfin"]
    # TERM.TAB:1235 E_natwfin: natwfin = natpfin + natxfin (%-change) ⇒ levels
    # natwfin = natpfin·natxfin, all bmk=1 indices. No ΣPUR_CS factor (that was a
    # compensating artefact of the xfin-as-flow mistranslation; see E_xfinb!).
    @constraint(m, [fi=1:4], natwfin[fi] == natpfin[fi] * natxfin[fi])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 26 — Commodity tax revenues
#
# xint (own benchmark = PUR(c,s,i,d)) is normalized on the purchaser-price
# basis, while USE is on the basic/user-price basis — ppur and puse are both
# bmk=1 but sit on different absolute-dollar bases, so converting the CES
# output to USE-basis dollars needs the benchmark USE/PUR rescaling factor.
# Multiplying through by PUR avoids dividing by it. Verified by Taylor-
# expanding to first order: reproduces TERM.TAB's own linearized formula
# 0.01*TAX*(xint+puse) + 0.01*PUR*tuser exactly (tuser_sud is absorbed
# automatically since ppur already incorporates it via E_ppur!).
# ═══════════════════════════════════════════════════════════════════════════
function E_delTAXint!(m, vars, na, nr, ns, params)
    delTAXint = vars["delTAXint"]; xint = vars["xint"]; ppur = vars["ppur"]; puse = vars["puse"]
    for c in 1:na, s in 1:ns, i in 1:na, d in 1:nr
        PUR = get(PUR_idx, (c,s,i,d), 0.0)
        if PUR > 1e-10
            USE = get(USE_idx, (c,s,i,d), 0.0); TAX = get(TAX_idx, (c,s,i,d), 0.0)
            @constraint(m, delTAXint[c,s,i,d] * PUR ==
                xint[c,s,i,d] * (ppur[c,s,i,d] * PUR - USE * puse[c,s,d]) - TAX * PUR)
        else
            # No purchase flow → no tax base → the tax-revenue change is
            # structurally 0 for any shock. GEMPACK omits the cell; the levels
            # port pins it, else delTAXint(c,s,i,d) is a free variable with no
            # defining equation (under-determines the square system).
            fix(delTAXint[c,s,i,d], 0.0; force=true)
        end
    end
end

const TAX_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
const PUR_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
const USE_idx = Dict{Tuple{Int,Int,Int,Int}, Float64}()
function TAX_PUR_setup!(na, nr, ns, nu, TAX_data, PUR_data, USE_data)
    empty!(TAX_idx); empty!(PUR_idx); empty!(USE_idx)
    for c in 1:na, s in 1:ns, u in 1:nu, d in 1:nr
        TAX_idx[(c,s,u,d)] = TAX_data[c,s,u,d]
        PUR_idx[(c,s,u,d)] = PUR_data[c,s,u,d]
        USE_idx[(c,s,u,d)] = USE_data[c,s,u,d]
    end
end

function E_delTAXhou!(m, vars, na, nr, ns, params)
    delTAXhou = vars["delTAXhou"]; xhou = vars["xhou"]; ppur = vars["ppur"]; puse = vars["puse"]
    u_hou = na + 1
    for c in 1:na, s in 1:ns, d in 1:nr
        PUR = get(PUR_idx, (c,s,u_hou,d), 0.0)
        if PUR > 1e-10
            USE = get(USE_idx, (c,s,u_hou,d), 0.0); TAX = get(TAX_idx, (c,s,u_hou,d), 0.0)
            @constraint(m, delTAXhou[c,s,d] * PUR ==
                xhou[c,s,d] * (ppur[c,s,u_hou,d] * PUR - USE * puse[c,s,d]) - TAX * PUR)
        else
            fix(delTAXhou[c,s,d], 0.0; force=true)   # no tax base → 0 revenue change (see E_delTAXint!)
        end
    end
end

function E_delTAXinv!(m, vars, na, nr, ns, params)
    delTAXinv = vars["delTAXinv"]; xinv = vars["xinv"]; ppur = vars["ppur"]; puse = vars["puse"]
    u_inv = na + 2
    for c in 1:na, s in 1:ns, d in 1:nr
        PUR = get(PUR_idx, (c,s,u_inv,d), 0.0)
        if PUR > 1e-10
            USE = get(USE_idx, (c,s,u_inv,d), 0.0); TAX = get(TAX_idx, (c,s,u_inv,d), 0.0)
            @constraint(m, delTAXinv[c,s,d] * PUR ==
                xinv[c,s,d] * (ppur[c,s,u_inv,d] * PUR - USE * puse[c,s,d]) - TAX * PUR)
        else
            fix(delTAXinv[c,s,d], 0.0; force=true)   # no tax base → 0 revenue change (see E_delTAXint!)
        end
    end
end

function E_delTAXgov!(m, vars, na, nr, ns, params)
    delTAXgov = vars["delTAXgov"]; xgov = vars["xgov"]; ppur = vars["ppur"]; puse = vars["puse"]
    u_gov = na + 3
    for c in 1:na, s in 1:ns, d in 1:nr
        PUR = get(PUR_idx, (c,s,u_gov,d), 0.0)
        if PUR > 1e-10
            USE = get(USE_idx, (c,s,u_gov,d), 0.0); TAX = get(TAX_idx, (c,s,u_gov,d), 0.0)
            @constraint(m, delTAXgov[c,s,d] * PUR ==
                xgov[c,s,d] * (ppur[c,s,u_gov,d] * PUR - USE * puse[c,s,d]) - TAX * PUR)
        else
            fix(delTAXgov[c,s,d], 0.0; force=true)   # no tax base → 0 revenue change (see E_delTAXint!)
        end
    end
end

function E_delTAXexp!(m, vars, na, nr, ns, params)
    delTAXexp = vars["delTAXexp"]; xexp = vars["xexp"]; ppur = vars["ppur"]; puse = vars["puse"]
    u_exp = na + 4
    for c in 1:na, s in 1:ns, d in 1:nr
        PUR = get(PUR_idx, (c,s,u_exp,d), 0.0)
        if PUR > 1e-10
            USE = get(USE_idx, (c,s,u_exp,d), 0.0); TAX = get(TAX_idx, (c,s,u_exp,d), 0.0)
            @constraint(m, delTAXexp[c,s,d] * PUR ==
                xexp[c,s,d] * (ppur[c,s,u_exp,d] * PUR - USE * puse[c,s,d]) - TAX * PUR)
        else
            fix(delTAXexp[c,s,d], 0.0; force=true)   # no tax base → 0 revenue change (see E_delTAXint!)
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 27 — Primary Factor Aggregates
# ═══════════════════════════════════════════════════════════════════════════
# wlnd_i/wcap_i pair a genuine quantity (xlnd/xcap) with a bmk=1 price
# (plnd/pcap) multiplicatively (Principle B) — the old code additively
# combined them as if both were still %-changes.
function E_wlnd_i!(m, vars, na, nr, params)
    wlnd_i = vars["wlnd_i"]; plnd = vars["plnd"]; xlnd = vars["xlnd"]
    LND_I = parent(params["LND_I"])
    for d in 1:nr
        LND_I[d] > 1e-10 || continue
        @constraint(m, LND_I[d] * wlnd_i[d] == sum(plnd[i,d] * xlnd[i,d] for i in 1:na))
    end
end

function E_wcap_i!(m, vars, na, nr, params)
    wcap_i = vars["wcap_i"]; pcap = vars["pcap"]; xcap = vars["xcap"]
    CAP_I = parent(params["CAP_I"])
    for d in 1:nr
        CAP_I[d] > 1e-10 || continue
        @constraint(m, CAP_I[d] * wcap_i[d] == sum(pcap[i,d] * xcap[i,d] for i in 1:na))
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
        # xlnd already carries the LND flow (xlnd = LND·χ, χ=sectoral quantity
        # index), so the value-weighted quantity index is CAP_I·xcap_i = Σ xcap,
        # i.e. RHS weight = 1 (NOT LND[i,d] again). Plain-sum-over-i normalized
        # by LND_I gives the bmk=1 index. cf. E_xint_i/E_xcom (Principle A).
        @constraint(m, LND_I[d] * xlnd_i[d] == sum(xlnd[i,d] for i in 1:na))
    end
end

function E_xcap_i!(m, vars, na, nr, params)
    xcap_i = vars["xcap_i"]; xcap = vars["xcap"]
    CAP_I = parent(params["CAP_I"]); CAP = parent(params["CAP"])
    for d in 1:nr
        CAP_I[d] > 1e-10 || continue
        # xcap already carries the CAP flow, so RHS weight = 1 (see E_xlnd_i!).
        @constraint(m, CAP_I[d] * xcap_i[d] == sum(xcap[i,d] for i in 1:na))
    end
end

# realwage == plab/pfin(HOU) is a ratio of two bmk=1 variables (still bmk=1);
# pins realwage, leaving flab/flab_i/... (Excerpt 38) to solve for the wage
# shifters via E_realwage! — see the decomposition-residual note there.
function E_plab!(m, vars, na, nr, no)
    realwage = vars["realwage"]; plab = vars["plab"]; pfin = vars["pfin"]
    @constraint(m, [i=1:na, o=1:no, d=1:nr], realwage[i,o,d] * pfin[1,d] == plab[i,o,d])
end

function E_plab_i!(m, vars, na, nr, no, params)
    plab_i = vars["plab_i"]; plab = vars["plab"]
    LAB = parent(params["LAB"])  # raw V1LAB[i,o,d] — no separate _idx cache needed
    for o in 1:no, d in 1:nr
        lab_sum = sum(LAB[i,o,d] for i in 1:na)
        lab_sum > 0 || continue
        @constraint(m, lab_sum * plab_i[o,d] == sum(LAB[i,o,d] * plab[i,o,d] for i in 1:na))
    end
end

function E_realwage_i!(m, vars, na, nr, no, params)
    realwage_i = vars["realwage_i"]; realwage = vars["realwage"]
    LAB = parent(params["LAB"])
    for o in 1:no, d in 1:nr
        rw_sum = sum(LAB[i,o,d] for i in 1:na)
        rw_sum > 0 || continue
        @constraint(m, rw_sum * realwage_i[o,d] == sum(LAB[i,o,d] * realwage[i,o,d] for i in 1:na))
    end
end

function E_xlab_i!(m, vars, na, nr, no, params)
    xlab_i = vars["xlab_i"]; xlab = vars["xlab"]
    LAB = parent(params["LAB"])
    for o in 1:no, d in 1:nr
        lab_sum = sum(LAB[i,o,d] for i in 1:na)
        lab_sum > 0 || continue
        # xlab already carries the LAB flow, so RHS weight = 1 (see E_xlnd_i!).
        @constraint(m, lab_sum * xlab_i[o,d] == sum(xlab[i,o,d] for i in 1:na))
    end
end

function E_wlab_i!(m, vars, na, nr, no)
    wlab_i = vars["wlab_i"]; xlab_i = vars["xlab_i"]; plab_i = vars["plab_i"]
    @constraint(m, [o=1:no, d=1:nr], wlab_i[o,d] == xlab_i[o,d] * plab_i[o,d])
end

function E_rlab_i!(m, vars, na, nr, no)
    rlab_i = vars["rlab_i"]; xlab_i = vars["xlab_i"]; realwage_i = vars["realwage_i"]
    @constraint(m, [o=1:no, d=1:nr], rlab_i[o,d] == xlab_i[o,d] * realwage_i[o,d])
end

function E_plab_id!(m, vars, na, nr, no, params)
    plab_id = vars["plab_id"]; plab_i = vars["plab_i"]
    SLAB_ID = parent(params["SLAB_ID"])
    @constraint(m, [o=1:no], plab_id[o] == sum(SLAB_ID[o,d] * plab_i[o,d] for d in 1:nr))
end

function E_realwage_id!(m, vars, na, nr, no, params)
    realwage_id = vars["realwage_id"]; realwage_i = vars["realwage_i"]
    SLAB_ID = parent(params["SLAB_ID"])
    @constraint(m, [o=1:no], realwage_id[o] == sum(SLAB_ID[o,d] * realwage_i[o,d] for d in 1:nr))
end

function E_xlab_id!(m, vars, na, nr, no, params)
    xlab_id = vars["xlab_id"]; xlab_i = vars["xlab_i"]
    SLAB_ID = parent(params["SLAB_ID"])
    @constraint(m, [o=1:no], xlab_id[o] == sum(SLAB_ID[o,d] * xlab_i[o,d] for d in 1:nr))
end

function E_wlab_id!(m, vars, na, nr, no)
    wlab_id = vars["wlab_id"]; xlab_id = vars["xlab_id"]; plab_id = vars["plab_id"]
    @constraint(m, [o=1:no], wlab_id[o] == xlab_id[o] * plab_id[o])
end

function E_rlab_id!(m, vars, na, nr, no)
    rlab_id = vars["rlab_id"]; xlab_id = vars["xlab_id"]; realwage_id = vars["realwage_id"]
    @constraint(m, [o=1:no], rlab_id[o] == xlab_id[o] * realwage_id[o])
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
    @constraint(m, [d=1:nr], wlab_io[d] == xlab_io[d] * plab_io[d])
end

function E_rlab_io!(m, vars, na, nr)
    rlab_io = vars["rlab_io"]; xlab_io = vars["xlab_io"]; realwage_io = vars["realwage_io"]
    @constraint(m, [d=1:nr], rlab_io[d] == xlab_io[d] * realwage_io[d])
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 28 — Income-side GDP
# ═══════════════════════════════════════════════════════════════════════════
# delGDPINC(*) is a GEMPACK (change) variable — a genuine absolute-dollar
# change from benchmark. wlnd_i/wcap_i/wlab_io are bmk=1 nominal-value
# ratios, so (ratio-1)*BENCHMARK gives exactly that dollar change.
function E_delGDPINCa!(m, vars, nr, params)
    delGDPINC = vars["delGDPINC"]; wlnd_i = vars["wlnd_i"]
    LND_I = parent(params["LND_I"])
    @constraint(m, [d=1:nr], delGDPINC[d,1] == LND_I[d] * (wlnd_i[d] - 1))
end

function E_delGDPINCb!(m, vars, nr, params)
    delGDPINC = vars["delGDPINC"]; wcap_i = vars["wcap_i"]
    CAP_I = parent(params["CAP_I"])
    # Slot 3, NOT 2. TERM.TAB indexes delGDPINC by CATEGORY NAME, and
    # GDPINCCAT = (Land, Labour, Capital, ProdTax, ComTax) — so E_delGDPINCb,
    # despite being the *second* equation, writes to "Capital" = slot 3
    # (TERM.TAB:1418-1419). Translating the equation letter positionally
    # (a→1, b→2, c→3) silently swapped Labour and Capital against GDPINCSUM,
    # which prepare_parameters.jl fills in the correct category order. Totals
    # were unaffected (wgdpinc sums over all k), so benchmark replication could
    # never catch it — only per-component GDP reporting would show it.
    @constraint(m, [d=1:nr], delGDPINC[d,3] == CAP_I[d] * (wcap_i[d] - 1))
end

function E_delGDPINCc!(m, vars, nr, params)
    delGDPINC = vars["delGDPINC"]; wlab_io = vars["wlab_io"]
    LAB_IO = parent(params["LAB_IO"])
    # Slot 2 ("Labour"), NOT 3 — see the note on E_delGDPINCb above
    # (TERM.TAB:1420-1421).
    @constraint(m, [d=1:nr], delGDPINC[d,2] == LAB_IO[d] * (wlab_io[d] - 1))
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
        @constraint(m, gdp_inc * wgdpinc[d] == gdp_inc + sum(delGDPINC[d,k] for k in 1:5))
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 29 — Expenditure-side GDP
# ═══════════════════════════════════════════════════════════════════════════
# GDPEXPCAT columns (matches prepare_parameters.jl / TERM.TAB's own Set order):
# 1=HOU 2=INV 3=GOV 4=STOCKS 5=EXP 6=Imports 7=RExports 8=RImports 9=NetMar
#
# delXGDPEXP/delPGDPEXP are (change) variables: genuine absolute-dollar
# changes from benchmark. For quantity components this is always
# (new genuine-quantity component total) - GDPEXPSUM(benchmark component
# total), with a sign flip for the two categories (Imports, RImports) whose
# GDPEXPSUM contribution is itself negative. For price components it is
# always sum(BENCHMARK_QUANTITY * (price_ratio - 1)) — automatically zero
# at benchmark, no separate subtraction needed.
function E_delXGDPEXPa_setup!(m, vars, na, nr, params)
    delXGDPEXP = vars["delXGDPEXP"]; xfin = vars["xfin"]
    PUR_CS = parent(params["PUR_CS"])
    col_list = [1, 2, 3, 5]  # HOU, INV, GOV, EXP (STOCKS=4 handled separately below)
    u_list = [na+1, na+2, na+3, na+4]
    for d in 1:nr, (fi, (col, u)) in enumerate(zip(col_list, u_list))
        # dollar change = benchmark value × (bmk=1 quantity index − 1), exactly
        # mirroring E_delPGDPEXPa's PUR_CS·(pfin−1). TERM.TAB:1464 uses
        # PUR_CS(u,d)·xfin(u,d) for this contribution; the earlier port's
        # xfin − PUR_CS assumed xfin was a flow (see E_xfinb! header).
        @constraint(m, delXGDPEXP[d,col] == PUR_CS[u,d] * (xfin[fi,d] - 1))
    end
end

# TERM.TAB:1465-1466 E_delXGDPEXPb: delXGDPEXP(d,"Stocks") = 0.01·Σ_i STOCKS·xstocks%.
# delXGDPEXP is an ordinary (change) dollar delta (bmk=0); STOCKS is the benchmark
# inventory value (weight); xstocks% is the volume %-change. Our xstocks is the LOG
# volume index (E_xstocks: xstocks=log(xtot)+fxstocks, bmk=0), so the volume ratio is
# exp(xstocks) and the fractional change is exp(xstocks)−1 — mirroring the price-side
# E_delPGDPEXPb's Σ STOCKS·(ptot−1). At the benchmark xstocks=0 ⇒ contribution 0.
# (The earlier port's Σ(xstocks−STOCKS) treated xstocks as a stock-value flow and
# never cancelled the benchmark Σ STOCKS, leaving a residual = benchmark inventory.)
function E_delXGDPEXPb!(m, vars, na, nr, params)
    delXGDPEXP = vars["delXGDPEXP"]; xstocks = vars["xstocks"]
    for d in 1:nr
        @constraint(m, delXGDPEXP[d,4] ==
            sum(get(STOCKS_idx, (i,d), 0.0) * (exp(xstocks[i,d]) - 1.0) for i in 1:na))
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
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for q in 1:nr
        @constraint(m, delXGDPEXP[q,6] == -sum(xtrad_d[c,2,q] for c in 1:na) - GDPEXPSUM[q,6])
    end
end

function E_delXGDPEXPd!(m, vars, na, nr, nm, params)
    delXGDPEXP = vars["delXGDPEXP"]; xsuppmar_d = vars["xsuppmar_d"]; xsuppmar = vars["xsuppmar"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for q in 1:nr
        rhs = sum(xsuppmar_d[m_,r,q] - sum(xsuppmar[m_,r,q,p] for p in 1:nr) for m_ in 1:nm, r in 1:nr)
        @constraint(m, delXGDPEXP[q,9] == rhs - GDPEXPSUM[q,9])
    end
end

function E_delXGDPEXPe!(m, vars, na, nr, ns, params)
    delXGDPEXP = vars["delXGDPEXP"]; xtrad_d = vars["xtrad_d"]; xtrad = vars["xtrad"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for q in 1:nr
        rhs = sum(xtrad_d[c,s,q] - xtrad[c,s,q,q] for c in 1:na, s in 1:ns)
        @constraint(m, delXGDPEXP[q,7] == rhs - GDPEXPSUM[q,7])
    end
end

function E_delXGDPEXPf!(m, vars, na, nr, ns, params)
    delXGDPEXP = vars["delXGDPEXP"]; xtrad_r = vars["xtrad_r"]; xtrad = vars["xtrad"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for q in 1:nr
        rhs = sum(xtrad_r[c,s,q] - xtrad[c,s,q,q] for c in 1:na, s in 1:ns)
        @constraint(m, delXGDPEXP[q,8] == -rhs - GDPEXPSUM[q,8])
    end
end

function E_xgdpexp!(m, vars, nr, params)
    xgdpexp = vars["xgdpexp"]; delXGDPEXP = vars["delXGDPEXP"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for d in 1:nr
        gdp_exp = sum(GDPEXPSUM[d,k] for k in 1:9)
        gdp_exp > 1e-10 || continue
        @constraint(m, gdp_exp * xgdpexp[d] == gdp_exp + sum(delXGDPEXP[d,k] for k in 1:9))
    end
end

function E_delPGDPEXPa!(m, vars, na, nr, nu, params)
    delPGDPEXP = vars["delPGDPEXP"]; pfin = vars["pfin"]
    PUR_CS = parent(params["PUR_CS"])
    col_list = [1, 2, 3, 5]
    u_list = [na+1, na+2, na+3, na+4]
    for d in 1:nr, (fi, (col, u)) in enumerate(zip(col_list, u_list))
        @constraint(m, delPGDPEXP[d,col] == PUR_CS[u,d] * (pfin[fi,d] - 1))
    end
end

function E_delPGDPEXPb!(m, vars, na, nr, params)
    delPGDPEXP = vars["delPGDPEXP"]; ptot = vars["ptot"]
    for d in 1:nr
        @constraint(m, delPGDPEXP[d,4] == sum(get(STOCKS_idx, (i,d), 0.0) * (ptot[i,d] - 1) for i in 1:na))
    end
end

function E_delPGDPEXPc!(m, vars, na, nr, ns, params)
    delPGDPEXP = vars["delPGDPEXP"]; pimp = vars["pimp"]
    TRADE_D = parent(params["TRADE_D"])
    for q in 1:nr
        @constraint(m, delPGDPEXP[q,6] == -sum(TRADE_D[c,2,q] * (pimp[c,q] - 1) for c in 1:na))
    end
end

function E_delPGDPEXPd!(m, vars, na, nr, nm, params)
    delPGDPEXP = vars["delPGDPEXP"]; pdom = vars["pdom"]
    SUPPMAR_RD = parent(params["SUPPMAR_RD"]); SUPPMAR_R = parent(params["SUPPMAR_R"])
    for q in 1:nr
        rhs = sum(SUPPMAR_RD[m_,q] * (pdom[m_,q] - 1) for m_ in 1:nm) -
              sum(SUPPMAR_R[m_,q,p] * (pdom[m_,p] - 1) for m_ in 1:nm, p in 1:nr)
        @constraint(m, delPGDPEXP[q,9] == rhs)
    end
end

function E_delPGDPEXPe!(m, vars, na, nr, ns, params)
    delPGDPEXP = vars["delPGDPEXP"]; pbasic = vars["pbasic"]
    TRADE_D = parent(params["TRADE_D"]); TRDIAG = parent(params["TRDIAG"])
    for q in 1:nr
        rhs = sum((TRADE_D[c,s,q] - TRDIAG[c,s,q]) * (pbasic[c,s,q] - 1) for c in 1:na, s in 1:ns)
        @constraint(m, delPGDPEXP[q,7] == rhs)
    end
end

function E_delPGDPEXPf!(m, vars, na, nr, ns, params)
    delPGDPEXP = vars["delPGDPEXP"]; pbasic = vars["pbasic"]
    TRDIAG = parent(params["TRDIAG"]); TRADE = parent(params["TRADE"])
    for q in 1:nr
        rhs = sum(sum(TRADE[c,s,r,q] * (pbasic[c,s,r] - 1) for r in 1:nr) - TRDIAG[c,s,q] * (pbasic[c,s,q] - 1)
                  for c in 1:na, s in 1:ns)
        @constraint(m, delPGDPEXP[q,8] == -rhs)
    end
end

function E_pgdpexp!(m, vars, nr, params)
    pgdpexp = vars["pgdpexp"]; delPGDPEXP = vars["delPGDPEXP"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for d in 1:nr
        gdp_exp = sum(GDPEXPSUM[d,k] for k in 1:9)
        gdp_exp > 1e-10 || continue
        @constraint(m, gdp_exp * pgdpexp[d] == gdp_exp + sum(delPGDPEXP[d,k] for k in 1:9))
    end
end

function E_wgdpexp!(m, vars, nr)
    wgdpexp = vars["wgdpexp"]; xgdpexp = vars["xgdpexp"]; pgdpexp = vars["pgdpexp"]
    @constraint(m, [d=1:nr], wgdpexp[d] == xgdpexp[d] * pgdpexp[d])
end

function E_wgdpdiff!(m, vars, nr)
    wgdpdiff = vars["wgdpdiff"]; wgdpinc = vars["wgdpinc"]; wgdpexp = vars["wgdpexp"]
    @constraint(m, [d=1:nr], wgdpdiff[d] == wgdpinc[d] - wgdpexp[d])
end

# GNECAT = {HOU,INV,GOV,STOCKS} = columns {1,2,3,4} (excludes EXP, unlike GDPEXPCAT)
function E_xgne!(m, vars, nr, params)
    xgne = vars["xgne"]; delXGDPEXP = vars["delXGDPEXP"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for d in 1:nr
        gne_gdp = sum(GDPEXPSUM[d,k] for k in [1,2,3,4])
        gne_gdp > 1e-10 || continue
        @constraint(m, gne_gdp * xgne[d] == gne_gdp + sum(delXGDPEXP[d,k] for k in [1,2,3,4]))
    end
end

function E_pgne!(m, vars, nr, params)
    pgne = vars["pgne"]; delPGDPEXP = vars["delPGDPEXP"]
    GDPEXPSUM = parent(params["GDPEXPSUM"])
    for d in 1:nr
        gne_gdp = sum(GDPEXPSUM[d,k] for k in [1,2,3,4])
        gne_gdp > 1e-10 || continue
        @constraint(m, gne_gdp * pgne[d] == gne_gdp + sum(delPGDPEXP[d,k] for k in [1,2,3,4]))
    end
end

# National (nr+1) slot: xgne/pgne/wgne collapse to a GDPEXPSUM-weighted
# average across regions (Category 2 — xgne(d)/pgne(d) are bmk=1 ratios).
function E_xgneB!(m, vars, nr, params)
    xgne = vars["xgne"]; GDPEXPSUM = parent(params["GDPEXPSUM"])
    GNE_d = [sum(GDPEXPSUM[d,k] for k in [1,2,3,4]) for d in 1:nr]
    tot = sum(GNE_d)
    tot > 1e-10 || return
    @constraint(m, tot * xgne[nr+1] == sum(GNE_d[d] * xgne[d] for d in 1:nr))
end

function E_pgneB!(m, vars, nr, params)
    pgne = vars["pgne"]; GDPEXPSUM = parent(params["GDPEXPSUM"])
    GNE_d = [sum(GDPEXPSUM[d,k] for k in [1,2,3,4]) for d in 1:nr]
    tot = sum(GNE_d)
    tot > 1e-10 || return
    @constraint(m, tot * pgne[nr+1] == sum(GNE_d[d] * pgne[d] for d in 1:nr))
end

function E_wgne!(m, vars, nr)
    wgne = vars["wgne"]; xgne = vars["xgne"]; pgne = vars["pgne"]
    @constraint(m, [d=1:nr+1], wgne[d] == xgne[d] * pgne[d])
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
# realwage_id/realwage_i/realwage are bmk=1 ratios; xlab_id/xlab_i are bmk=1
# quantity indices. Shifters (flabsup_id/flabsupA/labslack/flab/flab_i/
# flab_io/flab_id/flab_iod) carry NO 0.01 scaling — this codebase redefines
# them in raw log-points (see the Excerpt 15/16 gret/ggro/finv2 precedent),
# not GEMPACK's native percent-points.
function E_labslack!(m, vars, no)
    realwage_id = vars["realwage_id"]; xlab_id = vars["xlab_id"]; flabsup_id = vars["flabsup_id"]
    @constraint(m, [o=1:no],
        log(realwage_id[o] + _TINY) == 2.0 * log(xlab_id[o] + _TINY) + flabsup_id[o])
end

function E_flab_i!(m, vars, na, nr, no)
    realwage_i = vars["realwage_i"]; xlab_i = vars["xlab_i"]
    flabsupA = vars["flabsupA"]; labslack = vars["labslack"]
    @constraint(m, [o=1:no, d=1:nr],
        log(realwage_i[o,d] + _TINY) == (1.0/1.1) * log(xlab_i[o,d] + _TINY) + flabsupA[o,d] + labslack[o])
end

# realwage is pinned by E_plab! (Excerpt 27, from the factor-demand block).
# This equation instead solves for the otherwise-undetermined wage-shifter
# residuals flab/flab_i/flab_io/flab_id/flab_iod given that already-fixed
# realwage — the standard ORANI/TERM decomposition-residual pattern (see
# TERM.TAB's "Exogenous flab;" declaration and closure comments).
function E_realwage!(m, vars, na, nr, no)
    realwage = vars["realwage"]; flab = vars["flab"]
    flab_i = vars["flab_i"]; flab_io = vars["flab_io"]
    flab_id = vars["flab_id"]; flab_iod = vars["flab_iod"]
    @constraint(m, [i=1:na, o=1:no, d=1:nr],
        log(realwage[i,o,d] + _TINY) == flab[i,o,d] + flab_i[o,d] + flab_io[d] + flab_id[o] + flab_iod)
end

# ═══════════════════════════════════════════════════════════════════════════
# EXCERPT 39 — Household closure
# ═══════════════════════════════════════════════════════════════════════════
# whouhtot is a genuine VALUE (E_whouhtot!, Excerpt "wlnd_i"-style block),
# so it must be normalized by its own benchmark (HOUPUR_C_v) before
# comparing against the bmk=1 ratios wlab_io/wgdpexp.
function E_fhou!(m, vars, na, nr, ns, nu, params)
    whouhtot = vars["whouhtot"]; wlab_io = vars["wlab_io"]
    fhou = vars["fhou"]; houslack = vars["houslack"]
    HOUPUR_C_v = parent(params["HOUPUR_C"])
    @constraint(m, [d=1:nr],
        log(whouhtot[d] + _TINY) - log(HOUPUR_C_v[1,d] + _TINY) == log(wlab_io[d] + _TINY) + fhou[d] + houslack)
end

function E_fhou2!(m, vars, na, nr, ns, nu, params)
    whouhtot = vars["whouhtot"]; wgdpexp = vars["wgdpexp"]
    fhou2 = vars["fhou2"]; houslack = vars["houslack"]
    HOUPUR_C_v = parent(params["HOUPUR_C"])
    @constraint(m, [d=1:nr],
        log(whouhtot[d] + _TINY) - log(HOUPUR_C_v[1,d] + _TINY) == log(wgdpexp[d] + _TINY) + fhou2[d] + houslack)
end

# TERM.TAB:2052-2053
#   Variable natfhou # National ratio, nominal household consumption to GDP #;
#   Equation E_natfhou  NatMacro("NomHou") = natfhou + NatMacro("NomGDPexp");
#
# This was an empty stub until 2026-07-31, which left `natfhou` a free variable
# appearing in no constraint. `solve_newton!` pinned it as an "orphan" and the
# system squared up, so every closure that leaves natfhou endogenous ran
# correctly and the gap was invisible — its docstring even named natfhou as the
# one known orphan. But `origin/coalprice.CMF:83 swap houslack = natfhou` makes
# natfhou EXOGENOUS: the orphan disappears, nothing gets pinned, and the model
# comes out one variable short of square. A missing equation and a pinned orphan
# are indistinguishable until a closure swaps the orphan.
#
# natfhou is a bmk-0 log-additive shifter (it is in ZERO_START_NAMES), and both
# NatMacro entries are bmk-1 ratios, so TABLO's %-change sum becomes a log sum.
# At the benchmark this reads 0 == 0 + 0, so replication is untouched. It is also
# homogeneous of degree zero — NomHou and NomGDPexp scale together with the
# numeraire, leaving natfhou fixed — which matters because `:exrate` runs leave
# every nominal aggregate free.
#
# Must be called AFTER `build_macros!`, which is what declares NatMacro.
function E_natfhou!(m, vars, na, nr, ns, nu, params)
    haskey(vars, "NatMacro") || error(
        "E_natfhou! needs NatMacro, which build_macros! declares — call it after " *
        "the macro block, not with the rest of Excerpt 39.")
    NatMacro = vars["NatMacro"]; natfhou = vars["natfhou"]
    k_hou = findfirst(==("NomHou"), MAINMACROS)
    k_gdp = findfirst(==("NomGDPexp"), MAINMACROS)
    @constraint(m, log(NatMacro[k_hou]) == natfhou + log(NatMacro[k_gdp]))
end
