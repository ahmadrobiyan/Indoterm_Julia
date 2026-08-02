"""
    benchmark_levels(params::Dict) -> Dict{String,Any}

Benchmark ("no shock") LEVEL of every genuine quantity/value variable, keyed by
its JuMP variable name, for use as start values (and closure fix values) in
`initialize_model!`.

Why this exists: the model is calibrated in the GEMPACK "P = 1, Q = value-flow"
convention (see `prepare_parameters!`) — every CES/CET share is fitted so that
(all prices = 1, all quantities = their benchmark value-flow) is an exact
solution. `initialize_model!`'s generic default seeds every quantity variable at
`1.0`, which is only correct for *index* quantities (composites normalized to 1
at calibration, e.g. `xtot` — `ces_calibrate(MAKE[:,i,d], -SCET[i], 1.0)`). For
genuine *flow* quantities (`xmake`, `xint`, `xlab`, ...) the benchmark level is
the value flow, not 1, and starting them at 1 makes the value-level identities
and `delTAX` tax equations violate by roughly the flow magnitude (~1e10),
producing the massively-infeasible start point that stalled the Ipopt solve.

This map contains ONLY the flow variables — every variable NOT listed here keeps
`initialize_model!`'s `1.0` (ratio/index) or `0.0` (shifter) default, which is
correct for them. Each entry maps to an array whose shape matches the variable's
JuMP declaration exactly, so `initialize_model!` can index it element-wise.

The benchmark quantity of a flow variable equals its benchmark value flow (since
P = 1), so a value variable `w…` and its quantity twin `x…` share the same
benchmark array here.

Which variable is a flow vs. an index is taken from the calibration `total`
argument in `prepare_parameters!` (composite normalized to 1 ⇒ index ⇒ NOT
listed) and from the `delTAX`/value-identity equations in `build_equations.jl`
(quantity pinned to a flow ⇒ listed). Anything whose benchmark is still
uncertain is deliberately omitted and left to surface in the residual check
(`test/check_residual_6reg.jl`) rather than guessed.
"""
function benchmark_levels(params::Dict)
    p = params
    bmk = Dict{String,Any}()

    # dims
    PUR = p["PUR"]                    # na × ns × nu × nr   (purchaser value)
    na, ns, nu, nr = size(PUR)
    u_hou = na + 1; u_inv = na + 2; u_gov = na + 3; u_exp = na + 4

    # ── Demand-side source-specific quantities: benchmark = purchaser value PUR.
    #    Pinned by the delTAX tax equations (resid = TAX·(PUR − x); zero ⇔ x=PUR).
    bmk["xint"] = PUR[:, :, 1:na, :]           # na × ns × na × nr
    bmk["xhou"] = PUR[:, :, u_hou, :]          # na × ns × nr
    bmk["xinv"] = PUR[:, :, u_inv, :]
    bmk["xgov"] = PUR[:, :, u_gov, :]
    bmk["xexp"] = PUR[:, :, u_exp, :]

    # ── Source-composite (dom+imp) demand quantities: benchmark = PUR_S.
    #    Composite total in the Armington nest ces_calibrate(PUR[c,:,u,d], …, PUR_S).
    PUR_S = p["PUR_S"]                          # na × nu × nr
    bmk["xint_s"] = PUR_S[:, 1:na, :]           # na × na × nr
    bmk["xhou_s"] = PUR_S[:, u_hou, :]          # na × nr
    bmk["xinv_s"] = PUR_S[:, u_inv, :]
    bmk["xgov_s"] = PUR_S[:, u_gov, :]
    bmk["xexp_s"] = PUR_S[:, u_exp, :]

    # ── Factor demands & primary nest: benchmark = the calibration flows.
    #    ces_calibrate(V1LAB[i,:,d], …, LAB_O)  and
    #    ces_calibrate([LAB_O, V1CAP, V1LND], …, PRIM).
    bmk["xlab"]   = p["LAB"]                    # na × no × nr  (V1LAB)
    bmk["xlab_o"] = p["LAB_O"]                  # na × nr
    bmk["xcap"]   = p["CAP"]                    # na × nr   (also a closure var)
    bmk["xlnd"]   = p["LND"]                    # na × nr   (also a closure var)
    bmk["xprim"]  = p["PRIM"]                   # na × nr

    # ── Make/CET: components xmake benchmark = MAKE (composite xtot is an INDEX=1).
    bmk["xmake"] = p["MAKE"]                    # na × na × nr
    bmk["xcom"]  = p["MAKE_I"]                  # na × nr   (E_xcomA_B: xcom = Σ_i xmake = MAKE_I)

    # ── Regional sourcing / delivery (Excerpts 19-23): all genuine basic-value
    #    flows (E_xtrad/E_xtradmar/E_xsuppmar are log-form, exact at the data).
    bmk["xtrad"]   = p["TRADE"]                 # na × ns × nr × nr
    bmk["xtrad_d"] = p["TRADE_D"]              # na × ns × nr   (Σ_d xtrad)
    bmk["xtrad_r"] = p["TRADE_R"]              # na × ns × nr   (Σ_r xtrad)
    bmk["xtradmar"]    = p["TMAR"]            # na × ns × nm × nr × nr
    bmk["xsuppmar"]    = p["MARS"]           # nm × nr × nr × nr
    bmk["xsuppmar_p"]  = p["SUPPMAR_P"]      # nm × nr × nr   (Σ_{c,s} xtradmar)
    bmk["xsuppmar_d"]  = p["SUPPMAR_D"]      # nm × nr × nr   (Σ_d xsuppmar)
    bmk["xsuppmar_rd"] = p["SUPPMAR_RD"]     # nm × nr        (Σ_r xsuppmar_d)

    # ── xuse: total regional demand, seeded at the basic-value total TRADE_R so
    #    the dominant E_xtrad sourcing nest (na·ns·nr² eqs) is exact. The old
    #    valuation clash (E_xuse literal-summed purchaser-basis users, leaving a
    #    margin+tax wedge here) is GONE: E_xuse now uses the faithful basic-
    #    weighted form (USE_U/TRADE_R)·xuse = Σ_user (USE_user/PUR_user)·xuser,
    #    which is exact at this TRADE_R seed by the identity USE_U = Σ_u USE.
    bmk["xuse"]   = p["TRADE_R"]              # na × ns × nr
    #    xint_i = Σ_i xint (purchaser); makes E_xint_i exact.
    xint_i = dropdims(sum(PUR[:, :, 1:na, :], dims=3), dims=3)   # na × ns × nr
    bmk["xint_i"] = xint_i

    # ── Investment commodity demands: benchmark = INVEST (=V2PUR); Leontief in
    #    the industry index xinvitot (=1). E_xinvi: xinvi = INVEST·xinvitot.
    bmk["xinvi"] = p["INVEST"]                 # na × na × nr

    # ── Household Stone-Geary split: xlux/xsub benchmark = XLUX0/XSUB0
    #    (na × 1 × nr storage → slice the singleton middle axis). Their sum is
    #    xhou_s = PUR_S[hou] already seeded above.
    bmk["xlux"] = p["XLUX0"][:, 1, :]          # na × nr
    bmk["xsub"] = p["XSUB0"][:, 1, :]          # na × nr
    # wlux is a genuine nominal value (aggregate supernumerary expenditure),
    # bmk = WLUX0[d]; alux is the "shifter" exception whose true bmk is the
    # non-unity calibrated share ALUX0[c,d]=XLUX0/WLUX0 (E_xlux: xlux·phou =
    # wlux·alux ⇒ needs BOTH to reproduce XLUX0). See E_xlux!/E_alux! headers.
    bmk["wlux"] = p["WLUX0"][1, :]             # nr
    bmk["alux"] = p["ALUX0"][:, 1, :]          # na × nr

    # ── Exports: xexpd benchmark = XEXPD0 (E_xexpd log-form exact at data;
    #    E_xexp then pins xexp[·,1,·]=xexpd, xexp[·,2,·]=0, matched by PUR[exp]).
    bmk["xexpd"] = p["XEXPD0"]                 # na × nr

    # ── Inventories: xstocks is an ADDITIVE log-volume index (E_xstocks:
    #    xstocks = log(xtot)+fxstocks), so its benchmark is 0, not 1.
    bmk["xstocks"] = zeros(na, nr)             # na × nr

    # ── Total nominal household expenditure: benchmark = HOUPUR_C (E_whouhtot:
    #    whouhtot = HOUPUR_C·phouhtot·xhoutot, with phouhtot=xhoutot=1).
    bmk["whouhtot"] = p["HOUPUR_C"][1, :]      # nr

    # ── Final-demand aggregates (Excerpt 24) are all bmk=1 quantity INDICES
    #    (TERM.TAB:1209 "Final user quantity indices"), NOT flows: E_xfina/b/c/d
    #    divide the component sum by the category total PUR_CS(u,d), so xfin=1 at
    #    benchmark for every user — the generic 1.0 default is correct. Likewise
    #    wfin/natxfin/natwfin are bmk=1. (An earlier port dropped the PUR_CS
    #    divisor and seeded xfin/natxfin as flows here; both are removed now that
    #    E_xfinb!/E_natxfin! carry the divisor — see their headers.)

    return bmk
end
