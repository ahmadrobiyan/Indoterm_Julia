"""
Step 6 — Dynamic extension (TERM.TAB Excerpts 50-54, "following ORANIGRD.DOC").

Four blocks:
  * **Excerpt 50** capital accumulation — this year's capital grows by
    [last year's investment − depreciation].
  * **Excerpt 51** investment rule — investment/capital ratio responds to the gap
    between the expected and the normal rate of return.
  * **Excerpt 53** national dynamic reporting aggregates (`Backsolve` in TABLO).
  * **Excerpt 54** real-wage adjustment — an upward-sloping labour supply
    schedule that shifts out while actual employment exceeds trend.

**Levels convention.** TABLO writes these in %-change form, and the pattern
`0.01*COEFF*x̂` (COEFF × a percentage change) is a *level* change. This port
already represents such quantities as additive **log-gaps** benchmarked at 0 —
see `E_gret!` (`gret = log(pcap) − log(pinvitot)`) and `E_xinvitot!`
(`ggro = log(xinvitot) − log(xcap/CAP)`) in Excerpt 15. So `0.01*COEFF*x̂`
translates to `COEFF * (log-gap)` and every equation here is benchmark-consistent
with all gaps at 0. TABLO's `(change)` variables (`delgret`, `delgretexp`,
`delempratio`, `delwagerate`, `delfwage`, `delUnity`, `faccum`) are likewise
additive, benchmark 0.

**The block is a passive satellite under the base static closure**, exactly as
TERM.TAB intends ("to switch off, put faccum endogenous; delUnity is exogenous
but unshocked"). Given the core solve, it recursively computes
`gret → delgret → delgretexp → gretxp → mratio → gro → finv4`, with `faccum`,
`finv4` and `delfwage` absorbing. Variable/equation counts balance exactly
(8 equations against 8 endogenous + 2 exogenous per industry-region), so adding
it keeps the system square. `TERM.CMF`'s dynamic closure switches it on by
swapping `xcap`↔`faccum`, `finv1`↔`finv4` and `delfwage`↔`flabsup_id`.

Inter-period `Update` statements (CAPSTOK, RNORMAL, GROTREND, INVGRO, EMPRAT,
WAGERATE) belong to a recursive multi-period driver, not to a single solve;
`update_dynamics!` below applies them so a year-on-year loop can be built on top.
"""

"""
    build_dynamics!(m, vars, na, nr, no, params)

Declare the Excerpt 50-54 variables and equations on `m`, adding them to `vars`.
No-op unless `prepare_parameters!` derived the dynamic coefficients (i.e. the
`STOC`/`DPRC`/... headers survived aggregation).
"""
function build_dynamics!(m::JuMP.Model, vars::Dict{String,Any}, na, nr, no, params)
    haskey(params, "CAPSTOK") || return m, vars

    CAPSTOK      = parent(params["CAPSTOK"])
    CAPSTOK_OLDP = parent(params["CAPSTOK_OLDP"])
    CAPADD       = parent(params["CAPADD"])
    GROSSRET     = parent(params["GROSSRET"]);  GROSSRET0 = parent(params["GROSSRET0"])
    GRETEXP      = parent(params["GRETEXP"]);   GRETEXP0  = parent(params["GRETEXP0"])
    MCOEFF       = parent(params["MCOEFF"])
    RORADJ       = parent(params["RORADJ"])
    CAPv         = parent(params["CAP"])
    INVEST_C     = parent(params["INVEST_C"])
    CAPSTOK_D    = params["CAPSTOK_D"];  INVEST_CD = params["INVEST_CD"]

    xcap = vars["xcap"]; gret = vars["gret"]; ggro = vars["ggro"]
    xinvitot = vars["xinvitot"]; pinvitot = vars["pinvitot"]
    invslack = vars["invslack"]

    # ── variables ───────────────────────────────────────────────────────────
    # All additive log-gap / (change) variables: unbounded, benchmark 0.
    @variable(m, faccum[1:na, 1:nr])       # accumulation switch/shifter
    @variable(m, finv4[1:na, 1:nr])        # dynamic investment-rule shifter
    @variable(m, gro[1:na, 1:nr])          # planned investment/capital ratio
    @variable(m, rnorm[1:na, 1:nr])        # normal gross rate of return
    @variable(m, mratio[1:na, 1:nr])       # expected/normal rate-of-return ratio
    @variable(m, gretxp[1:na, 1:nr])       # expected rate of return
    @variable(m, delgret[1:na, 1:nr])      # (change) gross rate of return
    @variable(m, delgretexp[1:na, 1:nr])   # (change) expected rate of return
    @variable(m, frnorm[1:na, 1:nr])       # exogenous shifter (TERM.CMF)
    @variable(m, gtrend[1:na, 1:nr])       # exogenous trend I/K ratio (TERM.CMF)
    @variable(m, delUnity)                 # exogenous dummy, shocked to 1 to switch on
    @variable(m, frnorm_id)                # exogenous scalar shifter
    @variable(m, delfwage_o)               # exogenous scalar wage shifter
    @variable(m, emptrend[1:no])           # exogenous trend employment
    @variable(m, delempratio[1:no])        # (change) actual/trend employment
    @variable(m, delwagerate[1:no])        # (change) real wage index
    @variable(m, delfwage[1:no])           # real-wage-mechanism shifter
    # National reporting: quantity/price indices (benchmark 1) and a log-gap.
    @variable(m, natcapstok[1:na] >= 1e-6)
    @variable(m, natxinvitot[1:na] >= 1e-6)
    @variable(m, natpinvitot[1:na] >= 1e-6)
    @variable(m, natggro[1:na])

    for (nm, v) in (("faccum", faccum), ("finv4", finv4), ("gro", gro),
                    ("rnorm", rnorm), ("mratio", mratio), ("gretxp", gretxp),
                    ("delgret", delgret), ("delgretexp", delgretexp),
                    ("frnorm", frnorm), ("gtrend", gtrend),
                    ("delUnity", delUnity), ("frnorm_id", frnorm_id),
                    ("delfwage_o", delfwage_o), ("emptrend", emptrend),
                    ("delempratio", delempratio), ("delwagerate", delwagerate),
                    ("delfwage", delfwage), ("natcapstok", natcapstok),
                    ("natxinvitot", natxinvitot), ("natpinvitot", natpinvitot),
                    ("natggro", natggro))
        vars[nm] = v
    end

    # ── Excerpt 50: capital accumulation ────────────────────────────────────
    # TABLO: 0.01*CAPSTOK_OLDP*xcap = CAPADD*delUnity + faccum, where xcap is the
    # %-change in the capital stock. Here the capital quantity index is
    # xcap/CAP, so its log-gap is log(xcap) − log(CAP) (0 at the benchmark).
    for i in 1:na, d in 1:nr
        @constraint(m, CAPSTOK_OLDP[i,d] *
            (log(xcap[i,d] + _TINY) - log(CAPv[i,d] + _TINY)) ==
            CAPADD[i,d] * delUnity + faccum[i,d])
    end

    # ── Excerpt 51: investment rule ─────────────────────────────────────────
    @constraint(m, [i=1:na, d=1:nr], rnorm[i,d] == frnorm[i,d] + frnorm_id)
    @constraint(m, [i=1:na, d=1:nr], delgret[i,d] == GROSSRET[i,d] * gret[i,d])
    @constraint(m, [i=1:na, d=1:nr], ggro[i,d] == gro[i,d] + finv4[i,d])
    @constraint(m, [i=1:na, d=1:nr], mratio[i,d] == gretxp[i,d] - rnorm[i,d])
    @constraint(m, [i=1:na, d=1:nr], delgretexp[i,d] ==
        RORADJ[i,d] * ((GROSSRET0[i,d] - GRETEXP0[i,d]) * delUnity + delgret[i,d]))
    @constraint(m, [i=1:na, d=1:nr],
        gro[i,d] == gtrend[i,d] + invslack + MCOEFF[i,d] * mratio[i,d])
    # delgretexp = GRETEXP*gretxp. Where GRETEXP is 0 the equation carries no
    # information about gretxp (zero Jacobian column) — pin it instead, exactly
    # as the zero-flow guards elsewhere in build_equations.jl.
    for i in 1:na, d in 1:nr
        if GRETEXP[i,d] > 1e-10
            @constraint(m, delgretexp[i,d] == GRETEXP[i,d] * gretxp[i,d])
        else
            fix(gretxp[i,d], 0.0; force=true)
        end
    end

    # ── Excerpt 53: national dynamic reporting ──────────────────────────────
    for i in 1:na
        if CAPSTOK_D[i] > 1e-10
            @constraint(m, CAPSTOK_D[i] * natcapstok[i] ==
                sum(CAPSTOK[i,d] * (CAPv[i,d] > 1e-10 ? xcap[i,d] / CAPv[i,d] : 1.0)
                    for d in 1:nr))
        else
            fix(natcapstok[i], 1.0; force=true)
        end
        if INVEST_CD[i] > 1e-10
            @constraint(m, INVEST_CD[i] * natxinvitot[i] ==
                sum(INVEST_C[i,d] * xinvitot[i,d] for d in 1:nr))
            @constraint(m, INVEST_CD[i] * natpinvitot[i] ==
                sum(INVEST_C[i,d] * pinvitot[i,d] for d in 1:nr))
        else
            fix(natxinvitot[i], 1.0; force=true)
            fix(natpinvitot[i], 1.0; force=true)
        end
        @constraint(m, natggro[i] ==
            log(natxinvitot[i] + _TINY) - log(natcapstok[i] + _TINY))
    end

    # ── Excerpt 54: real-wage adjustment ────────────────────────────────────
    if haskey(params, "EMPRAT")
        EMPRAT = params["EMPRAT"]; EMPRAT0 = params["EMPRAT0"]
        ELASTWAGE = params["ELASTWAGE"]; WAGERATE = params["WAGERATE"]
        xlab_id = vars["xlab_id"]; realwage_id = vars["realwage_id"]
        for o in 1:no
            @constraint(m, delempratio[o] ==
                EMPRAT[o] * (log(xlab_id[o] + _TINY) - emptrend[o]))
            @constraint(m, delwagerate[o] == WAGERATE[o] * log(realwage_id[o] + _TINY))
            # LABMOMENTUM = ELASTWAGE*(EMPRAT0-1) enters via delUnity, so the
            # mechanism keeps pushing the wage up while employment exceeds trend.
            @constraint(m, delwagerate[o] == delfwage[o] + delfwage_o +
                ELASTWAGE[o] * ((EMPRAT0[o] - 1.0) * delUnity + delempratio[o]))
        end
    end

    return m, vars
end

"""
    update_dynamics!(params, sol) -> params

Apply TERM.TAB's inter-period `Update` statements so a recursive multi-period
run can advance the database between solves. `sol` is a `NewtonResult.values`
dict. Mutates and returns `params`.

Updates (Excerpts 50-51, 54): `CAPSTOK = xcap*pinvitot`, `CAPSTOK_OLDP = xcap`,
`RNORMAL = rnorm`, `GROTREND = gtrend`, `INVGRO = xinvitot`,
`EMPRAT += delempratio`, `WAGERATE += delwagerate`. Derived coefficients
(`GROSSRET`, `GROSSGRO`, `GROMAX`, `GRETEXP`, `MCOEFF`, `CAPADD`) are recomputed
from the updated stocks; the `*0` initial-value coefficients are deliberately
NOT updated, matching TABLO's `(parameter)`/`(initial)` declarations.
"""
function update_dynamics!(params::Dict{String,Any}, sol::Dict{String,Any})
    haskey(params, "CAPSTOK") || return params
    na, nr = size(parent(params["CAPSTOK"]))
    CAPSTOK = parent(params["CAPSTOK"])
    CAPv    = parent(params["CAP"])
    xcap    = sol["xcap"]; pinvitot = sol["pinvitot"]

    for i in 1:na, d in 1:nr
        # CAPSTOK measured in current asset prices; CAPSTOK_OLDP in last year's.
        capidx = CAPv[i,d] > 1e-10 ? xcap[i,d] / CAPv[i,d] : 1.0
        params["CAPSTOK_OLDP"][i,d] = CAPSTOK[i,d] * capidx
        CAPSTOK[i,d] = CAPSTOK[i,d] * capidx * pinvitot[i,d]
    end
    haskey(sol, "rnorm")  && (params["RNORMAL"]  .*= exp.(sol["rnorm"]))
    haskey(sol, "gtrend") && (params["GROTREND"] .*= exp.(sol["gtrend"]))
    if haskey(params, "EMPRAT") && haskey(sol, "delempratio")
        params["EMPRAT"]   .+= sol["delempratio"]
        params["WAGERATE"] .+= sol["delwagerate"]
    end

    # Recompute the derived coefficients from the updated stocks.
    INVEST_C = parent(params["INVEST_C"]); DPRC = parent(params["DPRC"])
    QRATIO = parent(params["QRATIO"]); GROTREND = parent(params["GROTREND"])
    ALPHA_D = parent(params["ALPHA_DYN"]); RNORMAL = parent(params["RNORMAL"])
    for i in 1:na, d in 1:nr
        cs = CAPSTOK[i,d]
        params["GROSSRET"][i,d] = cs > 1e-10 ? CAPv[i,d] / cs : 0.0
        params["GROSSGRO"][i,d] = cs > 1e-10 ? INVEST_C[i,d] / cs : 0.0
        params["GROMAX"][i,d]   = QRATIO[i,d] * GROTREND[i,d]
        gg = params["GROSSGRO"][i,d]; gmx = params["GROMAX"][i,d]
        den = gg > 1e-10 ? (gmx / gg) - 1.0 : 0.001
        den <= 0 && (den = 0.001)
        params["GRETEXP"][i,d] = ALPHA_D[i,d] > 1e-10 ?
            RNORMAL[i,d] * ((QRATIO[i,d] - 1.0) / den)^(1.0 / ALPHA_D[i,d]) : 0.0
        params["MCOEFF"][i,d] = gmx > 1e-10 ? ALPHA_D[i,d] * (1.0 - gg / gmx) : 0.0
        params["CAPADD"][i,d] = INVEST_C[i,d] - DPRC[i,d] * cs
    end
    return params
end
