"""
Excerpts 31-32 — regional macro reporting variables and their aggregation to
national macros (`MainMacro` / `NatMacro`).

This is the standard TERM reporting layer, and it is more than cosmetic: it
defines `NatMacro("GDPPI")`, which is what `TERM.CMF`'s active numeraire swap
(`swap phi = Natmacro("GDPPI");`) points at. Without it the port had to leave
`phi` as the numeraire — see `initialize_model!.jl`.

**Levels conventions.** Every `MainMacro`/`NatMacro` entry equals some existing
benchmark=1 index variable (`xfin`, `xgdpexp`, `pfin`, …), so both are bmk=1
ratios and the defining equations are plain equalities. Where TABLO combines
%-changes additively — the contribution decomposition and the balance-of-trade
shares — this port uses **log-gaps** of those bmk=1 indices, matching the
`gret`/`ggro` convention used elsewhere, so each such equation is 0 = 0 at the
benchmark. `contMainMacro`, `shrBoT`, `shrBoTnom`, `shrBoTnom2` and `finvgdp` are
therefore additive variables benchmarked at 0.

**Squareness.** Every declared variable has exactly one defining equation:
`MainMacro[:, 1:nr]` from E_MainMacroA-U, `MainMacro[:, National]` from
E_MainMacroV, `NatMacro` from E_NatMacro, and so on — so the block keeps the
system square, like `build_dynamics!.jl`.
"""

# Order is load-bearing — it is TERM.TAB's own MAINMACROS set order
# (TERM.TAB:1656-1659). Index by these names, never by position guessed from an
# equation's letter suffix (see the E_delGDPINC slot bug in calculate_gdp.jl).
const MAINMACROS = [
    "RealHou", "RealInv", "RealGov", "ExpVol", "ImpVolUsed", "ImpsLanded",
    "RealGDP", "RealGNE", "AggEmploy", "realwage_io", "plab_io", "AggCapStock",
    "GDPPI", "GNEPI", "CPI", "ExportPI", "ImpsLandedPI", "Population",
    "NomHou", "NomGDPexp", "NomGDPinc",
]

const SELMACROS = [
    "RealHou", "RealInv", "RealGov", "ExpVol", "ImpVolUsed", "RealGNE",
    "RealGDP", "AggEmploy", "realwage_io", "CPI", "AggCapStock",
]

_mac(name) = findfirst(==(name), MAINMACROS)

"""
    build_macros!(m, vars, na, nr, ns, params)

Declare and constrain Excerpt 31-32's reporting variables on `m`. `MainMacro` is
indexed `[macro, region]` with `region = nr+1` meaning "National" (TABLO's
`REGPLUS`).
"""
function build_macros!(m::JuMP.Model, vars::Dict{String,Any}, na, nr, ns, params)
    TRADE   = parent(params["TRADE"])
    TRADE_R = parent(params["TRADE_R"]); TRADE_D = parent(params["TRADE_D"])
    IMPUSED_C   = params["IMPUSED_C"]; IMPLANDED_C = params["IMPLANDED_C"]
    PUR_CS  = parent(params["PUR_CS"]); GDPEXP = params["GDPEXP"]
    GDPINC  = params["GDPINC"];         GDPEXPSUM = parent(params["GDPEXPSUM"])
    LAB_IO  = parent(params["LAB_IO"]); CAP_I = parent(params["CAP_I"])
    POP     = parent(params["PO01"])
    nmac = length(MAINMACROS); nsel = length(SELMACROS)
    u_hou, u_inv, u_gov, u_exp = na+1, na+2, na+3, na+4

    xtrad_r = vars["xtrad_r"]; xtrad_d = vars["xtrad_d"]; pimp = vars["pimp"]
    xfin = vars["xfin"]; pfin = vars["pfin"]; wfin = vars["wfin"]
    xgdpexp = vars["xgdpexp"]; pgdpexp = vars["pgdpexp"]
    wgdpexp = vars["wgdpexp"]; wgdpinc = vars["wgdpinc"]
    xgne = vars["xgne"]; pgne = vars["pgne"]
    xlab_io = vars["xlab_io"]; realwage_io = vars["realwage_io"]
    plab_io = vars["plab_io"]; xcap_i = vars["xcap_i"]; nhou = vars["nhou"]

    # ── Excerpt 31: import volume/price indices ─────────────────────────────
    @variable(m, ximps[1:na, 1:nr] >= 1e-6)
    @variable(m, ximpused[1:nr] >= 1e-6)
    @variable(m, pimpused[1:nr] >= 1e-6)
    @variable(m, ximplanded[1:nr] >= 1e-6)
    @variable(m, pimplanded[1:nr] >= 1e-6)

    # ximps is a bmk=1 INDEX while xtrad_r is a value FLOW (bmk TRADE_R), so the
    # translation of TABLO's `ximps = xtrad_r` is the ratio, not the identity.
    # Zero-flow cells carry no information about the index (zero Jacobian
    # column) — pin them, as elsewhere in this port.
    for c in 1:na, d in 1:nr
        if TRADE_R[c,2,d] > 1e-10
            @constraint(m, ximps[c,d] * TRADE_R[c,2,d] == xtrad_r[c,2,d])
        else
            fix(ximps[c,d], 1.0; force=true)
        end
    end
    for d in 1:nr
        if IMPUSED_C[d] > 1e-10
            @constraint(m, IMPUSED_C[d] * ximpused[d] ==
                sum(TRADE_R[c,2,d] * ximps[c,d] for c in 1:na))
            @constraint(m, IMPUSED_C[d] * pimpused[d] ==
                sum(TRADE[c,2,r,d] * pimp[c,r] for c in 1:na, r in 1:nr))
        else
            fix(ximpused[d], 1.0; force=true); fix(pimpused[d], 1.0; force=true)
        end
        if IMPLANDED_C[d] > 1e-10
            @constraint(m, IMPLANDED_C[d] * ximplanded[d] ==
                sum(TRADE_D[c,2,d] > 1e-10 ? xtrad_d[c,2,d] : 0.0 for c in 1:na))
            @constraint(m, IMPLANDED_C[d] * pimplanded[d] ==
                sum(TRADE_D[c,2,d] * pimp[c,d] for c in 1:na))
        else
            fix(ximplanded[d], 1.0; force=true); fix(pimplanded[d], 1.0; force=true)
        end
    end

    # ── Excerpt 31: MainMacro (column nr+1 = "National", TABLO's REGPLUS) ────
    @variable(m, MainMacro[1:nmac, 1:(nr+1)] >= 1e-6)
    @variable(m, NatMacro[1:nmac] >= 1e-6)
    @variable(m, SelMacro[1:nsel, 1:(nr+1)] >= 1e-6)

    src_of = Dict(
        "RealHou"      => d -> xfin[1,d],   "RealInv"    => d -> xfin[2,d],
        "RealGov"      => d -> xfin[3,d],   "ExpVol"     => d -> xfin[4,d],
        "ImpVolUsed"   => d -> ximpused[d], "ImpsLanded" => d -> ximplanded[d],
        "RealGDP"      => d -> xgdpexp[d],  "RealGNE"    => d -> xgne[d],
        "AggEmploy"    => d -> xlab_io[d],  "realwage_io"=> d -> realwage_io[d],
        "plab_io"      => d -> plab_io[d],  "AggCapStock"=> d -> xcap_i[d],
        "CPI"          => d -> pfin[1,d],   "GDPPI"      => d -> pgdpexp[d],
        "GNEPI"        => d -> pgne[d],     "ExportPI"   => d -> pfin[4,d],
        "ImpsLandedPI" => d -> pimplanded[d], "Population"=> d -> nhou[d],
        "NomHou"       => d -> wfin[1,d],   "NomGDPexp"  => d -> wgdpexp[d],
        "NomGDPinc"    => d -> wgdpinc[d],
    )
    for (name, f) in src_of, d in 1:nr
        @constraint(m, MainMacro[_mac(name), d] == f(d))
    end
    # E_MainMacroV: the "National" column IS NatMacro.
    @constraint(m, [k=1:nmac], MainMacro[k, nr+1] == NatMacro[k])
    @constraint(m, [j=1:nsel, q=1:(nr+1)],
        SelMacro[j,q] == MainMacro[_mac(SELMACROS[j]), q])

    # ── Excerpt 32: weights and national aggregation ────────────────────────
    GNE = [sum(GDPEXPSUM[d,k] for k in 1:4) for d in 1:nr]
    wsrc = Dict(
        "RealHou" => d -> PUR_CS[u_hou,d], "RealInv" => d -> PUR_CS[u_inv,d],
        "RealGov" => d -> PUR_CS[u_gov,d], "ExpVol"  => d -> PUR_CS[u_exp,d],
        "ImpVolUsed" => d -> IMPUSED_C[d], "ImpsLanded" => d -> IMPLANDED_C[d],
        "RealGDP" => d -> GDPEXP[d],       "RealGNE" => d -> GNE[d],
        "AggEmploy" => d -> LAB_IO[d],     "realwage_io" => d -> LAB_IO[d],
        "plab_io" => d -> LAB_IO[d],       "AggCapStock" => d -> CAP_I[d],
        "CPI" => d -> PUR_CS[u_hou,d],     "GDPPI" => d -> GDPEXP[d],
        "GNEPI" => d -> GNE[d],            "ExportPI" => d -> PUR_CS[u_exp,d],
        "ImpsLandedPI" => d -> IMPLANDED_C[d], "Population" => d -> POP[d],
        "NomHou" => d -> PUR_CS[u_hou,d],  "NomGDPexp" => d -> GDPEXP[d],
        "NomGDPinc" => d -> GDPINC[d],
    )
    WMAIN = zeros(Float64, nmac, nr)
    for (name, f) in wsrc, d in 1:nr
        WMAIN[_mac(name), d] = f(d)
    end
    WNAT = [sum(WMAIN[k,d] for d in 1:nr) for k in 1:nmac]

    for k in 1:nmac
        if WNAT[k] > 1e-10
            @constraint(m, WNAT[k] * NatMacro[k] ==
                sum(WMAIN[k,d] * MainMacro[k,d] for d in 1:nr))
        else
            fix(NatMacro[k], 1.0; force=true)
        end
    end

    # contMainMacro is a TABLO `(change)` contribution: in %-change form it reads
    # RATIOMMACRO*WMAIN*MainMacro, with MainMacro a percentage change. Here
    # MainMacro is a bmk=1 index, so the log-gap is the faithful counterpart and
    # the equation is 0 = 0 at the benchmark. RATIOMMACRO is 1 initially (it is
    # `Update`d between periods, which only a multi-period driver would do).
    @variable(m, contMainMacro[1:nmac, 1:nr])
    for k in 1:nmac, d in 1:nr
        if WNAT[k] > 1e-10
            @constraint(m, WNAT[k] * contMainMacro[k,d] ==
                WMAIN[k,d] * log(MainMacro[k,d] + _TINY))
        else
            fix(contMainMacro[k,d], 0.0; force=true)
        end
    end

    # ── Balance-of-trade shares and investment/GDP (Excerpt 32 tail) ────────
    @variable(m, shrBoT); @variable(m, shrBoTnom)
    @variable(m, shrBoTnom2); @variable(m, finvgdp)
    kX, kM = _mac("ExpVol"), _mac("ImpsLanded")
    kG, kNG = _mac("RealGDP"), _mac("NomGDPexp")
    kXP, kMP = _mac("ExportPI"), _mac("ImpsLandedPI")
    lg(v) = log(v + _TINY)

    @constraint(m, (WNAT[kX] - WNAT[kM]) * lg(NatMacro[kG]) + WNAT[kG] * shrBoT ==
        WNAT[kX] * lg(NatMacro[kX]) - WNAT[kM] * lg(NatMacro[kM]))
    @constraint(m, (WNAT[kX] - WNAT[kM]) * lg(NatMacro[kNG]) + WNAT[kG] * shrBoTnom ==
        WNAT[kX] * (lg(NatMacro[kX]) + lg(NatMacro[kXP])) -
        WNAT[kM] * (lg(NatMacro[kM]) + lg(NatMacro[kMP])))
    @constraint(m, (WNAT[kX] - WNAT[kM]) * (shrBoTnom2 + lg(NatMacro[kNG])) ==
        WNAT[kX] * (lg(NatMacro[kX]) + lg(NatMacro[kXP])) -
        WNAT[kM] * (lg(NatMacro[kM]) + lg(NatMacro[kMP])))
    @constraint(m, finvgdp == lg(NatMacro[_mac("RealInv")]) - lg(NatMacro[kG]))

    for (nm, v) in (("ximps", ximps), ("ximpused", ximpused), ("pimpused", pimpused),
                    ("ximplanded", ximplanded), ("pimplanded", pimplanded),
                    ("MainMacro", MainMacro), ("NatMacro", NatMacro),
                    ("SelMacro", SelMacro), ("contMainMacro", contMainMacro),
                    ("shrBoT", shrBoT), ("shrBoTnom", shrBoTnom),
                    ("shrBoTnom2", shrBoTnom2), ("finvgdp", finvgdp))
        vars[nm] = v
    end
    return m, vars
end
