# `build_equations.jl` — Equation & Variable Translation

Part-by-part translation of `src/build_equations.jl` and the build, closure, dynamics and macro files of IndotermJulia — a levels-form TERM-family CGE model of Indonesia — from `origin/TERM.TAB` and `origin/TERM.CMF` into JuMP constraints. Each entry names the TERM.TAB excerpt, the JuMP constraint, the variable shapes and the benchmark convention, cross-referenced to `src/build_model!.jl` (variable declarations), `src/prepare_parameters.jl` (shares/elasticities), `src/ces_helper.jl` (CES/CET primitive), `src/initialize_model!.jl` (closure fix/free) and `src/benchmark_levels.jl` (start values).

**The port's defining constraint** is that this is a **levels formulation** — every GEMPACK `%_change` linear equation becomes its exact nonlinear identity whose first-order log-differential around the benchmark reproduces the TABLO source. There is no Johansen / linearised solve path. The model is square `F(x) = 0` and solved with a sparse Newton on the MOI nonlinear evaluator (`src/solve_newton!.jl`).

---

## Part 1 — Project architecture and set conventions

### Set dimensions

All equation functions receive these scalars, computed in `build_model!` / `build_model_full!` from `agg["MAKE"]`:

| Symbol | Source | Value at 25×6 |
|---|---|---|
| `na` | `size(MAKE, 1)` | 25 commodities |
| `nr` | `size(MAKE, 3)` | 6 island regions |
| `ns` | hardcoded | 2 sources (domestic, import) |
| `nu` | `na + 4` | 29 users (25 ind + HOU/INV/GOV/EXP) |
| `no` | `size(V1LAB, 2)` | occupation count |
| `nm` | `size(TMAR, 3)` | 9 margin commodities |

### User-index constants

```julia
u_hou = na + 1   # 26 — household final demand
u_inv = na + 2   # 27 — investment
u_gov = na + 3   # 28 — government
u_exp = na + 4   # 29 — exports
```

Users 1:na are intermediate (industries driving intermediate demand).

### The three equation classes

Every constraint in the model falls into one of three strategies for converting GEMPACK's `%_change` linear form into a genuine levels equation:

**Class A — Principle-A plain sums.** Quantity components share the same benchmark unit (value-flows) so their aggregation is a literal sum: `xcom = Σ xmake`, `xtrad_d = Σ_d xtrad`. One-weight-per-component, exact at the benchmark. Used for: `E_xcomA_B!`, `E_xtrad_d!`, `E_xtrad_r!`, `E_xint_i!`, `E_xsuppmar_p!`, `E_xsuppmar_d!`, `E_xsuppmar_rd!`, `E_xhouh_s_agg!`, `E_xgov_s!`, `E_xexp_s!`.

**Class B — Principle-B value identities (Euler's theorem).** Aggregate price × aggregate current quantity = Σ component price × component current quantity. Exact for any CES/CET nest at any point. Used for: `E_ppur_s!`, `E_pprim!`, `E_plab_o!`, `E_xtotA_B!`.

**Class C — Log-linear ratio identities.** GEMPACK's `x% := a·p% + b·q%` is the first-order expansion of `log(X/X₀) = a·log(P/P₀) + b·log(Q/Q₀)` or `X/X₀ = (P/P₀)^a · (Q/Q₀)^b`. Every export demand, government demand, margin demand, sourcing and tax-revenue equation uses this form. Used for: `E_xexpd!`, `E_xgov!`, `E_xtrad!`, `E_xsuppmar!`, `E_xtradmar_na!`, `E_xstocks!`, `E_ggro!`, `E_fgovtot2!`, `E_fgovtot3!`, `E_delTAX*!`.

### Benchmark conventions

| Variable type | Benchmark value | Declared bounds | Identified by |
|---|---|---|---|
| Price index (pdom, pimp, phi, pbasic, …) | 1.0 | `>= 1e-6` | `_is_ratio_type` |
| Quantity index (xtot, xfin, natxfin, …) | 1.0 | `>= 1e-6` | `_is_ratio_type` |
| Flow quantity (xint, xhou, xmake, …) | own data value | `>= 0` | listed in `benchmark_levels` |
| Nominal value (whouhtot, wlux) | own data value | `>= 0` | listed in `benchmark_levels` |
| Shifter/`(change)` delta (delPTX, delTAX*, …) | 0.0 | unbounded / `>= -0.99` | in `ZERO_START_NAMES` |
| Additive log-gap (gret, ggro, pfexp, …) | 0.0 | unbounded | in `ZERO_START_NAMES` |

### _TINY and _MIN_PRICED_FLOW

```julia
const _TINY = 1e-9
const _MIN_PRICED_FLOW = 1e5 * _TINY   # = 1e-4
```

`_TINY` guards `log(p*X + _TINY)` inside value identities so the argument never reaches exactly zero. `_MIN_PRICED_FLOW` is the minimum benchmark flow for which a `log(price*flow + _TINY)` constraint determines its price — below this threshold the price column's numerical weight drops below ~1e4 and Newton amplifies noise. 26 of 4350 `PUR_S` cells fall below, all carrying ≤ 1e-4 Rp bn — economically nil, numerically fatal (a single such cell throttled a `blabnat` step to α ≈ 1.2e-3, before NaNs). For these cells the value identity is replaced by a share-weighted price average, homogeneous of degree 1 and O(1)-scaled.

---

## Part 2 — Lines 81–108: Excerpt 6-7, Basic prices and purchaser prices

### E_pimp! — import price = foreign-currency price × exchange rate

```julia
function E_pimp!(m, vars, na, nr, params)
    pimp = vars["pimp"]; pfimp = vars["pfimp"]; phi = vars["phi"]
    @constraint(m, [c=1:na, r=1:nr], log(pimp[c,r]) == log(pfimp[c]) + log(phi))
end
```

| Variable | Shape | Bmk | Source | Meaning |
|---|---|---|---|---|
| `pimp[c,r]` | na × nr | 1.0 | `>= 1e-6` | Import price in domestic currency |
| `pfimp[c]` | na | 1.0 | `>= 1e-6` | Foreign-currency import price (exogenous under small-country assumption) |
| `phi` | scalar | 1.0 | `>= 1e-6` | Exchange rate (numeraire under default closure) |

**Concept — multiplicative price pass-through.** GEMPACK's `pimp% = pfimp% + phi%` is the log-differential of `pimp = pfimp × phi`. In levels: `log(pimp) = log(pfimp) + log(phi)`.

### E_pbasic! — basic price = domestic or import price by source

```julia
function E_pbasic!(m, vars, na, nr, ns)
    pbasic = vars["pbasic"]; pdom = vars["pdom"]; pimp = vars["pimp"]
    for s_dom in 1:ns
        p = s_dom == 1 ? pdom : pimp
        @constraint(m, [c=1:na, r=1:nr], pbasic[c,s_dom,r] == p[c,r])
    end
end
```

| Variable | Shape | Bmk | Meaning |
|---|---|---|---|
| `pbasic[c,s,r]` | na × ns × nr | 1.0 | Basic (producer) price by source-region |
| `pdom[c,r]` | na × nr | 1.0 | Domestic output price |
| `pimp[c,r]` | na × nr | 1.0 | Import price |

`pbasic[c,1,r]` (domestic) ≡ `pdom[c,r]`; `pbasic[c,2,r]` (import) ≡ `pimp[c,r]`.

### E_tuser! — tax-power decomposition by user type and commodity

```julia
function E_tuser!(m, vars, na, nr, ns, nu)
    tuser = vars["tuser"]; tuser_ud = vars["tuser_ud"]; tuser_su = vars["tuser_su"]
    @constraint(m, [c=1:na, s=1:ns, u=1:nu, d=1:nr],
        log(tuser[c,s,u,d]) == log(tuser_ud[c,s]) + log(tuser_su[c,d]))
end
```

| Variable | Shape | Meaning |
|---|---|---|
| `tuser[c,s,u,d]` | na × ns × nu × nr | Tax "power" (1 + ad-valorem rate) by commodity/source/user/region |
| `tuser_ud[c,s]` | na × ns | User-dimension component |
| `tuser_su[c,d]` | na × nr | Source-dimension component |
| `tuser_sud[c]` | na | National residual |

### E_ppur! — purchaser price = basic price × tax × national residual

```julia
function E_ppur!(m, vars, na, nr, ns, nu)
    ppur = vars["ppur"]; puse = vars["puse"]
    tuser = vars["tuser"]; tuser_sud = vars["tuser_sud"]
    @constraint(m, [c=1:na, s=1:ns, u=1:nu, d=1:nr],
        log(ppur[c,s,u,d]) == log(puse[c,s,d]) + log(tuser[c,s,u,d]) + log(tuser_sud[c]))
end
```

| Variable | Shape | Bmk | Meaning |
|---|---|---|---|
| `ppur[c,s,u,d]` | na × ns × nu × nr | 1.0 | Purchaser (tax-inclusive) price by commodity/source/user/region |
| `puse[c,s,d]` | na × ns × nr | 1.0 | User price (basic+delivery cost), the "delivered" price from Excerpt 20 |
| `tuser_sud[c]` | na | 1.0 | National-level tax-power residual |

**Concept — the three-layer price system.** `tuser` decomposes the wedge between user price and purchaser price multiplicatively: `ppur = puse × tuser × tuser_sud`. Under a zero-shock benchmark all prices and tax powers are 1.

---

## Part 3 — Lines 111–233: Excerpt 8, Armington CES substitution

### E_ppur_s! — source-composite price identity (Principle B, with dust-flow guard)

```julia
function E_ppur_s!(m, vars, na, nr, ns, nu, params)
    ppur_s = vars["ppur_s"]; ppur = vars["ppur"]
    PUR_S = parent(params["PUR_S"]); SRCSHR = parent(params["SRCSHR"])
    for c in 1:na, i in 1:na, d in 1:nr
        if PUR_S[c,i,d] > _MIN_PRICED_FLOW
            @constraint(m, log(ppur_s[c,i,d] * xint_s[c,i,d] + _TINY) ==
                log(sum(xint[c,s,i,d] * ppur[c,s,i,d] for s in 1:ns) + _TINY))
        else
            @constraint(m, ppur_s[c,i,d] ==
                sum(SRCSHR[c,s,i,d] * ppur[c,s,i,d] for s in 1:ns))
        end
    end
    # … identical pattern for hou, inv, gov, exp
end
```

| Variable | Shape | Bmk | Meaning |
|---|---|---|---|
| `ppur_s[c,u,d]` | na × nu × nr | 1.0 | Source-composite (Armington aggregate) price |
| `xint_s[c,i,d]` | na × na × nr | PUR_S[c,i,d] | Armington-intermediate quantity |
| `SRCSHR[c,s,u,d]` | na × ns × nu × nr | value share | Benchmark domestic/import share |

**Concept — the two-branch identity.** The `PUR_S > _MIN_PRICED_FLOW` branch writes the genuine value identity `ppur_s·X_s = Σ_s ppur·X` as a log equation — exact at the benchmark, homogeneous of degree 1, and its derivative w.r.t. `ppur_s` is `X_s/(ppur_s·X_s + _TINY) ≈ 1` for flows above threshold. The dust-flow branch replaces this with a fixed-share-weighted price average `ppur_s = Σ_s SRCSHR·ppur` — the exact degenerate case of the identity, homogeneous of degree 1, and yielding 1 at benchmark.

**Prior bug: homogeneity-breaker.** The original dust-flow branch read `fix(ppur_s, 1.0)`, which breaks degree-zero homogeneity: scaling the numeraire by λ leaves these 26 cells at 1.0 while every other price becomes λ. The share-weighted average fixes this — `pcap(5,5)` and `pcap(5,6)` dropped from 3.7e-4 / 4.1e-4 error to 2.4e-15 / 2.2e-15.

### E_phou! — household Armington price = source-composite price for user HOU

```julia
function E_phou!(m, vars, na, nr)
    phou = vars["phou"]; ppur_s = vars["ppur_s"]
    u_hou = na + 1
    @constraint(m, [c=1:na, d=1:nr], phou[c,d] == ppur_s[c,u_hou,d])
end
```

### E_xint! / E_xhou! / E_xinv! — Armington CES demand

```julia
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
```

**Concept — genuine CES via `ces()`.** `ces(y, p, α, σ, γ)` returns the vector of cost-minimizing component demands for aggregate `y` at component prices `p`. For Armington, `σ = sigmadomimp[c]` (from `P015`, default 5.0) — the elasticity of substitution between domestic and imported varieties. `ces_calibrate` (in `ces_helper.jl`) fits the share parameters `ALPHA/GAMMA` so that at all-prices=1, `ces(x_s, ones(ns), α, σ, γ) = PUR[c,s,u,d]` exactly. `E_xhou!` and `E_xinv!` are identical in structure, differing only in the user index (HOU/INV).

---

## Part 4 — Lines 237–264: Excerpt 9, Intermediate demands (Leontief)

### E_aint_s! — intermediate technical-change shifters by source and from-all-industries

```julia
function E_aint_s!(m, vars, na, nr)
    aint_s = vars["aint_s"]; bint_scd = vars["bint_scd"]; bint_s = vars["bint_s"]
    @constraint(m, [c=1:na, i=1:na, d=1:nr],
        log(aint_s[c,i,d]) == log(bint_scd[i]) + log(bint_s[c,i,d]))
end
```

### E_xint_s! — intermediate demand = benchmark × atot × aint_s × xtot

```julia
function E_xint_s!(m, vars, na, nr, params)
    xint_s = vars["xint_s"]; atot = vars["atot"]
    aint_s = vars["aint_s"]; xtot = vars["xtot"]
    PUR_S = parent(params["PUR_S"])
    @constraint(m, [c=1:na, i=1:na, d=1:nr],
        xint_s[c,i,d] == PUR_S[c,i,d] * atot[i,d] * aint_s[c,i,d] * xtot[i,d])
end
```

**Concept — Leontief multiplier.** Every intermediate `(c,i,d)` is used in fixed proportion to industry `i`'s output `xtot[i,d]`, with sector-level `atot` and commodity-level `aint_s` shifters. No Armington substitution at the intermediate layer — that is already in `E_xint!` upstream.

### E_pint! — intermediate price index = benchmark-share-weighted average

```julia
function E_pint!(m, vars, na, nr, params)
    pint = vars["pint"]; ppur_s = vars["ppur_s"]; aint_s = vars["aint_s"]
    PUR_S = parent(params["PUR_S"]); PUR_CS = parent(params["PUR_CS"])
    for i in 1:na, d in 1:nr
        PUR_CS[i,d] > 1e-10 || continue
        rhs = sum(PUR_S[c,i,d] * ppur_s[c,i,d] * aint_s[c,i,d] for c in 1:na)
        @constraint(m, PUR_CS[i,d] * pint[i,d] == rhs)
    end
end
```

Guard: `|| continue` skips cells with no intermediate demand — avoids a `0 = 0` dead row.

---

## Part 5 — Lines 268–398: Excerpt 10-11, Labour and factor demands (CES)

### E_xlab! — labour by occupation = CES over skill types

```julia
function E_xlab!(m, vars, na, nr, no, sigmalab, params)
    dem = ces(xlab_o[i,d], [plab[i,o,d] for o in 1:no], αv, sigmalab[i], GAMMA[i,d])
    @constraint(m, [o=1:no], xlab[i,o,d] == dem[o])
end
```

Elasticity: `sigmalab[i]` from `SLAB` (default 0.5). `ALPHA_LAB/GAMMA_LAB` calibrated at `ces_calibrate(V1LAB[i,:,d], sigma, LAB_O[i,d])`.

### E_plab_o! — composite labour price identity (Principle B)

```julia
@constraint(m, [i=1:na, d=1:nr],
    log(plab_o[i,d] * xlab_o[i,d] + _TINY) ==
    log(sum(plab[i,o,d] * xlab[i,o,d] for o in 1:no) + _TINY))
```

### E_wlab_o! — labour value = price × quantity (normalized by benchmark LAB_O)

```julia
@constraint(m, LAB_O[i,d] * wlab_o[i,d] == plab_o[i,d] * xlab_o[i,d])
```

### E_xlab_o! / E_pcap! / E_plnd! — factor demands from the primary CES nest

These three together form the bottom layer of the primary-factor CES nest (Excerpt 11). Nest order `[LAB_O, CAP, LND]` matches `ALPHA_FAC/GAMMA_FAC` column convention:

```julia
# E_xlab_o!:  xlab_o = alab_o * ces(xprim, [plab_o*alab_o, pcap*acap, plnd*alnd], ...)[1]
# E_pcap!:    xcap   = acap   * ces(xprim, [plab_o*alab_o, pcap*acap, plnd*alnd], ...)[2]
# E_plnd!:    xlnd   = alnd   * ces(xprim, [plab_o*alab_o, pcap*acap, plnd*alnd], ...)[3]
```

Elasticity: `sigmaprim[i]` from `P028` (default 0.5). Each factor has a technical-change augmenting shifter `(alab_o/acap/alnd)` multiplying both its effective quantity and price — the effective price facing the CES nest is `p × a`, so a technical improvement (`a`↑) lowers the effective cost of that factor.

**E_plnd! zero-flow guard.** Sectors with no benchmark land use (`LND[i,d] ≤ 1e-10`) get `plnd[i,d]` pinned at 1.0 rather than a meaningless `0 = 0` constraint — 15 of 150 cells, the 15 of 637 null directions.

### E_pprim! — factor-composite value identity (Principle B)

```julia
@constraint(m, [i=1:na, d=1:nr],
    log(pprim[i,d] * xprim[i,d] + _TINY) == log(
        plab_o[i,d]*xlab_o[i,d] + pcap[i,d]*xcap[i,d] + plnd[i,d]*xlnd[i,d] + _TINY))
```

The `a`-augmentation shifters cancel exactly: `E_pprim!` is Euler's theorem for the CES nest and is exact for every point, not just the benchmark.

### E_xprim! — factor composite quantity = benchmark × xtot × atot × aprim

```julia
xprim[i,d] == PRIM[i,d] * xtot[i,d] * atot[i,d] * aprim[i,d]
```

### E_aprim! / E_alab_o! — technical-change decomposition

```julia
log(aprim[i,d]) == log(bprimnat) + log(bprim_d[i]) + log(bprim[i,d])
log(alab_o[i,d]) == log(blabnat) + log(blab_d[i]) + log(blab[i,d])
```

Both decompose into national (`bprimnat`/`blabnat`), sector-dimension (`bprim_d`/`blab_d`) and sector×region (`bprim`/`blab`) multiplicative components. `blabnat` is TERM.CMF's primary shock instrument: `Shock blabnat = -3`.

### E_wprim! — factor-cost value index (normalized by fixed benchmark PRIM)

```julia
PRIM[i,d] * wprim[i,d] == pcap[i,d]*xcap[i,d] + plnd[i,d]*xlnd[i,d] + Σ_o plab[i,o,d]*xlab[i,o,d]
```

---

## Part 6 — Lines 402–461: Excerpt 12, Output prices and production tax

### E_pvar! — variable-cost price (labour + intermediate)

```julia
pvar[i,d] * VARCST[i,d] == atot[i,d] * (LAB_O[i,d]*plab_o[i,d]*alab_o[i,d] + PUR_CS[i,d]*pint[i,d])
```

`VARCST[i,d]` is benchmark variable cost (labour + intermediate). At benchmark `atot=alab_o=plab_o=pint=pvar=1`.

### E_pcst! — ex-tax cost price (variable + fixed-factor)

```julia
pcst[i,d] * VCST[i,d] == atot[i,d] * (PRIM[i,d]*aprim[i,d]*pprim[i,d] + PUR_CS[i,d]*pint[i,d])
```

### E_delPTX! — production-tax revenue change (not a level)

```julia
delPTX[i,d] == PRODTAX[i,d] * (pcst[i,d]*xtot[i,d] - 1.0)
             + VCST[i,d] * pcst[i,d] * xtot[i,d] * delPTXRATE[i,d]
```

`delPTX` is a `(change)` variable — benchmark 0, not a level. The formula: (benchmark revenue) × (ex-tax cost deviation) + (current ex-tax value) × (rate change). At benchmark `pcst=xtot=1, delPTXRATE=0` ⇒ `delPTX=0`. This is the exact nonlinear version of TERM.TAB:547-548's linearization.

### E_ptot! — tax-inclusive output value

```julia
ptot[i,d] * xtot[i,d] * VTOT[i,d] == pcst[i,d] * xtot[i,d] * VCST[i,d] + PRODTAX[i,d] + delPTX[i,d]
```

At benchmark `ptot=xtot=pcst=1, delPTX=0, VTOT=VCST+PRODTAX` ⇒ holds exactly.

---

## Part 7 — Lines 466–566: Excerpt 13, Household demands (ELES / Stone-Geary)

INDOTERM's household block uses the Extended Linear Expenditure System (Lluch 1973). Benchmark household purchases `HOUPUR(c,d)` are split into subsistence `(1-BLUX)·HOUPUR` and luxury/supernumerary `BLUX·HOUPUR`, where `BLUX/SLUX` are the marginal-budget / supernumerary-budget shares from `prepare_parameters.jl`.

### E_asub! — subsistence shifter (geometric-mean-normalized)

```julia
geomean = Σ_k BUDGSHR[k,1,d] * log(ahou_s[k,d] + _TINY)
asub[c,d] == ahou_s[c,d] / exp(geomean)
```

`ahou_s` is TABLO-`Omit`ted (`origin/TERM.TAB:2575`) — permanently fixed at its benchmark 1.0. `asub` is its budget-share-geometric-mean-normalized version, which preserves the budget-share-weighted mean of `log(asub)` at exactly log(1) = 0.

### E_alux! — luxury shifter (calibrated non-unity benchmark)

```julia
geomean = Σ_k SLUX[k,1,d] * log(asub[k,d] + _TINY)
alux[c,d] == ALUX0[c,1,d] * asub[c,d] / exp(geomean)
```

`ALUX0` is the calibrated per-commodity supernumerary share `XLUX0(c,d)/WLUX0(d)` — `alux` has benchmark `ALUX0`, not 1.0.

### E_xsub! — subsistence quantity

```julia
xsub[c,d] == XSUB0[c,1,d] * nhou[d] * asub[c,d]
```

`nhou[d]` is the number of households (population index), benchmark 1.

### E_xlux! — luxury quantity = supernumerary expenditure / price × share

```julia
xlux[c,d] * phou[c,d] == wlux[d] * alux[c,d]
```

`wlux[d]` is nominal aggregate supernumerary expenditure (benchmark `WLUX0[d]`). This is the standard Stone-Geary supernumerary-demand equation: `x = (w·α)/p`.

### E_xhouh_s_agg! — total household composite = subsistence + luxury

```julia
xhou_s[c,d] == xlux[c,d] + xsub[c,d]
```

### E_wlux! — household volume index (historical name for E_xhouhtot)

```julia
xhoutot[d] == Σ_c xhou_s[c,d] / HOUPUR_C[1,d]
```

All composite quantities share the same benchmark flow unit (`HOUPUR`), so the sum is a Divisia index normalized by total benchmark consumption.

### E_phouhtot! — consumer price index (Laspeyres-style)

```julia
phouhtot[d] == Σ_c BUDGSHR[c,1,d] * phou[c,d]
```

### E_whouhtot! — total nominal household expenditure

```julia
whouhtot[d] == HOUPUR_C[1,d] * phouhtot[d] * xhoutot[d]
```

`whouhtot` is a genuine nominal value, benchmark `HOUPUR_C[d]`.

### E_xhoutot! / E_phoutot! — national household aggregates (single-type collapse)

```julia
xhoutot[d] == xhouhtot[d]
phoutot[d] == phouhtot[d]
```

With only one household type these are identity mappings.

---

## Part 8 — Lines 570–673: Excerpt 14-15, Investment demands and the flexible accelerator

### E_xinvi! — investment by industry (Leontief)

```julia
xinvi[c,i,d] == INVEST[c,i,d] * xinvitot[i,d]
```

`INVEST` is the raw `V2PUR` matrix — industry `i`'s benchmark investment composition by commodity. `xinvitot` is a real investment-volume index (benchmark 1).

### E_pinvest! / E_pinvitot! — investment price indices

`pinvest[c,d]` = `ppur_s[c,u_inv,d]` — the investment composite price equals the Armington-sourcing price for the investment user. `pinvitot[i,d]` is the `INVEST_C`-weighted average of `pinvest[c,d]` over commodities.

### E_xinv_s! — investment composite = ratio identity across data sources

```julia
xinv_s[c,d] / PUR_S[c,u_inv,d] == Σ_i xinvi[c,i,d] / INVEST_I[c,d]
```

**Prior bug — data-source mismatch.** A literal `xinv_s = Σ xinvi` assumes both benchmarks match, but `PUR_S[inv]` (from the purchaser matrix) and `INVEST_I = Σ_i INVEST` (from `2PUR`) come from different source data and differ by ~33570 in the worst cell. The ratio form `xinv_s/bmk_s = Σ xinvi/bmk_inv` is exact at the benchmark and linearises back to TERM.TAB's E_xinv_s.

### E_gret! / E_xinvitot! / E_ggro! / E_fgret! / E_finv2! — the flexible accelerator

These five equations form Excerpt 15's investment-rule block. All operate in logs:

| Equation | Julia constraint | Meaning |
|---|---|---|
| `E_gret!` | `gret = log(pcap) - log(pinvitot)` | Rate-of-return gap (rental vs. investment-good price) |
| `E_xinvitot!` | `log(xinvitot) - log(xcap/CAP) = ggro` | Investment growth = capital growth + gap |
| `E_ggro!` | `ggro = finv1 + 0.33·(2·gret - invslack)` | Investment responds with 0.66 elasticity to rate-of-return gap |
| `E_fgret!` | `gret = fgret + capslack` | Rate-of-return decomposition |
| `E_finv2!` | `log(xinvitot) = finv2 + log(xgdpexp)` | Investment follows GDP |

`gret, ggro, fgret, finv2` are additive log-gaps (benchmark 0, in `ZERO_START_NAMES`).

---

## Part 9 — Lines 676–743: Excerpt 16, Government, exports, inventories

### E_xgov! — government demand = benchmark × multiplicative shifters

```julia
xgov[c,s,d] == XGOV0[c,s,d] * exp(fgovtot[d] + fgov[c,s,d] + fgov_s[c,d] + fgovgen)
```

`fgov*` shifters are additive (benchmark 0). `exp(sum)` gives the multiplier from benchmark.

### E_fgovtot2! / E_fgovtot3! — government demand follows consumption or GDP

```julia
fgovtot[d] == fgovtot2[d] + log(xhoutot[d])
fgovtot[d] == fgovtot3[d] + log(xgdpexp[d])
```

Under `TERM.CMF`'s closure: `swap fgovtot = fgovtot3` makes regional real government spending follow regional real GDP.

### E_pfexp! — export price in foreign currency

```julia
pfexp[c,d] == log(ppur[c,1,u_exp,d]) - log(phi)
```

`pfexp` is a log-based export price (benchmark 0). The domestic export price `ppur[dom,exp]` divided by `phi` gives the foreign-currency price.

### E_xexpd! — export demand with constant-elasticity price response

```julia
log(xexpd[c,d]) - log(XEXPD0[c,d]) ==
    natfqexp + fqexp[c,d] + fqexp_d[c] -
    exp_elast[c] * (pfexp[c,d] - fpexp[c,d] - fpexp_d[c] - natfpexp)
```

`exp_elast[c]` from `P018` (default 2.0) is the export-demand elasticity. `fqexp`/`fpexp` shifters shift the world demand/price curve outward — a coal-price shock is `fpexp_d("Coal") = logpct(50)`.

### E_xexp! / E_xexp_s! — export allocation and composite

`xexp[c,1,d] == xexpd[c,d]` (all export is domestic-source), `xexp[c,2,d] == 0`, and `xexp_s[c,d] = Σ_s xexp[c,s,d]`.

### E_xstocks! — inventories as log-additive deviation from output

```julia
xstocks[i,d] == log(xtot[i,d]) + fxstocks[i,d]
```

`fxstocks` is TABLO-`Omit`ted (pinned at 0), so `xstocks` purely tracks `log(xtot)`.

---

## Part 10 — Lines 746–880: Excerpt 17, 19-20, Total demand / delivery / sourcing

### E_xint_i! — total intermediate use by commodity-source

```julia
xint_i[c,s,d] == sum(xint[c,s,i,d] for i in 1:na)
```

### E_xuse! — basic-value-weighted aggregation of all users

```julia
(USE_U/TRADE_R) * xuse[c,s,d] ==
    (USE_I/PUR_int_i) * xint_i[c,s,d] +
    (USE[u_hou]/PUR[u_hou]) * xhou[c,s,d] + …
```

**Prior bug — the "literal-equality" class.** A plain `xuse = Σ user` wrongly forces basic and purchaser benchmarks to match. The faithful levels form scales each variable by its own benchmark and cross-weights by TERM.TAB's `USE_U/USE_I/USE[fin]` basic-value coefficients. At the benchmark, LHS = `USE_U` and RHS = `USE_I + Σ_fin USE[fin]` = `USE_U` exactly.

### E_pdelivrd! — delivered price = basic + margins

```julia
log(pdelivrd[c,s,r,d]) ==
    BASSHR[c,s,r,d] * log(pbasic[c,s,r]) +
    Σ_m MARSHR[c,s,m,r,d] * (log(psuppmar_p[m,r,d]) + atradmar[c,s,m,r,d])
```

`BASSHR + Σ_m MARSHR = 1` by construction — a Divisia-style price index of the basic good and the margin services.

### E_xtradmar_na! — margin quantity = trade quantity × (1 + atradmar shock)

```julia
log(xtradmar[c,s,m,r,d]/TMAR[c,s,m,r,d]) ==
    log(xtrad[c,s,r,d]/TRADE[c,s,r,d]) + atradmar[c,s,m,r,d]
```

Leontief in the margin dimension: margin flow moves proportionally with trade flow.

### E_puse! — user price = DELIVRD-weighted average of delivered prices

```julia
DELIVRD_R[c,s,d] * log(puse[c,s,d]) ==
    Σ_r DELIVRD[c,s,r,d] * (log(pdelivrd[c,s,r,d]) + atrad[c,s,r,d])
```

### E_xtrad! — bilateral trade = CES over supply regions

```julia
log(xtrad[c,s,r,d]/TRADE[c,s,r,d]) - atrad[c,s,r,d] ==
    log(xuse[c,s,d]/TRADE_R[c,s,d]) + srctwist[c,s,r,d] - avesrctwist[c,s,d] -
    SGDD[c] * (log(pdelivrd[c,s,r,d]) + atrad[c,s,r,d] - log(puse[c,s,d]))
```

`SGDD[c]` (σ, the CES elasticity of substitution across regional sources) from `SMAR`.

### E_avesrctwist! — source-twist average

```julia
DELIVRD_R[c,s,d] * avesrctwist[c,s,d] == Σ_r DELIVRD[c,s,r,d] * srctwist[c,s,r,d]
```

---

## Part 11 — Lines 884–942: Excerpt 21, Margin supply

### E_xsuppmar_p! — total margin flow by route

```julia
xsuppmar_p[m,r,d] == sum(xtradmar[c,s,m,r,d] for c in 1:na, s in 1:ns)
```

### E_psuppmar_p! — margin composite price

```julia
SUPPMAR_P[m,r,d] * log(psuppmar_p[m,r,d]) ==
    Σ_p MARS[m,r,d,p] * (log(pdom[m,p]) + asuppmar[m,r,d,p])
```

Zero-flow guard: where `SUPPMAR_P[m,r,d] ≤ 1e-10`, pin `psuppmar_p` at 1.0 (34 of the 637 null directions).

### E_xsuppmar! — margin sourcing CES

```julia
log(xsuppmar[m,r,d,p] / MARS[m,r,d,p]) - asuppmar[m,r,d,p] ==
    log(xsuppmar_p[m,r,d] / SUPPMAR_P[m,r,d]) -
    SMAR[m] * (log(pdom[m,p]) + asuppmar[m,r,d,p] - log(psuppmar_p[m,r,d]))
```

`SMAR[m]` is the margin-sourcing elasticity (CES across production regions `p`).

### E_xsuppmar_d! / E_xsuppmar_rd! — margin aggregation

`xsuppmar_d[m,r,p] = Σ_d xsuppmar`, `xsuppmar_rd[m,p] = Σ_r xsuppmar_d`.

---

## Part 12 — Lines 946–1081: Excerpt 22-23, MAKE/CET and market clearing

### E_xmake! — multi-product CET allocation

```julia
dem = ces(xtot[i,d], [pmake[c,i,d] for c in 1:na], αv, -SCET[i], GAMMA[i,d])
xmake[c,i,d] == dem[c]   # where MAKE[c,i,d] > 1e-10
```

CET uses `-SCET[i]` (negative sigma): revenue-maximizing allocation of industry `i`'s activity `xtot` across commodities. Calibrated via `ces_calibrate(MAKE[:,i,d], -SCET[i], MAKE_I[i,d])`. Zero-flow cells pinned at 0.

### E_xtotA_B! — make revenue identity (Principle B)

```julia
ptot[i,d] * xtot[i,d] * MAKE_C[i,d] == Σ_c pmake[c,i,d] * xmake[c,i,d]
```

### E_xcomA_B! — commodity output = sum over industries

```julia
xcom[c,d] == sum(xmake[c,i,d] for i in 1:na)
```

### E_pmake! — make price = output price with ad hoc supply-curve slope

```julia
log(pmake[c,i,d]) == log(pdom[c,d]) -
    0.05 * (log(xmake[c,i,d]/MAKE[c,i,d]) - log(xcom[c,d]/MAKE_I[c,d]))
```

The 0.05 coefficient is TERM's own reduced-form supply-curve slope (not a CET parameter — the true CET response is in `E_xmake!`).

### E_xtrad_d! / E_xtrad_r! — trade aggregation

```julia
xtrad_d[c,s,r] = Σ_d xtrad[c,s,r,d]     # total supply from region r
xtrad_r[c,s,d] = Σ_r xtrad[c,s,r,d]     # total demand absorbed in d
```

### E_pdomA_sum! — commodity market clearing (two forms)

This one constraint implements both of TERM.TAB's split forms:

**E_pdomA — non-margins (ratio form):**
```julia
xcom[c,r] / MAKE_I[c,r] == xtrad_d[c,1,r] / TRADE_D[c,1,r]
```
The ratio form handles the `MAKE_I ≠ TRADE_D` imbalance at the benchmark (e.g. Crops has surplus output exported to RoW that the domestic trade matrix does not carry).

**E_pdomB — margins (additive-constant form):**
```julia
K = MAKE_I[c,r] - TRADE_D[c,1,r] - Σ SUPPMAR_RD[m,r]
xcom[c,r] == xtrad_d[c,1,r] + Σ xsuppmar_rd[m,r] + K
```
Margin commodities (Trade, Transport, InfoComm — aggregated-commodity indices 20, 21, 22) balance domestic demand plus the margin services they supply, with a benchmark-imbalance constant `K` so the benchmark `xcom=MAKE_I, xtrad_d=TRADE_D, xsuppmar_rd=SUPPMAR_RD` is exact.

---

## Part 13 — Lines 1084–1190: Excerpt 24, Final demand aggregates

### E_pfin! — final-demand price indices (by category)

```julia
PUR_CS[u,d] * pfin[fi,d] == Σ_c PUR_S[c,u,d] * ppur_s[c,u,d]
```

`pfin[1:4,d]` = {HOU, INV, GOV, EXP} prices, benchmark 1, `PUR_CS`-weighted average of the corresponding Armington prices.

### E_xfina! — household volume = xhoutot

```julia
xfin[1,d] == xhoutot[d]
```

### E_xfinb! / E_xfinc! / E_xfind! — investment, government, export volume indices

```julia
PUR_CS[u,d] * xfin[fi,d] == Σ_c flow[c,d]
```

These are bmk=1 QUANTITY INDICES, not flows. The `PUR_CS` divisor is what makes each index benchmark 1 — dropping it was a prior bug that cascaded compensating `PUR_CS` factors into four downstream equations.

### E_wfin! — final-demand value = quantity × price

```julia
wfin[fi,d] == xfin[fi,d] * pfin[fi,d]
```

### E_natpfin! / E_natxfin! / E_natwfin! — national final-demand aggregates

`PUR_CS`-weighted averages across regions of the regional `pfin/xfin/wfin` indices.

---

## Part 14 — Lines 1194–1292: Excerpt 26, Commodity tax revenues

### E_delTAX{int,hou,inv,gov,exp}! — tax-revenue change

```julia
delTAX[c,s,u,d] * PUR[c,s,u,d] ==
    x[c,s,u,d] * (ppur[c,s,u,d] * PUR[c,s,u,d] - USE[c,s,u,d] * puse[c,s,d]) - TAX[c,s,u,d] * PUR[c,s,u,d]
```

All five user categories share the identical formula. The core idea: `TAX = PUR - USE` is the benchmark tax per unit × flow, `ppur*PUR` is the current purchaser value, and `USE*puse` is the current basic value. Their difference minus benchmark tax gives the change.

`PUR`, `USE` and `TAX` are drawn from Dict caches (`PUR_idx`, `USE_idx`, `TAX_idx`) populated by `TAX_PUR_setup!`. Zero-flow cells (`PUR ≤ 1e-10`) pin `delTAX = 0` — no purchaser flow means no tax base, and the behaviour under any shock is structurally zero.

---

## Part 15 — Lines 1296–1462: Excerpt 27, Primary factor aggregates

### E_wlnd_i! / E_wcap_i! / E_wprim_i! — regional factor values

```julia
LND_I[d] * wlnd_i[d] == Σ_i plnd[i,d] * xlnd[i,d]
CAP_I[d] * wcap_i[d] == Σ_i pcap[i,d] * xcap[i,d]
PRIM_I[d] * wprim_i[d] == Σ_i PRIM[i,d] * wprim[i,d]
```

### E_xlnd_i! / E_xcap_i! — regional factor quantity indices

```julia
LND_I[d] * xlnd_i[d] == Σ_i xlnd[i,d]
CAP_I[d] * xcap_i[d] == Σ_i xcap[i,d]
```

Since `xlnd/xcap` already carry the `LND/CAP` flow, the RHS weight is 1 (plain sum), not `LND/CAP` again.

### E_plab! — real wage = nominal wage / household price index

```julia
realwage[i,o,d] * pfin[1,d] == plab[i,o,d]
```

### National and occupation-level aggregations (plab_i, xlab_i, plab_id, …)

A ladder of `LAB`-weighted averages from sector-level `(i,o,d)` to occupation-level `(o,d)` to national `(o)` to all-occupation `(d)`, with corresponding price, quantity, value and real-wage pairs.

---

## Part 16 — Lines 1465–1760: Excerpt 28-29, Income-side and expenditure-side GDP

### E_delGDPINCa-e! — income-side GDP components

`delGDPINC[d,k]` for k = {Land, Labour, Capital, ProdTax, ComTax}:

| k | Category | Equation |
|---|---|---|
| 1 | Land | `LND_I[d] * (wlnd_i[d] - 1)` |
| 2 | Labour | `LAB_IO[d] * (wlab_io[d] - 1)` |
| 3 | Capital | `CAP_I[d] * (wcap_i[d] - 1)` |
| 4 | ProdTax | `Σ_i delPTX[i,d]` |
| 5 | ComTax | `Σ delTAXint + delTAXhou + delTAXinv + delTAXgov + delTAXexp` |

**Slot bug — GDPINCCAT ordering.** TERM.TAB indexes by category name, not letter: `(Land, Labour, Capital, ProdTax, ComTax)`. The earlier port mapped equation letters positionally `(a→1, b→2, c→3)`, silently swapping Labour and Capital. Totals were unaffected (summation over all k), so benchmark replication was blind to it.

### E_wgdpinc! — nominal income-side GDP index

```julia
gdp_inc * wgdpinc[d] == gdp_inc + Σ_k delGDPINC[d,k]
```

### E_delXGDPEXP*! / E_delPGDPEXP*! — expenditure-side GDP decomposition

9 categories: HOU, INV, GOV, STOCKS, EXP, Imports(−), RExports(+), RImports(−), NetMar.

For quantity components: `delXGDPEXP = Σ current_flow - GDPEXPSUM` (sign-flipped for imports). For price components: `delPGDPEXP = Σ benchmark_flow × (price - 1)`. Stocks uses `exp(xstocks) - 1` for the volume ratio.

### E_xgdpexp! / E_pgdpexp! / E_wgdpexp! — GDP expenditure aggregates

`GDPEXPSUM`-weighted indices mirroring the income-side form.

### E_wgdpdiff! — GDP income/expenditure statistical discrepancy

```julia
wgdpdiff[d] == wgdpinc[d] - wgdpexp[d]
```

At benchmark both are 1 ⇒ discrepancy 0. Under closure the model should keep this near zero.

### GNE and budget aggregates

`xgne/pgne/wgne`: GNE = HOU+INV+GOV+STOCKS (excludes trade balance). `E_delINDTAX!`: indirect tax = ProdTax + ComTax change. `E_delBUDG1!/E_delBUDG2!`: budget balance decompositions.

---

## Part 17 — Lines 1762–1847: Excerpt 38-39, Labour market and household closure

### E_labslack! — national labour supply (real wage vs employment)

```julia
log(realwage_id[o]) == 2.0 * log(xlab_id[o]) + flabsup_id[o]
```

Elasticity 0.5: a 1% real-wage increase raises labour supply by 0.5%.

### E_flab_i! — regional labour supply

```julia
log(realwage_i[o,d]) == (1.0/1.1) * log(xlab_i[o,d]) + flabsupA[o,d] + labslack[o]
```

Elasticity 1.1: a more elastic regional response.

### E_realwage! — wage-shifter decomposition residual

```julia
log(realwage[i,o,d]) == flab[i,o,d] + flab_i[o,d] + flab_io[d] + flab_id[o] + flab_iod
```

`realwage` is already pinned by `E_plab!` from the factor-demand block. This equation **solves for the otherwise-undetermined wage-shifter residuals** — the standard ORANI/TERM decomposition-residual pattern.

### E_fhou! / E_fhou2! — household consumption closure

```julia
log(whouhtot[d] / HOUPUR_C[1,d]) == log(wlab_io[d]) + fhou[d] + houslack
log(whouhtot[d] / HOUPUR_C[1,d]) == log(wgdpexp[d]) + fhou2[d] + houslack
```

`whouhtot` is a nominal value (not an index), so it is normalized by its own benchmark `HOUPUR_C` before comparing against the bmk=1 wage/GDP indices. Under `TERM.CMF`'s closure: `swap xhouhtot = fhou` makes regional consumption follow wage income.

### E_natfhou! — national household consumption ratio

```julia
log(NatMacro[k_hou]) == natfhou + log(NatMacro[k_gdp])
```

`k_hou` = `findfirst(==("NomHou"), MAINMACROS)`, `k_gdp` = `findfirst(==("NomGDPexp"), MAINMACROS)`. Both are bmk=1 ratios from `build_macros!`, and `natfhou` is a bmk=0 log-additive shifter.

**Critical prior bug.** This was an empty stub until 2026-07-31, leaving `natfhou` a free orphan. `solve_newton!`'s defensive orphan-pinning hid this: pinning an orphan and implementing its equation produce the same answer under every closure that keeps it endogenous. But `coalprice.CMF:83 swap houslack = natfhou` makes `natfhou` EXOGENOUS — the orphan vanishes, no pin fires, and the model comes out one variable short of square. An orphan and a missing equation are only distinguishable by a closure that swaps that variable. This function must be called AFTER `build_macros!`, which is what declares `NatMacro`.

---

## Part 18 — Supporting files: dynamics, macros, CES helper, closures

### src/ces_helper.jl — the CES/CET primitive

```julia
ces(y, p, α, σ, γ)   →   vector of component demands
ces_calibrate(q, σ, y)   →   (α, γ) such that ces(y, ones, α, σ, γ) == q
```

`σ > 0` is CES (substitution), `σ < 0` is CET (transformation — used in `E_xmake!` with `-SCET[i]`). Handles the Cobb-Douglas singularity `σ = 1` analytically.

### src/build_macros!.jl — Excerpts 31-32 (regional/national macro reporting)

Declares `MainMacro[1:nmac, 1:nr+1]` and `NatMacro[1:nmac]`, where every row is an equality to an existing benchmark-1 index:

| MACRO name | Equation |
|---|---|
| `RealHou` | `xfin[1,d]` |
| `RealInv` | `xfin[2,d]` |
| `RealGov` | `xfin[3,d]` |
| `ExpVol` | `xfin[4,d]` |
| `ImpVolUsed` | `ximpused[d]` |
| `RealGDP` | `xgdpexp[d]` |
| `GDPPI` | `pgdpexp[d]` — **this is what `TERM.CMF`'s numeraire swap pins** |
| `CPI` | `phouhtot[d]` |
| `AggEmploy` | `xlab_io[d]` |
| `AggCapStock` | `xcap_i[d]` |
| `NomHou` | `whouhtot[d] / HOUPUR_C[1,d]` |
| `NomGDPexp` | `wgdpexp[d]` |
| `NomGDPinc` | `wgdpinc[d]` |

National-level `NatMacro[k]` is the `WMAIN`-weighted average of `MainMacro[k, 1:nr]`.

### src/build_dynamics!.jl — Excerpts 50-54 (dynamic extension)

Declares the dynamic variables (`delUnity`, `frnorm`, `frnorm_id`, `emptrend`, `gtrend`, `delfwage`/`delfwage_o`, `faccum`/`finv4`, `gro`/`rnorm`/`mratio`, `gretxp`/`delgret`/`delgretexp`, `delempratio`/`delwagerate`, `natggro`/`natcapstok`/`natxinvitot`/`natpinvitot`) and their defining equations as a passive satellite block. Under the base static closure all dynamic shifters are fixed at benchmark 0, and the block is benchmark-consistent.

### src/initialize_model!.jl — closure fix/free + OMIT_NAMES + numeraire

Three lists encode the closure:

- `BASE_CLOSURE_SCALARS` (13 names) and `BASE_CLOSURE_ARRAYS` (39 names) = TERM.CMF's 48 `Exogenous` entries, minus `phi` (the numeraire variable).
- `OMIT_NAMES = [ahou_s, asuppmar, bint_s, fxstocks]` — TABLO-`Omit`ted variables, not shockable.
- `ZERO_START_NAMES` — free shifters/`(change)` deltas that benchmark at 0, not 1.

The numeraire swap `phi ↔ NatMacro("GDPPI")` frees `phi` and pins the national GDP price index at 1.0.

### src/benchmark_levels.jl — start-value seeding

Returns a Dict of true benchmark flow levels for genuine quantity/value variables (`xint→PUR`, `xmake→MAKE`, `xlab→V1LAB`, etc.) — by `initialize_model!` alongside benchmark=1 for indices and benchmark=0 for shifters. Without it the start point violates value identities by ~1e10.

---

## Cross-reference: equation → TERM.TAB source

| Julia function | TERM.TAB excerpt | TERM.TAB equation name | Lines (approx) |
|---|---|---|---|
| `E_pimp!` | Excerpt 6 | `E_pimp` | 240-241 |
| `E_pbasic!` | Excerpt 6 | `E_pbasic` | 243-249 |
| `E_tuser!` | Excerpt 7 | `E_tuser` | 275-276 |
| `E_ppur!` | Excerpt 7 | `E_ppur` | 286-288 |
| `E_ppur_s!` | Excerpt 8 | `E_ppur_s` | 305-312 |
| `E_phou!` | Excerpt 8 | `E_phou` | 315 |
| `E_xint!` | Excerpt 8 | `E_xint` | 317-321 |
| `E_xhou!` | Excerpt 8 | `E_xhou` | 323-327 |
| `E_xinv!` | Excerpt 8 | `E_xinv` | 329-332 |
| `E_aint_s!` | Excerpt 9 | `E_aint_s` | 353-354 |
| `E_xint_s!` | Excerpt 9 | `E_xint_s` | 358 |
| `E_pint!` | Excerpt 9 | `E_pint` | 376 |
| `E_xlab!` | Excerpt 10 | `E_xlab` | 398-402 |
| `E_plab_o!` | Excerpt 10 | `E_plab_o` | 414 |
| `E_wlab_o!` | Excerpt 10 | `E_wlab_o` | 420 |
| `E_xlab_o!` | Excerpt 11 | `E_xlab_o` | 440-444 |
| `E_pcap!` | Excerpt 11 | `E_pcap` | 446-450 |
| `E_plnd!` | Excerpt 11 | `E_plnd` | 453 |
| `E_pprim!` | Excerpt 11 | `E_pprim` | 468 |
| `E_xprim!` | Excerpt 11 | `E_xprim` | 476 |
| `E_aprim!` | Excerpt 11 | `E_aprim` | 484-486 |
| `E_alab_o!` | Excerpt 11 | `E_alab_o` | 478-480 |
| `E_wprim!` | Excerpt 11 | `E_wprim` | 497-498 |
| `E_pvar!` | Excerpt 12 | `E_pvar` | 513-514 |
| `E_pcst!` | Excerpt 12 | `E_pcst` | 537 |
| `E_delPTX!` | Excerpt 12 | `E_delPTX` | 547-548 |
| `E_ptot!` | Excerpt 12 | `E_ptot` | 549-551 |
| `E_asub!` | Excerpt 13 | `E_asub` | 581-583 |
| `E_alux!` | Excerpt 13 | `E_alux` | 591-593 |
| `E_xsub!` | Excerpt 13 | `E_xsub` | 626-627 |
| `E_xlux!` | Excerpt 13 | `E_xlux` | 632-633 |
| `E_xhouh_s_agg!` | Excerpt 13 | – (aggregation) | – |
| `E_wlux!` | Excerpt 13 | `E_wlux` | 659-660 |
| `E_phouhtot!` | Excerpt 13 | `E_phouhtot` | 644-645 |
| `E_whouhtot!` | Excerpt 13 | `E_whouhtot` | 646-647 |
| `E_xinvi!` | Excerpt 14 | `E_xinvi` | 717-719 |
| `E_pinvest!` | Excerpt 14 | `E_pinvest` | 722 |
| `E_pinvitot!` | Excerpt 14 | `E_pinvitot` | 730-731 |
| `E_xinv_s!` | Excerpt 14 | `E_xinv_s` | 738-740 |
| `E_gret!` | Excerpt 15 | `E_gret` | 771-772 |
| `E_xinvitot!` | Excerpt 15 | `E_xinvitot` | 774 |
| `E_ggro!` | Excerpt 15 | `E_ggro` | 777 |
| `E_fgret!` | Excerpt 15 | `E_fgret` | 782 |
| `E_finv2!` | Excerpt 15 | `E_finv2` | 787 |
| `E_xgov!` | Excerpt 16 | `E_xgov` | 809-810 |
| `E_fgovtot2!` | Excerpt 16 | `E_fgovtot2` | 819 |
| `E_fgovtot3!` | Excerpt 16 | `E_fgovtot3` | 821 |
| `E_pfexp!` | Excerpt 16 | `E_pfexp` | 827 |
| `E_xexpd!` | Excerpt 16 | `E_xexpd` | 833-835 |
| `E_xstocks!` | Excerpt 16 | `E_xstocks` | 839-840 |
| `E_xint_i!` | Excerpt 17 | `E_xint_i` | 857 |
| `E_xuse!` | Excerpt 17 | `E_xuse_s` | 855-861 |
| `E_pdelivrd!` | Excerpt 19 | `E_pdelivrd` | 932-934 |
| `E_xtradmar_na!` | Excerpt 19 | `E_xtradmar` | 922-924 |
| `E_puse!` | Excerpt 20 | `E_puse` | 979-980 |
| `E_xtrad!` | Excerpt 20 | `E_xtrad` | 987-989 |
| `E_avesrctwist!` | Excerpt 20 | `E_avesrctwist` | 997-998 |
| `E_xsuppmar_p!` | Excerpt 21 | `E_xsuppmar_p` | 1033 |
| `E_psuppmar_p!` | Excerpt 21 | `E_psuppmar_p` | 1039-1041 |
| `E_xsuppmar!` | Excerpt 21 | `E_xsuppmar` | 1045-1047 |
| `E_xmake!` | Excerpt 22 | `E_xmake` | 1081-1085 |
| `E_xtotA_B!` | Excerpt 22 | `E_xtotA` | 1104 |
| `E_xcomA_B!` | Excerpt 22 | `E_xcomA` | 1106 |
| `E_pmake!` | Excerpt 22 | `E_pmake` | 1108-1109 |
| `E_xtrad_d!` | Excerpt 23 | `E_xtrad_d` | 1190 |
| `E_xtrad_r!` | Excerpt 23 | `E_xtrad_r` | 1193 |
| `E_pdomA_sum!` | Excerpt 23 | `E_pdomA` / `E_pdomB` | 1196-1203 |
| `E_pfin!` | Excerpt 24 | `E_pfin` | 1207-1208 |
| `E_xfina-d!` | Excerpt 24 | `E_xfina-d` | 1215-1224 |
| `E_wfin!` | Excerpt 24 | `E_wfin` | 1226 |
| `E_nat*pfin!` | Excerpt 24 | `E_natpfin/natxfin/natwfin` | 1230-1235 |
| `E_delTAX*!` | Excerpt 26 | `E_delTAXint…exp` | 1301-1325 |
| `E_wlnd_i!` etc. | Excerpt 27 | Factor aggregates | 1391-1411 |
| `E_plab!` … `E_rlab_io!` | Excerpt 27 | Labour aggregates | 1355-1391 |
| `E_delGDPINC*!` | Excerpt 28 | Income-side GDP | 1416-1425 |
| `E_delXGDPEXP*!` / `E_delPGDPEXP*!` | Excerpt 29 | Expenditure-side GDP | 1461-1497 |
| `E_labslack!` / `E_flab_i!` / `E_realwage!` | Excerpt 38 | Labour closure | 1955-1987 |
| `E_fhou!` / `E_fhou2!` / `E_natfhou!` | Excerpt 39 | Household closure | 2008-2053 |