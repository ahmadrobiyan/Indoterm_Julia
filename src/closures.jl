"""
Closure swaps — GEMPACK `swap` applied to the levels port.

`initialize_model!` fixes TERM.CMF's 48 `Exogenous` names. That list is TABmate's
*automatic* closure: the state before any modelling choice is made. Every closure
file in the source then applies swaps on top of it, and **the port was applying
only one of them** (`phi = NatMacro("GDPPI")`, the numeraire). Running an
economy-wide shock against the bare automatic closure is what produced the fold
at `blabnat = 0.99975` documented in PLAN.md §J: with `xcap` AND `xlnd` frozen in
all 25 sectors, fixed-factor *prices* carry the whole adjustment, at ~200×
amplification in the cell with the smallest fixed-factor share.

TERM.CMF's own scenario is `Shock blabnat = -3` — 120× beyond where that closure
folds. It reaches it because `TERM.CMF:88` swaps `xcap` endogenous.

**`swap` is symmetric.** GEMPACK requires exactly one of the pair to be
exogenous; the swap exchanges their status. Do NOT read the left-hand side as
"the one that becomes endogenous": `TERM.CMF:110 swap delfwage = flabsup_id`
contradicts its own `!old exog = new exog` header (line 74), because `delfwage`
is *not* in the Exogenous list while `flabsup_id` is. Implementing the symmetric
rule gets that pair right and every other pair too.

Fixing is done at each variable's CURRENT start value, which under
`initialize_model!` is its benchmark seed. That is what keeps the benchmark an
exact solution across a swap — the regression gate every closure change must
pass before its results mean anything.
"""

"""
TERM.CMF's active (uncommented) swaps, minus the numeraire swap that
`initialize_model!` already applies internally.

Source lines, in file order:
- `:79`  `xhouhtot = fhou`        — regional consumption follows wage income
- `:82`  `houslack = shrBoTnom2`  — fix nominal BOT/GDP
- `:86`  `fgovtot = fgovtot3`     — regional real government follows regional real GDP
- `:88`  `xcap = faccum`          — **switch on capital accumulation** (frees `xcap`)
- `:89`  `finv1 = finv4`          — dynamic investment rule
- `:110` `delfwage = flabsup_id`  — switch on the national labour adjustment mechanism

`xlnd` has no swap in ANY source closure file, so holding land fixed is faithful;
holding capital fixed is not.
"""
const TERM_CMF_SWAPS = [
    ("xhouhtot", "fhou"),
    ("houslack", "shrBoTnom2"),
    ("fgovtot",  "fgovtot3"),
    ("xcap",     "faccum"),
    ("finv1",    "finv4"),
    ("delfwage", "flabsup_id"),
]

"""
`termsr.cmf`'s short-run swaps (`:66–69`), for reference and for a genuine
short-run closure. Note it keeps `xcap`/`xlnd` exogenous — appropriate for its
own scenario, a targeted investment shock in one cell, NOT for an economy-wide
labour shock.
"""
const TERM_SR_SWAPS = [
    ("invslack", ("NatMacro", "RealInv")),   # termsr.cmf:66 — do NOT hold national investment fixed
    ("xhouhtot", "fhou"),
    ("flab_i",   "flabsupA"),
    ("labslack", "flabsup_id"),
]

"""
`origin/coalprice.CMF`'s swaps — the closure that produced the published draft
report, and therefore the only closure in this project with an EXTERNAL reference
result to compare against (`test/reference/draftreport_coalprice.md`).

Source lines, in file order:
- `:79`  `xhouhtot = fhou`         — regional consumption follows wage income
- `:83`  `houslack = natfhou`      — fix the national propensity to consume
- `:108` `flab_i = flabsupA`       — switch off regional wage differentials
- `:109` `labslack = flabsup_id`   — switch off the national labour supply mechanism

TWO THINGS THIS CLOSURE DOES **NOT** DO, both by commented-out lines in the
source, and both of which change results if we get them wrong:

1. `:72` `!  swap phi = Natmacro("GDPPI");` is COMMENTED OUT. The numeraire is the
   exchange rate `phi`, NOT the GDP price index. `initialize_model!` applies the
   GDPPI swap unconditionally, so this closure must be built with
   `numeraire=:exrate` or every nominal result — the CPI column above all — is
   measured against the wrong yardstick.
2. `:113–114` `!Shock blabnat = -3;` and `!Shock delunity=1;` are BOTH commented
   out. This is a pure comparative-static run with no capital accumulation and no
   labour-productivity trend — unlike the TERM.CMF reference scenario, which has
   both. Do not carry `TERM_CMF_SWAPS` habits into it.

The shock itself is `:116` `shock fpexp_d("coal") = 50;` and lives in
`scenarios.jl`, not here — a closure says which variables are free, not what
happens to them.
"""
const COALPRICE_SWAPS = [
    ("xhouhtot", "fhou"),
    ("houslack", "natfhou"),
    ("flab_i",   "flabsupA"),
    ("labslack", "flabsup_id"),
]

"""
`TERM_LR_SWAPS_MATCHED` — same capital/investment/consumption/government swaps as
`TERM_CMF_SWAPS`, but with `TERM.CMF`'s own labour swap (`:110 delfwage =
flabsup_id`, a third variant that is neither textbook long-run nor short-run —
see the module docstring) replaced by `termlr.cmf`'s own textbook long-run
labour swap:

- `swap flabsup_id = xlab_id;` (`origin/termlr.cmf`, "Swaps for long-run
  closure") — fixes national employment; the real wage is what carries the
  adjustment. `flabsup_id` is exogenous in the automatic closure and `xlab_id`
  is endogenous, so this swap applies cleanly with no other changes needed.

Exists ONLY to give V6 (`test/verify_closure_ordering.jl`) a long-run leg that
is genuinely matched against `TERM_SR_SWAPS` on the labour axis (short-run:
`flab_i=flabsupA` + `labslack=flabsup_id`, regional wage differentials AND the
national labour-supply mechanism both switched off) as well as the capital
axis, per `origin/TERM.TAB:1971-2036` Excerpt 38's own definition of what
short-run/long-run means. Do NOT substitute this for `TERM_CMF_SWAPS` in any
gate that claims to reproduce `TERM.CMF`'s own scenario (V1, the coalprice
reference, etc.) — `TERM_CMF_SWAPS` stays the faithful implementation of the
source file; this is a deliberately different (more textbook) closure built to
make V6's comparison valid.
"""
const TERM_LR_SWAPS_MATCHED = [
    ("xhouhtot",   "fhou"),
    ("houslack",   "shrBoTnom2"),
    ("fgovtot",    "fgovtot3"),
    ("xcap",       "faccum"),
    ("finv1",      "finv4"),
    ("flabsup_id", "xlab_id"),
]

_swap_refs(v::VariableRef) = [v]
_swap_refs(v::AbstractArray) = vec(collect(v))

# A swap side is either a variable name, or a ("NatMacro", "GDPPI")-style pair
# naming one element of an indexed macro vector — TERM.CMF and termsr.cmf both
# swap against single macros (`NatMacro("RealInv")`, `NatMacro("GDPPI")`), which a
# bare name lookup cannot express.
_side_name(s::AbstractString) = String(s)
_side_name(s::Tuple) = "$(s[1])(\"$(s[2])\")"

function _side_refs(vars::Dict{String,Any}, s::AbstractString)
    haskey(vars, s) || error("swap side \"$s\" is not a declared variable")
    return _swap_refs(vars[s])
end

function _side_refs(vars::Dict{String,Any}, s::Tuple)
    nm, key = String(s[1]), String(s[2])
    haskey(vars, nm) || error("swap side $(_side_name(s)): \"$nm\" is not declared")
    labels = nm == "NatMacro" ? MAINMACROS :
             nm == "SelMacro" ? SELMACROS  :
             error("swap side $(_side_name(s)): don't know the label set for \"$nm\"")
    k = findfirst(==(key), labels)
    k === nothing && error("swap side $(_side_name(s)): \"$key\" is not in $nm's label set")
    return [vars[nm][k]]
end

"""
    apply_swaps!(m, vars, swaps; verbose=true) -> Int

Apply GEMPACK `swap` statements to an already-closed model, and return the number
applied.

Each entry is a `(name_a, name_b)` pair. Exactly one side must currently be
exogenous, elementwise; that side is freed (started at its fixed value) and the
other is fixed at its current start value. A pair whose sides have different
lengths, or where both/neither side is exogenous, raises — silently skipping a
swap would leave the system non-square and turn a closure bug into a mysterious
solver failure.
"""
function apply_swaps!(m::JuMP.Model, vars::Dict{String,Any},
                      swaps::AbstractVector; verbose::Bool=true)
    log(msg) = verbose && println(msg)
    napplied = 0

    for (sa, sb) in swaps
        na_ = _side_name(sa); nb_ = _side_name(sb)
        A = _side_refs(vars, sa); B = _side_refs(vars, sb)
        length(A) == length(B) || error(
            "swap $na_ = $nb_: shape mismatch — $(length(A)) vs $(length(B)) elements. " *
            "A GEMPACK swap pairs variables over the same set.")

        fa = count(JuMP.is_fixed, A); fb = count(JuMP.is_fixed, B)
        # Demand a clean all-or-nothing status on each side. A partially fixed
        # side means an earlier swap half-applied, and continuing would silently
        # unbalance the square system.
        (fa == 0 || fa == length(A)) || error(
            "swap $na_ = $nb_: \"$na_\" is only partly exogenous ($fa/$(length(A)) fixed)")
        (fb == 0 || fb == length(B)) || error(
            "swap $na_ = $nb_: \"$nb_\" is only partly exogenous ($fb/$(length(B)) fixed)")

        exog, endo, exname, enname =
            if fa == length(A) && fb == 0
                A, B, na_, nb_
            elseif fb == length(B) && fa == 0
                B, A, nb_, na_
            else
                error("swap $na_ = $nb_: exactly one side must be exogenous, but " *
                      "\"$na_\" is $(fa == 0 ? "endogenous" : "exogenous") and " *
                      "\"$nb_\" is $(fb == 0 ? "endogenous" : "exogenous")")
            end

        for (a, b) in zip(exog, endo)
            # Free the old exogenous variable, starting it where it was pinned.
            fv = JuMP.fix_value(a)
            JuMP.unfix(a)
            JuMP.set_start_value(a, fv)
            # Pin the new exogenous variable at its benchmark seed, so the
            # benchmark remains an exact solution of the swapped system.
            sv = JuMP.start_value(b)
            JuMP.fix(b, sv === nothing ? 0.0 : sv; force=true)
        end
        napplied += 1
        log("  swap $(rpad(exname, 12)) → endogenous,  $(rpad(enname, 12)) → exogenous " *
            "($(length(exog)) element$(length(exog) == 1 ? "" : "s"))")
    end

    log("  applied $napplied swap(s)")
    return napplied
end
