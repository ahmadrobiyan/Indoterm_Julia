"""
Scenario definitions — one declarative place per source `.CMF` file.

**Why this file exists.** A scenario is a closure PLUS a shock list, and until
now only the closure had a home (`closures.jl`); the shock list was written
inline at each call site. That is how a real wrong result got produced and
believed:

> `test/verify_arclength.jl` traced `blabnat` alone, found a limit point at
> λ\\* = 0.9908083, and reported "the branch genuinely folds — `blabnat = 0.97`
> does not exist in this closure". It does exist. The run had silently omitted
> TERM.CMF's **second** shock, `delUnity = 1` (`TERM.CMF:113`). With both shocks
> applied the branch passes straight through and reaches 0.97 in 22 steps with 0
> rejections (PLAN.md §N.19e). `delUnity = 0` is not a neutral default — it hard-
> pins all 150 `xcap` cells to their benchmark values, and removing that
> adjustment margin is what manufactured the fold.

The mechanism of the error is worth naming, because it will recur otherwise:
`continuation_solve!` and `arclength_solve!` both take **one** `VariableRef` as
the continuation parameter. Their signature can only express a one-instrument
experiment, so a two-shock scenario had nowhere to put its second shock and
nothing anywhere raised. A missing shock was not a visible absence — it was an
API that could not represent the thing being left out.

The fix is to name each scenario once, transcribed from its source file with the
source line number on every entry, and to make the drivers take a `Scenario`
rather than a bare list. Then omitting a shock requires deleting a line that
carries the citation it is violating, and `describe` prints the whole thing next
to the source so a reader can diff it by eye.

**Transcription rule.** Every field carries the file and line it came from. If a
line in the source has no counterpart here, say so explicitly in the docstring
with the reason — silence is what caused the problem this file exists to prevent.
"""

"""
A named scenario: which closure, which shocks, and where both came from.

- `name`     short label, used in logs and result objects.
- `source`   the file (and line range) this was transcribed from.
- `swaps`    closure swaps, in source order (see `closures.jl`).
- `shocks`   `spec => levels_target`, where `spec` is a variable name or a
  `(name, indices...)` tuple. Targets are LEVELS values, so a GEMPACK
  `Shock x = -3` on a `bmk = 1` variable is `pct(-3)`.
- `pct_shocks` `spec => percent`, applied to each variable's OWN base value.
  Use only when the benchmark level is not 1.
- `numeraire` `:gdppi` (pin NatMacro("GDPPI"), free phi) or `:exrate` (leave phi
  exogenous). This belongs to the scenario, not to the call site: the `.CMF`
  files disagree with each other about it — TERM.CMF activates the swap at `:72`,
  coalprice.CMF has the same line commented out — and getting it wrong rebases
  every nominal result in the run against the wrong yardstick.
- `notes`    anything a reader must know before quoting a number from this run.

Construct with keywords; every field except `name` and `source` has a default, so
an empty scenario is a benchmark re-solve under the default closure.
"""
struct Scenario
    name::String
    source::String
    swaps::Vector{Tuple}
    shocks::Vector{Pair}
    pct_shocks::Vector{Pair}
    numeraire::Symbol
    notes::String
end

function Scenario(; name::AbstractString, source::AbstractString,
                    swaps = TERM_CMF_SWAPS, shocks = (), pct_shocks = (),
                    numeraire::Symbol = :gdppi, notes::AbstractString = "")
    numeraire in (:gdppi, :exrate, :cpi) ||
        error("Scenario $name: numeraire must be :gdppi, :exrate, or :cpi, got :$numeraire")
    Scenario(String(name), String(source),
             Tuple[Tuple(s) for s in swaps],
             Pair[p for p in shocks],
             Pair[p for p in pct_shocks],
             numeraire,
             String(notes))
end

"""
    describe(sc::Scenario; io = stdout)

Print a scenario in full — closure, every shock, and the notes — so a run log
records what was actually asked for. Call this from any script that runs a
scenario; it is the cheapest available guard against the omission described in
this file's header, because a missing shock becomes a missing line in the log
rather than a number nobody can reconstruct.
"""
function describe(sc::Scenario; io::IO = stdout)
    println(io, "─"^72)
    println(io, "scenario: $(sc.name)")
    println(io, "source:   $(sc.source)")
    println(io, "─"^72)
    println(io, "numeraire: " * (sc.numeraire === :gdppi ?
        ":gdppi — swap phi = NatMacro(\"GDPPI\") IS applied" :
        sc.numeraire === :cpi ?
        ":cpi — swap phi = NatMacro(\"CPI\") IS applied" :
        ":exrate — the GDPPI swap is NOT applied; phi stays exogenous"))
    println(io, "closure swaps ($(length(sc.swaps))):")
    if isempty(sc.swaps)
        println(io, "  (none — bare automatic closure; xcap/xlnd stay exogenous)")
    else
        for s in sc.swaps
            println(io, "  swap $(_side_name(s[1])) = $(_side_name(s[2]))")
        end
    end
    n = length(sc.shocks) + length(sc.pct_shocks)
    println(io, "shocks ($n):")
    n == 0 && println(io, "  (none — this is a benchmark re-solve)")
    for (spec, tgt) in sc.shocks
        println(io, "  $(rpad(_spec_label(spec), 22)) → $tgt   (levels)")
    end
    for (spec, p) in sc.pct_shocks
        println(io, "  $(rpad(_spec_label(spec), 22)) → $(p)% of its own base")
    end
    isempty(sc.notes) || (println(io, "notes:"); println(io, sc.notes))
    println(io, "─"^72)
    flush(io)
end

_spec_label(spec) = spec isa Tuple ? "$(spec[1])$(spec[2:end])" : String(spec)

"""
    run_model!(agg, params, sc::Scenario; kwargs...) -> ScenarioResult

Run a named [`Scenario`](@ref): its closure, its complete shock list, its notes.
This is the form every script should use — the keyword form is the primitive, and
spelling a shock list out at a call site is exactly how `TERM.CMF:113` went
missing (see this file's header).

`describe(sc)` is printed first, so the run log always records the full scenario
next to its source citation. Keyword arguments are forwarded to the primitive for
step control (`h0`, `tol`, `report`, …); `name`, `swaps`, `shocks` and
`pct_shocks` come from the scenario and passing them here is an error rather than
a silent override — a scenario that can be edited at the call site is not a
single declarative place.
"""
function run_model!(agg, params::Dict{String,Any}, sc::Scenario;
                    verbose::Bool = true, kwargs...)
    for k in (:name, :swaps, :shocks, :pct_shocks, :numeraire)
        haskey(kwargs, k) && error(
            "run_model!(…, sc::Scenario) does not accept `$k` — it is part of the " *
            "scenario. Overriding it at the call site defeats the point of naming " *
            "the scenario: edit `$(sc.name)` in src/scenarios.jl, or build an " *
            "explicitly-named Scenario for the variant you want.")
    end
    verbose && describe(sc)
    run_model!(agg, params; name=sc.name, swaps=sc.swaps,
               shocks=sc.shocks, pct_shocks=sc.pct_shocks,
               numeraire=sc.numeraire, verbose=verbose, kwargs...)
end

"""
TERM.CMF's own reference scenario, transcribed from `origin/TERM.CMF`.

**Shocks — the complete list from the source file.** TERM.CMF contains exactly
two `Shock` statements and BOTH are here:

- `:112` `Shock blabnat = -3;`  — the comment on that line reads
  `! 3% increase labour productivity`. `blabnat` drives `alab_o`, *labour-
  augmenting technical change* (`TERM.TAB:460`, `:468`), so this is a
  productivity **gain** and it is expansionary. It is NOT a labour-supply cut;
  describing it that way once produced a test that asserted real GDP must fall.
- `:113` `Shock delunity = 1;`  — `! move forward 1 year`. `delUnity` is a
  GEMPACK `(change)` dummy declared "always exogenously set to one"
  (`TERM.CMF:36`). Its unambiguous effect here is to release the capital
  block: at `delUnity = 0` the accumulation equations weld all 150 `xcap` cells
  to their benchmark values, which is precisely the pin that produced the
  spurious fold at λ\\* = 0.9908083.

Omitting either one changes the experiment. The second was omitted once and the
result was reported as a property of the model.

Closure is `TERM_CMF_SWAPS` — the six active swaps at `:79, :82, :86, :88, :89,
:110`. The seventh, `:72 swap phi = Natmacro("GDPPI")`, is applied inside
`initialize_model!` as the numeraire choice and so is deliberately absent from
that list; see `closures.jl`.
"""
const TERM_CMF_REFERENCE = Scenario(
    name   = "TERM.CMF reference — blabnat = -3%, delUnity = 1",
    source = "origin/TERM.CMF:72–113",
    swaps  = TERM_CMF_SWAPS,
    shocks = [
        "blabnat"  => pct(-3),   # TERM.CMF:112  ! 3% increase labour productivity
        "delUnity" => 1.0,       # TERM.CMF:113  ! move forward 1 year
    ],
    notes = """
  `blabnat` is labour-AUGMENTING TECHNICAL CHANGE, not labour supply, so -3 is a
  3% productivity GAIN and real GDP must RISE. `delUnity = 1` is load-bearing:
  at 0 it hard-pins all 150 xcap cells and the branch folds at λ ≈ 0.9908083.
  Quote no number from this scenario without naming the closure — the base static
  closure reads the same cell 43× differently.""",
)

"""
The same productivity shock with `delUnity` left at 0 — i.e. the INCOMPLETE
scenario that produced the spurious fold.

This is kept deliberately, and named for what it is, so the comparison can be
re-run on demand and so nobody reconstructs it by accident thinking it is the
reference. It is not a variant of TERM.CMF; it is TERM.CMF with a line missing.
"""
const TERM_CMF_NO_DELUNITY = Scenario(
    name   = "INCOMPLETE — blabnat = -3% with delUnity = 0",
    source = "origin/TERM.CMF:112 only (:113 deliberately omitted)",
    swaps  = TERM_CMF_SWAPS,
    shocks = ["blabnat" => pct(-3)],
    notes = """
  ⚠️ NOT the reference scenario. TERM.CMF:113 `Shock delunity=1` is missing on
  purpose. Expect a limit point at λ ≈ 0.9908083; that fold is an artifact of the
  omission, not a property of the model. Use TERM_CMF_REFERENCE for any result.""",
)

"""
`origin/coalprice.CMF`, transcribed in full. **This is the only scenario in the
project with a published external result to check against** — the draft report's
Tables 2 and 3, recovered into `test/reference/draftreport_coalprice.md`. Every
other verification we have is internal (the model against itself); this one is the
model against a number someone else computed with GEMPACK.

**Shock — the complete list.** coalprice.CMF contains exactly three `Shock` lines
and only ONE is active:

- `:116` `shock fpexp_d("coal") = 50;` — a 50% rise in the world price of coal
  exports. Coal is sector 5 (`aggregation_data.jl:42`).
- `:113` `!Shock blabnat = -3;`   — COMMENTED OUT.
- `:114` `!Shock delunity=1;`     — COMMENTED OUT.

The last two are the two shocks that DO fire in `TERM_CMF_REFERENCE`, so this is
a pure comparative static: no capital accumulation, no productivity trend. Do not
reach for `TERM_CMF_SWAPS` or `delUnity = 1` out of habit here.

**Why `logpct(50)` and not `pct(50)`.** `fpexp_d` is declared unbounded
(`build_model!.jl:214`), so it is benchmark 0, and `E_xexpd` uses it additively
beside `pfexp = log(ppur) − log(phi)`. The term is in natural logs, so +50% is
`log(1.5)`. `pct(50) = 1.5` would apply `e^1.5 − 1` ≈ **+348%** and still solve
without complaint. See `logpct`.

**Numeraire is `:exrate`.** `:72 !  swap phi = Natmacro("GDPPI");` is commented
out in this file. The report's CPI column is measured against a fixed exchange
rate, not a fixed GDP price index.

**Before quoting any comparison**, read the caveats in the reference file: the
report is an 8-step Euler approximation and this port solves the exact levels
system, so agreement is expected to be close but not exact; and Table 3's
Imports/Exports columns could not be recovered from the PDF and must not be
reconstructed.
"""
const COALPRICE_REFERENCE = Scenario(
    name   = "coalprice.CMF — world coal export price +50%",
    source = "origin/coalprice.CMF:79–116",
    swaps  = COALPRICE_SWAPS,
    shocks = [("fpexp_d", 5) => logpct(50)],   # coalprice.CMF:116, coal = sector 5
    numeraire = :exrate,                        # coalprice.CMF:72 is COMMENTED OUT
    notes = """
  Comparative static: blabnat and delUnity are both commented out in the source,
  unlike TERM_CMF_REFERENCE. Numeraire is the exchange rate, NOT GDPPI.
  fpexp_d is a LOG-space additive shifter — the shock is log(1.5), not 1.5.
  Reference results: test/reference/draftreport_coalprice.md (national row of
  Table 2: Real GDP -0.09, CPI 1.54, Export vol -3.41, Employment -0.27).""",
)

# ═══════════════════════════════════════════════════════════════════════════
# Industrial / GVC policy scenarios — "hilirisasi" (downstreaming)
#
# These are NOT transcribed from a `.CMF` file. They are constructed here, and
# `source` says so. The transcription rule in this file's header applies to
# scenarios that have a GEMPACK original; for a constructed scenario the
# equivalent obligation is to state the POLICY MAPPING — which real instrument
# each shock is standing in for, and what that mapping does not capture. That
# is what the `notes` field carries below, and no number from these runs should
# be quoted without it.
#
# Closure is `COALPRICE_SWAPS` with `numeraire = :exrate` throughout — the only
# closure in this project with an external reference result, and a genuine
# comparative static (no capital accumulation, no productivity trend). Capital
# is therefore EXOGENOUS, which is what makes an explicit `xcap` injection the
# right way to represent a smelter build: the policy puts the capital there, the
# model does not have to decide to.
#
# Sector indices (aggregation_data.jl): 7 = OthMining, 14 = BasMetals.
# Region indices (REG6): 1 Sumatra, 2 Java, 3 Kalimantan, 4 Sulawesi,
# 5 BaliNusa, 6 MalukuPapua.
# ═══════════════════════════════════════════════════════════════════════════

"""
Ore export restriction — the *ban* leg of hilirisasi, with no downstream build.

**Instrument mapping.** `fqexp_d[7] = logpct(-30)` shifts the foreign demand
schedule for `OthMining` (metal ores and concentrates) down 30%. `fqexp_d` is
the national, log-additive quantity shifter in `E_xexpd`
(`build_equations.jl:714`), so a downward shift means the same domestic price
now clears a smaller export volume — ore that used to leave the country is
pushed back onto the domestic market and the domestic ore price falls. That
price fall IS the mechanism hilirisasi relies on: cheap feedstock for smelters.

**What this mapping does NOT capture.** A real export ban is a quantity
prohibition, and a real export tax raises revenue. A demand-schedule shift
raises neither revenue nor a licence rent, so the fiscal leg of the policy is
absent and the welfare result here is, if anything, generous to the ban's
critics on revenue and silent on rent capture. It also treats the restriction
as permanent and fully anticipated.

**Aggregation caveat, and it is the big one.** In the 2016 benchmark the
`OthMining` export base is concentrated in **MalukuPapua (35,093 of 48,733 =
72%)** and **BaliNusa (8,803 = 18%)** — i.e. copper concentrate from Papua and
Sumbawa, not Sulawesi nickel, which by 2016 was already being smelted in-region
(Sulawesi `BasMetals` output 19,290, exports 11,145). So this scenario is an
experiment about the COPPER concentrate chain in the two poorest island groups.
Reading it as "the 2020 nickel ore ban" is a category error the sector labels
invite and the data refuses.
"""
const HILIRISASI_BAN = Scenario(
    name   = "hilirisasi — ore export restriction (OthMining demand −30%)",
    source = "constructed; analysis/gvc_policy.jl",
    swaps  = COALPRICE_SWAPS,
    shocks = [("fqexp_d", 7) => logpct(-30)],
    numeraire = :exrate,
    notes = """
  fqexp_d is LOG-additive — the shock is log(0.7), not 0.7. No export-tax
  revenue and no licence rent is modelled. The OthMining export base is 72%
  MalukuPapua / 18% BaliNusa: this is the copper-concentrate chain, NOT the
  2020 nickel ban.""",
)

"""
Smelter build sited in **Java** — `xcap[14, 2]` +20%, no export restriction.

`xcap` is exogenous under `COALPRICE_SWAPS` (no `xcap = faccum` swap), so this
is a direct capital injection into `BasMetals` in Java, standing in for a
greenfield smelter financed from outside the model. `pct_shocks` is used, not
`shocks`, because `xcap`'s benchmark is the CAP data value, not 1.

Java's BasMetals capital base is 41,742, so +20% = **8,348 units** injected.
Compare against `HILIRISASI_SMELTER_EAST`, which injects only 1,048, and
normalise before drawing a siting conclusion — see `analysis/gvc_policy.jl`.
"""
const HILIRISASI_SMELTER_JAVA = Scenario(
    name   = "hilirisasi — smelter capital +20% in Java BasMetals",
    source = "constructed; analysis/gvc_policy.jl",
    swaps  = COALPRICE_SWAPS,
    pct_shocks = [("xcap", 14, 2) => 20.0],
    numeraire = :exrate,
    notes = """
  Capital appears from outside the model — no crowding-out of the investment
  that financed it, and no construction phase. Absolute injection 8,348 units
  (20% of a 41,742 base), ~8x the eastern variant: normalise per unit of
  capital before comparing sites.""",
)

"""
The same build sited **at the resource** — `xcap[14, 4]` (Sulawesi) and
`xcap[14, 6]` (MalukuPapua) each +20%.

Absolute injection is 20% × (5,110 + 129) = **1,048 units**, about one eighth of
the Java variant, because eastern BasMetals capital is tiny in the benchmark.
That asymmetry is the point of the pair, not a flaw in it: it is what "siting
downstream capacity where the ore is" actually costs relative to what it moves.
"""
const HILIRISASI_SMELTER_EAST = Scenario(
    name   = "hilirisasi — smelter capital +20% in Sulawesi + MalukuPapua BasMetals",
    source = "constructed; analysis/gvc_policy.jl",
    swaps  = COALPRICE_SWAPS,
    pct_shocks = [("xcap", 14, 4) => 20.0, ("xcap", 14, 6) => 20.0],
    numeraire = :exrate,
    notes = """
  Absolute injection 1,048 units vs 8,348 for the Java variant — a like-for-like
  siting comparison requires normalising by capital injected, not comparing the
  raw percentage results.""",
)

"""
The full policy: restriction **and** an at-the-resource build, run jointly so the
interaction can be recovered as `S4 - S1 - S3` against the two single-leg runs.

The interaction is the whole empirical question. Hilirisasi's case is that the
ban and the build are complements — the ban is what makes the smelter viable. If
the joint result is no better than the sum of the parts, the ban is buying
nothing that the build does not already deliver on its own.
"""
const HILIRISASI_FULL = Scenario(
    name   = "hilirisasi — ore restriction + eastern smelter build",
    source = "constructed; analysis/gvc_policy.jl",
    swaps  = COALPRICE_SWAPS,
    shocks = [("fqexp_d", 7) => logpct(-30)],
    pct_shocks = [("xcap", 14, 4) => 20.0, ("xcap", 14, 6) => 20.0],
    numeraire = :exrate,
    notes = """
  Three shocks, so `run_model!`'s arclength fallback is unavailable (it traces a
  single VariableRef). If the homotopy stalls, stage the legs rather than
  reporting a partial t as a result.""",
)
