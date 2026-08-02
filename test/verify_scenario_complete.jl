"""
Does `TERM_CMF_REFERENCE` still contain every `Shock` line in `origin/TERM.CMF`?

This test exists because the answer was once **no**, and nothing noticed.
`test/verify_arclength.jl` traced `blabnat` alone, hit a limit point at
λ\\* = 0.9908083, and reported "the branch genuinely folds — `blabnat = 0.97`
does not exist in this closure". It does exist. The run had omitted
`TERM.CMF:113 Shock delunity=1`, whose effect is to release a hard pin on all 150
`xcap` cells; with it, the branch reaches the target in 22 steps with 0
rejections (PLAN.md §N.19e). A wrong claim about the model was published because
a scenario was one line short of its source.

A docstring asserting "these are the only two Shock statements" is exactly the
kind of claim that rots the moment someone edits the source. So this test does
not trust the docstring — it **parses `origin/TERM.CMF` and counts**, and fails
if the file and `src/scenarios.jl` disagree about how many shocks the scenario
has or about which variables they name.

Note this is a completeness check, not a value check: it verifies that every
shocked variable in the source has a counterpart in the Scenario, not that the
targets were transcribed with the right numbers. The magnitudes are asserted by
the run gates (`test/run_scenario_6reg.jl`).

Run:  julia --project=IndotermJulia IndotermJulia/test/verify_scenario_complete.jl
"""

using IndotermJulia

const CMF = joinpath(@__DIR__, "..", "origin", "TERM.CMF")

# GEMPACK is case-insensitive: the file has `Shock blabnat` and `Shock delunity`
# on adjacent lines with different capitalisation, which is part of why a
# case-sensitive skim missed one. Strip comments first — `!...` is a GEMPACK
# comment and a commented-out shock is not a shock.
function source_shocks(path)
    found = Tuple{Int,String}[]
    for (n, raw) in enumerate(eachline(path))
        line = split(raw, '!')[1]                    # drop trailing comment
        m = match(r"^\s*shock\s+([A-Za-z_][A-Za-z0-9_]*)"i, line)
        m === nothing || push!(found, (n, String(m.captures[1])))
    end
    found
end

isfile(CMF) || error("cannot find $CMF — this test reads the shipped source, not a copy")

src = source_shocks(CMF)
println("`Shock` statements found in origin/TERM.CMF:")
for (n, v) in src
    println("  :$n  $v")
end

sc = TERM_CMF_REFERENCE
declared = String[]
for (spec, _) in sc.shocks
    push!(declared, spec isa Tuple ? String(spec[1]) : String(spec))
end
for (spec, _) in sc.pct_shocks
    push!(declared, spec isa Tuple ? String(spec[1]) : String(spec))
end

println("\nvariables shocked by TERM_CMF_REFERENCE:")
foreach(v -> println("  $v"), declared)

# Compare case-insensitively: TERM.CMF:113 writes `delunity`, the model variable
# is `delUnity`. GEMPACK does not distinguish them and neither should this check.
srcset  = Set(lowercase(v) for (_, v) in src)
declset = Set(lowercase(v) for v in declared)

missing_here  = sort(collect(setdiff(srcset, declset)))
missing_there = sort(collect(setdiff(declset, srcset)))

println()
if !isempty(missing_here)
    error("""
    TERM_CMF_REFERENCE is MISSING shock(s) that origin/TERM.CMF applies: $missing_here

    This is the exact failure that produced the spurious fold at λ* = 0.9908083.
    A scenario short one shock is a DIFFERENT experiment, and every number from it
    describes that different experiment. Add the missing shock to src/scenarios.jl
    with its source line number, or — if it is being left out deliberately — build a
    separate Scenario named INCOMPLETE, as TERM_CMF_NO_DELUNITY is.""")
end
if !isempty(missing_there)
    error("""
    TERM_CMF_REFERENCE shocks variable(s) that origin/TERM.CMF does not: $missing_there

    The reference scenario must be a transcription, not an extension. If this is a
    deliberate variant it needs its own named Scenario, so nobody quotes its results
    as the authors' reference simulation.""")
end

@assert length(src) == length(declared) """
    count mismatch: origin/TERM.CMF has $(length(src)) Shock statement(s), \
    TERM_CMF_REFERENCE declares $(length(declared)). The variable NAMES agree, so \
    one of them shocks the same variable twice — check for a duplicate."""

println("✅ TERM_CMF_REFERENCE matches origin/TERM.CMF: $(length(src)) shock(s), " *
        "same variables, parsed from the source file.")

# The incomplete variant must stay incomplete, and must stay labelled. If someone
# "fixes" it by adding delUnity, the §N.19c comparison silently stops reproducing
# the fold and the historical record becomes unverifiable.
@assert length(TERM_CMF_NO_DELUNITY.shocks) == 1 """
    TERM_CMF_NO_DELUNITY is supposed to be the INCOMPLETE scenario (blabnat only).
    It now has $(length(TERM_CMF_NO_DELUNITY.shocks)) shocks. If it was 'fixed', the
    §N.19c fold comparison can no longer be reproduced — revert it; it is a fixture,
    not a scenario anyone should run for results."""
@assert occursin("INCOMPLETE", TERM_CMF_NO_DELUNITY.name) """
    TERM_CMF_NO_DELUNITY lost its INCOMPLETE label. That label is the only thing
    stopping its results being read as the authors' reference simulation."""

println("✅ TERM_CMF_NO_DELUNITY is still the labelled INCOMPLETE fixture (1 shock).")
println("\nDone.")
