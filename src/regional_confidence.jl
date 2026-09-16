"""
Regional confidence warning — a reusable early-warning check for anyone reading region-level
results out of a solved `ScenarioResult`, for ANY scenario (not just the ones V4 happened to
test).

**What this is answering.** VV_PLAN.md §V4 ran a 6-region solve and a 12-region solve of the
SAME shock, aggregated the 12-region result down to 6, and compared. After five diagnostic
rounds (ruling out aggregation-composition bugs, benchmark-base thinness, and region-weight
imbalance as the driver — see VV_PLAN.md), the one clean empirical pattern that survived every
round was: **every confirmed sign disagreement between the two resolutions had a region-level
deviation-from-benchmark under 0.3%** (`v4_run12_gate_b.log`, the 29 region-aggregate flips
under option (b), largest was `xsuppmar(6,3,3)` at 0.298%). Nothing with a materially large
regional move (1%+) ever flipped sign between resolutions. This isn't a coincidence: a region's
merged response can only straddle zero (and so be fragile to how it's resolved/weighted) if it's
already close to zero — a large, one-directional response can't flip from a small change in how
it's aggregated.

**What this function does — and does NOT do.** It does NOT rerun the model at a second
resolution (that costs a full second solve, ~15-30 min). It applies the empirical threshold
above as a cheap, always-available heuristic: any region-level flow deviation under
`REGIONAL_FRAGILITY_THRESHOLD` is flagged as a CANDIDATE for resolution-sensitivity — a
near-zero result that should be read as "no material regional effect" rather than cited by
sign/direction, unless independently re-checked (e.g. by an actual 12-region rerun of that
specific scenario, following `test/verify_aggregation_consistency.jl`'s method). This will flag
some cells that are in fact perfectly resolution-robust (false positives are the intended,
safe failure mode for a warning); it is not a proof of fragility, only a screen for it.

**Threshold choice.** 0.5% — set with margin above the largest confirmed V4 flip (0.298%), so
every historically-observed flip is caught, while staying well below the smallest observed
"real" (sign-stable) regional move (0.6%+ in the cell-level material-flip list).

**Noise floor.** A zero-shock (or barely-shocked) re-solve leaves Newton/homotopy residual
noise of order 1e-9 in every variable, which reads as a "deviation" of the same tiny order —
technically nonzero, not economically anything. `noise_floor` (default `1e-6`, two orders above
the ~1e-9 solver residual floor documented throughout this project) excludes these before the
fragility threshold is applied, so a cell only gets flagged if the shock actually moved it a
*little* (not moved it by nothing and been read as noise). Found by the first smoke test of this
function: an un-shocked benchmark re-solve flagged `xtradmar` at "-0.0%" in every region before
this floor was added — solver noise, not a result.
"""
const REGIONAL_FRAGILITY_THRESHOLD = 0.005
const REGIONAL_NOISE_FLOOR = 1e-6

"""
    regional_confidence_report(vals, bmk, nr=6; threshold=REGIONAL_FRAGILITY_THRESHOLD,
                                noise_floor=REGIONAL_NOISE_FLOOR, regnames=REG6)

For every flow variable present in both `vals` (a solved `ScenarioResult.values`) and `bmk`
(from `benchmark_levels`), collapse it to its region axis/axes (same rule as V4's
`regional_totals`: sum over every non-region axis, then over any extra region axis beyond the
first — e.g. `xtrad[c,s,r,d]` reduces to one total per origin region `r`), compute
`%dev = value/base - 1` per region, and flag every cell with `noise_floor <= |dev| < threshold`.

Returns a `Vector{NamedTuple}` sorted by ascending `|dev|` (closest to zero, most fragile,
first). Each entry has fields `variable`, `region`, `dev`.
"""
function regional_confidence_report(vals::Dict, bmk::Dict, nr::Int = 6;
                                     threshold::Float64 = REGIONAL_FRAGILITY_THRESHOLD,
                                     noise_floor::Float64 = REGIONAL_NOISE_FLOOR,
                                     regnames::Vector{String} = REG6)
    flagged = NamedTuple[]
    for (nm, b) in bmk
        b isa AbstractArray || continue
        v = get(vals, nm, nothing)
        (v isa AbstractArray && size(v) == size(b)) || continue

        rax = [d for d in 1:ndims(v) if size(v, d) == nr]
        isempty(rax) && continue
        otherax = [d for d in 1:ndims(v) if !(d in rax)]
        vr = isempty(otherax) ? v : dropdims(sum(v; dims = otherax); dims = Tuple(otherax))
        br = isempty(otherax) ? b : dropdims(sum(b; dims = otherax); dims = Tuple(otherax))
        while ndims(vr) > 1
            vr = dropdims(sum(vr; dims = 2); dims = 2)
            br = dropdims(sum(br; dims = 2); dims = 2)
        end

        for i in eachindex(vr)
            br[i] > 1e-8 || continue
            dev = vr[i] / br[i] - 1.0
            if noise_floor <= abs(dev) < threshold
                label = (nr == length(regnames) && i <= length(regnames)) ? regnames[i] : "region$i"
                push!(flagged, (variable = nm, region = label, dev = dev))
            end
        end
    end
    sort!(flagged; by = x -> abs(x.dev))
    flagged
end

"""
    print_regional_confidence_report(vals, bmk, nr=6; kwargs...)

Prints `regional_confidence_report`'s result as a warning list, meant to be run right after any
`run_model!` call whose output will be used to report region-level findings. Returns the same
flagged vector (also usable programmatically, e.g. to gate a report-generation script).
"""
function print_regional_confidence_report(vals::Dict, bmk::Dict, nr::Int = 6; kwargs...)
    flagged = regional_confidence_report(vals, bmk, nr; kwargs...)
    if isempty(flagged)
        println("✅ regional confidence check: no region-level flow sits in the resolution-fragile near-zero band (< $(round(REGIONAL_FRAGILITY_THRESHOLD*100, digits=2))%).")
        return flagged
    end
    println("⚠  REGIONAL CONFIDENCE WARNING — $(length(flagged)) region-level result(s) under $(round(REGIONAL_FRAGILITY_THRESHOLD*100, digits=2))% deviation from benchmark.")
    println("   Per V4 (VV_PLAN.md): every confirmed resolution-dependent sign flip found across five")
    println("   diagnostic rounds had |dev| < 0.3%. A value this small should be read as \"no material")
    println("   regional effect\" — do NOT cite its sign/direction as a finding without an independent")
    println("   check (e.g. a 12-region rerun of this scenario via test/verify_aggregation_consistency.jl).")
    for f in flagged
        println("     $(f.variable) [$(f.region)]: $(round(f.dev * 100, digits = 4))%")
    end
    flagged
end
