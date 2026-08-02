"""
Phase-1 region collapse: reduce the 25×34 aggregated model (`agg` dict produced
by `aggregate_model!`) to 25×`n_new` by summing/averaging over region axes.

Motivation (see PLAN.md, Gate 5c contingency): every stage downstream of
`aggregate_model!` derives `nr` from `size(agg["MAKE"], 3)`, so a single pass
that collapses the region axes of the `agg` dict — followed by re-running the
existing `prepare_parameters!` — rescales the *entire* model (parameter
derivation, equation building, initialization, solve) to `n_new` regions with
zero equation-code changes. At 6 island groups this yields ~80k variables,
inside the regime where Ipopt+MUMPS already factors cleanly, letting the levels
equations be validated against a benchmark solve before tackling the full-scale
2.4M-variable factorization problem.

Correctness note: `prepare_parameters!` *re-derives* every CES/CET calibration
share from the summed flows, so only the raw value flows must aggregate exactly
(sum). Intensive region-indexed parameters (distances, elasticities, the Frisch
parameter) are collapsed by an unweighted mean — imperfect weighting there
cannot break benchmark replication (calibration absorbs benchmark levels); it
only perturbs shock magnitudes, which is acceptable for an equation-correctness
gate.

Region axes are detected automatically as any array axis whose length equals
`nr_old` (= 34). This is unambiguous for INDOTERM: no other set dimension
(na=25, nu=29, nm=9, no=4, ns=2) collides with 34.
"""

# Reduce a single axis of `arr` from `nr_old` groups down to `n_new`, using
# `rmap` (old-index → new-group). `:sum` accumulates (value flows); `:avg`
# takes the unweighted mean over the provinces in each new group (intensive
# quantities).
function _reduce_axis(arr::AbstractArray{T}, axis::Int, rmap::Vector{Int},
                      n_new::Int, mode::Symbol) where {T}
    sz = collect(size(arr))
    @assert sz[axis] == length(rmap) "axis $axis size $(sz[axis]) != rmap length $(length(rmap))"
    outsz = copy(sz)
    outsz[axis] = n_new
    out = zeros(T, outsz...)
    nd = ndims(arr)
    for I in CartesianIndices(arr)
        t = Tuple(I)
        g = rmap[t[axis]]
        Iout = ntuple(k -> k == axis ? g : t[k], nd)
        out[Iout...] += arr[I]
    end
    if mode == :avg
        counts = zeros(Int, n_new)
        for g in rmap
            counts[g] += 1
        end
        for I in CartesianIndices(out)
            g = Tuple(I)[axis]
            if counts[g] > 0
                out[I] /= counts[g]
            end
        end
    end
    out
end

# Reduce every region axis (length == nr_old) of `arr`.
function _reduce_region_axes(arr::AbstractArray, rmap::Vector{Int}, n_new::Int,
                             nr_old::Int, mode::Symbol)
    out = arr
    ax = 1
    while ax <= ndims(out)
        if size(out, ax) == nr_old
            out = _reduce_axis(out, ax, rmap, n_new, mode)
        end
        ax += 1
    end
    out
end

# Value flows: sum-aggregate over every region axis.
const _REGION_FLOW_KEYS = Set([
    "MAKE", "TRAD", "TMAR", "MARS",
    "1LAB", "1CAP", "1LND", "1PTX", "STOK",
    "STOC",   # CAPSTOK (Excerpt 50) — a capital-stock VALUE flow, not a rate
    "BSMR", "UTAX", "2PUR",
    "PO01",   # population is extensive → sum
])

# Intensive region-indexed parameters: unweighted mean over the region axis.
const _REGION_AVG_KEYS = Set([
    "DIST", "XPEL", "P021",
    "ALFA", "RADJ", "REXP",   # Step-6 behavioural params — not flow ratios
    "EMPR", "ELWG",
])

# Capital-block RATES: weighted by the capital stock, not an unweighted mean.
# These are ratios whose denominator is CAPSTOK — DPRC = deprec/K, TARG
# (= RNORMAL) = CAP/K, TFRO (= GROTREND) = trend investment/K, and QRAT as a
# ratio of two such ratios. Averaging them unweighted while summing "STOC" broke
# the steady-state identity GROSSRET == RNORMAL that `build_premod!.jl:131-139`
# constructs the stock to satisfy; see `test/diag_rnormal_consistency.jl`, which
# measured the identity failing for 15 of 25 sectors before this was fixed.
const _REGION_CAPWTD_KEYS = Set(["DPRC", "TARG", "TFRO", "QRAT"])

"""
    aggregate_regions!(agg, rmap=REG_MAP_34_to_6, n_new=length(REG6)) -> Dict

Collapse the region axes of an `agg` dict (from `aggregate_model!`) to `n_new`
groups. Returns a fresh `Dict{String,Any}`; the input is not mutated despite the
`!` (kept for naming symmetry with `aggregate_model!`). Feed the result straight
into `prepare_parameters!`.
"""
function aggregate_regions!(agg::Dict{String,Any},
                            rmap::Vector{Int}=REG_MAP_34_to_6,
                            n_new::Int=length(REG6))
    # `rmap` and `n_new` are two independent arguments describing ONE partition, so
    # nothing structural stops a caller passing a map that does not cover `1:n_new`.
    # A short map silently leaves the tail groups identically zero and the result
    # still looks like a valid database — `test/pipeline_cache.jl` did exactly that
    # for every nr ∉ {6,34} until 2026-08-01. Surjectivity is the cheap invariant
    # that makes the mistake impossible rather than merely undocumented.
    sort(unique(rmap)) == collect(1:n_new) || error(
        "aggregate_regions!: rmap must be onto 1:$n_new, got groups " *
        "$(sort(unique(rmap))) — every target region must receive at least one source region")
    nr_old = length(rmap)
    out = Dict{String,Any}()

    # Capital-stock weights, read BEFORE the loop so they are the pre-reduction
    # 34-region values. QRAT is a ratio of investment/capital ratios, so its exact
    # weight is trend investment = GROTREND*CAPSTOK rather than the stock itself —
    # same rule as `aggregate_model!`.
    _stoc = haskey(agg, "STOC") && agg["STOC"] !== nothing ?
        (agg["STOC"] isa NamedArray ? parent(agg["STOC"]) : agg["STOC"]) : nothing
    _tfro = haskey(agg, "TFRO") && agg["TFRO"] !== nothing ?
        (agg["TFRO"] isa NamedArray ? parent(agg["TFRO"]) : agg["TFRO"]) : nothing
    _capwt_of(k) = (k == "QRAT" && _stoc !== nothing && _tfro !== nothing &&
                    size(_tfro) == size(_stoc)) ? _tfro .* _stoc : _stoc

    for (k, v) in agg
        if v === nothing
            out[k] = v
            continue
        end
        arr = v isa NamedArray ? parent(v) : v
        if !(arr isa AbstractArray) || ndims(arr) == 0
            out[k] = arr            # scalars / non-arrays pass through
            continue
        end
        if !any(d -> size(arr, d) == nr_old, 1:ndims(arr))
            out[k] = arr            # no region axis → region-free, pass through
            continue
        end
        # Capital-block rates get the same treatment as in `aggregate_model!`: an
        # unweighted mean over merged provinces would break the steady-state
        # identity GROSSRET == RNORMAL, because CAPSTOK ("STOC") is summed on the
        # line above while the rate would be averaged. Weighting by the stock
        # makes sum(CAP)/sum(CAPSTOK) come out exactly right instead.
        _capwt = k in _REGION_CAPWTD_KEYS ? _capwt_of(k) : nothing
        if _capwt !== nothing && size(arr) == size(_capwt)
            num = _reduce_region_axes(arr .* _capwt, rmap, n_new, nr_old, :sum)
            den = _reduce_region_axes(_capwt,        rmap, n_new, nr_old, :sum)
            res = similar(num)
            for I in CartesianIndices(num)
                res[I] = den[I] > 0 ? num[I] / den[I] : 0.0
            end
            out[k] = res
            continue
        end
        mode = if k in _REGION_FLOW_KEYS
            :sum
        elseif k in _REGION_AVG_KEYS
            :avg
        elseif k in _REGION_CAPWTD_KEYS
            # Reached only if the stock weights are missing or mis-shaped. Fall
            # back to the old unweighted mean — inexact, but a rate must never
            # fall through to the :sum default below, which would multiply it by
            # the number of merged provinces.
            @warn "aggregate_regions!: no usable capital-stock weight for \"$k\"; " *
                  "falling back to an unweighted mean (steady-state identity will be approximate)"
            :avg
        else
            @warn "aggregate_regions!: key \"$k\" has a region axis but no reduction rule; defaulting to :sum (value-flow)"
            :sum
        end
        out[k] = _reduce_region_axes(arr, rmap, n_new, nr_old, mode)
    end
    out
end
