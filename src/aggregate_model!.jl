Base.@kwdef struct AggResult
    agg::Dict{String,Any}
    mapping::Dict{String,Any}
end

"Sum-aggregate over dimension 1 (COM/IND) into `nnew` groups"
function _agg_first(arr::AbstractArray{T}, map_old::Vector{Int}, nnew::Int) where {T}
    sz = size(arr)
    rest_axes = axes(arr)[2:end]
    out = zeros(T, nnew, sz[2:end]...)
    for i in 1:length(map_old)
        out[map_old[i], rest_axes...] .+= arr[i, rest_axes...]
    end
    out
end

"Sum-aggregate MAKE array: both COM (dim 1) and IND (dim 2) → AGG"
function _agg_make(MAKE_raw::AbstractArray{T}, map::Vector{Int}, na::Int) where {T}
    out = zeros(T, na, na, size(MAKE_raw,3))
    for c in 1:length(map), i in 1:length(map), d in axes(MAKE_raw,3)
        out[map[c], map[i], d] += MAKE_raw[c,i,d]
    end
    out
end

"Aggregate USE/TAX array: COM×SRC×USR×DST → AGG×SRC×USR_AGG×DST.
 USR = [IND(old 185) + FINDEM(Hou,Inv,Gov,Exp)] → [AGG(25) + FINDEM(4)]"
function _agg_use(USE_raw::AbstractArray{T}, map::Vector{Int}, na::Int, nr::Int) where {T}
    ns = size(USE_raw, 2)
    nu_old = size(USE_raw, 3)
    nu_new = na + 4
    out = zeros(T, na, ns, nu_new, nr)
    for c_old in 1:length(map)
        cn = map[c_old]
        for s in 1:ns
            for u_old in 1:nu_old
                un = u_old <= 185 ? map[u_old] : (na + (u_old - 185))
                @assert 1 <= un <= nu_new "un=$un out of range [1,$nu_new]"
                for d in 1:nr
                    out[cn, s, un, d] += USE_raw[c_old, s, u_old, d]
                end
            end
        end
    end
    out
end

"Recursively unwrap NamedArray to plain array"
_unwrap(x) = x
_unwrap(x::NamedArray) = parent(x)

function aggregate_model!(premod::Dict{String,Any})
    na = length(AGGCOM)
    nr = length(REG)
    nm = length(MAR)
    mp = SEC_MAP_185_to_25

    agg = Dict{String,Any}()

    # MAKE: COM×IND×DST → AGG×AGG×DST
    agg["MAKE"] = _agg_make(_unwrap(premod["MAKE"]), mp, na)

    # TRAD: COM×SRC×ORG×DST → AGG×SRC×ORG×DST
    agg["TRAD"] = _agg_first(_unwrap(premod["TRAD"]), mp, na)

    # TMAR: COM×SRC×MAR×ORG×DST → AGG×SRC×MAR×ORG×DST
    if haskey(premod, "TMAR") && premod["TMAR"] !== nothing
        agg["TMAR"] = _agg_first(_unwrap(premod["TMAR"]), mp, na)
    end

    # MARS: MAR×ORG×DST×PRD → identity (no COM/IND dims)
    if haskey(premod, "MARS") && premod["MARS"] !== nothing
        agg["MARS"] = copy(_unwrap(premod["MARS"]))
    end

    # DIST: ORG×DST — no COM/IND dimension, pass through
    if haskey(premod, "DIST") && premod["DIST"] !== nothing
        agg["DIST"] = copy(_unwrap(premod["DIST"]))
    end

    # IND-indexed flows: sum-aggregate over IND
    for k in ["1LAB", "1CAP", "1LND"]
        if haskey(premod, k)
            agg[k] = _agg_first(_unwrap(premod[k]), mp, na)
        end
    end

    # IND-indexed parameters: weighted-average over IND
    # First build IND value weights from MAKE
    MAKE_raw = _unwrap(premod["MAKE"])
    MAKE_i = [sum(MAKE_raw[:,i,:]) for i in 1:size(MAKE_raw,2)]
    MAKE_i_agg = zeros(Float64, na)
    for i in 1:length(mp)
        MAKE_i_agg[mp[i]] += MAKE_i[i]
    end

    # COM-indexed parameters: weighted-average over COM
    # Build COM weights from TRAD
    TRAD_raw = _unwrap(premod["TRAD"])
    TRAD_c = [sum(TRAD_raw[c,:,:,:]) for c in 1:size(TRAD_raw,1)]
    TRAD_c_agg = zeros(Float64, na)
    for c in 1:length(mp)
        TRAD_c_agg[mp[c]] += TRAD_c[c]
    end

    function _wavg_1d(vec, map, wts_old, wts_new)
        out = zeros(Float64, na)
        for i in 1:length(map)
            if wts_old[i] > 0
                out[map[i]] += vec[i] * wts_old[i]
            end
        end
        for j in 1:na
            if wts_new[j] > 0
                out[j] /= wts_new[j]
            end
        end
        out
    end

    function _wavg_nd(arr, map, wts_old, wts_new)
        out = zeros(Float64, na, size(arr)[2:end]...)
        for i in 1:length(map), rest in CartesianIndices(size(arr)[2:end])
            if wts_old[i] > 0
                out[map[i], rest] += arr[i, rest] * wts_old[i]
            end
        end
        for j in 1:na, rest in CartesianIndices(size(arr)[2:end])
            if wts_new[j] > 0
                out[j, rest] /= wts_new[j]
            end
        end
        out
    end

    # IND-indexed parameters (weight by MAKE output)
    for k in ["SLAB", "P028", "SCET", "P018"]
        if haskey(premod, k)
            data = _unwrap(premod[k])
            if ndims(data) == 1 && size(data,1) == 185
                agg[k] = _wavg_1d(vec(data), mp, MAKE_i, MAKE_i_agg)
            elseif ndims(data) >= 2 && size(data,1) == 185
                agg[k] = _wavg_nd(data, mp, MAKE_i, MAKE_i_agg)
            else
                agg[k] = data
            end
        end
    end

    # IND×REG dynamic parameters (weight by MAKE output)
    for k in ["DPRC", "TARG", "TFRO", "QRAT", "ALFA", "RADJ", "REXP"]
        if haskey(premod, k)
            agg[k] = _wavg_nd(_unwrap(premod[k]), mp, MAKE_i, MAKE_i_agg)
        end
    end

    # COM-indexed parameters (weight by TRAD)
    for k in ["SGDD"]
        if haskey(premod, k)
            data = _unwrap(premod[k])
            if ndims(data) == 1 && size(data,1) == 185
                agg[k] = _wavg_1d(vec(data), mp, TRAD_c, TRAD_c_agg)
            elseif ndims(data) >= 2 && size(data,1) == 185
                agg[k] = _wavg_nd(data, mp, TRAD_c, TRAD_c_agg)
            else
                agg[k] = data
            end
        end
    end

    # COM-indexed: XPEL (COM×REG), LCOM (COM)
    if haskey(premod, "XPEL")
        agg["XPEL"] = _wavg_nd(_unwrap(premod["XPEL"]), mp, TRAD_c, TRAD_c_agg)
    end
    if haskey(premod, "LCOM") && premod["LCOM"] !== nothing
        lcom = _unwrap(premod["LCOM"])
        if length(lcom) == 185
            agg_lcom = zeros(Float64, na)
            for c in 1:185
                agg_lcom[mp[c]] += lcom[c]
            end
            agg["LCOM"] = agg_lcom
        else
            agg["LCOM"] = lcom
        end
    end

    # BSMR, UTAX: COM×SRC×USR×DST → AGG×SRC×USR_AGG×DST
    #   USR = IND(185) + FINDEM(4) → AGG(25) + FINDEM(4)
    for k in ["BSMR", "UTAX"]
        if haskey(premod, k) && premod[k] !== nothing
            data = _unwrap(premod[k])
            @assert size(data,1) == 185 "Expected $k dim 1 = 185, got $(size(data,1))"
            agg[k] = _agg_use(data, mp, na, nr)
        end
    end

    # 2PUR: COM×IND×DST → AGG×AGG×DST (like MAKE)
    if haskey(premod, "2PUR") && premod["2PUR"] !== nothing
        agg["2PUR"] = _agg_make(_unwrap(premod["2PUR"]), mp, na)
    end

    # STOK, 1PTX: IND×DST → AGG×DST
    for k in ["STOK", "1PTX"]
        if haskey(premod, k) && premod[k] !== nothing
            data = _unwrap(premod[k])
            if ndims(data) == 2 && size(data,1) == 185
                agg[k] = _agg_first(data, mp, na)
            else
                agg[k] = data
            end
        end
    end

    # REG-only and MAR-only (no COM/IND dims): pass through
    for k in ["PO01", "P021", "EMPR", "ELWG", "SMAR"]
        if haskey(premod, k) && premod[k] !== nothing
            agg[k] = _unwrap(premod[k])
        end
    end

    mapping = Dict{String,Any}(
        "MAPCOM" => mp,
        "MAPIND" => mp,
        "AGGCOM" => AGGCOM,
        "AGGIND" => AGGCOM,
    )

    AggResult(agg=agg, mapping=mapping)
end
