"""
national.har / regsupp.har / DISTGONE.HAR ingestion.

`HeaderArrayFile.jl` v0.2.0 fails on three national.har headers (`P21H`, `XPLH`, `3PUR`) with
"Data dimensions do not match metadata" -- confirmed (see PLAN.md Step 1) to be a package
limitation around headers with a size-1 dimension (`HOU` has only 1 element, `AllHou`), not an
"Old Lahey" format issue (this file is the modern format; `harpy` reads it correctly, including
those three headers, e.g. `3PUR`'s shape is `(185, 1)` over `[COM, HOU]`). Rather than patch a
vendored copy of that package, this project exports every header from each `.har` file once via
`harpy` (`data/export_har_csv.py`, already run -- see `data/*_data.csv`) into a generic long-format
CSV (one `__header__:NAME,dim1,...,Value` block per header, in HAR file order, with set-element
string labels, dense -- zeros included), analogous to WayangJulia's `har2csv.exe` conversion but
using the already-validated `harpy` reader as the extraction tool. IndotermJulia has no runtime
GEMPACK/harpy dependency; only these checked-in CSVs are read at runtime.
"""

const DATA_DIR = joinpath(@__DIR__, "..", "data")

"""
    _read_har_csv(path) -> Dict{String, Any}

Parses one generic long-format HAR-export CSV into a `Dict` keyed by HAR header code (e.g.
`"3PUR"`), each value a `NamedArray` (dimension names/labels are anonymous positional indices
here -- Step 1 does not yet map header codes to TABLO coefficient names/set semantics, that is
Step 4's job) or a bare `Float64` for scalar headers.
"""
function _read_har_csv(path::String)
    data = Dict{String,Any}()
    header_code = ""
    ndim = 0
    col_labels = Vector{Vector{String}}()   # per-dim ordered-unique labels, first-seen order
    col_index = Vector{Dict{String,Int}}()
    rows = Vector{Tuple{Vector{String},Float64}}()

    function flush!()
        isempty(header_code) && return
        if ndim == 0
            length(rows) == 1 || error("Expected exactly one row for scalar header $header_code")
            data[header_code] = rows[1][2]
        else
            dims = length.(col_labels)
            arr = zeros(Float64, dims...)
            for (labels, val) in rows
                idx = ntuple(k -> col_index[k][labels[k]], ndim)
                arr[idx...] = val
            end
            data[header_code] = length(dims) == 1 ? NamedArray(arr, col_labels[1]) :
                                 NamedArray(arr, Tuple(col_labels))
        end
    end

    for line in eachline(path)
        isempty(line) && continue
        fields = split(line, ',')
        if startswith(fields[1], "__header__:") || fields[1] == "__scalar__"
            flush!()
            header_code = startswith(fields[1], "__header__:") ? String(fields[1][12:end]) : ""
            if fields[1] == "__scalar__"
                header_code = ""  # filled in from the next data line's first field below
                ndim = 0
                rows = Tuple{Vector{String},Float64}[]
                col_labels = Vector{Vector{String}}()
                col_index = Vector{Dict{String,Int}}()
                # scalar blocks store the header name as dim1 of the single data row; handled below
                ndim = -1
                continue
            end
            ndim = length(fields) - 2   # minus header-name field and trailing "Value"
            col_labels = [String[] for _ in 1:ndim]
            col_index = [Dict{String,Int}() for _ in 1:ndim]
            rows = Tuple{Vector{String},Float64}[]
            continue
        end

        if ndim == -1
            # scalar block: "HEADERNAME,value"
            header_code = String(fields[1])
            data[header_code] = parse(Float64, fields[2])
            ndim = 0
            header_code = ""
            continue
        end

        labels = String.(fields[1:end-1])
        val = parse(Float64, fields[end])
        for k in 1:ndim
            if !haskey(col_index[k], labels[k])
                push!(col_labels[k], labels[k])
                col_index[k][labels[k]] = length(col_labels[k])
            end
        end
        push!(rows, (labels, val))
    end
    flush!()
    return data
end

"""
    read_national_data() -> Dict{String, Any}

Raw ingestion of `national.har` (via its CSV export), keyed by HAR header code.
"""
read_national_data() = _read_har_csv(joinpath(DATA_DIR, "national_data.csv"))

"""
    read_regsupp_data() -> Dict{String, Any}

Raw ingestion of `regsupp.har` (via its CSV export), keyed by HAR header code.
"""
read_regsupp_data() = _read_har_csv(joinpath(DATA_DIR, "regsupp_data.csv"))

"""
    read_distgone_data() -> Dict{String, Any}

Raw ingestion of `DISTGONE.HAR` (via its CSV export), keyed by HAR header code.
"""
read_distgone_data() = _read_har_csv(joinpath(DATA_DIR, "distgone_data.csv"))
