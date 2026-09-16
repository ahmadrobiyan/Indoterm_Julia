"""
IO 2020 gap audit — Step 1 of the 2020 rebase plan.

**Read-only.** Writes nothing, solves nothing, touches no cache. Safe to run any time.
Cost: ~20-60 s, almost all of it streaming the 71 MB `data/national_data.csv`.

## Why this exists

The model's whole database is 2016 vintage. BPS has published the official 2020 Input-Output
table (185 products, domestic transactions, basic prices). Before spending any effort on a
rebase, this script answers three questions with evidence rather than assertion:

1. **Is the 2020 classification still the model's classification?** The project's concordances
   (`SEC_MAP_185_to_25` in `src/aggregation_data.jl`, `COM`/`IND` in `src/prepare_sets.jl`) are
   *positional integer arrays*. Nothing in the codebase has ever validated that position N still
   names the same product. If BPS reordered or re-scoped anything, every downstream aggregate
   would be silently wrong. §1 turns that into a hard assertion.
2. **How much of a `national.har` does the 2020 table actually supply?** §2 walks the 30
   `nat["..."]` headers `src/build_reg0!.jl` consumes and classifies each as supplied / partial /
   absent, with the dimension the pipeline expects. This is the sourcing shopping list.
3. **Is the rebase worth it?** §3 aggregates both vintages' domestic intermediate matrices to the
   25 model sectors and reports where the structure actually moved. If nothing moved, the 2016
   baseline is fine and the rebase is not worth the data-construction effort.

§4 checks the 2020 table balances against its own published control totals before anything is
built on top of it.

## What it does NOT do

- Does not build a 2020 database, write to `data/`, or modify the 2016 vintage in any way.
- Does not compare *levels* between vintages. The 2016 HAR and the 2020 CSV are in different
  units (the CSV is Juta Rupiah; the HAR's `1BAS` is in the authors' own scale). Every
  cross-vintage number below is a **share**, never a level. Do not read §3 as growth rates.
- Does not touch `origin/`.

Run with:
```bash
julia --project=IndotermJulia IndotermJulia/test/audit_io2020.jl
```
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using Printf

# ── locating the 2020 table ───────────────────────────────────────────────────
# Preferred home is `data/io2020/` (per the rebase plan — 2020 inputs are kept out of `data/`
# proper so the validated 2016 database is never shadowed). Falls back to the original download
# location, and to an INDOTERM_IO2020 env override.
const IO2020_CANDIDATES = [
    get(ENV, "INDOTERM_IO2020", ""),
    joinpath(@__DIR__, "..", "data", "io2020",
             "io2020_domestic_basic_185.csv"),
    raw"C:\Users\ahmad\OneDrive\Documents\13. KERJA 2026\Tabel Input-Output Indonesia Transaksi Domestik Atas Dasar Harga Dasar (185 Produk), 2020 (Juta Rupiah).csv",
]

function locate_io2020()
    for p in IO2020_CANDIDATES
        isempty(p) && continue
        isfile(p) && return p
    end
    error("""
    2020 IO table not found. Looked in:
    $(join(["  - " * p for p in IO2020_CANDIDATES if !isempty(p)], "\n"))
    Set INDOTERM_IO2020 to the CSV path, or place it at data/io2020/io2020_domestic_basic_185.csv.""")
end

# ── minimal quoted-CSV reader ─────────────────────────────────────────────────
# Product names contain commas ("Jasa Kesenian, Hiburan dan Rekreasi"), so a plain `split(l,',')`
# corrupts the row. This project deliberately carries no CSV.jl dependency (`src/read_data.jl`
# hand-parses too), so parse quotes here rather than adding one for a single audit.
function split_csv_line(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    inquote = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if inquote
            if c == '"'
                if i < lastindex(line) && line[nextind(line, i)] == '"'
                    write(buf, '"'); i = nextind(line, i)   # escaped ""
                else
                    inquote = false
                end
            else
                write(buf, c)
            end
        elseif c == '"'
            inquote = true
        elseif c == ','
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    fields
end

"""
    bps_number(s) -> Float64

BPS cell formats seen in this release: `" 21,245,051 "`, `" - "` (zero), and negatives as either
`"-1,234"` or `"(1,234)"`. A bare `"-"` is BPS's zero, NOT a minus sign — getting that backwards
would poison every subsequent total, so it is handled explicitly before any sign logic.
"""
function bps_number(s::AbstractString)
    t = strip(s)
    (isempty(t) || t == "-" || t == "–") && return 0.0
    neg = false
    if startswith(t, "(") && endswith(t, ")")
        neg = true; t = t[2:end-1]
    end
    t = replace(t, "," => "", " " => "")
    isempty(t) && return 0.0
    v = tryparse(Float64, t)
    v === nothing && error("unparseable BPS cell: $(repr(s))")
    neg ? -v : v
end

"""
    read_io2020(path) -> NamedTuple

Layout of the BPS 2020 domestic/basic release, 0-indexed rows as the file ships:
row 3 = column codes, row 4 = column names, rows 5-189 = products 1-185, row 190 = i.f.
adjustment, rows 191-198 = aggregate rows (1900, 1950, 2000, 2010, 2020, 2030, 2090, 2100).
Columns: 0 blank, 1 `Kode`, 2 `Produk`, 3-187 = industries 1-185, 188-210 = the 23 aggregate
columns (1800, 3011…3100, 4011…4019, 5011…5019, 6090, 7000, 8000).
"""
function read_io2020(path::String)
    rows = [split_csv_line(l) for l in eachline(path)]
    length(rows) >= 199 || error("2020 CSV has only $(length(rows)) rows; expected >= 199")

    colcodes = strip.(rows[4])          # 1-indexed Julia view of 0-indexed row 3
    colnames = strip.(rows[5])
    prodrows = rows[6:190]              # the 185 product rows
    length(prodrows) == 185 || error("expected 185 product rows, got $(length(prodrows))")

    prodnames = [strip(r[3]) for r in prodrows]
    prodcodes = [strip(r[2]) for r in prodrows]

    # industry j lives at column index 3+j (1-indexed)
    Z = [bps_number(prodrows[i][3+j]) for i in 1:185, j in 1:185]

    # aggregate columns, addressed by their BPS code
    colof = Dict(c => k for (k, c) in enumerate(colcodes) if !isempty(c))
    getcol(code) = [bps_number(prodrows[i][colof[code]]) for i in 1:185]

    # aggregate rows, addressed by BPS code, over the 185 industry columns
    rowof = Dict(strip(r[2]) => r for r in rows if length(r) > 2 && !isempty(strip(r[2])))
    getrow(code) = [bps_number(rowof[code][3+j]) for j in 1:185]

    (; path, prodnames, prodcodes, Z, getcol, getrow, colof, colnames, rows)
end

# ── §1  classification assertion ──────────────────────────────────────────────
"""
Positional spot-checks: the model's `COM[i]` against the 2020 table's product i. These pairs were
verified by hand when the table was first inspected; encoding them here means a future BPS
re-release that reorders products fails loudly instead of silently mis-mapping every aggregate.
"""
const ANCHOR_CHECKS = [
    (1,   "Rice",         "Padi"),
    (26,  "Livestock",    "Ternak"),
    (37,  "CoalLignite",  "Batubara"),
    (38,  "CrudeOil",     "Minyak Bumi"),
    (50,  "CoarseSalt",   "Garam"),
    (96,  "BasChemicals", "Kimia Dasar"),
    (111, "Glass",        "Kaca"),
    (147, "WaterSupply",  "Pengadaan Air"),
    (178, "GovEducation", "Jasa Pendidikan Pemerintah"),
    (185, "OthSvc",       "Jasa Lainnya"),
]

function section1_classification(io)
    println("\n", "="^78)
    println("§1  CLASSIFICATION — 2020 product list vs. the model's COM array")
    println("="^78)

    @assert length(io.prodnames) == 185 "2020 table must have exactly 185 products"
    @assert length(COM) == 185

    # BPS codes must be 1..185 in order — this is what makes the positional maps legitimate
    badcode = findfirst(i -> tryparse(Int, io.prodcodes[i]) != i, 1:185)
    if badcode !== nothing
        error("2020 product codes are not 1..185 in order (position $badcode has code " *
              "$(repr(io.prodcodes[badcode]))). SEC_MAP_185_to_25 is a POSITIONAL map and is " *
              "invalid against this file — stop here.")
    end

    fails = String[]
    for (i, com, indo) in ANCHOR_CHECKS
        ok = occursin(lowercase(indo), lowercase(io.prodnames[i]))
        ok || push!(fails, "  pos $i: COM=$com expected Indonesian name containing " *
                           "$(repr(indo)), got $(repr(io.prodnames[i]))")
        @printf("  %3d  %-14s  %-8s  %s\n", i, com, ok ? "OK" : "MISMATCH", io.prodnames[i])
    end

    if isempty(fails)
        println("\n  ✅ 185 products, codes 1-185 in order, all $(length(ANCHOR_CHECKS)) anchors match.")
        println("     COM / IND / SEC_MAP_185_to_25 carry over to 2020 UNCHANGED.")
    else
        println("\n  ❌ CLASSIFICATION MISMATCH — the positional concordance is NOT valid:")
        foreach(println, fails)
        error("classification check failed; a rebase on this file would silently mis-map sectors")
    end
    nothing
end

# ── §2  header coverage ───────────────────────────────────────────────────────
# The 30 headers `src/build_reg0!.jl` reads, with the dimension the pipeline expects (confirmed
# against row counts in data/national_data.csv, e.g. 1BAS = 68450 = 185×2×185) and what the 2020
# domestic/basic release can supply.
const HEADER_COVERAGE = [
    ("1BAS", "COM×SRC×IND  185×2×185", :partial,
     "domestic half is observed (rows 5-189 × cols 3-187); src=imp absent"),
    ("2BAS", "COM×SRC×IND  investment", :partial,
     "only the 3030 PMTB column total by product; no investing-industry or src split"),
    ("3BAS", "COM×SRC×HOU  household", :partial,
     "only the 3011 column total by product; no src split"),
    ("4BAS", "COM×EXP      exports", :partial,
     "3050 + 3060 column totals by product; f.o.b. so already domestic-only"),
    ("5BAS", "COM×SRC      government", :partial,
     "only the 3020 column total by product; no src split"),
    ("6BAS", "COM×SRC      inventories", :partial,
     "only the 3040 column total by product; no src split"),
    ("1MAR", "COM×SRC×IND×MAR", :absent,
     "only product-level margin totals 5011/5012/5013; no by-user, no by-margin-commodity"),
    ("2MAR", "investment margins", :absent, "same — no by-user breakdown exists in this release"),
    ("3MAR", "household margins", :absent, "same"),
    ("4MAR", "export margins", :absent, "same"),
    ("5MAR", "government margins", :absent, "same"),
    ("1TAX", "COM×SRC×IND  taxes", :absent,
     "only row 1950 (pajak dikurang subsidi atas produk) per industry; no per-cell detail"),
    ("2TAX", "investment taxes", :absent, "same"),
    ("3TAX", "household taxes", :absent, "same"),
    ("4TAX", "export taxes", :absent, "same"),
    ("5TAX", "government taxes", :absent, "same"),
    ("1PTX", "IND  production tax", :supplied,
     "row 2030 (pajak dikurang subsidi lainnya atas produksi) per industry"),
    ("0TAR", "COM  tariffs", :absent, "requires the import table; not in a domestic release"),
    ("1LAB", "IND×OCC  185×4", :partial,
     "row 2010 gives the industry total only; the 4-way OCC split needs Sakernas"),
    ("1CAP", "IND  capital rentals", :partial,
     "lumped with land and mixed income in row 2020 (surplus usaha bruto)"),
    ("1LND", "IND  land rentals", :partial, "same — row 2020 is not decomposed"),
    ("1OCT", "IND  other cost tickets", :absent, "no counterpart in the published table"),
    ("MAKE", "COM×IND  185×185", :absent,
     "this is a USE table; the make/supply matrix is the separate SUT publication"),
    ("LCOM", "COM  local-market flag", :absent, "ORANI-G construct, not an IO-table field"),
    ("1ARM", "COM  Armington sigma", :absent, "behavioural parameter — 2016 carry-over by design"),
    ("R001", "COM×REG  185×34", :absent, "national table; needs IRIO 2020 or provincial PDRB"),
    ("R002", "COM×REG  185×34", :absent, "same"),
    ("R003", "COM×REG  185×34", :absent, "same"),
    ("R004", "COM×REG  185×34", :absent, "same"),
    ("R005", "COM×REG  185×34", :absent, "same"),
]

function section2_coverage(io)
    println("\n", "="^78)
    println("§2  HEADER COVERAGE — what build_reg0! needs vs. what 2020 domestic/basic supplies")
    println("="^78)

    mark = Dict(:supplied => "✅ supplied", :partial => "🟡 partial ", :absent => "❌ absent  ")
    @printf("  %-6s %-26s %-12s %s\n", "header", "expected dimension", "2020", "note")
    println("  " * "-"^100)
    for (h, dim, status, note) in HEADER_COVERAGE
        @printf("  %-6s %-26s %-12s %s\n", h, dim, mark[status], note)
    end

    n = length(HEADER_COVERAGE)
    ns = count(t -> t[3] === :supplied, HEADER_COVERAGE)
    np = count(t -> t[3] === :partial,  HEADER_COVERAGE)
    na = count(t -> t[3] === :absent,   HEADER_COVERAGE)
    @assert n == 30 "coverage table must account for all 30 nat[...] headers build_reg0! reads"
    println("\n  $n headers: $ns fully supplied, $np partial, $na absent.")
    println("""
      Reading: "partial" means the 2020 table gives a CONTROL TOTAL but not the interior split.
      Those are RAS targets — the 2016 interior proportions get rebalanced onto 2020 totals.
      "absent" splits into two very different cases: behavioural parameters (1ARM and the other
      elasticities) are 2016 carry-over BY DESIGN and are not a gap; MAKE, 0TAR, the margin and
      tax matrices, and R001-R005 are genuine missing DATA that must be sourced.""")
    nothing
end

# ── §3  structural change, 2016 vs 2020 ───────────────────────────────────────
"""
    stream_har_header(path, code; ndim) -> Vector{Tuple{Vector{String},Float64}}

Pull ONE header out of the 71 MB long-format HAR export without materialising the other 56.
`read_national_data()` (src/read_data.jl) parses the whole file into ~57 NamedArrays, which is
right for the pipeline and wasteful for an audit that wants four of them.
"""
function stream_har_header(path::String, code::String)
    want = "__header__:" * code * ","
    out = Tuple{Vector{String},Float64}[]
    active = false
    for line in eachline(path)
        if startswith(line, "__header__:")
            active && break            # header blocks are contiguous — stop at the next one
            active = startswith(line, want)
            continue
        end
        active || continue
        f = split(line, ',')
        push!(out, (String.(f[1:end-1]), parse(Float64, f[end])))
    end
    isempty(out) && error("header $code not found in $path")
    out
end

function section3_structure(io)
    println("\n", "="^78)
    println("§3  STRUCTURAL CHANGE — 2016 vs 2020, domestic intermediate matrix at 25 sectors")
    println("="^78)
    println("  (SHARES only — the two vintages are in different units and are not comparable in levels.)")

    natpath = joinpath(@__DIR__, "..", "data", "national_data.csv")
    isfile(natpath) || error("2016 database not found at $natpath")
    comidx = Dict(c => i for (i, c) in enumerate(COM))
    S = SEC_MAP_185_to_25
    na = length(AGGCOM)

    # 2016: domestic half of 1BAS (COM × SRC × IND), aggregated straight to 25×25
    print("  streaming 1BAS from the 2016 database … ")
    B16 = zeros(na, na)
    for (labels, v) in stream_har_header(natpath, "1BAS")
        labels[2] == "dom" || continue
        B16[S[comidx[labels[1]]], S[comidx[labels[3]]]] += v
    end
    println("done")

    # ── margin-accounting correction — do NOT skip this ───────────────────────
    # ORANI-G strips trade/transport margins OUT of `1BAS` into a separate `1MAR`
    # (COM×SRC×IND×MAR) matrix. The BPS basic-price use table does not: each cell is valued at
    # basic price and the margin appears as a direct purchase of trade/transport services by the
    # using industry, which is standard SNA convention. Comparing raw `1BAS` against the 2020
    # quadrant therefore reads the margin convention as a structural change — measured here, that
    # alone moved Trade's supplier share by ~7.3 pp, which is larger than every genuine sectoral
    # shift combined. `1MAR` is folded back in at [margin commodity, using industry] so the two
    # vintages are on the same accounting basis before any share is differenced.
    print("  streaming 1MAR (margin correction) … ")
    M16 = zeros(na, na)
    for (labels, v) in stream_har_header(natpath, "1MAR")
        labels[2] == "dom" || continue
        M16[S[comidx[labels[4]]], S[comidx[labels[3]]]] += v   # dim4 = margin com, dim3 = user ind
    end
    println("done")
    Z16 = B16 .+ M16
    @printf("    1MAR is %.1f%% of 1BAS; Trade supplier share %.2f%% → %.2f%% once folded in.\n",
            100 * sum(M16) / sum(B16),
            100 * sum(B16[20, :]) / sum(B16), 100 * sum(Z16[20, :]) / sum(Z16))

    # 2020: the observed domestic intermediate quadrant
    Z20 = zeros(na, na)
    for i in 1:185, j in 1:185
        Z20[S[i], S[j]] += io.Z[i, j]
    end

    t16, t20 = sum(Z16), sum(Z20)

    println("\n  Sector shares of total domestic intermediate flow (%):")
    @printf("  %-14s %8s %8s %8s   %8s %8s %8s\n",
            "sector", "row'16", "row'20", "Δrow", "col'16", "col'20", "Δcol")
    println("  " * "-"^70)
    drow = zeros(na); dcol = zeros(na)
    for a in 1:na
        r16 = 100 * sum(Z16[a, :]) / t16; r20 = 100 * sum(Z20[a, :]) / t20
        c16 = 100 * sum(Z16[:, a]) / t16; c20 = 100 * sum(Z20[:, a]) / t20
        drow[a] = r20 - r16; dcol[a] = c20 - c16
        @printf("  %-14s %8.2f %8.2f %+8.2f   %8.2f %8.2f %+8.2f\n",
                AGGCOM[a], r16, r20, drow[a], c16, c20, dcol[a])
    end

    println("\n  Largest movers (row = sector as SUPPLIER of intermediates):")
    for a in sortperm(abs.(drow); rev=true)[1:5]
        @printf("    %-14s %+6.2f pp\n", AGGCOM[a], drow[a])
    end
    println("  Largest movers (col = sector as USER of intermediates):")
    for a in sortperm(abs.(dcol); rev=true)[1:5]
        @printf("    %-14s %+6.2f pp\n", AGGCOM[a], dcol[a])
    end

    println("\n  Largest cell-level share shifts (share of total intermediate flow):")
    D = 100 .* (Z20 ./ t20 .- Z16 ./ t16)
    ord = sortperm(vec(abs.(D)); rev=true)[1:12]
    @printf("    %-14s → %-14s %8s %8s %8s\n", "supplier", "user", "'16 %", "'20 %", "Δpp")
    for k in ord
        a, b = Tuple(CartesianIndices(D)[k])
        @printf("    %-14s → %-14s %8.3f %8.3f %+8.3f\n",
                AGGCOM[a], AGGCOM[b], 100*Z16[a,b]/t16, 100*Z20[a,b]/t20, D[a,b])
    end

    # primary-factor shares — small headers, worth streaming for the labour-share comparison
    println("\n  Primary factor shares of total primary input (%):")
    lab16 = sum(v for (_, v) in stream_har_header(natpath, "1LAB"))
    cap16 = sum(v for (_, v) in stream_har_header(natpath, "1CAP"))
    lnd16 = sum(v for (_, v) in stream_har_header(natpath, "1LND"))
    tot16 = lab16 + cap16 + lnd16
    lab20 = sum(io.getrow("2010"))
    gos20 = sum(io.getrow("2020"))          # capital + land + mixed income, not decomposed
    tot20 = lab20 + gos20
    @printf("    2016  labour %5.1f   capital %5.1f   land %5.1f\n",
            100*lab16/tot16, 100*cap16/tot16, 100*lnd16/tot16)
    @printf("    2020  labour %5.1f   gross operating surplus %5.1f  (capital+land NOT split)\n",
            100*lab20/tot20, 100*gos20/tot20)
    println("""
      ⚠ These two labour shares are NOT comparable and no conclusion should be drawn from the
      gap. BPS row 2020 (surplus usaha bruto) contains mixed income of the self-employed, which
      is large in Indonesia and which ORANI-G's 1LAB/1CAP split allocates differently. The gap
      measures a definitional difference of unknown size, not a change in factor shares.
      Decomposing row 2020 into 1CAP / 1LND / mixed income is item 11 on the plan's sourcing
      list; until it is sourced, 1CAP/1LND must inherit the 2016 proportions.""")

    maxmove = maximum(abs.(vcat(drow, dcol)))
    println("\n  Verdict: largest sector-share move is $(round(maxmove; digits=2)) pp.")
    println(maxmove < 0.5 ?
        "    → Structure is close to unchanged. A rebase would move little; the 2016 baseline\n" *
        "      is defensible and the sourcing effort is hard to justify on these grounds." :
        "    → Structure has moved materially. The 2016 baseline is stale in ways that bear\n" *
        "      directly on the intermediate block this model solves; the rebase is justified.")
    (; Z16, Z20, drow, dcol)
end

# ── §4  internal balance of the 2020 table ────────────────────────────────────
function section4_balance(io)
    println("\n", "="^78)
    println("§4  INTERNAL BALANCE of the 2020 table (before anything is built on it)")
    println("="^78)

    reldiff(a, b) = maximum(abs.(a .- b) ./ max.(1.0, abs.(b)))

    # (a) per industry: column sum of the domestic quadrant vs published 1900 Total Input Antara
    colsum = vec(sum(io.Z; dims=1))
    pub_int = io.getrow("1900")
    d1 = reldiff(colsum, pub_int)
    @printf("  intermediate column sums vs row 1900   max rel diff %.3e  %s\n",
            d1, d1 < 1e-6 ? "✅" : "❌")

    # (b) per industry: 2090 Total Input Primer + 1900 + 1950 + 2000 vs 2100 Total Input
    lhs = io.getrow("1900") .+ io.getrow("1950") .+ io.getrow("2000") .+ io.getrow("2090")
    d2 = reldiff(lhs, io.getrow("2100"))
    @printf("  1900+1950+2000+2090 vs row 2100        max rel diff %.3e  %s\n",
            d2, d2 < 1e-6 ? "✅" : "❌")

    # (c) per product: row sum of the domestic quadrant vs published 1800 Total Permintaan Antara
    rowsum = vec(sum(io.Z; dims=2))
    d3 = reldiff(rowsum, io.getcol("1800"))
    @printf("  intermediate row sums vs col 1800      max rel diff %.3e  %s\n",
            d3, d3 < 1e-6 ? "✅" : "❌")

    # (d) per product: 1800 + 3090 Total Permintaan Akhir vs 3100 Total Permintaan
    d4 = reldiff(io.getcol("1800") .+ io.getcol("3090"), io.getcol("3100"))
    @printf("  1800+3090 vs col 3100                  max rel diff %.3e  %s\n",
            d4, d4 < 1e-6 ? "✅" : "❌")

    # (e) per product: supply side 8000 vs demand side 3100
    d5 = reldiff(io.getcol("8000"), io.getcol("3100"))
    @printf("  col 8000 (supply) vs col 3100 (demand) max rel diff %.3e  %s\n",
            d5, d5 < 1e-6 ? "✅" : "❌")

    worst = maximum((d1, d2, d3, d4, d5))
    println(worst < 1e-6 ?
        "\n  ✅ The published 2020 table is internally consistent. Safe to build on." :
        "\n  ⚠ The 2020 table does NOT close to its own control totals (worst $(worst)).\n" *
        "    Find the cause BEFORE propagating it — a rounding artifact in the published\n" *
        "    release and a parsing bug in read_io2020 look identical at this stage.")
    (; d1, d2, d3, d4, d5)
end

# ── driver ────────────────────────────────────────────────────────────────────
function run_audit()
    path = locate_io2020()
    println("="^78)
    println("IO 2020 GAP AUDIT — read-only, writes nothing")
    println("="^78)
    println("  2020 source : $path")
    println("  2016 source : ", joinpath(@__DIR__, "..", "data", "national_data.csv"))

    io = read_io2020(path)
    section1_classification(io)
    section2_coverage(io)
    bal = section4_balance(io)        # balance first — §3 is meaningless on an unbalanced table
    st  = section3_structure(io)

    println("\n", "="^78)
    println("NEXT STEP")
    println("="^78)
    println("""
    The absent-DATA rows in §2 are the sourcing shopping list. In priority order:
      1. IO 2020 Transaksi Total (harga dasar)   → imports by subtraction; unblocks 1BAS/2BAS/3BAS/5BAS
      2. IO 2020 Transaksi Total (harga pembeli) → margin+tax wedge; unblocks 1MAR…5MAR, 1TAX…5TAX
      3. SUT Indonesia 2020, supply side         → MAKE
      4. IRIO 2020 (or provincial PDRB 2020)     → R001-R005
      5. Sakernas 2020 industry×occupation       → 1LAB's 4-way split
    Nothing downstream of build_reg0! changes — that is the authors' pipeline and it stays as is.""")
    (; io, bal, st)
end

# Runs when this file is executed directly, same convention as the test/verify_*.jl gates.
run_audit()
