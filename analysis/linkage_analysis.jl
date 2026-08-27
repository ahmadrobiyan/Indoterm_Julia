"""
Backward/forward linkage (Rasmussen–Hirschman) analysis of the 2016 Indonesia national
input-output table underlying INDOTERM, run entirely from the Julia port (no GEMPACK/harpy
needed at runtime; `origin/sec.agg` is legacy and not used here — the 185→25 sector
aggregation comes from `AGGCOM`/`SEC_MAP_185_to_25` in `src/aggregation_data.jl`).

Focus: interconnectedness of the "Metals & Machinery & Transport-Equipment" manufacturing
complex — the three adjacent aggregated sectors `BasMetals`, `MetalMach`, `TranspEquip` —
with each other and with the rest of the economy.

Method
------
1. Read `1BAS` (intermediate commodity flows: commodity × source × using-industry, basic
   values) and `MAKE` (commodity × industry output) from `data/national_data.csv` via
   `read_national_data()`.
2. Keep only the *domestic*-sourced slice of `1BAS` (imports are a demand leakage, not part
   of the domestic production network) to build the 185×185 domestic intermediate-flow
   matrix Z.
3. Aggregate Z and the industry-output vector X (column sums of MAKE) from 185 to the 25
   `AGGCOM` sectors using `SEC_MAP_185_to_25`.
4. Technical coefficients A[i,j] = Z[i,j] / X[j] (demand side) and allocation coefficients
   B[i,j] = Z[i,j] / X[i] (supply side); Leontief inverse L = (I-A)^-1, Ghosh inverse
   G = (I-B)^-1.
5. Rasmussen–Hirschman indices, normalised so the 25-sector average is 1.0:
     BL_j        = colsum_j(L) / mean(colsums(L))   -- "power of dispersion" (backward linkage)
     FL_leontief = rowsum_i(L) / mean(rowsums(L))   -- demand-side forward-linkage proxy
     FL_ghosh    = rowsum_i(G) / mean(rowsums(G))   -- supply-side forward linkage (preferred)
   BL/FL > 1 = above-average; sectors with both > 1 are classified "key sectors".

`data/national_data.csv` is gitignored (regenerated from the checked-in
`data/national_data.zip`); unzip it once before running:
    unzip -o data/national_data.zip -d data/

Run with:  julia --project=. analysis/linkage_analysis.jl
Writes:    analysis/output/linkage_25sector.csv, analysis/output/triad_technical_coeffs.csv
"""

using IndotermJulia
using LinearAlgebra
using Printf

const OUTDIR = joinpath(@__DIR__, "output")
mkpath(OUTDIR)

# --- 1. Load base-year national flows -----------------------------------------------------
data = read_national_data()
bas1 = Array(data["1BAS"])   # (COM, SRC, IND) = (185, 2, 185), basic values
make = Array(data["MAKE"])   # (COM, IND)      = (185, 185), output values

@assert size(bas1) == (length(COM), length(SRC), length(IND))
@assert size(make) == (length(COM), length(IND))
dom_idx = findfirst(==("dom"), SRC)

Z185 = bas1[:, dom_idx, :]          # domestic intermediate flow, commodity i -> industry j
X185 = vec(sum(make; dims = 1))     # gross output by industry (column sum of MAKE)

# --- 2. Aggregate 185 -> 25 sectors (AGGCOM / SEC_MAP_185_to_25; sec.agg is legacy) --------
n = length(AGGCOM)
Zagg = zeros(n, n)
for i in 1:length(COM), j in 1:length(IND)
    Zagg[SEC_MAP_185_to_25[i], SEC_MAP_185_to_25[j]] += Z185[i, j]
end
Xagg = zeros(n)
for i in 1:length(COM)
    Xagg[SEC_MAP_185_to_25[i]] += X185[i]
end

# --- 3. Technical coefficients, Leontief & Ghosh inverses ----------------------------------
A = Zagg ./ reshape(Xagg, 1, n)   # a[i,j] = Z[i,j] / X[j]
B = Zagg ./ reshape(Xagg, n, 1)   # b[i,j] = Z[i,j] / X[i]
L = inv(I - A)                    # Leontief inverse (demand-driven)
G = inv(I - B)                    # Ghosh inverse (supply-driven)

col_sums_L = vec(sum(L; dims = 1))
row_sums_L = vec(sum(L; dims = 2))
row_sums_G = vec(sum(G; dims = 2))

BL = col_sums_L ./ (sum(col_sums_L) / n)
FL_leontief = row_sums_L ./ (sum(row_sums_L) / n)
FL_ghosh = row_sums_G ./ (sum(row_sums_G) / n)

sector_type(bl, fl) = bl > 1 && fl > 1 ? "Key sector" :
                       bl > 1 && fl <= 1 ? "Backward-oriented" :
                       bl <= 1 && fl > 1 ? "Forward-oriented" : "Weakly linked"

open(joinpath(OUTDIR, "linkage_25sector.csv"), "w") do io
    println(io, "sector,output_Xj,BL_backward,FL_leontief,FL_ghosh,type")
    for k in sortperm(BL; rev = true)
        @printf(io, "%s,%.6f,%.6f,%.6f,%.6f,%s\n",
            AGGCOM[k], Xagg[k], BL[k], FL_leontief[k], FL_ghosh[k],
            sector_type(BL[k], FL_ghosh[k]))
    end
end

println("Backward/forward linkages, 25 aggregated sectors (sorted by backward linkage):")
println(rpad("sector", 14), rpad("output_Xj", 16), rpad("BL", 8), rpad("FL_leo", 8), rpad("FL_ghosh", 10), "type")
for k in sortperm(BL; rev = true)
    @printf("%-14s%16.0f%8.3f%8.3f%10.3f  %s\n",
        AGGCOM[k], Xagg[k], BL[k], FL_leontief[k], FL_ghosh[k], sector_type(BL[k], FL_ghosh[k]))
end

# --- 4. Focus: Metals / Machinery&Electrical / Transport-Equipment triad -------------------
target = ["BasMetals", "MetalMach", "TranspEquip"]
tidx = [findfirst(==(t), AGGCOM) for t in target]

println("\nTarget-sector summary:")
for (t, k) in zip(target, tidx)
    @printf("  %-12s BL=%.3f (rank %d/%d)  FL_ghosh=%.3f (rank %d/%d)  -> %s\n",
        t, BL[k], count(>(BL[k]), BL) + 1, n,
        FL_ghosh[k], count(>(FL_ghosh[k]), FL_ghosh) + 1, n,
        sector_type(BL[k], FL_ghosh[k]))
end

println("\nDirect technical-coefficient sub-matrix among the triad (a[i,j] = share of j's total input cost sourced from i):")
println(rpad("from \\ to", 14), join(rpad.(target, 14)))
for i in tidx
    print(rpad(AGGCOM[i], 14))
    for j in tidx
        @printf("%-14.4f", A[i, j])
    end
    println()
end

open(joinpath(OUTDIR, "triad_technical_coeffs.csv"), "w") do io
    println(io, "from,to,a_ij")
    for i in tidx, j in tidx
        @printf(io, "%s,%s,%.6f\n", AGGCOM[i], AGGCOM[j], A[i, j])
    end
end

within_inputs = sum(Zagg[tidx, tidx])
total_inputs_to_triad = sum(Zagg[:, tidx])
within_sales = within_inputs
total_sales_from_triad = sum(Zagg[tidx, :])

@printf("\nShare of the triad's combined intermediate INPUTS sourced from within the triad: %.1f%%\n",
    100 * within_inputs / total_inputs_to_triad)
@printf("Share of the triad's combined intermediate SALES going to within the triad:       %.1f%%\n",
    100 * within_sales / total_sales_from_triad)

println("\nTop-5 input suppliers (by a[i,j]) for each target sector:")
for (t, j) in zip(target, tidx)
    order = sortperm(A[:, j]; rev = true)[1:5]
    @printf("  %s (total intermediate input coefficient = %.3f):\n", t, sum(A[:, j]))
    for i in order
        @printf("      %-14s %.4f\n", AGGCOM[i], A[i, j])
    end
end

println("\nTop-5 output customers for each target sector (share of the sector's own domestic intermediate sales):")
for (t, i) in zip(target, tidx)
    total_sales = sum(Zagg[i, :])
    row_share = Zagg[i, :] ./ total_sales
    order = sortperm(row_share; rev = true)[1:5]
    @printf("  %s (total domestic intermediate sales = %.0f):\n", t, total_sales)
    for j in order
        @printf("      %-14s %.4f\n", AGGCOM[j], row_share[j])
    end
end
