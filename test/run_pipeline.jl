push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia

println("="^60)
println("INDOTERM PIPELINE (STEPS 0-3) SMOKE TEST")
println("="^60)

# ── Load raw CSV data ────────────────────────────────────────────────────────
println("\n--- Loading raw data ---")
nat = read_national_data()
regsup = read_regsupp_data()
distgone = read_distgone_data()

println("National headers: $(length(nat))")
println("REGSUPP headers: $(length(regsup))")
println("DISTGONE: keys=$(keys(distgone))")

# ── Step 2a: reg0 ────────────────────────────────────────────────────────────
println("\n--- build_reg0! ---")
reg0_r = build_reg0!(nat, regsup)
println("reg0 keys: $(sort!(collect(keys(reg0_r.reg0))))")
println("elast keys: $(sort!(collect(keys(reg0_r.elast))))")
@assert haskey(reg0_r.reg0, "FACT") "FACT missing from reg0"
println("FACT: size=$(size(reg0_r.reg0["FACT"]))")

# ── Step 2b: reg1 ────────────────────────────────────────────────────────────
println("\n--- build_reg1! ---")
reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
println("reg1 keys: $(sort!(collect(keys(reg1_r.reg1))))")
@assert haskey(reg1_r.reg1, "FACT") "FACT missing from reg1"
println("FACT: size=$(size(reg1_r.reg1["FACT"]))")

# ── Step 2c: reg2 ────────────────────────────────────────────────────────────
println("\n--- build_reg2! ---")
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
println("reg2 keys: $(sort!(collect(keys(reg2_r.reg2))))")
@assert haskey(reg2_r.reg2, "TRAD") "TRAD missing from reg2"
println("TRAD: size=$(size(reg2_r.reg2["TRAD"]))")

# ── Step 2d: RAS ─────────────────────────────────────────────────────────────
println("\n--- ras_balance! ---")
ras_r = ras_balance!(reg2_r.reg2)
println("ras keys: $(sort!(collect(keys(ras_r.ras))))")
@assert haskey(ras_r.ras, "TRAD") "TRAD missing from ras output"

# ── Step 2e: pstras ──────────────────────────────────────────────────────────
println("\n--- build_pstras! ---")
pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone)
println("pstras keys: $(sort!(collect(keys(pstra_r.pstras))))")
@assert haskey(pstra_r.pstras, "MAKE") "MAKE missing from pstras"
println("Converged: $(pstra_r.converged)")

# ── Step 2f: premod ──────────────────────────────────────────────────────────
println("\n--- build_premod! ---")
# Need to check premod signature: (pstras, regsupp, elast)
prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast)
println("premod keys: $(sort!(collect(keys(prem_r.premod))))")
@assert haskey(prem_r.premod, "MAKE") "MAKE missing from premod"

# ── Step 3: Aggregation 185→25 ──────────────────────────────────────────────
println("\n--- aggregate_model! (185→25) ---")
agg_r = aggregate_model!(prem_r.premod)
println("agg keys: $(length(agg_r.agg))")
println("MAKE: size=$(size(agg_r.agg["MAKE"]))")
println("TRAD: size=$(size(agg_r.agg["TRAD"]))")
println("TMAR: size=$(size(agg_r.agg["TMAR"]))")
@assert size(agg_r.agg["MAKE"]) == (25,25,34) "MAKE should be (25,25,34)"
@assert size(agg_r.agg["TRAD"]) == (25,2,34,34) "TRAD should be (25,2,34,34)"
@assert size(agg_r.agg["TMAR"]) == (25,2,9,34,34) "TMAR should be (25,2,9,34,34)"
println("Aggregation shapes verified ✓")

println("\n"^2)
println("="^60)
println("ALL STAGES (0-3) COMPLETED SUCCESSFULLY")
println("="^60)
