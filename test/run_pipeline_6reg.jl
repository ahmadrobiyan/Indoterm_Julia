push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia

println("="^60)
println("INDOTERM PIPELINE (STEPS 0-4) @ 25×6 REGION-REDUCED")
println("="^60)

# ── Load raw CSV data ────────────────────────────────────────────────────────
println("\n--- Loading raw data ---")
nat = read_national_data()
regsup = read_regsupp_data()
distgone = read_distgone_data()

# ── Steps 2a–2f: regional build ──────────────────────────────────────────────
println("\n--- build_reg0! … build_premod! ---")

let t = time(); global reg0_r = build_reg0!(nat, regsup); println("  build_reg0!  : $(round(time()-t,digits=2))s"); end
let t = time(); global reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone); println("  build_reg1!  : $(round(time()-t,digits=2))s"); end
let t = time(); global reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag); println("  build_reg2!  : $(round(time()-t,digits=2))s"); end
let t = time(); global ras_r = ras_balance!(reg2_r.reg2); println("  ras_balance! : $(round(time()-t,digits=2))s"); end
let t = time(); global pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone); println("  build_pstras!: $(round(time()-t,digits=2))s"); end
let t = time(); global prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast); println("  build_premod!: $(round(time()-t,digits=2))s"); end
@assert haskey(prem_r.premod, "MAKE") "MAKE missing from premod"
println("premod OK")

# ── Step 3: Aggregation 185→25 (regions still 34) ────────────────────────────
println("\n--- aggregate_model! (185→25 sectors) ---")
agg_r = aggregate_model!(prem_r.premod)
@assert size(agg_r.agg["MAKE"]) == (25,25,34) "MAKE should be (25,25,34)"
@assert size(agg_r.agg["TRAD"]) == (25,2,34,34) "TRAD should be (25,2,34,34)"
@assert size(agg_r.agg["TMAR"]) == (25,2,9,34,34) "TMAR should be (25,2,9,34,34)"
println("34-region agg shapes verified ✓")

# ── Step 3b: Region collapse 34→6 (Phase 1) ──────────────────────────────────
println("\n--- aggregate_regions! (34→6 island groups) ---")
println("REG6 = $(REG6)")
@assert length(REG_MAP_34_to_6) == 34 "REG_MAP must cover all 34 provinces"
@assert sort(unique(REG_MAP_34_to_6)) == collect(1:6) "REG_MAP must map onto 1:6"
agg6 = aggregate_regions!(agg_r.agg, REG_MAP_34_to_6, 6)

println("MAKE: size=$(size(agg6["MAKE"]))")
println("TRAD: size=$(size(agg6["TRAD"]))")
println("TMAR: size=$(size(agg6["TMAR"]))")
println("MARS: size=$(size(agg6["MARS"]))")
@assert size(agg6["MAKE"]) == (25,25,6) "MAKE should be (25,25,6)"
@assert size(agg6["TRAD"]) == (25,2,6,6) "TRAD should be (25,2,6,6)"
@assert size(agg6["TMAR"]) == (25,2,9,6,6) "TMAR should be (25,2,9,6,6)"
@assert size(agg6["MARS"]) == (9,6,6,6) "MARS should be (9,6,6,6)"
@assert size(agg6["1LAB"]) == (25,size(agg6["1LAB"],2),6) "1LAB region axis should be 6"
@assert size(agg6["1CAP"]) == (25,6) "1CAP should be (25,6)"
@assert size(agg6["BSMR"]) == (25,2,29,6) "BSMR should be (25,2,29,6)"
@assert size(agg6["2PUR"]) == (25,25,6) "2PUR should be (25,25,6)"
println("6-region agg shapes verified ✓")

# ── Conservation check: national totals must be preserved by the collapse ─────
for k in ("MAKE", "TRAD", "TMAR", "MARS", "1LAB", "1CAP", "1LND", "BSMR", "UTAX", "2PUR")
    haskey(agg_r.agg, k) || continue
    s34 = sum(parent(agg_r.agg[k]))
    s6  = sum(agg6[k])
    rel = abs(s34) > 0 ? abs(s6 - s34) / abs(s34) : abs(s6 - s34)
    @assert rel < 1e-8 "$k grand total not conserved: 34reg=$s34 6reg=$s6 (rel=$rel)"
end
println("Flow grand-totals conserved across 34→6 collapse ✓")

# ── Step 4: Derived parameters at 25×6 ───────────────────────────────────────
println("\n--- prepare_parameters! @ 25×6 ---")
params = prepare_parameters!(agg6)
println("Derived $(length(params)) parameters")
@assert size(params["DELIVRD"]) == (25,2,6,6) "DELIVRD should be (25,2,6,6)"
@assert size(params["PUR_CS"]) == (29,6) "PUR_CS should be (29,6)"
@assert size(params["SRCSHR"]) == (25,2,29,6) "SRCSHR should be (25,2,29,6)"
@assert size(params["COSTMAT"]) == (25,7,6) "COSTMAT should be (25,7,6)"
@assert size(params["GDPEXPSUM"]) == (6,9) "GDPEXPSUM should be (6,9)"
println("25×6 derived parameter shapes verified ✓")

println("\n"^2)
println("="^60)
println("PIPELINE @ 25×6 COMPLETED SUCCESSFULLY")
println("="^60)
