push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP

nat = read_national_data()
regsup = read_regsupp_data()
distgone = read_distgone_data()
reg0_r = build_reg0!(nat, regsup)
reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
ras_r  = ras_balance!(reg2_r.reg2)
pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone)
prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast)
agg_r  = aggregate_model!(prem_r.premod)
params = prepare_parameters!(agg_r.agg)

m, vars = build_model_full!(agg_r.agg, params)

println("--- lookup dict population checks (post build_model_full!) ---")
for nm in [:V1LAB_idx, :TRADMAR_idx, :SUPPMAR_idx, :SUPPMAR_D_idx, :TRADE_idx,
           :PUR_src_idx, :TAX_idx, :PUR_idx, :INVEST_C_idx, :USE_IS_idx, :USE_usc_idx, :STOCKS_idx]
    d = getfield(IndotermJulia, nm)
    nz = count(x -> abs(x) > 1e-10, values(d))
    println("  $nm: length=$(length(d)) nonzero=$nz sum=$(sum(values(d)))")
end

println("\n--- STOK / STOCKS pass-through check ---")
println("agg has STOK: ", haskey(agg_r.agg, "STOK"))
if haskey(agg_r.agg, "STOK")
    println("agg[STOK] sum: ", sum(agg_r.agg["STOK"]))
end
println("params[STOCKS] sum: ", sum(params["STOCKS"]))

initialize_model!(m, vars)

println("\n--- plab_o / wlab_o / wprim constraint RHS spot-check ---")
LAB_O = parent(params["LAB_O"])
na, nr = size(LAB_O)
checked = 0
for i in 1:na, d in 1:nr
    global checked
    if LAB_O[i,d] > 1e-10
        checked += 1
        if checked <= 5
            println("  i=$i d=$d LAB_O=$(LAB_O[i,d])")
        end
    end
end
println("Total (i,d) with LAB_O>1e-10: $checked")
println("Done.")
