push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP, Statistics

function run_shock()
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
agg6 = aggregate_regions!(agg_r.agg, REG_MAP_34_to_6, 6)
params = prepare_parameters!(agg6)

println("Building model...")
t0 = time(); m, vars = build_model_full!(agg6, params)
tbuild = time() - t0
println("Build: ", round(tbuild; digits=1), "s")

bmk = benchmark_levels(params)
initialize_model!(m, vars; bmk_levels=bmk)

# Helper: % change from benchmark
function pct_change(arr, nm)
    bl = get(bmk, nm, nothing)
    if bl === nothing
        return (arr .- 1.0) * 100
    else
        return (arr ./ bl .- 1.0) * 100
    end
end

# Run continuation: blabnat from 1.0 to 0.97 (3% wage-bill reduction)
blabnat = vars["blabnat"]
println("\n-- Euler continuation: blabnat 1.0 -> 0.97 (50 steps) --")
t0 = time()
res = euler_continuation!(m, vars, blabnat, 0.97, 50; tol=1e-7, verbose=true)
t_cont = time() - t0
println("Continuation: ", round(t_cont; digits=1), "s @ ", round(t_cont/51; digits=2), "s/step")
println("Final scaled ||F||_inf = ", round(res.residual; digits=6))
println("Raw ||F||_inf = ", round(res.raw_residual; digits=6))

vals = res.values

println("\n-- Economic results (blabnat = -3%) --")
for (nm, desc) in [("alab_o","Labour augment. tech"), ("xlab_o","Labour demand"),
       ("plab_o","Labour composite price"), ("xprim","Primary factor demand"),
       ("pprim","Primary factor price"), ("xcap","Capital demand (FIXED)"),
       ("pcap","Capital price"), ("xlnd","Land demand (FIXED)"),
       ("plnd","Land rental price"), ("xmake","Output index"),
       ("xint","Intermediate demand"), ("ppur","Purchase price")]
    haskey(vals, nm) || continue
    v = vals[nm]
    pct = pct_change(v, nm)
    if v isa AbstractArray
        println("  $(rpad(nm,10)) $(rpad(desc,25))  min=$(round(minimum(pct);digits=3))%  " *
                "max=$(round(maximum(pct);digits=3))%  avg|=$(round(mean(abs.(pct));digits=3))%")
    else
        println("  $(rpad(nm,10)) $(rpad(desc,25))  value=$(round(pct;digits=3))%")
    end
end

# Count significant movements
n1 = 0; n5 = 0; total = 0
for (nm, v) in vals
    v isa AbstractArray || continue
    bl = get(bmk, nm, nothing)
    for idx in eachindex(v)
        total += 1
        ref = bl === nothing ? 1.0 : (bl isa AbstractArray ? bl[idx] : bl)
        chg = abs(v[idx] / ref - 1.0) * 100
        chg > 1.0 && (n1 += 1)
        chg > 5.0 && (n5 += 1)
    end
end
println("\n-- Summary --")
println("Variables moved >1%: $n1 / $total")
println("Variables moved >5%: $n5 / $total")
println("Total runtime: ", round(tbuild + t_cont; digits=1), "s")

return vals
end

result = run_shock()
