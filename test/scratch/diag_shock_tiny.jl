push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP, Statistics

const LOG = open(joinpath(@__DIR__, "..", "logs", "diag_shock_tiny.log"), "w")
logp(msg) = (println(LOG, msg); flush(LOG); println(msg); flush(stdout))

function main()
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

    logp("Building...")
    t0 = time(); m, vars = build_model_full!(agg6, params)
    logp("Build: $(round(time()-t0; digits=1))s")
    bmk = benchmark_levels(params)
    initialize_model!(m, vars; bmk_levels=bmk)

    logp("\n-- Benchmark --")
    t0 = time()
    res0 = solve_newton!(m, vars; maxit=1, tol=1e-8, verbose=true)
    logp("status=$(res0.status) ||F||=$(res0.residual) t=$(round(time()-t0;digits=1))s")

    target = 0.9997
    JuMP.fix(vars["blabnat"], target; force=true)

    alab = vars["alab_o"]; blab_d = vars["blab_d"]; blab = vars["blab"]
    for i in axes(alab, 1), d in axes(alab, 2)
        bd = JuMP.is_fixed(blab_d[i]) ? JuMP.fix_value(blab_d[i]) : JuMP.start_value(blab_d[i])
        b  = JuMP.is_fixed(blab[i,d]) ? JuMP.fix_value(blab[i,d]) : JuMP.start_value(blab[i,d])
        JuMP.set_start_value(alab[i,d], target * bd * b)
    end
    logp("Seeded alab_o = $target * blab_d * blab")

    logp("\n-- Direct Newton blabnat=$target --")
    t0 = time()
    @time res1 = solve_newton!(m, vars; maxit=10, tol=1e-7, verbose=true)
    logp("status=$(res1.status) ||F||=$(res1.residual) t=$(round(time()-t0;digits=1))s")
    a = res1.values["alab_o"]
    logp("alab_o mean=$(mean(a)) min=$(minimum(a)) max=$(maximum(a))")
    logp("phi=$(res1.values["phi"])")

    if res1.status == :converged
        logp("\n-- Euler to 0.97 (20 steps) --")
        t0 = time()
        res2 = euler_continuation!(m, vars, vars["blabnat"], 0.97, 20; tol=1e-6, verbose=true)
        logp("status=$(res2.status) ||F||=$(res2.residual) t=$(round(time()-t0;digits=1))s")
        for nm in ("xlab_o", "plab_o", "xprim", "pcap", "xmake", "phi")
            haskey(res2.values, nm) || continue
            v = res2.values[nm]; bl = get(bmk, nm, nothing)
            pct = bl === nothing ? (v ./ 1.0 .- 1) .* 100 : (v ./ bl .- 1) .* 100
            logp("  $nm avg|%|=$(round(mean(abs.(pct));digits=4)) min=$(round(minimum(pct);digits=4)) max=$(round(maximum(pct);digits=4))")
        end
    end
    close(LOG)
end

main()
