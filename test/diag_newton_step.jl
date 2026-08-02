# One-shot Newton step diagnostic: build, shock, one LU, report d stats.
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP, SparseArrays, LinearAlgebra, Statistics

function main()
    println("pipeline..."); flush(stdout)
    nat = read_national_data(); regsup = read_regsupp_data(); distgone = read_distgone_data()
    reg0_r = build_reg0!(nat, regsup)
    reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
    reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
    ras_r  = ras_balance!(reg2_r.reg2)
    pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone)
    prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast)
    agg6 = aggregate_regions!(aggregate_model!(prem_r.premod).agg, REG_MAP_34_to_6, 6)
    params = prepare_parameters!(agg6)

    println("build..."); flush(stdout)
    m, vars = build_model_full!(agg6, params)
    bmk = benchmark_levels(params)
    initialize_model!(m, vars; bmk_levels=bmk)

    # Pin orphan
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        # nothing
    end
    res0 = solve_newton!(m, vars; maxit=1, tol=1e-8, verbose=false)
    println("bmk ||F||=$(res0.residual)"); flush(stdout)

    JuMP.fix(vars["blabnat"], 0.9997; force=true)
    # seed alab_o
    alab = vars["alab_o"]
    for idx in eachindex(alab)
        JuMP.set_start_value(alab[idx], 0.9997)
    end

    # Manual one-step via solve_newton verbose
    println("newton..."); flush(stdout)
    t0 = time()
    # Use undamped path only — temporarily simpler solve
    res = solve_newton!(m, vars; maxit=5, tol=1e-7, verbose=true)
    println("done t=$(round(time()-t0;digits=1))s status=$(res.status) ||F||=$(res.residual)")
    println("alab_o mean=$(mean(res.values["alab_o"]))")
    println("phi=$(res.values["phi"])")
end
main()
