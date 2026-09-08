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
    res0 = solve_newton!(m, vars; maxit=1, tol=1e-8, verbose=false)
    println("bmk ||F||=$(res0.residual)"); flush(stdout)

    # ── shock blabnat ────────────────────────────────────────────────────────
    target = 0.9997
    JuMP.fix(vars["blabnat"], target; force=true)

    # Properly seed alab_o via the decomposition: alab_o = blabnat * blab_d * blab
    alab = vars["alab_o"]; blab_d = vars["blab_d"]; blab = vars["blab"]
    for i in axes(alab,1), d in axes(alab,2)
        bd = JuMP.is_fixed(blab_d[i]) ? JuMP.fix_value(blab_d[i]) : JuMP.start_value(blab_d[i])
        b  = JuMP.is_fixed(blab[i,d]) ? JuMP.fix_value(blab[i,d]) : JuMP.start_value(blab[i,d])
        JuMP.set_start_value(alab[i,d], target * bd * b)
    end
    println("seeded alab_o via decomposition"); flush(stdout)

    # ── build the nonlinear evaluator manually for one-step diagnostic ───────
    freeid = Dict{Int64,Int}(); freeref = VariableRef[]
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) && continue
        iv = vr.index.value; haskey(freeid, iv) && continue
        push!(freeref, vr); freeid[iv] = length(freeref)
    end
    seen = falses(length(freeref))
    deadcons = Any[]
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            buf = Int[]
            f = JuMP.constraint_object(con).func
            if f isa VariableRef
                id = get(freeid, f.index.value, 0); id != 0 && push!(buf, id)
            elseif f isa GenericAffExpr
                for (v, _) in f.terms
                    id = get(freeid, v.index.value, 0); id != 0 && push!(buf, id)
                end
            elseif f isa GenericQuadExpr
                for (v, _) in f.aff.terms
                    id = get(freeid, v.index.value, 0); id != 0 && push!(buf, id)
                end
                for (p, _) in f.terms
                    id = get(freeid, p.a.index.value, 0); id != 0 && push!(buf, id)
                    id = get(freeid, p.b.index.value, 0); id != 0 && push!(buf, id)
                end
            elseif f isa GenericNonlinearExpr
                function coll(nf)
                    if nf isa VariableRef
                        id = get(freeid, nf.index.value, 0); id != 0 && push!(buf, id)
                    elseif nf isa GenericNonlinearExpr
                        for a in nf.args; coll(a); end
                    else
                    end
                end; coll(f)
            end
            if isempty(buf)
                push!(deadcons, con)
            else
                for id in buf; seen[id] = true; end
            end
        end
    end
    orphans = [freeref[i] for i in eachindex(freeref) if !seen[i]]
    for vr in orphans
        sv = JuMP.start_value(vr)
        JuMP.fix(vr, sv === nothing ? 0.0 : sv; force=true)
    end
    for con in deadcons; JuMP.delete(m, con); end
    println("  pinned $(length(orphans)) orphans, deleted $(length(deadcons)) dead")

    nowfree = VariableRef[]
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) || push!(nowfree, vr)
    end
    N = length(nowfree)
    allv = JuMP.all_variables(m)
    allidx = [JuMP.index(v) for v in allv]
    pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))
    freecols = [pos[JuMP.index(v).value] for v in nowfree]

    nlmodel = MOI.Nonlinear.Model()
    rhs = Float64[]
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            co = JuMP.constraint_object(con)
            r = co.set isa MOI.EqualTo ? co.set.value :
                co.set isa MOI.LessThan ? co.set.upper :
                co.set isa MOI.GreaterThan ? co.set.lower : 0.0
            MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r))
            push!(rhs, r)
        end
    end
    NC = length(rhs)
    println("  system: $NC eqs x $N vars")
    evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
    MOI.initialize(evaluator, [:Grad, :Jac])

    x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
         (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]

    # Evaluate F(x)
    g = zeros(NC)
    MOI.eval_constraint(evaluator, g, x)
    g .-= rhs
    g_inf = norm(g, Inf)
    println("  raw residual ||F||_inf = $g_inf")

    # Find which constraint has the max residual
    maxk = argmax(abs.(g))
    # Print the constraint (hard to get the name, but we can try)
    ci = 0
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            ci += 1
            if ci == maxk
                println("  max residual at constraint #$maxk ($(name(con))) = $(g[maxk])")
            end
        end
    end

    # Evaluate Jacobian
    Jval = zeros(length(MOI.jacobian_structure(evaluator)))
    MOI.eval_constraint_jacobian(evaluator, Jval, x)
    st = MOI.jacobian_structure(evaluator)
    jrows_all = getindex.(st, 1); jcols_all = getindex.(st, 2)
    col2free = zeros(Int, length(allv))
    for (i, c) in enumerate(freecols); col2free[c] = i; end
    keep = [col2free[c] != 0 for c in jcols_all]
    jr = jrows_all[keep]; jc = [col2free[c] for c in jcols_all[keep]]

    # Row equilibration
    rowscale = ones(NC)
    let v0 = Jval[keep], mx = zeros(NC)
        for k in eachindex(jr)
            a = abs(v0[k]); a > mx[jr[k]] && (mx[jr[k]] = a)
        end
        for i in 1:NC; rowscale[i] = mx[i] > 1e-12 ? mx[i] : 1.0; end
    end
    gs = g ./ rowscale
    println("  scaled ||F||_inf = $(norm(gs, Inf))")

    # Build column-scaled J
    vals = Jval[keep]
    J = sparse(jr, jc, [vals[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, N)
    colnrm = [norm(J[:, j]) for j in 1:N]
    println("  column norms: min=$(minimum(colnrm)) max=$(maximum(colnrm))")

    # LU solve without LM
    println("\n── LU solve (no damping) ──"); flush(stdout)
    d0 = try -(lu(J) \ gs) catch; zeros(N); end
    if any(!isfinite, d0); d0 = zeros(N); end
    maxdi, maxd = argmax(abs.(d0)), maximum(abs, d0)
    println("  max|d| = $maxd at free var #$maxdi: $(JuMP.name(nowfree[maxdi]))")
    # Top 10
    top10 = sortperm(abs.(d0), rev=true)[1:min(10, N)]
    for ti in top10
        println("    $(rpad(JuMP.name(nowfree[ti]), 25)) d=$(d0[ti])")
    end

    # LU solve with LM=1e-2
    println("\n── LU solve (LM μ=1e-2) ──"); flush(stdout)
    colnrm_s = copy(colnrm); colnrm_s[colnrm.<1e-10] .= 1.0
    Jsc = J * spdiagm(1.0 ./ colnrm_s)
    d2 = try -(lu(Jsc + 1e-2*I) \ gs) catch; -(qr(Jsc + 1e-2*I) \ gs); end
    d2 ./= colnrm_s
    for i in 1:N; colnrm[i] < 1e-10 && (d2[i] = 0.0); end
    maxd2i, maxd2 = argmax(abs.(d2)), maximum(abs, d2)
    println("  max|d| = $maxd2 at free var #$maxd2i: $(JuMP.name(nowfree[maxd2i]))")
    top10 = sortperm(abs.(d2), rev=true)[1:min(10, N)]
    for ti in top10
        println("    $(rpad(JuMP.name(nowfree[ti]), 25)) d=$(d2[ti])")
    end

    # CGNR solve
    println("\n── CGNR solve ──"); flush(stdout)
    function cgnr_solve(J, b; maxit=200, tol=1e-10)
        n = size(J, 1); d = zeros(n); r = copy(b)
        Jt = J'; z = Jt * r; p = copy(z)
        γ = dot(z, z); γ0 = γ
        for it in 1:maxit
            Jp = J * p
            α = γ / dot(Jp, Jp)
            for i in 1:n; d[i] += α * p[i]; r[i] -= α * Jp[i]; end
            z .= Jt * r
            γnext = dot(z, z)
            sqrt(γnext) < tol * sqrt(γ0) && break
            β = γnext / γ
            for i in 1:n; p[i] = z[i] + β * p[i]; end
            γ = γnext
        end
        d, sqrt(γ) / sqrt(γ0)
    end

    dcgnr, relres = cgnr_solve(J, -gs)
    maxdci, maxdc = argmax(abs.(dcgnr)), maximum(abs, dcgnr)
    println("  CGNR rel_res = $relres")
    println("  max|d| = $maxdc at free var #$maxdci: $(JuMP.name(nowfree[maxdci]))")
    top10 = sortperm(abs.(dcgnr), rev=true)[1:min(10, N)]
    for ti in top10
        println("    $(rpad(JuMP.name(nowfree[ti]), 25)) d=$(dcgnr[ti])")
    end

    # Check if CGNR + α=1 is a descent direction
    Jdcgnr = J * dcgnr
    pred_reduction = -dot(gs, Jdcgnr) / (norm(gs)^2 + 1e-30)
    println("  predicted linear reduction factor = $pred_reduction (should be ≈1)")

    # Try α=1 step from CGNR
    x1 = copy(x)
    for i in 1:N; x1[freecols[i]] = x[freecols[i]] + dcgnr[i]; end
    g1 = zeros(NC)
    MOI.eval_constraint(evaluator, g1, x1)
    g1 .-= rhs
    g1s = g1 ./ rowscale
    println("  CGNR actual ||F|| after α=1 step = $(norm(g1s, Inf))  (vs $(norm(gs, Inf)))")

    # Try α=1 step from LU with LM
    x2 = copy(x)
    for i in 1:N; x2[freecols[i]] = x[freecols[i]] + d2[i]; end
    g2 = zeros(NC)
    MOI.eval_constraint(evaluator, g2, x2)
    g2 .-= rhs
    g2s = g2 ./ rowscale
    println("  LM(LU) actual ||F|| after α=1 step = $(norm(g2s, Inf))  (vs $(norm(gs, Inf)))")

    println("\n── done ──"); flush(stdout)
end
main()
