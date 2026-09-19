# Restricted condition estimate around the stall rows (H5 gate, scoped).
#
# Context: Options A and B both dead (tautology / yardstick-breaker). Scaling
# is the flashlight, not the disease. The remaining H5 question — scoped, not
# whole-system: is the ~10-row submatrix around rows 37426/37276/36784 (the
# [2,5] factor-nest plateau) itself ill-conditioned, or does the stall come
# from the coupling to the other 83k rows?
#
# Method: at the staged stall point, extract the dense sub-Jacobian on the
# plateau rows x their touching free columns, SVD it, report σmax/σmin.
# Small (<= ~15 cols) so dense SVD is trivial. Compares against the same
# construction at the benchmark (control): if the stall submatrix is orders
# more ill-conditioned than the control, H5-confirmed-local; if both are
# mild, H5-falsified and the stall is coupling/global, not local.
#
# Run: julia --project=. test/scratch/_probe_h5_local.jl (verbose=ON)
# Log: logs/probe_h5_local_2026-09-19.log
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra
import JuMP.MOI as MOI

agg6, params = cached_pipeline(6)
COAL = findfirst(==("Coal"), AGGCOM)
(m, vars) = build_model_full!(agg6, params)
initialize_model!(m, vars; bmk_levels=benchmark_levels(params), numeraire=:exrate)
apply_swaps!(m, vars, collect(COALPRICE_SWAPS); verbose=false)
r0 = solve_newton!(m, vars; maxit=5, tol=1e-8, verbose=false)
@assert r0.residual <= 1e-8

function build_system()
    allv = JuMP.all_variables(m)
    allidx = [JuMP.index(v) for v in allv]
    pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))
    freeid = Dict{Int64,Int}()
    freeref = JuMP.VariableRef[]
    for (vnm, v) in vars, vr2 in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr2) && continue
        iv = vr2.index.value
        haskey(freeid, iv) && continue
        push!(freeref, vr2)
        freeid[iv] = length(freeref)
    end
    nlmodel = MOI.Nonlinear.Model()
    rhs = Float64[]
    for (Ftype, S) in JuMP.list_of_constraint_types(m)
        Ftype <: JuMP.VariableRef && continue
        for con in JuMP.all_constraints(m, Ftype, S)
            co = JuMP.constraint_object(con)
            r2 = co.set isa MOI.EqualTo ? co.set.value : 0.0
            MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r2))
            push!(rhs, r2)
        end
    end
    evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
    MOI.initialize(evaluator, [:Grad, :Jac])
    st = MOI.jacobian_structure(evaluator)
    srows = getindex.(st, 1)
    scols = getindex.(st, 2)
    Jval = zeros(length(srows))
    freecols = [pos[JuMP.index(v).value] for v in freeref]
    col2free = zeros(Int, length(allv))
    for (i, c) in enumerate(freecols)
        col2free[c] = i
    end
    keep = [col2free[c] != 0 for c in scols]
    jrr = srows[keep]
    jcc = [col2free[c] for c in scols[keep]]
    startval(v) = JuMP.is_fixed(v) ? JuMP.fix_value(v) : begin
        sv = JuMP.start_value(v)
        sv === nothing ? 1.0 : sv
    end
    xs = Float64[startval(v) for v in allv]
    MOI.eval_constraint_jacobian(evaluator, Jval, xs)
    v0 = Jval[keep]
    mx = zeros(Float64, length(rhs))
    @inbounds for k in eachindex(jrr)
        a = abs(v0[k])
        a > mx[jrr[k]] && (mx[jrr[k]] = a)
    end
    rs = [mx[i] > 1e-12 ? mx[i] : 1.0 for i in 1:length(rhs)]
    Jb = sparse(jrr, jcc, [v0[k] / rs[jrr[k]] for k in eachindex(jrr)],
                length(rhs), length(freeref))
    return Jb, freeref, rs
end

function subcond(Jb, rows, label)
    # touching columns of the given rows (single O(nnz) pass):
    touch = Set{Int}()
    rowset = Set(rows)
    for j in 1:size(Jb, 2)
        for k in Jb.colptr[j]:(Jb.colptr[j+1]-1)
            if Jb.rowval[k] in rowset
                push!(touch, j)
                break
            end
        end
    end
    tc = sort(collect(touch))
    S = Matrix(Jb[rows, tc])
    sv = svdvals(S)
    println("$label: $(length(rows)) rows x $(length(tc)) cols  " *
        "σmax=$(round(maximum(sv); sigdigits=4))  σmin=$(round(minimum(sv); sigdigits=4))  " *
        "κ=$(round(maximum(sv)/max(minimum(sv),eps()); sigdigits=4))")
    flush(stdout)
end

# control at benchmark
Jb0, fr0, _ = build_system()
println("CONTROL (benchmark point):")
flush(stdout)
# locate E_xprim[2,5]-area rows by touching xprim[2,5] column
xp0 = findfirst(v -> JuMP.name(v) == "xprim[2,5]", fr0)
rows0 = Int[]
if xp0 !== nothing
    for j in [xp0]
        for k in Jb0.colptr[j]:(Jb0.colptr[j+1]-1)
            push!(rows0, Jb0.rowval[k])
        end
    end
    rows0 = unique(rows0)
end
subcond(Jb0, rows0, "control")

# stage to stall
nm, idx, vr = IndotermJulia._shock_ref(vars, ("fpexp_d", COAL))
b = JuMP.fix_value(vr)
tg = logpct(-20.0)
for t_ in (0.015625, 0.03125, 0.046875)
    JuMP.fix(vr, b + t_ * (tg - b); force=true)
    r = solve_newton!(m, vars; maxit=30, tol=1e-8, verbose=false)
    @assert r.residual <= 1e-8 "staging failed at $t_"
end
JuMP.fix(vr, b + 0.0625 * (tg - b); force=true)
r = solve_newton!(m, vars; maxit=12, tol=1e-8, verbose=false)
println("stall-point: status=$(r.status) residual=$(r.residual)")
flush(stdout)

Jb1, fr1, _ = build_system()
xp1 = findfirst(v -> JuMP.name(v) == "xprim[2,5]", fr1)
rows1 = Int[]
if xp1 !== nothing
    for k in Jb1.colptr[xp1]:(Jb1.colptr[xp1+1]-1)
        push!(rows1, Jb1.rowval[k])
    end
    rows1 = unique(rows1)
end
println("STALL (t=0.0625 point):")
flush(stdout)
subcond(Jb1, rows1, "stall")
println("verdict: κ_stall >> κ_control → H5-confirmed-local; " *
    "κ similar → H5-falsified, stall is coupling/global.")
println("Done.")
