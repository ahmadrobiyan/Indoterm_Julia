# Step-2 probe: E_xprim! rowscale vs xprim[2,5] column entries at the stall.
#
# Step-1 finding (file read, no solve): xprim's own column entries come from
# rows where xprim is the differentiated variable — biggest is E_xprim!'s
# identity (dF/dxprim = 1). Shares (ALPHA_FAC) shape the DEMAND columns, not
# xprim's own. So a quiet xprim column (0.0005) with raw entry 1 implies a
# rowscale ~2000 on its owning row — the frozen ruler muting it, not economics.
# This probe prints, for every row touching free xprim[2,5]: row id, JuMP name,
# rowscale, raw entry, scaled entry. Confirm-or-kill in one staged run.
#
# Run: julia --project=. test/scratch/_probe_xprim_rs.jl (verbose=ON)
# Log: logs/probe_xprim_rs_2026-09-19.log
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
println("stall residual: $(r.residual)")
flush(stdout)
println("PRIM[2,5]=$(parent(params["PRIM"])[2,5])")
flush(stdout)

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
connames = String[]
for (Ftype, S) in JuMP.list_of_constraint_types(m)
    Ftype <: JuMP.VariableRef && continue
    for con in JuMP.all_constraints(m, Ftype, S)
        co = JuMP.constraint_object(con)
        r2 = co.set isa MOI.EqualTo ? co.set.value : 0.0
        MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r2))
        push!(rhs, r2)
        push!(connames, try JuMP.name(con) catch; "?" end)
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
xp_col = findfirst(v -> JuMP.name(v) == "xprim[2,5]", freeref)
println("xprim[2,5] free col: $xp_col")
flush(stdout)
for i in 1:length(rhs)
    touch = [jcc[k] for k in eachindex(jrr) if jrr[k] == i]
    if xp_col !== nothing && xp_col in touch
        raw = [v0[k] for k in eachindex(jrr) if jrr[k] == i && jcc[k] == xp_col]
        println("row $i $(connames[i]) rs=$(rs[i]) raw=$(raw) scaled=$(raw[1] / rs[i])")
        flush(stdout)
    end
end
println("Done.")
