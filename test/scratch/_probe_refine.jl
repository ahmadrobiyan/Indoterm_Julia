# Refinement probe: can one extra polish step rescue the exact LU direction?
#
# What we know: at t=0.0625 the direct LU solve starts sharp (1e-10) then
# softens to ~1e-6 and gets rejected by the 1e-6 gate. Every later step then
# falls to the CGNR backup, which never converges (200/200 cap-hits).
# H8 dead (fresh ruler same answer). H7 dead (no single row).
#
# Idea to test: reuse the SAME LU factorization for one cleanup pass —
# measure the leftover error, solve for the correction, add it back.
# If the error drops back below 1e-6, the exact path reopens with no new
# machinery. If it stays stuck, the matrix itself is too sick and we go to
# the full condition check.
#
# Run: julia --project=. test/scratch/_probe_refine.jl (verbose=ON)
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra
import JuMP.MOI as MOI

agg6, params = cached_pipeline(6)
COAL = findfirst(==("Coal"), AGGCOM)

(m, vars) = build_model_full!(agg6, params)
initialize_model!(m, vars; bmk_levels=benchmark_levels(params), numeraire=:exrate)
apply_swaps!(m, vars, collect(COALPRICE_SWAPS); verbose=true)
r0 = solve_newton!(m, vars; maxit=5, tol=1e-8, verbose=true)
@assert r0.residual <= 1e-8

nm, idx, vr = IndotermJulia._shock_ref(vars, ("fpexp_d", COAL))
b = JuMP.fix_value(vr); tg = logpct(-20.0)
for t_ in (0.015625, 0.03125, 0.046875)
    JuMP.fix(vr, b + t_ * (tg - b); force=true)
    r = solve_newton!(m, vars; maxit=30, tol=1e-8, verbose=true)
    @assert r.residual <= 1e-8 "staging failed at $t_"
end
JuMP.fix(vr, b + 0.0625 * (tg - b); force=true)

# Build the scaled system at the stuck point, exactly like the solver does.
allv = JuMP.all_variables(m)
allidx = [JuMP.index(v) for v in allv]
pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))
freeid = Dict{Int64,Int}(); freeref = JuMP.VariableRef[]
for (vnm, v) in vars, vr2 in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr2) && continue
    iv = vr2.index.value; haskey(freeid, iv) && continue
    push!(freeref, vr2); freeid[iv] = length(freeref)
end
nlmodel = MOI.Nonlinear.Model()
rhs = Float64[]
for (Ftype, S) in JuMP.list_of_constraint_types(m)
    Ftype <: JuMP.VariableRef && continue
    for con in JuMP.all_constraints(m, Ftype, S)
        co = JuMP.constraint_object(con)
        r2 = co.set isa MOI.EqualTo ? co.set.value :
             co.set isa MOI.LessThan ? co.set.upper :
             co.set isa MOI.GreaterThan ? co.set.lower : 0.0
        MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r2))
        push!(rhs, r2)
    end
end
evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
MOI.initialize(evaluator, [:Grad, :Jac])
st = MOI.jacobian_structure(evaluator)
jrows_all = getindex.(st, 1); jcols_all = getindex.(st, 2)
Jval = zeros(length(jrows_all))
freecols = [pos[JuMP.index(v).value] for v in freeref]
col2free = zeros(Int, length(allv))
for (i, c) in enumerate(freecols); col2free[c] = i; end
keep = [col2free[c] != 0 for c in jcols_all]
jr = jrows_all[keep]; jc = [col2free[c] for c in jcols_all[keep]]
xs = Float64[(JuMP.is_fixed(v) ? JuMP.fix_value(v) :
    (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v))) for v in allv]
MOI.eval_constraint_jacobian(evaluator, Jval, xs)
v0 = Jval[keep]
mx = zeros(Float64, length(rhs))
@inbounds for k in eachindex(jr)
    a = abs(v0[k]); a > mx[jr[k]] && (mx[jr[k]] = a)
end
rowscale = [mx[i] > 1e-12 ? mx[i] : 1.0 for i in 1:length(rhs)]
Jbase = sparse(jr, jc, [v0[k] / rowscale[jr[k]] for k in eachindex(jr)], length(rhs), length(freeref))
g = zeros(length(rhs))
MOI.eval_constraint(evaluator, g, xs)
@inbounds for i in eachindex(g); g[i] -= rhs[i]; end
gs = g ./ rowscale
# One LU factorization, then measure + one cleanup pass.
println("factoring J ($(size(Jbase,1))x$(size(Jbase,2)), nnz=$(nnz(Jbase)))...")
flush(stdout)
F = lu(Jbase)
dy = F \ (-gs)
lin1 = norm(Jbase * dy .+ gs) / max(norm(gs), eps())
println("plain LU error = $(lin1)")
flush(stdout)
# Cleanup pass: solve for the leftover and add it back.
resid = Jbase * dy .+ gs
corr = F \ (-resid)
dy2 = dy .+ corr
lin2 = norm(Jbase * dy2 .+ gs) / max(norm(gs), eps())
println("after one cleanup pass error = $(lin2)")
flush(stdout)
println("gate is 1e-6 → " * (lin2 < 1e-6 ? "REOPENS exact path" : "still rejected"))
println("Done.")
