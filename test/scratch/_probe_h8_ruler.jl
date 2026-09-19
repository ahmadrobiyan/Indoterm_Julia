# H8 probe: is the plateau a stale-ruler artifact?
#
# The frozen rowscale is fixed at the benchmark start and never updated.
# At t=0.0625 off-benchmark, rows whose true max entry drifted orders are
# measured against a stale yardstick: scaled residual understates (or
# overstates) the raw violation by the drift factor. Row 83070 (rs=1.08e8)
# is the extreme case.
#
# Method: same staging + stall point as the H7 probes, then score the SAME
# iterate twice — once with the benchmark-frozen rowscale (solver convention),
# once with a freshly recomputed rowscale at the stall point — and report both
# ||F||∞ + argmax rows. Verdict:
#   fresh-||F||∞ < tol  → stall was measurement artifact, not a stall.
#   fresh-||F||∞ ~ same on a DIFFERENT row → genuine plateau, correctly
#     measured; back to the H5 condition-estimate gate.
#   fresh-||F||∞ ~ same on the SAME rows → ruler is fine; plateau is real.
#
# Run: julia --project=. test/scratch/_probe_h8_ruler.jl
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

# frozen rowscale = replicate solver convention: computed at the BENCHMARK
# start vector (all vars at their fixed/benchmark values), before any shock.
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
NC = length(rhs)
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

frozen_x = Float64[(JuMP.is_fixed(v) ? JuMP.fix_value(v) :
    (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v))) for v in allv]
MOI.eval_constraint_jacobian(evaluator, Jval, frozen_x)
v0 = Jval[keep]
mx = zeros(Float64, NC)
@inbounds for k in eachindex(jr)
    a = abs(v0[k]); a > mx[jr[k]] && (mx[jr[k]] = a)
end
rowscale_frozen = [mx[i] > 1e-12 ? mx[i] : 1.0 for i in 1:NC]

# ── stage + stall ─────────────────────────────────────────────────────────
nm, idx, vr = IndotermJulia._shock_ref(vars, ("fpexp_d", COAL))
b = JuMP.fix_value(vr); tg = logpct(-20.0)
for t_ in (0.015625, 0.03125, 0.046875)
    JuMP.fix(vr, b + t_ * (tg - b); force=true)
    r = solve_newton!(m, vars; maxit=30, tol=1e-8, verbose=false)
    @assert r.residual <= 1e-8 "staging failed at $t_"
end
JuMP.fix(vr, b + 0.0625 * (tg - b); force=true)
r = solve_newton!(m, vars; maxit=12, tol=1e-8, verbose=false)
println("stall-point: status=$(r.status) residual=$(r.residual)")
flush(stdout)

xs = Float64[(JuMP.is_fixed(v) ? JuMP.fix_value(v) :
    (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v))) for v in allv]
# fresh rowscale at the STALLED point
MOI.eval_constraint_jacobian(evaluator, Jval, xs)
v1 = Jval[keep]
mx1 = zeros(Float64, NC)
@inbounds for k in eachindex(jr)
    a = abs(v1[k]); a > mx1[jr[k]] && (mx1[jr[k]] = a)
end
rowscale_fresh = [mx1[i] > 1e-12 ? mx1[i] : 1.0 for i in 1:NC]
g = zeros(NC)
MOI.eval_constraint(evaluator, g, xs)
@inbounds for i in eachindex(g); g[i] -= rhs[i]; end

gs_frozen = g ./ rowscale_frozen
gs_fresh = g ./ rowscale_fresh
of = sortperm(abs.(gs_frozen); rev=true)
on = sortperm(abs.(gs_fresh); rev=true)
println("FROZEN ruler: ||F||∞=$(maximum(abs, gs_frozen)) merit=$(0.5*sum(abs2, gs_frozen)) argmax=row $(of[1])")
println("FRESH  ruler: ||F||∞=$(maximum(abs, gs_fresh)) merit=$(0.5*sum(abs2, gs_fresh)) argmax=row $(on[1])")
println("top-5 frozen rows: $([(of[k], round(abs(gs_frozen[of[k]]); sigdigits=4)) for k in 1:5])")
println("top-5 fresh  rows: $([(on[k], round(abs(gs_fresh[on[k]]); sigdigits=4)) for k in 1:5])")
# ruler drift on the plateau rows
println("ruler drift rs_fresh/rs_frozen on top-5 frozen rows:")
for k in 1:5
    i = of[k]
    println("  row $i  frozen=$(round(rowscale_frozen[i]; sigdigits=3))  fresh=$(round(rowscale_fresh[i]; sigdigits=3))  ratio=$(round(rowscale_fresh[i]/max(rowscale_frozen[i],eps()); sigdigits=3))")
end
println("tol=1e-8 → fresh ruler verdict: " *
    (maximum(abs, gs_fresh) <= 1e-8 ? "CONVERGED (stall was artifact)" : "STILL STALLED"))
println("Done.")
