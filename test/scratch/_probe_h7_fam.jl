# H7 family-map probe: which EQUATION FAMILIES own the ||F||∞ plateau?
#
# The row-index probe falsified single-row H7 (top-10 band 4.18e-7→2.62e-7,
# diffuse) but row 37426 named aprim[2,5]/E_aprim. This probe maps the top-50
# scaled-residual rows back to E_* families via insertion-order bookkeeping:
# for each E_* call in build_model_full!'s order, record the constraint count
# added, so row i → family by interval lookup. Same evaluator + frozen
# rowscale construction as _probe_h7_row.jl.
#
# Decides: tech-shifter block (aprim/alab_o/acap/atot + their price duals) →
# targeted polish of the shifter subsystem; heterogeneous band → back to the
# whole-system condition estimate for the H5 gate.
#
# Run: julia --project=. test/scratch/_probe_h7_fam.jl
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

# ── family intervals: count constraints per E_* block in CALL ORDER ──────
# build_model! equation blocks (base) run first, then build_model_full!
# Excerpts 6-54. We recover intervals empirically: walk list_of_constraint_types
# in insertion order is NOT grouped by call — instead, tag each constraint by
# the free-variable family majority of its Jacobian row (structural attribution).
allv = JuMP.all_variables(m)
allidx = [JuMP.index(v) for v in allv]
pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))
freeid = Dict{Int64,Int}(); freeref = JuMP.VariableRef[]
fam_of = String[]
for (vnm, v) in vars, vr2 in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr2) && continue
    iv = vr2.index.value; haskey(freeid, iv) && continue
    push!(freeref, vr2); push!(fam_of, String(vnm)); freeid[iv] = length(freeref)
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
xv = Float64[(JuMP.is_fixed(v) ? JuMP.fix_value(v) :
    (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v))) for v in allv]
MOI.eval_constraint_jacobian(evaluator, Jval, xv)
v0 = Jval[keep]
mx = zeros(Float64, NC)
@inbounds for k in eachindex(jr)
    a = abs(v0[k]); a > mx[jr[k]] && (mx[jr[k]] = a)
end
rowscale = [mx[i] > 1e-12 ? mx[i] : 1.0 for i in 1:NC]
g = zeros(NC)
MOI.eval_constraint(evaluator, g, xv)
@inbounds for i in eachindex(g); g[i] -= rhs[i]; end
gs = g ./ rowscale
ord = sortperm(abs.(gs); rev=true)
println("top-30 rows: scaled|F|, dominant free-var family in row, rowscale:")
for k in 1:min(30, NC)
    i = ord[k]
    touched = unique([fam_of[jc[q]] for q in eachindex(jr) if jr[q] == i])
    dom = isempty(touched) ? "(no free vars)" : join(touched[1:min(3, length(touched))], ",")
    println("  row $i  |F|=$(round(abs(gs[i]); sigdigits=4))  rs=$(round(rowscale[i]; sigdigits=3))  fam=[$dom]")
end
println("||F||∞=$(maximum(abs, gs))  merit=$(0.5*sum(abs2, gs))")
# family histogram over top-50
hist = Dict{String,Int}()
for k in 1:min(50, NC)
    i = ord[k]
    for q in eachindex(jr)
        jr[q] == i || continue
        f = fam_of[jc[q]]
        hist[f] = get(hist, f, 0) + 1
    end
end
println("family histogram over top-50 rows (touches):")
for (f, c) in sort(collect(hist), by=x->-x[2])[1:min(20, length(hist))]
    println("  $(rpad(f,14)) $c")
end
println("Done.")
