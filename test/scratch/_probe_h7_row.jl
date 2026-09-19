# H7 row-localization probe: which equation owns the ||F||∞ plateau?
#
# Context: tight repro stalls at t=0.0625 with merit ~1.1e-12 (2-norm
# machine-zero) while ||F||∞ sits at 4.1e-7 (40x above tol). If one row (or a
# few) owns the ∞-norm while the 2-norm merit sees nothing, the fix is row
# targeting — not more CGNR budget, not a global condition estimate.
#
# Method: replicate solve_newton!'s own evaluator construction (MOI Nonlinear,
# SparseReverseMode over ALL model vars), build the FROZEN rowscale at the
# benchmark-start vector exactly as the solver does, then evaluate the SCALED
# residual at the stalled iterate and report argmax rows with equation family.
# Row→equation identity comes from the same insertion order solve_newton! uses
# (list_of_constraint_types/all_constraints loop), so row i here IS row i there.
#
# Checks at the stall point:
#   (a) argmax(|F_i|) — equation family + index (H7 localization)
#   (b) colnrm of the columns that row touches (weak-column? H4 individual)
#   (c) _TINY/_MIN_PRICED_FLOW proximity of the touched variables (H6)
#
# Run: julia --project=. test/scratch/_probe_h7_row.jl
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

# ── replicate the solver's evaluator + frozen rowscale ───────────────────
allv = JuMP.all_variables(m)
allidx = [JuMP.index(v) for v in allv]
pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))
freeid = Dict{Int64,Int}(); freeref = JuMP.VariableRef[]
for (_, v) in vars, vr2 in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr2) && continue
    iv = vr2.index.value; haskey(freeid, iv) && continue
    push!(freeref, vr2); freeid[iv] = length(freeref)
end
nlmodel = MOI.Nonlinear.Model()
rhs = Float64[]; rowdesc = String[]
for (Ftype, S) in JuMP.list_of_constraint_types(m)
    Ftype <: JuMP.VariableRef && continue
    for con in JuMP.all_constraints(m, Ftype, S)
        co = JuMP.constraint_object(con)
        r2 = co.set isa MOI.EqualTo ? co.set.value :
             co.set isa MOI.LessThan ? co.set.upper :
             co.set isa MOI.GreaterThan ? co.set.lower : 0.0
        MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r2))
        push!(rhs, r2)
        # equation identity: owning variable family via the constraint's function
        push!(rowdesc, "$(Ftype)")
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
# frozen rowscale at the CURRENT (stalled) point — matches solver convention
# (theirs is frozen at start; recomputing here only changes the yardstick
# scale, not the argmax ordering, unless a row's max entry moved orders).
MOI.eval_constraint_jacobian(evaluator, Jval, Float64[
    JuMP.is_fixed(v) ? JuMP.fix_value(v) :
    (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv])
v0 = Jval[keep]
mx = zeros(Float64, NC)
@inbounds for k in eachindex(jr)
    a = abs(v0[k]); a > mx[jr[k]] && (mx[jr[k]] = a)
end
rowscale = [mx[i] > 1e-12 ? mx[i] : 1.0 for i in 1:NC]
# scaled residual at stalled iterate
x = Float64[(JuMP.is_fixed(v) ? JuMP.fix_value(v) :
    (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v))) for v in allv]
g = zeros(NC)
MOI.eval_constraint(evaluator, g, x)
@inbounds for i in eachindex(g); g[i] -= rhs[i]; end
gs = g ./ rowscale
ord = sortperm(abs.(gs); rev=true)
println("top-10 scaled |F_i| rows at stall point:")
for k in 1:min(10, NC)
    i = ord[k]
    println("  row $i  |F|=$(gs[i])  type=$(rowdesc[i])")
end
println("||F||∞=$(maximum(abs, gs))  merit=$(0.5*sum(abs2, gs))")
flush(stdout)

# ── (b) column norms of the J_s columns touching the top row ─────────────
top = ord[1]
touch = unique([jc[k] for k in eachindex(jr) if jr[k] == top])
# column 2-norms in the row-scaled (pre-Ruiz) Jacobian
Jb = sparse(jr, jc, [v0[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, length(freeref))
cn = [norm(Jb[:, j]) for j in 1:length(freeref)]
println("top row touches $(length(touch)) free cols; weakest 5:")
for j in sort(touch, by=j->cn[j])[1:min(5, length(touch))]
    println("  col $j  $(JuMP.name(freeref[j]))  colnrm=$(cn[j])")
end
flush(stdout)

# ── (c) floor proximity of touched variables ──────────────────────────────
println("touched vars within 10x of _TINY (|x|<1e-8):")
for j in touch
    vr2 = freeref[j]
    val = JuMP.is_fixed(vr2) ? JuMP.fix_value(vr2) :
          (JuMP.start_value(vr2) === nothing ? NaN : JuMP.start_value(vr2))
    if isfinite(val) && abs(val) < 1e-8
        println("  $(JuMP.name(vr2)) = $val")
    end
end
println("Done.")
