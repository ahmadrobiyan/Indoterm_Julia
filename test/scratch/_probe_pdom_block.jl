# pdom block diagnostic at the t=0.0625 stall point.
#
# Context: normalized hit-rate analysis over the plateau band showed pdom at
# 39% (58/150 rows) vs pmake at background (51/3750, 1.4%). The [2,5]
# factor-nest story covers 3 of the top 5 rows; the price block dominates the
# wider band. This probe inspects the pdom family directly: count in band,
# residual distribution, rowscale spread, touching column norms.
#
# Staging mirrors the tight repro (0.015625/0.03125/0.046875 converged), then
# lands on the plateau with maxit=20 at t=0.0625. connames come from the
# solver's own row→equation identity (solve_newton!.jl); family match is
# substring "pdom" on the constraint name, falling back to touching-family
# attribution when names are anonymous ("?").
#
# Run: julia --project=. test/scratch/_probe_pdom_block.jl (verbose=ON)
# Log: logs/probe_pdom_block_2026-09-18.log ; heartbeat: logs/heartbeat.jsonl
include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra
import JuMP.MOI as MOI
using Dates

HEARTBEAT = joinpath(@__DIR__, "..", "..", "logs", "heartbeat.jsonl")
hb(stage, status, extra="") = open(HEARTBEAT, "a") do io
    write(io, "{\"stage\": \"$stage\", \"timestamp\": \"$(now())\", \"status\": \"$status\"$extra}\n")
    flush(io)
end

hb("pdom_probe", "running", ", \"phase\": \"pipeline\"")
agg6, params = cached_pipeline(6)
COAL = findfirst(==("Coal"), AGGCOM)

hb("pdom_probe", "running", ", \"phase\": \"build\"")
(m, vars) = build_model_full!(agg6, params)
initialize_model!(m, vars; bmk_levels=benchmark_levels(params), numeraire=:exrate)
apply_swaps!(m, vars, collect(COALPRICE_SWAPS); verbose=true)
r0 = solve_newton!(m, vars; maxit=5, tol=1e-8, verbose=true)
@assert r0.residual <= 1e-8 "benchmark must reproduce under this closure"

nm, idx, vr = IndotermJulia._shock_ref(vars, ("fpexp_d", COAL))
b = JuMP.fix_value(vr); tg = logpct(-20.0)
println("shock base=$b target=$tg t_stall=0.0625")
flush(stdout)
hb("pdom_probe", "running", ", \"phase\": \"staging\"")
for t_ in (0.015625, 0.03125, 0.046875)
    JuMP.fix(vr, b + t_ * (tg - b); force=true)
    r = solve_newton!(m, vars; maxit=30, tol=1e-8, verbose=true)
    @assert r.residual <= 1e-8 "staging failed at $t_"
end
hb("pdom_probe", "running", ", \"phase\": \"land\"")
JuMP.fix(vr, b + 0.0625 * (tg - b); force=true)
r = solve_newton!(m, vars; maxit=20, tol=1e-8, verbose=true)
println("landed: status=$(r.status) residual=$(r.residual)")
flush(stdout)

# ── pdom subspace extraction ─────────────────────────────────────────────
hb("pdom_probe", "running", ", \"phase\": \"extract\"")
allv = JuMP.all_variables(m)
allidx = [JuMP.index(v) for v in allv]
pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))
freeid = Dict{Int64,Int}(); freeref = JuMP.VariableRef[]; fam_of = String[]
for (vnm, v) in vars, vr2 in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr2) && continue
    iv = vr2.index.value; haskey(freeid, iv) && continue
    push!(freeref, vr2); push!(fam_of, String(vnm)); freeid[iv] = length(freeref)
end
nlmodel = MOI.Nonlinear.Model()
rhs = Float64[]; connames = String[]
for (Ftype, S) in JuMP.list_of_constraint_types(m)
    Ftype <: JuMP.VariableRef && continue
    for con in JuMP.all_constraints(m, Ftype, S)
        co = JuMP.constraint_object(con)
        r2 = co.set isa MOI.EqualTo ? co.set.value :
             co.set isa MOI.LessThan ? co.set.upper :
             co.set isa MOI.GreaterThan ? co.set.lower : 0.0
        MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r2))
        push!(rhs, r2)
        push!(connames, try JuMP.name(con) catch; "?" end)
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
xs = Float64[(JuMP.is_fixed(v) ? JuMP.fix_value(v) :
    (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v))) for v in allv]
MOI.eval_constraint_jacobian(evaluator, Jval, xs)
v0 = Jval[keep]
mx = zeros(Float64, NC)
@inbounds for k in eachindex(jr)
    a = abs(v0[k]); a > mx[jr[k]] && (mx[jr[k]] = a)
end
rowscale = [mx[i] > 1e-12 ? mx[i] : 1.0 for i in 1:NC]
g = zeros(NC)
MOI.eval_constraint(evaluator, g, xs)
@inbounds for i in eachindex(g); g[i] -= rhs[i]; end
gs = g ./ rowscale
Jb = sparse(jr, jc, [v0[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, length(freeref))
# Column norms in O(nnz) — explicit accumulator, NOT [norm(Jb[:,j]) for j].
# The slicing form allocates a SparseVector per column and (with the per-row
# family scans below) made the first run of this probe sit in 'extract' for
# 7+ min without output. See AGENTS.md on generator/accumulator patterns.
cn = zeros(Float64, length(freeref))
@inbounds for k in eachindex(jc)
    v = Jb.nzval[k]; cn[jc[k]] += v * v
end
cn .= sqrt.(cn)

# Row → touching-family sets, built ONCE in O(nnz). The per-row
# `any(... for q in eachindex(jr) if jr[q] == i)` form is O(NC·nnz) ≈ 2.6e10
# and never finishes.
rowfams = [Set{String}() for _ in 1:NC]
@inbounds for k in eachindex(jr)
    push!(rowfams[jr[k]], fam_of[jc[k]])
end

is_pdom_row(i) = occursin("pdom", connames[i]) || ("pdom" in rowfams[i])
pdom_rows = [i for i in 1:NC if is_pdom_row(i)]
println("pdom equations in model: $(length(pdom_rows)) / $NC")
inband = [i for i in pdom_rows if abs(gs[i]) > 1e-7]
println("pdom rows in plateau band (|F|>1e-7): $(length(inband)) ($(round(100*length(inband)/max(length(pdom_rows),1); digits=1))%)")
vals = sort(abs.(gs[pdom_rows]))
println("pdom |F| dist: min=$(round(minimum(vals); sigdigits=4)) median=$(round(vals[clamp(length(vals)÷2,1,end)]; sigdigits=4)) max=$(round(maximum(vals); sigdigits=4))")
rs_pdom = rowscale[pdom_rows]
println("pdom rowscale: min=$(round(minimum(rs_pdom); sigdigits=3)) median=$(round(sort(rs_pdom)[length(rs_pdom)÷2]; sigdigits=3)) max=$(round(maximum(rs_pdom); sigdigits=3))")
# touching columns of the worst 5 pdom rows
ord = sort(pdom_rows, by=i->abs(gs[i]); rev=true)[1:min(5, length(pdom_rows))]
for i in ord
    touched = unique([jc[q] for q in eachindex(jr) if jr[q] == i])
    println("  pdom row $i $(connames[i]) |F|=$(round(gs[i]; sigdigits=4)) rs=$(round(rowscale[i]; sigdigits=3))")
    for j in touched[1:min(4, length(touched))]
        println("      col $(JuMP.name(freeref[j])) colnrm=$(round(cn[j]; sigdigits=3))")
    end
end
flush(stdout)
hb("pdom_probe", "done", ", \"phase\": \"extract\"")
println("Done.")
