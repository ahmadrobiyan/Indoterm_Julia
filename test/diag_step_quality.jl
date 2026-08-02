"""
C4'b — is the (now correctly transformed) Newton step actually good?

`solve_newton!.jl:443` used to convert the CGNR direction back to original units with
`dy .* Dc_ruiz`, but `ruiz_scales` returns Dr/Dc as DIVISORS (J_s = Dr⁻¹·Jbase·Dc⁻¹,
rhs_s = gs./Dr), so the correct transform is `dy ./ Dc`. That is fixed. The shock now
descends instead of growing linearly in α — but it still crawls, for two candidate
reasons visible in the line search:

  * fraction-to-boundary throttles α to ~1.2e-3 (some variable sits near its lower bound);
  * the step cap `max(0.10*|x|, 10.0)` is meaningless for value-denominated additive
    variables such as delXGDPEXP (benchmark 0, natural scale 1e3–1e6 Rp bn), so capping
    turns the Newton direction into a non-Newton one.

This script separates those from a genuinely bad step. It takes the exact LU step at the
stalled point and evaluates ‖F‖∞ along it with NO cap and NO bound handling, starting at
α = 1. If ‖F‖∞ collapses at α ≈ 1, the direction is correct and the entire remaining
problem is the cap/bound heuristics. It also names:

  * the variable whose lower-bound gap sets the fraction-to-boundary limit;
  * every variable the cap would bind, with the cap it would impose vs the step it wants;
  * how many variables the full step drives below their lower bound, and by how much.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_step_quality.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra
const MOI = JuMP.MOI
const IJ = IndotermJulia

canon(s::AbstractString) = replace(s, r"\[[0-9, ]+\]" => "[.]")

SHOCK_VAR = get(ENV, "STALL_SHOCK_VAR", "blabnat")
SHOCK_VAL = parse(Float64, get(ENV, "STALL_SHOCK_VAL", "0.9997"))

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

ITERS = parse(Int, get(ENV, "STALL_ITERS", "12"))

r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("benchmark: status=$(r0.status)  ‖F‖∞=$(r0.residual)")

# Snapshot the benchmark point so drift can be measured. With STALL_ITERS=0 the
# analysis below runs at the SHOCK START POINT (= benchmark), which separates
# "the first Newton step is already runaway" from "the iterates drift there".
bmk = Dict(JuMP.name(v) => (JuMP.is_fixed(v) ? JuMP.fix_value(v) :
           (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)))
           for v in JuMP.all_variables(m))

sv = vars[SHOCK_VAR]
JuMP.fix(sv isa AbstractArray ? first(sv) : sv, SHOCK_VAL; force=true)
if ITERS > 0
    r1 = solve_newton!(m, vars; maxit=ITERS, verbose=false)
    println("shock $SHOCK_VAR=$SHOCK_VAL: status=$(r1.status) iters=$(r1.iters) ‖F‖∞=$(r1.residual)")
else
    println("shock $SHOCK_VAR=$SHOCK_VAL: no iterations taken (analysing the START point)")
end

# ── rebuild the evaluator at the solver's own x ─────────────────────────────
nowfree = VariableRef[]
for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr) || push!(nowfree, vr)
end
N = length(nowfree)
varname = [JuMP.name(v) for v in nowfree]
lb = [JuMP.has_lower_bound(v) ? JuMP.lower_bound(v) : -Inf for v in nowfree]
allv = JuMP.all_variables(m); allidx = [JuMP.index(v) for v in allv]
pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))

nlmodel = MOI.Nonlinear.Model(); rhs = Float64[]; rowsig = String[]
for (Ftype, S) in list_of_constraint_types(m)
    Ftype <: VariableRef && continue
    for con in all_constraints(m, Ftype, S)
        co = JuMP.constraint_object(con); r = IJ._rhs_of(co.set)
        MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r))
        push!(rhs, r); push!(rowsig, canon(string(co.func)))
    end
end
NC = length(rhs)
evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
MOI.initialize(evaluator, [:Grad, :Jac])

x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
     (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]
freecols = [pos[JuMP.index(v).value] for v in nowfree]

st = MOI.jacobian_structure(evaluator)
jrows_all = getindex.(st, 1); jcols_all = getindex.(st, 2)
Jval = zeros(length(jrows_all))
col2free = zeros(Int, length(allv)); for (i, c) in enumerate(freecols); col2free[c] = i; end
keep = [col2free[c] != 0 for c in jcols_all]
jr = jrows_all[keep]; jc = [col2free[c] for c in jcols_all[keep]]

MOI.eval_constraint_jacobian(evaluator, Jval, x); vals = Jval[keep]
rowscale = ones(NC)
let mx = zeros(NC)
    for k in eachindex(jr); a = abs(vals[k]); a > mx[jr[k]] && (mx[jr[k]] = a); end
    for i in 1:NC; rowscale[i] = mx[i] > 1e-12 ? mx[i] : 1.0; end
end
g = zeros(NC); MOI.eval_constraint(evaluator, g, x); g .-= rhs
gs = g ./ rowscale

function ruiz(J0::SparseMatrixCSC{Float64,Int}; maxiter=20, tol=1e-6)
    mm, nn = size(J0); Dr = ones(mm); Dc = ones(nn); Jw = copy(J0)
    for _ in 1:maxiter
        rowmax = zeros(mm)
        for k in 1:nnz(Jw); i = Jw.rowval[k]; a = abs(Jw.nzval[k]); a > rowmax[i] && (rowmax[i] = a); end
        for i in 1:mm; Dr[i] *= rowmax[i] > 1e-12 ? sqrt(rowmax[i]) : 1.0; end
        for k in 1:nnz(Jw); i = Jw.rowval[k]; Jw.nzval[k] /= (rowmax[i] > 1e-12 ? sqrt(rowmax[i]) : 1.0); end
        colmax = zeros(nn)
        for j in 1:nn, k in Jw.colptr[j]:(Jw.colptr[j+1]-1)
            a = abs(Jw.nzval[k]); a > colmax[j] && (colmax[j] = a)
        end
        for j in 1:nn; Dc[j] *= colmax[j] > 1e-12 ? sqrt(colmax[j]) : 1.0; end
        for j in 1:nn
            s = colmax[j] > 1e-12 ? sqrt(colmax[j]) : 1.0
            for k in Jw.colptr[j]:(Jw.colptr[j+1]-1); Jw.nzval[k] /= s; end
        end
        dev = 0.0
        for i in 1:mm; s = rowmax[i] > 1e-12 ? sqrt(rowmax[i]) : 1.0; dev = max(dev, abs(s-1)); end
        for j in 1:nn; s = colmax[j] > 1e-12 ? sqrt(colmax[j]) : 1.0; dev = max(dev, abs(s-1)); end
        dev < tol && break
    end
    return Dr, Dc
end

Jbase = sparse(jr, jc, [vals[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, N)
Dr, Dc = ruiz(Jbase)
J_s = sparse(jr, jc, [vals[k] / (rowscale[jr[k]] * Dr[jr[k]] * Dc[jc[k]]) for k in eachindex(jr)], NC, N)
b = -(gs ./ Dr)
println("\nJ_s $(size(J_s))  nnz=$(nnz(J_s))  ‖b‖₂=$(norm(b))   ‖F‖∞(scaled)=$(norm(gs,Inf))")

F = lu(J_s)
dy = F \ b
println("lu: ‖J_s·dy−b‖/‖b‖ = $(norm(J_s*dy .- b)/norm(b))")

# The step in ORIGINAL variable units. Dc is a divisor, so this is dy ./ Dc.
d = dy ./ Dc
println("‖dy‖∞ (Ruiz units) = $(maximum(abs, dy))     ‖d‖∞ (original units) = $(maximum(abs, d))")
# Sanity: Jbase*d should equal -gs.
println("check ‖Jbase·d + gs‖₂/‖gs‖₂ = $(norm(Jbase*d .+ gs)/norm(gs))")

# ── residual along the UNCONSTRAINED step ──────────────────────────────────
gbuf = zeros(NC)
function res_inf(alpha)
    xt = copy(x)
    @inbounds for i in 1:N; xt[freecols[i]] = x[freecols[i]] + alpha * d[i]; end
    try
        MOI.eval_constraint(evaluator, gbuf, xt)
        gt = (gbuf .- rhs) ./ rowscale
        return norm(gt, Inf), count(!isfinite, gt)
    catch e
        return NaN, -1
    end
end

function res_at(alpha)
    xt = copy(x)
    @inbounds for i in 1:N; xt[freecols[i]] = x[freecols[i]] + alpha * d[i]; end
    MOI.eval_constraint(evaluator, gbuf, xt)
    return (gbuf .- rhs) ./ rowscale
end

println("\n── ‖F‖∞ along the FULL Newton step (no cap, no fraction-to-boundary) ──")
println("   α = 0 gives ‖F‖∞ = $(norm(gs, Inf))")
println("   Newton predicts (1-α)·‖F‖∞ for every row where the linearisation holds;")
println("   the row named below is the one that does NOT, i.e. the one capping α.")
for a in (1.0, 0.9, 0.75, 0.5, 0.25, 0.1, 0.05, 0.01, 1e-3, 1e-4)
    gt = res_at(a)
    r, nbad = norm(gt, Inf), count(!isfinite, gt)
    i = argmax(map(v -> isfinite(v) ? abs(v) : -Inf, gt))
    println("   α=$(rpad(a,8))  ‖F‖∞=$(rpad(round(r, sigdigits=6), 12))  nonfin=$(rpad(nbad,4))" *
            "  argmax=row $(rpad(i,7)) $(first(rowsig[i], 90))")
end

# Which rows blow up worst relative to what Newton promised? A row obeying the
# linearisation lands at (1-α)F_i; the ratio below is 1 for such a row and large
# for the nonlinear rows that are actually forcing the backtrack.
println("\n── rows most violating the Newton prediction at α=1 (|F_i(x+d)| / |F_i(x)|) ──")
let g1 = res_at(1.0)
    sc = [(isfinite(g1[i]) && abs(gs[i]) > 1e-14) ? abs(g1[i]) / abs(gs[i]) : -Inf for i in 1:NC]
    for i in sortperm(sc; rev=true)[1:15]
        isfinite(sc[i]) || continue
        println("   $(rpad(round(sc[i], sigdigits=5),12))×   F0=$(rpad(round(gs[i],sigdigits=4),12)) " *
                "F1=$(rpad(round(g1[i],sigdigits=4),12)) row $(rpad(i,7)) $(first(rowsig[i], 80))")
    end
end

# The residual FLOOR: which equation families still carry residual at α=0?
println("\n── largest residuals at the stalled point (α=0), by row ──")
for i in sortperm(abs.(gs); rev=true)[1:15]
    println("   |F|=$(rpad(round(abs(gs[i]),sigdigits=5),12)) row $(rpad(i,7)) $(first(rowsig[i], 110))")
end
let fam = Dict{String,Float64}()
    for i in 1:NC; fam[rowsig[i]] = get(fam, rowsig[i], 0.0) + gs[i]^2; end
    tot = sum(values(fam))
    println("\n── residual mass by equation family (share of ‖F‖₂²) ──")
    for k in sort!(collect(keys(fam)); by = k -> -fam[k])[1:10]
        println("   $(rpad(round(100*fam[k]/tot, sigdigits=4),8))%  $(first(k, 110))")
    end
end

# ── what would throttle α in the solver ────────────────────────────────────
# ── which ROWS go non-finite, and at what α ────────────────────────────────
# This is the real blocker: the step is exact but leaves the domain of some
# log(x + _TINY). Name the rows that break and the variables driving them.
function bad_rows(alpha, k=12)
    xt = copy(x)
    @inbounds for i in 1:N; xt[freecols[i]] = x[freecols[i]] + alpha * d[i]; end
    MOI.eval_constraint(evaluator, gbuf, xt)
    gt = (gbuf .- rhs) ./ rowscale
    idx = findall(!isfinite, gt)
    println("\n── α=$alpha : $(length(idx)) non-finite rows ──")
    for i in first(idx, k); println("   row $i  $(first(rowsig[i], 150))"); end
    return idx
end
for a in (0.01, 0.25, 1.0); bad_rows(a); end

# Which free variables go negative (or below lb) along the step, earliest first?
function domain_exits(k=25)
    println("\n── variables driven NEGATIVE by the full step (earliest α first) ──")
    rows = Tuple{Float64,Int}[]
    for i in 1:N
        xi = x[freecols[i]]
        if d[i] < -1e-30 && xi > 0
            a = xi / (-d[i])          # α at which this variable hits zero
            a < 1.0 && push!(rows, (a, i))
        end
    end
    sort!(rows; by = t -> t[1])
    println("   $(length(rows)) variables cross zero before α=1")
    for (a, i) in first(rows, k)
        println("   α*=$(rpad(round(a, sigdigits=4),12)) $(rpad(varname[i],30)) x=$(rpad(round(x[freecols[i]],sigdigits=5),13)) lb=$(rpad(lb[i],9)) d=$(round(d[i],sigdigits=5))")
    end
end
domain_exits()

println("\n── how far have the iterates DRIFTED from the benchmark? (top 25 relative) ──")
println("   a $(round(100*(1-SHOCK_VAL), sigdigits=3))% shock should move ratio-type variables by O(that).")
let drift = [(abs(bmk[varname[i]]) > 1e-6 ?
              abs(x[freecols[i]] - bmk[varname[i]]) / abs(bmk[varname[i]]) : -Inf) for i in 1:N]
    for i in sortperm(drift; rev=true)[1:25]
        isfinite(drift[i]) || continue
        println("   $(rpad(round(100*drift[i], sigdigits=5),12))%  $(rpad(varname[i],28)) " *
                "bmk=$(rpad(round(bmk[varname[i]],sigdigits=6),14)) now=$(round(x[freecols[i]],sigdigits=6))")
    end
end

println("\n── top 25 |d| — what the step actually wants to move ──")
for i in sortperm(abs.(d); rev=true)[1:25]
    println("   d=$(rpad(round(d[i],sigdigits=6),14)) $(rpad(varname[i],30)) x=$(rpad(round(x[freecols[i]],sigdigits=6),14)) lb=$(lb[i])")
end

# The absolute-|d| ranking is dominated by value-denominated passengers (del*,
# xcom, xmake — all in Rp bn). What matters for the line search is the RELATIVE
# move, because that is what sets the curvature of every log(·) the variable
# enters. A 0.03% shock should produce O(0.03%) relative moves; anything far
# above that is the amplified mode.
println("\n── top 30 RELATIVE step |d_i / x_i| — the amplified mode ──")
let rel = [(abs(x[freecols[i]]) > 1e-8 ? abs(d[i]) / abs(x[freecols[i]]) : -Inf) for i in 1:N]
    ord = sortperm(rel; rev=true)
    for i in first(ord, 30)
        isfinite(rel[i]) || continue
        println("   $(rpad(round(100*rel[i], sigdigits=5),12))%  $(rpad(varname[i],30)) " *
                "x=$(rpad(round(x[freecols[i]],sigdigits=6),14)) d=$(round(d[i],sigdigits=6))")
    end
    # How many variables move by more than 10x / 100x the shock size?
    shock = abs(log(SHOCK_VAL))
    for f in (10, 100, 1000)
        println("   variables moving > $(f)x the shock ($(round(f*shock,sigdigits=3))): " *
                "$(count(v -> isfinite(v) && v > f*shock, rel))")
    end
end

println("\n── fraction-to-boundary: which variable binds? ──")
function ftb_limit()
    ba = 1.0; bi = 0
    for i in 1:N
        if isfinite(lb[i]) && d[i] < -1e-10
            gap = x[freecols[i]] - lb[i]
            if gap > 1e-8
                cap = 0.99 * gap / (-d[i])
                if cap > 0 && cap < ba; ba = cap; bi = i; end
            end
        end
    end
    return ba, bi
end
let (ba, bi) = ftb_limit()
    if bi > 0
        println("   a0 = $ba  set by $(varname[bi])")
        println("      x=$(x[freecols[bi]])  lb=$(lb[bi])  gap=$(x[freecols[bi]]-lb[bi])  d=$(d[bi])")
    end
end
println("   variables driven below their lower bound by the FULL step: " *
        "$(count(i -> isfinite(lb[i]) && x[freecols[i]] + d[i] < lb[i], 1:N)) of $N")

println("\n── top 10 fraction-to-boundary limiters ──")
ratios = [(isfinite(lb[i]) && d[i] < -1e-10 && (x[freecols[i]] - lb[i]) > 1e-8) ?
          0.99*(x[freecols[i]] - lb[i])/(-d[i]) : Inf for i in 1:N]
for i in sortperm(ratios)[1:10]
    isfinite(ratios[i]) || break
    println("   a=$(rpad(round(ratios[i], sigdigits=4),12)) $(rpad(varname[i],28)) x=$(rpad(round(x[freecols[i]],sigdigits=5),12)) lb=$(rpad(lb[i],10)) d=$(round(d[i],sigdigits=5))")
end

# ── which variables would the solver's cap bind, and how hard? ─────────────
println("\n── step cap `max(0.10|x|, 10.0)` (additive) / `0.10*max(|x|,1e-3)` (ratio) ──")
rel_cap = 0.10
binds = Tuple{Float64,String,Float64,Float64}[]
for i in 1:N
    absx = abs(x[freecols[i]])
    maxstep = (isfinite(lb[i]) && lb[i] >= 1e-6 - 1e-9) ? rel_cap*max(absx, 1e-3) :
                                                          max(rel_cap*absx, 10.0)
    if abs(d[i]) > maxstep
        push!(binds, (abs(d[i])/maxstep, varname[i], d[i], maxstep))
    end
end
sort!(binds; by = t -> -t[1])
println("   variables the cap binds: $(length(binds)) of $N")
for (ratio, nm, di, mx) in first(binds, 20)
    println("   $(rpad(round(ratio, sigdigits=5),12))× over cap   $(rpad(nm,30)) d=$(rpad(round(di,sigdigits=6),13)) cap=$(round(mx,sigdigits=6))")
end

println("\nDone.")
