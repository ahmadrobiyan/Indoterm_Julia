"""
Is the damping — not the step — what stops the shock converging?

Everything upstream now checks out. At the shock start point the LU step satisfies
`‖Jbase·d + gs‖₂/‖gs‖₂ = 3.2e-10`, no row goes non-finite at α=1, no variable
crosses zero, and the relative step is O(shock) for all but three variables. The
three exceptions — `pcap[2,5]`, `plnd[2,5]` (−4.65%) and `xinv[2,2,5]` (−3.44%) —
are not a defect: `xcap`/`xlnd` are exogenous in the base static closure, so the
sector-specific fixed factors of Livestock/BaliNusa absorb the whole adjustment in
their *price*. `pcap ∝ pprim` follows from inverting two CES demands with both
quantities pinned, which is why pcap and plnd move bit-identically.

The consequence is that row 38476, `gret[i,d] = log(pcap+τ) − log(pinvitot+τ)`,
picks up ½(0.046)² ≈ 1.1e-3 of pure log-curvature at α=1 — larger than the whole
starting residual of 3.0e-4. `solve_newton!` runs a MONOTONE ‖F‖∞ line search, so
that one row vetoes the full step; it settles on α=0.25, the residual falls by
exactly (1−α), and α then shrinks every iteration. That is linear convergence
with a shrinking rate — the observed stall.

A curvature overshoot at α=1 is exactly what Newton is supposed to absorb: the
next step corrects it at second order. This script tests that directly by running
UNDAMPED Newton (α=1, fraction-to-boundary as the only safeguard) and watching
‖F‖∞. If it converges quadratically, the equations and the step are both right
and the fix belongs in the line search (non-monotone / watchdog), not the model.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_full_newton.jl
      STALL_SHOCK_VAL=0.99 julia --project=IndotermJulia .../diag_full_newton.jl
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra
const MOI = JuMP.MOI
const IJ = IndotermJulia

SHOCK_VAR = get(ENV, "STALL_SHOCK_VAR", "blabnat")
SHOCK_VAL = parse(Float64, get(ENV, "STALL_SHOCK_VAL", "0.9997"))
MAXIT     = parse(Int, get(ENV, "STALL_MAXIT", "15"))

agg6, params = cached_pipeline(6)
t = @elapsed ((m, vars) = build_model_full!(agg6, params))
println("build_model_full!: $(round(t, digits=1))s")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))
r0 = solve_newton!(m, vars; maxit=5, verbose=false)
println("benchmark: status=$(r0.status)  ‖F‖∞=$(r0.residual)")

sv = vars[SHOCK_VAR]
JuMP.fix(sv isa AbstractArray ? first(sv) : sv, SHOCK_VAL; force=true)
println("shock $SHOCK_VAR = $SHOCK_VAL\n")

nowfree = VariableRef[]
for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr) || push!(nowfree, vr)
end
N = length(nowfree)
varname = [JuMP.name(v) for v in nowfree]
lb = [JuMP.has_lower_bound(v) ? JuMP.lower_bound(v) : -Inf for v in nowfree]
allv = JuMP.all_variables(m); allidx = [JuMP.index(v) for v in allv]
pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))

nlmodel = MOI.Nonlinear.Model(); rhs = Float64[]
for (Ftype, S) in list_of_constraint_types(m)
    Ftype <: VariableRef && continue
    for con in all_constraints(m, Ftype, S)
        co = JuMP.constraint_object(con); r = IJ._rhs_of(co.set)
        MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(r))
        push!(rhs, r)
    end
end
NC = length(rhs)
evaluator = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
MOI.initialize(evaluator, [:Grad, :Jac])
println("system: $NC equations, $N unknowns (square = $(NC == N))")

x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
     (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]
freecols = [pos[JuMP.index(v).value] for v in nowfree]

st = MOI.jacobian_structure(evaluator)
jrows_all = getindex.(st, 1); jcols_all = getindex.(st, 2)
Jval = zeros(length(jrows_all))
col2free = zeros(Int, length(allv)); for (i, c) in enumerate(freecols); col2free[c] = i; end
keep = [col2free[c] != 0 for c in jcols_all]
jr = jrows_all[keep]; jc = [col2free[c] for c in jcols_all[keep]]

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

gbuf = zeros(NC)
function resid(xv)
    MOI.eval_constraint(evaluator, gbuf, xv)
    return gbuf .- rhs
end

# A FROZEN reference scaling, computed once at the start point. `solve_newton!`
# recomputes its row scaling every iteration and then compares ‖F‖∞ ACROSS
# iterations — so its residuals are measured in units that change between
# readings. Reporting both here separates a genuine stall from that artifact.
rowscale0 = ones(NC)
let
    MOI.eval_constraint_jacobian(evaluator, Jval, x); v0 = Jval[keep]
    mx = zeros(NC)
    for k in eachindex(jr); a = abs(v0[k]); a > mx[jr[k]] && (mx[jr[k]] = a); end
    for i in 1:NC; rowscale0[i] = mx[i] > 1e-12 ? mx[i] : 1.0; end
end

# Label each row by the VARIABLE PREFIXES it touches, read off the Jacobian
# sparsity pattern. `string(co.func)` is the obvious alternative and is a trap:
# it renders the whole expression tree as text for all 83k constraints (some
# aggregation rows have thousands of terms) before any truncation, and never
# finished in 20 minutes. This costs one pass over the nonzeros, and the prefix
# set identifies the equation family exactly — E_gret's rows and nobody else's
# read "gret|pcap|pinvitot".
allprefix = [String(first(split(JuMP.name(v), "["))) for v in nowfree]
rowsig = let acc = [Set{String}() for _ in 1:NC]
    for k in eachindex(jr); push!(acc[jr[k]], allprefix[jc[k]]); end
    [join(sort!(collect(s)), "|") for s in acc]
end
famkey(s) = first(s, 90)

# For a single named row we also want the INDICES, so the culprit is identified
# as `gret[2,5]` (Livestock/BaliNusa) rather than just the family.
rowcols = let acc = [Int[] for _ in 1:NC]
    for k in eachindex(jr); push!(acc[jr[k]], jc[k]); end
    acc
end
rowvars(i, k=6) = join((varname[c] for c in first(rowcols[i], k)), " ")

println("\n── UNDAMPED Newton (α = 1, fraction-to-boundary only) ──")
println("iter  ‖F‖∞(periter)   ‖F‖∞(frozen)   ‖F‖₂(frozen)   α       max|d|      lin_rel")
prev = Inf
for it in 1:MAXIT
    MOI.eval_constraint_jacobian(evaluator, Jval, x); vals = Jval[keep]
    rowscale = ones(NC)
    let mx = zeros(NC)
        for k in eachindex(jr); a = abs(vals[k]); a > mx[jr[k]] && (mx[jr[k]] = a); end
        for i in 1:NC; rowscale[i] = mx[i] > 1e-12 ? mx[i] : 1.0; end
    end
    graw = copy(resid(x))
    gs = graw ./ rowscale
    nrm = norm(gs, Inf)

    Jbase = sparse(jr, jc, [vals[k] / rowscale[jr[k]] for k in eachindex(jr)], NC, N)
    Dr, Dc = ruiz(Jbase)
    J_s = sparse(jr, jc, [vals[k] / (rowscale[jr[k]] * Dr[jr[k]] * Dc[jc[k]]) for k in eachindex(jr)], NC, N)
    b = -(gs ./ Dr)
    dy = lu(J_s) \ b
    lin = norm(J_s * dy .- b) / max(norm(b), eps())
    d = dy ./ Dc                     # Dc is a DIVISOR: J_s = Dr⁻¹·Jbase·Dc⁻¹

    # Only safeguard: never leave the domain. No magnitude cap, no backtracking.
    a = 1.0
    for i in 1:N
        if isfinite(lb[i]) && d[i] < -1e-30
            gap = x[freecols[i]] - lb[i]
            gap > 1e-12 && (a = min(a, 0.99 * gap / (-d[i])))
        end
    end
    a = min(a, 1.0)

    # Measure in FROZEN units too, so successive iterations are comparable, and
    # report ‖F‖₂ — the merit function we are proposing to switch to.
    gf   = graw ./ rowscale0
    nrmf = norm(gf, Inf)
    n2f  = norm(gf, 2)
    println("$(rpad(it,4))  $(rpad(round(nrm, sigdigits=6),14))  $(rpad(round(nrmf, sigdigits=6),14))  " *
            "$(rpad(round(n2f, sigdigits=6),14))  $(rpad(round(a,sigdigits=4),7))  " *
            "$(rpad(round(maximum(abs,d),sigdigits=6),11))  $(round(lin, sigdigits=3))")

    # WHICH EQUATIONS carry the error. (The previous metric — largest relative
    # variable step — was meaningless: it flagged del* shift variables that
    # legitimately sit at 0, so |d/x| is huge for a step of no consequence.)
    let ip = argmax(abs.(gs)), jf = argmax(abs.(gf))
        println("        argmax per-iter units: |F|=$(round(abs(gs[ip]),sigdigits=4))  raw=$(round(graw[ip],sigdigits=4))  rs=$(round(rowscale[ip],sigdigits=4))  row $ip  $(rowsig[ip])")
        println("            vars: $(rowvars(ip))")
        if ip != jf
            println("        argmax frozen  units: |F|=$(round(abs(gf[jf]),sigdigits=4))  raw=$(round(graw[jf],sigdigits=4))  rs=$(round(rowscale0[jf],sigdigits=4))  row $jf  $(rowsig[jf])")
            println("            vars: $(rowvars(jf))")
        end
    end
    for i in first(sortperm(abs.(gf); rev=true), 5)
        println("        |F|=$(rpad(round(abs(gf[i]),sigdigits=5),12)) raw=$(rpad(round(graw[i],sigdigits=5),12)) $(first(rowsig[i],60))  [$(rowvars(i,3))]")
    end
    let fam = Dict{String,Float64}()
        for i in 1:NC; k = famkey(rowsig[i]); fam[k] = get(fam, k, 0.0) + gf[i]^2; end
        tot = max(sum(values(fam)), eps())
        println("        ── residual mass by family (share of ‖F‖₂², frozen units) ──")
        for k in first(sort!(collect(keys(fam)); by = k -> -fam[k]), 4)
            println("        $(rpad(round(100*fam[k]/tot, sigdigits=4),8))%  $k")
        end
    end

    if nrm < 1e-8
        println("\n✅ CONVERGED at iteration $it with ‖F‖∞ = $nrm")
        break
    end
    for i in 1:N; x[freecols[i]] += a * d[i]; end
    global prev = nrm
end

gs_end = resid(x)
println("\nfinal raw ‖F‖∞ = $(norm(gs_end, Inf))")

# Report the answer, so a converged solve can be sanity-checked economically.
println("\n── selected results (ratio-type variables, benchmark = 1) ──")
for nm in ("phi", "realwage", "pcap", "xtot", "xlab_o")
    haskey(vars, nm) || continue
    v = vars[nm]
    vals = [ (JuMP.is_fixed(vr) ? JuMP.fix_value(vr) : x[pos[JuMP.index(vr).value]])
             for vr in (v isa AbstractArray ? vec(v) : [v]) ]
    println("   $(rpad(nm,10)) n=$(rpad(length(vals),6)) min=$(rpad(round(minimum(vals),sigdigits=6),12)) " *
            "max=$(rpad(round(maximum(vals),sigdigits=6),12)) mean=$(round(sum(vals)/length(vals),sigdigits=6))")
end
println("\nDone.")
