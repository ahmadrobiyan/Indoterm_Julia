"""
Is `swap houslack = shrBoTnom2` an EXACT rank deficiency, and what kills ALL SIX?

`test/diag_swap_bisect.jl` isolated one culprit: of TERM.CMF's six active swaps,
five clear a 1e-5 probe in 2-3 Newton iterations, and `houslack = shrBoTnom2`
alone reproduces the whole ALL-SIX failure (maxit at 20, ‖F‖∞ ~ 3e-6) — while
leaving the benchmark bit-identical at 1.4260876923799515e-9.

`TERM.TAB:2047-2050` says why that combination is special:

    E_fhou  (all,h,HOU)(all,d,DST)  whouhtot(h,d) = wlab_io(d) + fhou(h,d)  + houslack;
    E_fhou2 (all,h,HOU)(all,d,DST)  whouhtot(h,d) = wgdpexp(d) + fhou2(h,d) + houslack;

`houslack` enters both blocks additively, with coefficient 1 in every row, right
alongside the per-region shifters. So its Jacobian column is *exactly*

    col(houslack) = Σ_d col(fhou[d]) + Σ_d col(fhou2[d])

whenever all three are free. Not approximately — the entries are literal 1s from
an additive term, so the cancellation is exact in floating point. While
`houslack` is exogenous this is harmless. Free it and the Jacobian is singular:
the benchmark remains a solution (it always was), but Newton has no direction to
move in. That is precisely the observed signature — bit-identical benchmark,
zero accepted continuation steps, residual falling only like h^1.5 and
NON-monotonically in h.

Two cases are measured, because the theory is clean for one and NOT obviously
right for the other:

  A. houslack swap ALONE. Predicted null vector: houslack=+1, every fhou[d]=-1,
     every fhou2[d]=-1, all else 0. If ‖Jv‖ is at round-off, the diagnosis is
     proven rather than argued.

  B. ALL SIX. Here `xhouhtot = fhou` pins `fhou`, so col(fhou) is gone from the
     free set and the exact dependency above should BREAK — yet ALL SIX failed
     with the same signature. The reduced candidate (houslack=+1, fhou2[d]=-1)
     is therefore predicted NOT to be null. Measuring it either exposes a second
     mechanism or shows the residual dependency is merely near-exact.

The test is a direct null check, not a rank algorithm: `‖Jv‖∞` compared against
`‖ |J|·|v| ‖∞`, the size the row entries would have without cancellation. A ratio
at 1e-16 means exact linear dependence; a ratio of order 1 means no dependence.
No solve is performed — the Jacobian is taken at the benchmark seed, so nothing
here depends on the solver's behaviour.

Run:  julia --project=IndotermJulia IndotermJulia/test/diag_houslack_null.jl
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra
const MOI = JuMP.MOI

agg6, params = cached_pipeline(6)
bmk = benchmark_levels(params)

"""
Assemble the constraint Jacobian over the FREE columns at the current point.

No solve is run first: `solve_newton!` pins orphan variables and deletes rows, so
solving would change the very column set under examination. The benchmark seed is
the point we care about anyway — the swap keeps it an exact solution.
"""
function free_jacobian(m::JuMP.Model, vars::Dict{String,Any})
    allv = JuMP.all_variables(m)
    allidx = [JuMP.index(v) for v in allv]
    pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))

    nowfree = VariableRef[]
    seenid = Set{Int64}()
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) && continue
        iv = vr.index.value
        iv in seenid && continue
        push!(seenid, iv); push!(nowfree, vr)
    end
    N = length(nowfree)
    col2free = zeros(Int, length(allv))
    for (i, vr) in enumerate(nowfree); col2free[pos[JuMP.index(vr).value]] = i; end

    # The RHS constant is irrelevant to a Jacobian, so every row is registered as
    # EqualTo(0.0) rather than reaching for the module-internal `_rhs_of`.
    nlmodel = MOI.Nonlinear.Model()
    NC = 0
    for (Ftype, S) in list_of_constraint_types(m)
        Ftype <: VariableRef && continue
        for con in all_constraints(m, Ftype, S)
            co = JuMP.constraint_object(con)
            MOI.Nonlinear.add_constraint(nlmodel, JuMP.moi_function(co.func), MOI.EqualTo(0.0))
            NC += 1
        end
    end
    ev = MOI.Nonlinear.Evaluator(nlmodel, MOI.Nonlinear.SparseReverseMode(), allidx)
    MOI.initialize(ev, [:Grad, :Jac])

    x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
         (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]
    st = MOI.jacobian_structure(ev)
    rr = getindex.(st, 1); cc = getindex.(st, 2)
    valbuf = zeros(length(rr))
    MOI.eval_constraint_jacobian(ev, valbuf, x)
    keep = [col2free[c] != 0 for c in cc]
    J = sparse(rr[keep], [col2free[c] for c in cc[keep]], valbuf[keep], NC, N)
    return J, nowfree
end

"""Index of each free variable by its JuMP name, for building a candidate vector."""
name_index(nowfree) = Dict(JuMP.name(vr) => i for (i, vr) in enumerate(nowfree))

function null_check(label, subset, want::Vector{String})
    println("\n", "="^72)
    println("══ $label ══")
    println("="^72)
    (m, vars) = build_model_full!(agg6, params)
    initialize_model!(m, vars; bmk_levels=bmk)
    isempty(subset) || apply_swaps!(m, vars, subset)

    J, nowfree = free_jacobian(m, vars)
    idx = name_index(nowfree)
    println("  Jacobian: $(size(J,1)) rows × $(size(J,2)) free columns, nnz=$(nnz(J))")

    v = zeros(size(J, 2))
    missing_names = String[]
    for w in want
        sgn, nm = w[1] == '-' ? (-1.0, w[2:end]) : (1.0, w)
        k = get(idx, nm, 0)
        k == 0 ? push!(missing_names, nm) : (v[k] = sgn)
    end
    if !isempty(missing_names)
        println("  ⚠️  not free (so not in the candidate): $(join(missing_names, ", "))")
    end
    nset = count(!=(0.0), v)
    println("  candidate null vector: $nset non-zero entries")
    if nset == 0
        println("  ⛔ nothing to test — every named variable is exogenous here.")
        return
    end

    r = J * v
    scale = abs.(J) * abs.(v)          # row magnitude WITHOUT cancellation
    rn, sn = norm(r, Inf), norm(scale, Inf)
    ratio = sn == 0 ? NaN : rn / sn
    println("  ‖J·v‖∞          = $rn")
    println("  ‖ |J|·|v| ‖∞    = $sn   (the size the entries have before cancelling)")
    println("  ratio           = $ratio")
    println(ratio < 1e-12 ?
        "  ✅ EXACT null direction — the Jacobian is rank-deficient along it." :
        ratio < 1e-4 ?
        "  ⚠️  near-null (ratio $ratio): not an exact dependency, but a very weak\n" *
        "     direction — enough to stall Newton without making the matrix singular." :
        "  ❌ NOT a null direction — this candidate does not explain the failure.")

    # Name the rows that survive, so a failed prediction is a lead and not a dead end.
    if ratio >= 1e-12
        top = sortperm(abs.(r); rev=true)[1:min(5, length(r))]
        println("  rows where v is NOT absorbed (row index, |residual|):")
        for i in top
            abs(r[i]) == 0 && continue
            println("    row $(i)   $(abs(r[i]))")
        end
    end
end

HOU_FHOU  = ["-fhou[$d]"  for d in 1:6]
HOU_FHOU2 = ["-fhou2[$d]" for d in 1:6]

# A: the singleton. Full predicted vector.
null_check("houslack = shrBoTnom2 ALONE  (predict: EXACT null)",
           Any[("houslack", "shrBoTnom2")],
           vcat(["houslack"], HOU_FHOU, HOU_FHOU2))

# B: all six. `fhou` is pinned by `xhouhtot = fhou`, so only the fhou2 half of the
# dependency can remain — predicted NOT to be null.
null_check("ALL SIX  (predict: NOT null — fhou is pinned by swap 1)",
           collect(TERM_CMF_SWAPS),
           vcat(["houslack"], HOU_FHOU2))

println("\nDone.")
