"""
Krylov-vs-condensation assessment, measurement 1: how much of the LU fill-in is the
ORDERING, and how dense are the densest rows/columns?

`test/scale_probe.jl` measured fill ratios of 15 → 22 → 55 → 70 at 6/10/14/20 regions with
Julia's default `lu` (UMFPACK, AMD/COLAMD ordering, auto strategy). Before deciding between a
Krylov rewrite (days) and Excerpt 49 condensation (weeks), this checks the zero-cost lever:
the same factorization under every ordering/strategy UMFPACK offers, plus the row/column
degree profile that tells whether a few aggregation rows are poisoning the elimination tree.

Run:  julia --project=. test/scratch/_probe_lu_ordering.jl [nregions...]   (default 6 10)
Log:  logs/lu_ordering_<date>.log
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra, Printf, Dates
const MOI = JuMP.MOI
const U = SparseArrays.UMFPACK

sizes = length(ARGS) >= 1 ? parse.(Int, ARGS) : [6, 10]
say(msg) = (println(msg); flush(stdout))

function region_map(n::Int)
    n == 6 && return REG_MAP_34_to_6
    n >= 34 && return collect(1:34)
    edges = [round(Int, i * 34 / n) for i in 0:n]
    map = zeros(Int, 34)
    for g in 1:n, p in (edges[g]+1):edges[g+1]; map[p] = g; end
    map
end

# Free-column Jacobian at the benchmark, row-scaled — same construction as scale_probe.jl.
# Also returns the variable-family name of every free column and the equation index → row.
function free_jacobian(m, vars)
    nowfree = VariableRef[]; fam = String[]
    for (nm, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) && continue
        push!(nowfree, vr); push!(fam, nm)
    end
    allv = JuMP.all_variables(m)
    allidx = [JuMP.index(v) for v in allv]
    pos = Dict(JuMP.index(v).value => k for (k, v) in enumerate(allv))
    nl = MOI.Nonlinear.Model(); rhs = Float64[]
    for (F, S) in list_of_constraint_types(m)
        F <: VariableRef && continue
        for con in all_constraints(m, F, S)
            co = JuMP.constraint_object(con); s = co.set
            r = s isa MOI.EqualTo ? s.value : 0.0
            MOI.Nonlinear.add_constraint(nl, JuMP.moi_function(co.func), MOI.EqualTo(r))
            push!(rhs, r)
        end
    end
    ev = MOI.Nonlinear.Evaluator(nl, MOI.Nonlinear.SparseReverseMode(), allidx)
    MOI.initialize(ev, [:Jac])
    x = [JuMP.is_fixed(v) ? JuMP.fix_value(v) :
         (JuMP.start_value(v) === nothing ? 1.0 : JuMP.start_value(v)) for v in allv]
    st = MOI.jacobian_structure(ev)
    jr_all = getindex.(st, 1); jc_all = getindex.(st, 2)
    Jv = zeros(length(jr_all)); MOI.eval_constraint_jacobian(ev, Jv, x)
    freecols = [pos[JuMP.index(v).value] for v in nowfree]
    c2f = zeros(Int, length(allv))
    for (i, c) in enumerate(freecols); c2f[c] = i; end
    keep = [c2f[c] != 0 for c in jc_all]
    jr = jr_all[keep]; jc = [c2f[c] for c in jc_all[keep]]; vals = Jv[keep]
    N = length(nowfree); NC = length(rhs)
    mx = zeros(NC)
    for k in eachindex(jr); a = abs(vals[k]); a > mx[jr[k]] && (mx[jr[k]] = a); end
    rs = [mx[i] > 1e-12 ? mx[i] : 1.0 for i in 1:NC]
    J = sparse(jr, jc, [vals[k]/rs[jr[k]] for k in eachindex(jr)], NC, N)
    dropzeros!(J)
    J, fam
end

function degree_profile(J, fam)
    NC, N = size(J)
    rowdeg = zeros(Int, NC); coldeg = zeros(Int, N)
    rv = rowvals(J)
    for j in 1:N, k in nzrange(J, j); rowdeg[rv[k]] += 1; coldeg[j] += 1; end
    say(@sprintf("   rows %d  cols %d  nnz %d  avg row deg %.1f  avg col deg %.1f",
                 NC, N, nnz(J), nnz(J)/NC, nnz(J)/N))
    for thr in (50, 200, 1000, 5000)
        nr_ = count(>=(thr), rowdeg); nc_ = count(>=(thr), coldeg)
        say(@sprintf("   deg ≥ %5d : %6d rows (%.2f%%)  %6d cols (%.2f%%)", thr, nr_, 100nr_/NC, nc_, 100nc_/N))
    end
    say("   top-10 row degrees: " * join(string.(sort(rowdeg; rev=true)[1:10]), " "))
    # densest columns, by variable family
    ord = sortperm(coldeg; rev=true)[1:10]
    say("   top-10 col degrees: " * join(["$(fam[j])=$(coldeg[j])" for j in ord], " "))
    # nnz share of the 1% densest rows
    srt = sort(rowdeg; rev=true); k = max(1, NC ÷ 100)
    say(@sprintf("   densest 1%% of rows (%d) hold %.1f%% of nnz", k, 100sum(srt[1:k])/nnz(J)))
end

# UMFPACK control indices (Julia 1-based = C 0-based + 1)
const STRATEGY_IDX = 6   # UMFPACK_STRATEGY (C index 5)
const ORDERING_IDX = U.JL_UMFPACK_ORDERING
# values from umfpack.h
const ORDERINGS  = [("CHOLMOD", 0.0), ("AMD", 1.0), ("METIS", 3.0), ("BEST", 4.0)]
const STRATEGIES = [("auto", 0.0), ("unsym", 1.0), ("sym", 3.0)]

function try_lu(J, ordering, strategy)
    c = U.get_umfpack_control(Float64, Int)
    c[ORDERING_IDX] = ordering; c[STRATEGY_IDX] = strategy
    GC.gc()
    t0 = time()
    F = lu(J; control = c)
    dt = time() - t0
    lunnz = nnz(F.L) + nnz(F.U)
    # residual of one solve as a sanity check that the factorization is usable
    b = J * ones(size(J, 2)); x = F \ b
    (lunnz, dt, norm(x .- 1, Inf))
end

function main(sizes)
    say("\n" * "="^78)
    say("LU ordering / degree probe   $(Dates.now())")
    say("="^78)
    agg34, _ = cached_pipeline(34)
    for n in sizes
        say("\n" * "─"^78); say("nregions = $n"); say("─"^78)
        agg = n == 34 ? agg34 : aggregate_regions!(agg34, region_map(n), n)
        params = prepare_parameters!(agg)
        t0 = time(); m, vars = build_model_full!(agg, params)
        initialize_model!(m, vars; bmk_levels = benchmark_levels(params))
        say(@sprintf("   built in %.0fs", time() - t0))
        J, fam = free_jacobian(m, vars)
        m = nothing; GC.gc()
        degree_profile(J, fam)
        say("")
        say(@sprintf("   %-8s %-6s %14s %8s %8s %10s", "ordering", "strat", "LU nnz", "fill", "sec", "solve err"))
        for (oname, o) in ORDERINGS, (sname, s) in STRATEGIES
            try
                lunnz, dt, err = try_lu(J, o, s)
                say(@sprintf("   %-8s %-6s %14d %8.2f %8.2f %10.1e", oname, sname, lunnz, lunnz/nnz(J), dt, err))
            catch e
                say(@sprintf("   %-8s %-6s   FAILED: %s", oname, sname, sprint(showerror, e)))
            end
            flush(stdout)
        end
        J = nothing; GC.gc()
    end
    say("\nDone  $(Dates.now())")
end

main(sizes)
