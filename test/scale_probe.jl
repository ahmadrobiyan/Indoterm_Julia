# Region-scaling probe: build + solve the model at several region counts and
# report the resource profile of each, to decide EMPIRICALLY whether Phase 3
# condensation is needed for the full 34-region model and which lever to use.
#
#   julia --project=IndotermJulia IndotermJulia/test/scale_probe.jl 6 10 14 20 34
#
# The Steps 0-4 data pipeline is region-independent (it produces the 34-region
# `agg`), costs ~750 s, and is therefore run ONCE and reused for every size —
# re-running it per size, as an earlier version did, dominated the measurement.
# Each size prints one `METRIC` line and flushes, so a crash at a large size
# still leaves the smaller data points on disk.
#
# `xtrad` scales as nr^2 and `asuppmar`/`MARS` as nr^3, so a naive 6→34
# extrapolation predicts ~20x variables. But the zero-flow fraction also grows
# with region count (small province pairs simply do not trade), and the
# zero-flow guards fix those cells out of the Jacobian. This measures which
# effect wins — that is the whole point.
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP
using SparseArrays, LinearAlgebra
const MOI = JuMP.MOI

sizes = length(ARGS) >= 1 ? parse.(Int, ARGS) : [6, 10, 14, 20]

"""
Contiguous partition of the 34 provinces into `n` groups. REG is already ordered
by island group, so contiguous buckets keep geography — and therefore the trade
sparsity pattern — realistic. `n == 6` returns the canonical island-group map so
the probe is directly comparable with the validation model.
"""
function region_map(n::Int)
    n == 6 && return REG_MAP_34_to_6
    n >= 34 && return collect(1:34)
    edges = [round(Int, i * 34 / n) for i in 0:n]
    map = zeros(Int, 34)
    for g in 1:n, p in (edges[g]+1):edges[g+1]
        map[p] = g
    end
    @assert all(>(0), map) "region map has unassigned provinces"
    return map
end

# NOTE: every accumulator here lives inside a function on purpose. Julia's
# `try`/`for` blocks introduce a new scope, so assigning to an outer variable
# from inside one silently creates a local — an earlier version of this probe
# reported jac_nnz = -1 for exactly that reason, with no error raised.
"""
Jacobian nonzeros, LU fill-in and factorization time of the equilibrated
Jacobian at the given point. Returns `(jac_nnz, lu_nnz, lu_seconds)`, or
`(-1, -1, -1.0)` if the probe could not run.
"""
function jacobian_profile(m, vars)
    try
        nowfree = VariableRef[]
        for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
            JuMP.is_fixed(vr) || push!(nowfree, vr)
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
        mx = zeros(Float64, NC)
        for k in eachindex(jr); a = abs(vals[k]); a > mx[jr[k]] && (mx[jr[k]] = a); end
        rs = [mx[i] > 1e-12 ? mx[i] : 1.0 for i in 1:NC]
        J = sparse(jr, jc, [vals[k]/rs[jr[k]] for k in eachindex(jr)], NC, N)
        t0 = time(); Fac = lu(J); lus = time() - t0
        return (nnz(J), nnz(Fac.L) + nnz(Fac.U), lus)
    catch e
        println("  [jacobian/LU probe failed: $(sprint(showerror, e))]")
        return (-1, -1, -1.0)
    end
end

"""Build, solve and profile the model at `n` regions. Prints one METRIC line."""
function probe(agg_full, n::Int)
    println("\n", "="^70); println("SCALE PROBE — $n regions"); println("="^70)
    agg = aggregate_regions!(agg_full, region_map(n), n)
    params = prepare_parameters!(agg)

    t0 = time(); m, vars = build_model_full!(agg, params); build_s = time() - t0
    nv = num_variables(m)
    nc = num_constraints(m; count_variable_in_set_constraints=false)
    initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

    t0 = time(); res = solve_newton!(m, vars; verbose=false); solve_s = time() - t0
    jnnz, lunnz, lus = jacobian_profile(m, vars)
    rss_mb = try Sys.maxrss() / 2^20 catch; -1.0 end

    println("METRIC nregions=$n vars=$nv cons=$nc free=$(res.n_free) eqs=$(res.n_eqs) " *
            "status=$(res.status) iters=$(res.iters) resid=$(res.residual) " *
            "build_s=$(round(build_s,digits=1)) solve_s=$(round(solve_s,digits=1)) " *
            "jac_nnz=$jnnz lu_nnz=$lunnz lu_s=$(round(lus,digits=2)) " *
            "fill_ratio=$(jnnz>0 ? round(lunnz/jnnz,digits=2) : -1) " *
            "maxrss_mb=$(round(rss_mb,digits=0))")
    flush(stdout)
    return nothing
end

function main(sizes)
    t_pipe = time()
    nat = read_national_data(); regsup = read_regsupp_data(); distgone = read_distgone_data()
    reg0_r = build_reg0!(nat, regsup)
    reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
    reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
    ras_r  = ras_balance!(reg2_r.reg2)
    pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone)
    prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast)
    agg_full = aggregate_model!(prem_r.premod).agg
    println("pipeline (Steps 0-4, region-independent): $(round(time()-t_pipe, digits=1))s")
    flush(stdout)

    for n in sizes
        try
            probe(agg_full, n)
        catch e
            println("METRIC nregions=$n FAILED: $(sprint(showerror, e))")
            flush(stdout)
        end
        GC.gc()
    end
end

main(sizes)
