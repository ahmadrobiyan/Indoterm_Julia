"""
Stage 1 probe: Jacobian-level condensation (exact Schur complement).

Eliminates the Excerpt-49 substitute set S at the linear-solver level:
  J = [D B; A K]  ->  K̃ = K - A·D⁻¹·B, factorize K̃ with default AMD lu.
Measures whether nnz(LU(K̃)) at 6/10/14 regions extrapolates to ≤250M at 34.

Row matching: Hopcroft-Karp maximum bipartite matching of S-columns to rows
on the stored Jacobian pattern, exploring each column's rows in
descending-|entry| order so matched rows prefer each variable's own defining
equation (shared aggregation rows are the stealing hazard for greedy).
Then a topological order of the S-dependency graph makes D block-triangular;
D⁻¹·B is a sparse triangular solve, never a factorization.
If matching is imperfect, S shrinks to the matched subset (reported per
family). If the S-graph has cycles, the cyclic nodes drop back to the core and
the acyclic remainder is measured (families reported); only a fully-cyclic S
is skipped.

Run:  julia --project=. test/scratch/_probe_condense.jl [nregions...]   (default 6)
Log:  logs/condense_<date>.log   (run one size per process; stdout is flushed)
No src/ changes. Scratch + logs only.
"""

include(joinpath(@__DIR__, "..", "pipeline_cache.jl"))
using JuMP, SparseArrays, LinearAlgebra, Printf, Dates

const MOI = JuMP.MOI

sizes = length(ARGS) >= 1 ? parse.(Int, ARGS) : [6]
say(msg) = (println(msg); flush(stdout))

function region_map(n::Int)
    n == 6 && return REG_MAP_34_to_6
    n >= 34 && return collect(1:34)
    edges = [round(Int, i * 34 / n) for i in 0:n]
    map = zeros(Int, 34)
    for g in 1:n, p in (edges[g]+1):edges[g+1]; map[p] = g; end
    map
end

# Free-column Jacobian at the benchmark, row-scaled — verbatim from _probe_lu_ordering.jl.
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

# ── Substitute sets (TERM.TAB Excerpt 49; Julia `vars` keys) ────────────────

const PRIMARY_S_VARS = ["xtradmar", "xtrad", "pdelivrd", "xsuppmar"]
const TRADENEST_S_VARS = ["xtradmar"]
const FULL_S_VARS = [
    "aint_s", "delTAXint", "pmake", "ppur", "ppur_s", "xint_s", "xinvi", "xmake", "xtradmar", "xint",
    "contCPI", "continccom", "contincind_d", "contMainMacro", "contnatxtot", "contxprim_i", "delPTX",
    "delTAXexp", "delTAXgov", "delTAXhou", "delTAXinv", "pcst", "pdelivrd", "pinvest", "plab",
    "plab_o", "pprim", "psuppmar_p", "realwage", "wlab_i", "xgov", "xhouh_s", "xhou_s", "xinv",
    "xinv_s", "xlab", "xlab_o", "xlux", "xprim", "xrowdem", "xrowdem_d", "xstocks", "xsub",
    "xsuppmar", "xsuppmar_d", "xtrad_d"
]

# ── Matching + ordering machinery ──────────────────────────────────────────

function scol_census(fam::Vector{String}, names::Vector{String})
    want = Set(names)
    Scols = Int[]
    per = Dict{String,Int}()
    for (i, nm) in enumerate(fam)
        if nm in want
            push!(Scols, i)
            per[nm] = get(per, nm, 0) + 1
        end
    end
    return Scols, per
end

# Adjacency: for each S-column, its stored rows sorted by descending |entry|.
function build_adj(J::SparseMatrixCSC{Float64,Int}, Scols::Vector{Int})
    rv = rowvals(J); nzv = nonzeros(J)
    adj = Vector{Vector{Int}}(undef, length(Scols))
    for (i, c) in enumerate(Scols)
        rng = nzrange(J, c)
        rows = rv[rng]; vs = nzv[rng]
        adj[i] = rows[sortperm(abs.(vs); rev=true)]
    end
    return adj
end

# Hopcroft-Karp maximum bipartite matching, left = S-columns, right = rows.
# Accepts an optional initial (valid) matching for tiered runs.
function hopcroft_karp(adj::Vector{Vector{Int}}, nrows::Int,
                       pairU::Vector{Int}=zeros(Int, length(adj)),
                       pairV::Vector{Int}=zeros(Int, nrows))
    nU = length(adj)
    dist = zeros(Int, nU)
    INF = typemax(Int)
    function bfs()
        q = Int[]
        for u in 1:nU
            if pairU[u] == 0
                dist[u] = 0; push!(q, u)
            else
                dist[u] = INF
            end
        end
        found = false; head = 1
        while head <= length(q)
            u = q[head]; head += 1
            for v in adj[u]
                pu = pairV[v]
                if pu == 0
                    found = true
                elseif dist[pu] == INF
                    dist[pu] = dist[u] + 1; push!(q, pu)
                end
            end
        end
        return found
    end
    function dfs(u)
        for v in adj[u]
            pu = pairV[v]
            if pu == 0 || (dist[pu] == dist[u] + 1 && dfs(pu))
                pairU[u] = v; pairV[v] = u
                return true
            end
        end
        dist[u] = INF
        return false
    end
    matching = count(!iszero, pairU)
    while bfs()
        for u in 1:nU
            if pairU[u] == 0 && dfs(u)
                matching += 1
            end
        end
    end
    return pairU, pairV, matching
end

# Topological order of the S-dependency graph of D (edge j->i when the
# defining row of S-col i also contains S-col j). Returns (order, cyclic).
function topo_order(D::SparseMatrixCSC{Float64,Int}; tol=1e-12)
    n = size(D, 1)
    indeg = zeros(Int, n)
    children = [Int[] for _ in 1:n]
    rv = rowvals(D); nzv = nonzeros(D)
    cutoff = tol * max(maximum(abs, nzv), 1e-300)
    for j in 1:n
        for k in nzrange(D, j)
            i = rv[k]
            if i != j && abs(nzv[k]) > cutoff
                push!(children[j], i); indeg[i] += 1
            end
        end
    end
    q = [i for i in 1:n if indeg[i] == 0]
    order = Int[]; head = 1
    while head <= length(q)
        j = q[head]; head += 1; push!(order, j)
        for i in children[j]
            indeg[i] -= 1
            indeg[i] == 0 && push!(q, i)
        end
    end
    if length(order) == n
        return order, Int[]
    else
        return order, [i for i in 1:n if indeg[i] > 0]
    end
end

# Strict-upper mass of a square sparse matrix (for triangularity check).
function strict_upper_mass(A::SparseMatrixCSC{Float64,Int})
    n = size(A, 1); rv = rowvals(A); nzv = nonzeros(A)
    cnt = 0; mx = 0.0
    for j in 1:n
        for k in nzrange(A, j)
            if rv[k] < j
                cnt += 1
                a = abs(nzv[k]); a > mx && (mx = a)
            end
        end
    end
    return cnt, mx
end

# Sparse forward substitution Dp*Y = Bp with lower-triangular Dp (stored
# pattern). Dense workspace + touched-list sparse accumulator; Y stays sparse.
function tri_solve_sparse(Dp::SparseMatrixCSC{Float64,Int}, Bp::SparseMatrixCSC{Float64,Int})
    nS = size(Dp, 1); nC = size(Bp, 2)
    @assert size(Dp, 2) == nS && size(Bp, 1) == nS
    Drv = rowvals(Dp); Dnz = nonzeros(Dp)
    diag = zeros(Float64, nS)
    for j in 1:nS
        for k in nzrange(Dp, j)
            if Drv[k] == j; diag[j] = Dnz[k]; break; end
        end
        abs(diag[j]) < 1e-14 && error("tri_solve_sparse: (near-)zero diagonal at $j")
    end
    Brv = rowvals(Bp); Bnz = nonzeros(Bp)
    y = zeros(Float64, nS); mark = zeros(Int, nS); epoch = 0
    I = Int[]; Jc = Int[]; V = Float64[]
    sizehint!(I, min(nnz(Bp) * 2, 10_000_000)); sizehint!(Jc, min(nnz(Bp) * 2, 10_000_000)); sizehint!(V, min(nnz(Bp) * 2, 10_000_000))
    tlist = Int[]
    for t in 1:nC
        epoch += 1
        empty!(tlist)
        for k in nzrange(Bp, t)
            i = Brv[k]
            if mark[i] != epoch; mark[i] = epoch; push!(tlist, i); end
            y[i] += Bnz[k]
        end
        sort!(tlist)
        idx = 1
        while idx <= length(tlist)
            i = tlist[idx]; idx += 1
            y[i] /= diag[i]
            iszero(y[i]) && continue
            yi = y[i]
            for k in nzrange(Dp, i)
                r = Drv[k]; r <= i && continue
                if mark[r] != epoch
                    mark[r] = epoch
                    pos = searchsortedfirst(tlist, r, idx, length(tlist), Base.Order.Forward)
                    insert!(tlist, pos, r)
                end
                y[r] -= Dnz[k] * yi
            end
        end
        for i in tlist
            v = y[i]; y[i] = 0.0
            if v != 0.0; push!(I, i); push!(Jc, t); push!(V, v); end
        end
    end
    return sparse(I, Jc, V, nS, nC)
end

# Dense-vector forward substitution with lower-triangular Dp.
function tri_solve_vec(Dp::SparseMatrixCSC{Float64,Int}, b::Vector{Float64})
    n = length(b); x = copy(b)
    Drv = rowvals(Dp); Dnz = nonzeros(Dp)
    for j in 1:n
        r = nzrange(Dp, j); d = 0.0
        for k in r; if Drv[k] == j; d = Dnz[k]; break; end; end
        abs(d) < 1e-14 && error("tri_solve_vec: (near-)zero diagonal at $j")
        x[j] /= d
        xj = x[j]
        iszero(xj) || for k in r
            rr = Drv[k]; rr > j && (x[rr] -= Dnz[k] * xj)
        end
    end
    return x
end

gb(x) = x / 1e9

function run_set(J::SparseMatrixCSC{Float64,Int}, fam::Vector{String},
                 setname::String, setvars::Vector{String}, Fref, nreg::Int)
    N = size(J, 2); NC = size(J, 1)
    say("\nSet: $setname")
    Scols_all, census = scol_census(fam, setvars)
    for nm in sort(collect(keys(census)))
        say(@sprintf("   census %-14s free cols %8d", nm, census[nm]))
    end
    say(@sprintf("   S total free cols: %d", length(Scols_all)))
    if isempty(Scols_all)
        say("   No free variables in this set. Skipping.")
        return
    end
    # S-degree of each row over this set's S-columns: true defining rows have
    # O(1) S-entries; shared aggregation rows have 625+ (sums over c,s) or nr
    # (sums over one region index). Tier 1 excludes everything wide; tier 2
    # relaxes to single-region-index sums for leftover columns.
    t0 = time()
    sdeg = zeros(Int, NC)
    rvJ = rowvals(J)
    for c in Scols_all
        for k in nzrange(J, c)
            sdeg[rvJ[k]] += 1
        end
    end
    adj_full = build_adj(J, Scols_all)
    K1 = 6; K2 = 4 * nreg
    adj1 = [filter(r -> sdeg[r] <= K1, rows) for rows in adj_full]
    pairU = zeros(Int, length(Scols_all)); pairV = zeros(Int, NC)
    pairU, pairV, m1 = hopcroft_karp(adj1, NC, pairU, pairV)
    say(@sprintf("   HK tier1 (sdeg<=%d): matched %d / %d (%.2fs)", K1, m1, length(Scols_all), time() - t0))
    m = m1
    if m1 < length(Scols_all)
        t0 = time()
        adj2 = [pairU[i] != 0 ? Int[] : filter(r -> sdeg[r] <= K2, adj_full[i]) for i in 1:length(Scols_all)]
        pairU, pairV, m = hopcroft_karp(adj2, NC, pairU, pairV)
        say(@sprintf("   HK tier2 (sdeg<=%d): matched %d / %d (%.2fs)", K2, m, length(Scols_all), time() - t0))
    end
    maxsdeg = 0; widematched = 0
    for i in 1:length(Scols_all)
        if pairU[i] != 0
            d = sdeg[pairU[i]]
            d > maxsdeg && (maxsdeg = d)
            d > K1 && (widematched += 1)
        end
    end
    say(@sprintf("   matched rows sdeg: max=%d  rows with sdeg>K1: %d", maxsdeg, widematched))
    # Pivot-dominance floor: a matched pair whose coefficient is ~0 in
    # row-scaled units is a shared/wrong row, not the defining equation —
    # keeping it puts a ~1e-14 pivot on D's diagonal, which then poisons
    # D⁻¹B (×1e14 amplification), K̃, and the verify. True defining rows
    # carry O(1) coefficients. Drop offenders back to the core, per family.
    DFLOOR = 1e-6
    dropped_dom = 0
    dropped_dom_per = Dict{String,Int}()
    for i in 1:length(Scols_all)
        r = pairU[i]
        if r != 0 && abs(J[r, Scols_all[i]]) < DFLOOR
            pairU[i] = 0; pairV[r] = 0; dropped_dom += 1
            nm = fam[Scols_all[i]]
            dropped_dom_per[nm] = get(dropped_dom_per, nm, 0) + 1
        end
    end
    say(@sprintf("   dominance floor %.0e: dropped %d matched pairs", DFLOOR, dropped_dom))
    for nm in sort(collect(keys(dropped_dom_per)))
        say(@sprintf("     %-14s dropped %d", nm, dropped_dom_per[nm]))
    end
    # per-family matched report
    matched_per = Dict{String,Int}()
    for (i, c) in enumerate(Scols_all)
        if pairU[i] != 0
            nm = fam[c]
            matched_per[nm] = get(matched_per, nm, 0) + 1
        end
    end
    for nm in sort(collect(keys(census)))
        got = get(matched_per, nm, 0)
        flag = got < census[nm] ? "   <-- DROPPED $(census[nm]-got)" : ""
        say(@sprintf("   matched %-14s %8d / %-8d%s", nm, got, census[nm], flag))
    end
    keep = findall(!iszero, pairU)
    S = Scols_all[keep]; R = pairU[keep]
    if length(S) < length(Scols_all)
        say(@sprintf("   S shrunk: %d -> %d cols", length(Scols_all), length(S)))
    end
    inS = falses(N); inS[S] .= true; C = findall(.!inS)
    inR = falses(NC); inR[R] .= true; restR = findall(.!inR)
    say(@sprintf("   core |C|=%d  rest rows=%d", length(C), length(restR)))
    D = J[R, S]; B = J[R, C]; A = J[restR, S]; K = J[restR, C]
    say(@sprintf("   blocks nnz: D=%d B=%d A=%d K=%d", nnz(D), nnz(B), nnz(A), nnz(K)))
    t0 = time()
    order, cyclic = topo_order(D)
    say(@sprintf("   topo sort: %.2fs", time() - t0))
    if !isempty(cyclic)
        cycfams = sort(unique(fam[S[cyclic]]))
        say(@sprintf("   CYCLIC S-graph: %d nodes; families: %s", length(cyclic), join(cycfams, ",")))
        cycset = Set(cyclic)
        keepidx = [i for i in 1:length(S) if !(i in cycset)]
        if isempty(keepidx)
            say("   Entire S cyclic. Skipping this set.")
            return
        end
        say(@sprintf("   dropping %d cyclic nodes to core; proceeding with %d acyclic S-cols",
                     length(cyclic), length(keepidx)))
        S = S[keepidx]; R = R[keepidx]
        inS = falses(N); inS[S] .= true; C = findall(.!inS)
        inR = falses(NC); inR[R] .= true; restR = findall(.!inR)
        say(@sprintf("   core |C|=%d  rest rows=%d", length(C), length(restR)))
        D = J[R, S]; B = J[R, C]; A = J[restR, S]; K = J[restR, C]
        say(@sprintf("   blocks nnz: D=%d B=%d A=%d K=%d", nnz(D), nnz(B), nnz(A), nnz(K)))
        t0 = time()
        order, cyclic = topo_order(D)
        say(@sprintf("   topo sort (acyclic subset): %.2fs", time() - t0))
        if !isempty(cyclic)
            say("   STILL cyclic after drop. Skipping this set (needs positional defining-row blocks).")
            return
        end
    end
    P = order; invP = invperm(P)
    Dp = D[P, P]
    cnt, mx = strict_upper_mass(Dp)
    dvals = [Dp[j, j] for j in 1:size(Dp, 1)]
    dmin = minimum(abs, dvals)
    dmax = maximum(abs, nonzeros(Dp))
    nDub6 = count(a -> abs(a) < 1e-6, dvals); nDub3 = count(a -> abs(a) < 1e-3, dvals)
    say(@sprintf("   D diag spectrum: <1e-6: %d  <1e-3: %d  of %d", nDub6, nDub3, length(dvals)))
    weakfam = Dict{String,Int}()
    for i in 1:length(dvals)
        if abs(dvals[i]) < 1e-3
            nm = fam[S[P[i]]]
            weakfam[nm] = get(weakfam, nm, 0) + 1
        end
    end
    for nm in sort(collect(keys(weakfam)))
        say(@sprintf("     weak-pivot %-14s %d", nm, weakfam[nm]))
    end
    say(@sprintf("   D permuted: strict-upper cnt=%d max=%.2e  diag min=%.2e max=%.2e", cnt, mx, dmin, dmax))
    if cnt > 0 && mx > 1e-12 * dmax
        say("   D not triangularisable by topo order. Skipping this set.")
        return
    end
    t0 = time()
    local Yp
    try
        Yp = tri_solve_sparse(Dp, B[P, :])
    catch e
        say("   tri solve failed: $(sprint(showerror, e)). Skipping this set.")
        return
    end
    say(@sprintf("   DinvB: %.2fs  nnz(Y)=%d (B nnz=%d)", time() - t0, nnz(Yp), nnz(B)))
    # dust cleanup on Y (relative), then form K̃ sparsely
    ymax = maximum(abs, nonzeros(Yp); init=0.0)
    if ymax > 0
        drop_tol = 1e-14 * ymax
        I2, J2, V2 = findnz(Yp)
        keep2 = abs.(V2) .> drop_tol
        Yp = sparse(I2[keep2], J2[keep2], V2[keep2], size(Yp, 1), size(Yp, 2))
        say(@sprintf("   Y after 1e-14 dust drop: nnz=%d", nnz(Yp)))
    end
    Y = Yp[invP, :]
    t0 = time()
    Kt = K - A * Y
    dropzeros!(Kt)
    say(@sprintf("   Ktilde formed: %.2fs  nnz(K)=%d nnz(Ktilde)=%d  Ktilde/K=%.2f", time() - t0, nnz(K), nnz(Kt), nnz(Kt) / max(nnz(K), 1)))
    GC.gc()
    t0 = time()
    # Singularity diagnostics first: exact singularity of K̃ with nonsingular J
    # and D is impossible, so a SingularException names a premise — find which.
    kt_r_empty = 0; kt_c_empty = 0
    for j in 1:size(Kt, 2); if nnz(Kt[:, j]) == 0; kt_c_empty += 1; end; end
    for i in 1:size(Kt, 1); if nnz(Kt[i, :]) == 0; kt_r_empty += 1; end; end
    say(@sprintf("   Ktilde struct: empty cols=%d empty rows=%d  diag min=%.2e",
                 kt_c_empty, kt_r_empty, minimum(abs, [Kt[j, j] for j in 1:size(Kt, 1)])))
    local Fk
    try
        Fk = lu(Kt)
    catch e
        say(@sprintf("   lu(Ktilde) threw %s — trying shifted lu(Ktilde+1e-12 I) to classify structural vs numeric",
                     typeof(e)))
        Fk = lu(Kt + 1e-12 * sparse(I, size(Kt, 1), size(Kt, 2)))
        say("   shifted factorization succeeded: zero-pivot is numeric, not structural")
    end
    dt = time() - t0
    lunnz = nnz(Fk.L) + nnz(Fk.U)
    say(@sprintf("   LU(Ktilde): nnz=%-12d fill=%-8.2f time=%.2fs  RSS=%.2f GB", lunnz, lunnz / max(nnz(Kt), 1), dt, gb(Sys.maxrss())))
    # verification: Schur path vs direct lu(J)\b
    b = Vector(J * ones(N))
    bR = b[R]; brest = b[restR]
    w = tri_solve_vec(Dp, bR[P]); w0 = w[invP]
    rhsC = brest - A * w0
    xC = Fk \ rhsC
    u = bR - B * xC
    xSp = tri_solve_vec(Dp, u[P]); xS = xSp[invP]
    xfull = zeros(N); xfull[S] .= xS; xfull[C] .= xC
    xd = Fref \ b
    err = norm(xfull .- xd, Inf) / norm(xd, Inf)
    say(@sprintf("   verify: max rel err Schur vs direct = %.2e", err))
end

function main(sizes)
    say("\n" * "="^80)
    say("Jacobian Condensation Probe (Stage 1)   $(Dates.now())")
    say("="^80)
    agg34, _ = cached_pipeline(34)
    for n in sizes
        say("\n" * "-"^80); say("nregions = $n"); say("-"^80)
        agg = n == 34 ? agg34 : aggregate_regions!(agg34, region_map(n), n)
        params = prepare_parameters!(agg)
        t0 = time()
        m, vars = build_model_full!(agg, params)
        initialize_model!(m, vars; bmk_levels=benchmark_levels(params))
        say(@sprintf("   built+initialized in %.0fs  RSS=%.2f GB", time() - t0, gb(Sys.maxrss())))
        J, fam = free_jacobian(m, vars)
        m = nothing; GC.gc()
        say(@sprintf("   J: %d x %d  nnz=%d  RSS=%.2f GB", size(J, 1), size(J, 2), nnz(J), gb(Sys.maxrss())))
        @assert size(J, 1) == size(J, 2) "non-square free Jacobian"
        t0 = time()
        Fref = lu(J)
        say(@sprintf("   LU(J) baseline: nnz=%d fill=%.2f time=%.2fs  RSS=%.2f GB",
                     nnz(Fref.L) + nnz(Fref.U), (nnz(Fref.L) + nnz(Fref.U)) / nnz(J),
                     time() - t0, gb(Sys.maxrss())))
        run_set(J, fam, "PRIMARY", PRIMARY_S_VARS, Fref, n)
        GC.gc()
        run_set(J, fam, "TRADENEST", TRADENEST_S_VARS, Fref, n)
        GC.gc()
        run_set(J, fam, "FULL", FULL_S_VARS, Fref, n)
        J = nothing; Fref = nothing; GC.gc()
    end
    say("\nDone  $(Dates.now())")
end

main(sizes)
