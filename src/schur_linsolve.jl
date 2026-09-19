"""
Schur-complement linear solve (Stage 2a, SCALE34_PLAN.md §4).

Exact Jacobian-level condensation of the Excerpt-49 PRIMARY substitute set
(`xtradmar`, `xtrad`, `pdelivrd`, `xsuppmar`): partition J = [D B; A K],
solve Y = D⁻¹B by sparse triangular substitution (D is block-triangular under
a topological order — verified, never assumed), factor K̃ = K − A·Y with
default AMD LU, recover the S-block. Plus iterative refinement on the Schur
path, which the Stage-1 probe showed is required (weak pivots cap the raw
Schur solve at ~1e-3; refinement with full-J residuals closes it).

Entry point: `schur_linsolve(J, rhs, fam; verbose) -> Vector or nothing`.
Returns `nothing` on ANY failure (bad match, cyclic S, tri-solve error,
singular K̃) — the caller (`solve_newton!`) then takes its existing CGNR/LU
fallback, so a Schur failure can never corrupt a solve, only cost time.

The plan is rebuilt from the current J on every call. Matching + topo cost
~0.1 s against ~2 s LU at 6 regions; rebuilding avoids every staleness
question (values change per Newton iteration, pattern does not — but
rebuilding is cheaper than proving that).

Provenance: ported from `test/scratch/_probe_condense.jl` (Stage-1 probe),
which measured fill 4.18–9.59 vs 10.97–69.62 baseline at 6–20 regions and
verified the Schur path against direct LU.
"""

const SCHUR_S_VARS = ["xtradmar", "xtrad", "pdelivrd", "xsuppmar"]
const SCHUR_DFLOOR = 1e-3   # floor sweep 2026-09-13: 1e-3 verifies 7.5e-7 (gate-clean) at fill 7.62 vs 10.97 baseline; 1e-6 verifies only 2.9e-3
const SCHUR_NREFINE = 2     # refinement steps on the Schur path

struct SchurPlan
    S::Vector{Int}
    R::Vector{Int}
    P::Vector{Int}
    invP::Vector{Int}
    C::Vector{Int}
    restR::Vector{Int}
end

function _schur_hopcroft_karp(adj::Vector{Vector{Int}}, nrows::Int,
                              pairU::Vector{Int}, pairV::Vector{Int})
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

function _schur_topo_order(D::SparseMatrixCSC{Float64,Int}; tol=1e-12)
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

function _schur_tri_sparse(Dp::SparseMatrixCSC{Float64,Int}, Bp::SparseMatrixCSC{Float64,Int})
    nS = size(Dp, 1); nC = size(Bp, 2)
    Drv = rowvals(Dp); Dnz = nonzeros(Dp)
    diag = zeros(Float64, nS)
    for j in 1:nS
        for k in nzrange(Dp, j)
            if Drv[k] == j; diag[j] = Dnz[k]; break; end
        end
        abs(diag[j]) < 1e-14 && error("schur: (near-)zero triangular diagonal at $j")
    end
    Brv = rowvals(Bp); Bnz = nonzeros(Bp)
    y = zeros(Float64, nS); mark = zeros(Int, nS); epoch = 0
    I = Int[]; Jc = Int[]; V = Float64[]
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

function _schur_tri_vec(Dp::SparseMatrixCSC{Float64,Int}, b::Vector{Float64})
    n = length(b); x = copy(b)
    Drv = rowvals(Dp); Dnz = nonzeros(Dp)
    for j in 1:n
        r = nzrange(Dp, j); d = 0.0
        for k in r; if Drv[k] == j; d = Dnz[k]; break; end; end
        abs(d) < 1e-14 && error("schur: (near-)zero diagonal at $j")
        x[j] /= d
        xj = x[j]
        iszero(xj) || for k in r
            rr = Drv[k]; rr > j && (x[rr] -= Dnz[k] * xj)
        end
    end
    return x
end

"""
    build_schur_plan(J, fam; verbose=false) -> SchurPlan or nothing

Tiered Hopcroft-Karp matching of PRIMARY-substitute columns to rows (narrow
rows first), 1e-6 pivot-dominance floor, cyclic-node drop to core. `fam[i]`
is the variable-family name of free column i. Returns `nothing` when no
usable S remains — never throws.
"""
function build_schur_plan(J::SparseMatrixCSC{Float64,Int}, fam::Vector{String};
                          verbose::Bool=false, dfloor::Float64=SCHUR_DFLOOR)
    log(msg) = verbose && println(msg)
    N = size(J, 2); NC = size(J, 1)
    want = Set(SCHUR_S_VARS)
    Scols_all = [i for i in 1:N if fam[i] in want]
    isempty(Scols_all) && return nothing
    rvJ = rowvals(J); nzv = nonzeros(J)
    sdeg = zeros(Int, NC)
    for c in Scols_all
        for k in nzrange(J, c)
            sdeg[rvJ[k]] += 1
        end
    end
    adj_full = Vector{Vector{Int}}(undef, length(Scols_all))
    for (i, c) in enumerate(Scols_all)
        rng = nzrange(J, c)
        rows = rvJ[rng]; vs = nzv[rng]
        adj_full[i] = rows[sortperm(abs.(vs); rev=true)]
    end
    K1 = 6; K2 = 4 * 6
    adj1 = [filter(r -> sdeg[r] <= K1, rows) for rows in adj_full]
    pairU = zeros(Int, length(Scols_all)); pairV = zeros(Int, NC)
    pairU, pairV, m1 = _schur_hopcroft_karp(adj1, NC, pairU, pairV)
    if m1 < length(Scols_all)
        adj2 = [pairU[i] != 0 ? Int[] : filter(r -> sdeg[r] <= K2, adj_full[i])
                for i in 1:length(Scols_all)]
        pairU, pairV, _ = _schur_hopcroft_karp(adj2, NC, pairU, pairV)
    end
    dropped = 0
    for i in 1:length(Scols_all)
        r = pairU[i]
        if r != 0 && abs(J[r, Scols_all[i]]) < dfloor
            pairU[i] = 0; pairV[r] = 0; dropped += 1
        end
    end
    keep = findall(!iszero, pairU)
    isempty(keep) && return nothing
    S = Scols_all[keep]; R = pairU[keep]
    inS = falses(N); inS[S] .= true; C = findall(.!inS)
    inR = falses(NC); inR[R] .= true; restR = findall(.!inR)
    D = J[R, S]
    order, cyclic = _schur_topo_order(D)
    if !isempty(cyclic)
        cycset = Set(cyclic)
        keepidx = [i for i in 1:length(S) if !(i in cycset)]
        isempty(keepidx) && return nothing
        S = S[keepidx]; R = R[keepidx]
        inS = falses(N); inS[S] .= true; C = findall(.!inS)
        inR = falses(NC); inR[R] .= true; restR = findall(.!inR)
        D = J[R, S]
        order, cyclic = _schur_topo_order(D)
        isempty(cyclic) || return nothing
    end
    P = order; invP = invperm(P)
    log("  schur plan: |S|=" * string(length(S)) * " core=" * string(length(C)) *
        " dropped_dom=" * string(dropped) * " dropped_cyc=" * string(length(cyclic)))
    return SchurPlan(S, R, P, invP, C, restR)
end

function _schur_correct(J::SparseMatrixCSC{Float64,Int}, plan::SchurPlan,
                        Y::SparseMatrixCSC{Float64,Int}, Fk, rhs::Vector{Float64})
    S, R, P, invP, C, restR = plan.S, plan.R, plan.P, plan.invP, plan.C, plan.restR
    B = J[R, C]; A = J[restR, S]
    bR = rhs[R]; brest = rhs[restR]
    Dp_rows = R[P]
    Dp = J[Dp_rows, S[P]]
    w = _schur_tri_vec(Dp, bR[P]); w0 = w[invP]
    rhsC = brest - A * w0
    xC = Fk \ rhsC
    u = bR - B * xC
    xSp = _schur_tri_vec(Dp, u[P]); xS = xSp[invP]
    x = zeros(Float64, size(J, 2)); x[S] .= xS; x[C] .= xC
    return x
end

"""
    schur_linsolve(J, rhs, fam; verbose=false) -> Vector or nothing

One condensed solve with iterative refinement. `nothing` on any failure —
caller falls back to its existing path.
"""
function schur_linsolve(J::SparseMatrixCSC{Float64,Int}, rhs::Vector{Float64},
                        fam::Vector{String}; verbose::Bool=false,
                        dfloor::Float64=SCHUR_DFLOOR)
    log(msg) = verbose && println(msg)
    fz = schur_factorize(J, fam; verbose=verbose, dfloor=dfloor)
    fz === nothing && return nothing
    plan, Y, Fk = fz
    local x
    try
        x = _schur_correct(J, plan, Y, Fk, rhs)
        for _ in 1:SCHUR_NREFINE
            r = rhs - J * x
            rel = norm(r, Inf) / max(norm(rhs, Inf), eps())
            rel < 1e-12 && break
            dx = _schur_correct(J, plan, Y, Fk, r)
            all(isfinite, dx) || break
            xn = x + dx
            rn = rhs - J * xn
            norm(rn, Inf) < norm(r, Inf) || break
            x = xn
        end
    catch e
        log("  schur: solve/refine failed — " * sprint(showerror, e))
        return nothing
    end
    return x
end

"""
    schur_factorize(J, fam; verbose, dfloor) -> (plan, Y, Fk) or nothing

Plan + triangular factor + core LU, shared by `schur_linsolve` and the
bordered solver below (one factorization serves both bordered RHS).
"""
function schur_factorize(J::SparseMatrixCSC{Float64,Int}, fam::Vector{String};
                         verbose::Bool=false, dfloor::Float64=SCHUR_DFLOOR)
    log(msg) = verbose && println(msg)
    # verbose is overloaded as the memory ledger switch by solve_newton! (memlog);
    # default false keeps the production path silent and behavior-identical.
    mem(tag) = verbose && (println("  [mem] schur " * tag * "  peakRSS=" *
                                   string(round(Sys.maxrss() / 2^30; digits=2)) * " GiB");
                            flush(stdout))
    local plan
    try
        plan = build_schur_plan(J, fam; verbose=verbose, dfloor=dfloor)
    catch e
        log("  schur: plan build failed — " * sprint(showerror, e))
        return nothing
    end
    plan === nothing && return nothing
    mem("plan built")
    S, R, P, invP, C, restR = plan.S, plan.R, plan.P, plan.invP, plan.C, plan.restR
    log("  schur: |S|=$(length(S)) |C|=$(length(C)) restR=$(length(restR))")
    try
        B = J[R, C]; A = J[restR, S]
        mem("blocks sliced: nnz(B)=$(nnz(B)) nnz(A)=$(nnz(A)) nnz(D)=$(nnz(J[R, S]))")
        Dp = J[R[P], S[P]]
        Yp = _schur_tri_sparse(Dp, B[P, :])
        mem("DinvB: nnz(Y)=$(nnz(Yp))")
        ymax = maximum(abs, nonzeros(Yp); init=0.0)
        if ymax > 0
            I2, J2, V2 = findnz(Yp)
            keep2 = abs.(V2) .> 1e-14 * ymax
            Yp = sparse(I2[keep2], J2[keep2], V2[keep2], size(Yp, 1), size(Yp, 2))
        end
        Y = Yp[invP, :]
        Kt = J[restR, C] - A * Y
        dropzeros!(Kt)
        mem("Ktilde formed: nnz(Ktilde)=$(nnz(Kt))")
        Fk = lu(Kt)
        mem("lu(Ktilde) done: nnz(L)+nnz(U)=$(nnz(Fk.L)+nnz(Fk.U))")
        return plan, Y, Fk
    catch e
        log("  schur: factorize failed — " * sprint(showerror, e))
        return nothing
    end
end

"""
    schur_bordered(A, b, N, fam; verbose, dfloor, tol=1e-9) -> Vector(N+1) or nothing

Solve the bordered system `A*z = b` where `A = [J Bcol; Crow d]` is
(NC+1)×(N+1) with NC == N: two condensed solves sharing one factorization,
plus the scalar border correction. Self-verifying against A to `tol`
(relative); `nothing` on any failure — caller falls back to `lu(A)`.
"""
function schur_bordered(A::SparseMatrixCSC{Float64,Int}, b::Vector{Float64}, N::Int,
                        fam::Vector{String}; verbose::Bool=false,
                        dfloor::Float64=SCHUR_DFLOOR, tol::Float64=1e-6)
    log(msg) = verbose && println(msg)
    NC = size(A, 1) - 1
    (size(A, 2) == N + 1 && NC == N) || return nothing
    J = A[1:NC, 1:N]
    bcol = Vector(A[1:NC, N+1])
    # Use a SparseVector for findnz compatibility
    row_vec = SparseVector(A[NC+1, :])
    # findnz returns (row_indices, col_ptr, values) for matrices, 
    # but for SparseVector it returns (indices, values).
    res = findnz(row_vec)
    if length(res) == 2
        Crow, Cval = res
    else
        Crow, _, Cval = res
    end
    d = 0.0
    c = zeros(Float64, N)
    for (i, v) in zip(Crow, Cval)
        if i == N + 1
            d = v
        else
            c[i] = v
        end
    end
    f = b[1:NC]; g = b[NC+1]
    fz = schur_factorize(J, fam; verbose=verbose, dfloor=dfloor)
    fz === nothing && return nothing
    plan, Y, Fk = fz
    local z
    try
        v1 = _schur_correct(J, plan, Y, Fk, Vector{Float64}(f))
        v2 = _schur_correct(J, plan, Y, Fk, bcol)
        for (vv, rr) in ((v1, Vector{Float64}(f)), (v2, bcol))
            for _ in 1:SCHUR_NREFINE
                r = rr - J * vv
                rel = norm(r, Inf) / max(norm(rr, Inf), eps())
                rel < 1e-12 && break
                dx = _schur_correct(J, plan, Y, Fk, r)
                all(isfinite, dx) || break
                vn = vv + dx
                norm(rr - J * vn, Inf) < norm(r, Inf) || break
                vv .= vn
            end
        end
        denom = d - dot(c, v2)
        abs(denom) <= 1e-14 * max(abs(d), 1.0) && return nothing
        y = (g - dot(c, v1)) / denom
        x = v1 - y * v2
        z = Vector{Float64}(undef, N + 1)
        z[1:N] .= x; z[N+1] = y
        rel = norm(A * z - b, Inf) / max(norm(b, Inf), eps())
        if !(rel <= tol && all(isfinite, z))
            log("  schur_bordered: self-verify failed rel=" * string(rel))
            return nothing
        end
    catch e
        log("  schur_bordered failed — " * sprint(showerror, e))
        return nothing
    end
    return z
end
