push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP

# ── Hopcroft-Karp maximum bipartite matching ─────────────────────────────
# left side = constraints (L), right side = variables (R)
# adj[l] = list of r that l connects to
function hopcroft_karp(L, R, adj)
    INF = typemax(Int)
    pairL = zeros(Int, L)   # which r is matched to l (0 = unmatched)
    pairR = zeros(Int, R)   # which l is matched to r (0 = unmatched)
    dist = zeros(Int, L)

    function bfs()
        q = Int[]
        for l in 1:L
            if pairL[l] == 0
                dist[l] = 0; push!(q, l)
            else
                dist[l] = INF
            end
        end
        found = false
        for cur in q
            for r in adj[cur]
                nxt = pairR[r]
                if nxt != 0 && dist[nxt] == INF
                    dist[nxt] = dist[cur] + 1; push!(q, nxt)
                elseif nxt == 0
                    found = true
                end
            end
        end
        return found
    end

    function dfs(l)
        for r in adj[l]
            nxt = pairR[r]
            if nxt == 0 || (dist[nxt] == dist[l] + 1 && dfs(nxt))
                pairL[l] = r; pairR[r] = l; return true
            end
        end
        dist[l] = INF; return false
    end

    matching = 0
    while bfs()
        for l in 1:L
            if pairL[l] == 0 && dfs(l)
                matching += 1
            end
        end
    end
    return matching, pairL, pairR
end

# ── Build model ──────────────────────────────────────────────────────────
println("="^64)
println("RANK-DEFICIENCY DIAGNOSTIC (bipartite matching @ 25x6)")
println("="^64)

println("--- Running data pipeline ---")
nat = read_national_data()
regsup = read_regsupp_data()
distgone = read_distgone_data()
reg0_r = build_reg0!(nat, regsup)
reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
ras_r  = ras_balance!(reg2_r.reg2)
pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone)
prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast)
agg_r  = aggregate_model!(prem_r.premod)

println("--- aggregate_regions! (34->6) ---")
agg6 = aggregate_regions!(agg_r.agg, REG_MAP_34_to_6, 6)
params = prepare_parameters!(agg6)

println("--- build_model_full! @ 25x6 (no shocks) ---")
t0 = time()
m, vars = build_model_full!(agg6, params)
println("Build: $(round(time()-t0, digits=1))s")

# ── Initialise at benchmark (no shocks) ──────────────────────────────────
bmk = benchmark_levels(params)
initialize_model!(m, vars; bmk_levels=bmk)

# ── Build incidence from the Jacobian structure ──────────────────────────
# This mirrors solve_newton!'s setup but without solving.
freeid = Dict{Int64,Int}(); freeref = VariableRef[]
for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
    JuMP.is_fixed(vr) && continue
    iv = vr.index.value
    haskey(freeid, iv) && continue
    push!(freeref, vr); freeid[iv] = length(freeref)
end
N = length(freeref)

# Build list_of_constraint_types
list_of_constraint_types = JuMP.list_of_constraint_types(m)
allidx = VariableRef[]
for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
    push!(allidx, vr)
end
pos = Dict{Int,Int}(JuMP.index(v).value => k for (k, v) in enumerate(allidx))

# Collect constraint rows and their referenced free vars
allv = JuMP.all_variables(m)
allvidx = [JuMP.index(v).value for v in allv]
buf = Int[]
cons_freevars = Vector{Int}[]
con_names = String[]
for (Ftype, S) in list_of_constraint_types
    Ftype <: VariableRef && continue
    for con in all_constraints(m, Ftype, S)
        co = JuMP.constraint_object(con)
        empty!(buf)
                IndotermJulia._collect_free!(buf, co.func, freeid)
        if isempty(buf)
            push!(con_names, "*DEAD* $(nameof(Ftype))")
            push!(cons_freevars, Int[])
        else
            push!(con_names, string(nameof(Ftype)))
            push!(cons_freevars, copy(buf))
        end
    end
end
NC = length(cons_freevars)
println("  $(NC) constraints, $(N) free variables")

# Build adjacency list for constraints -> variables
adj = [Int[] for _ in 1:NC]
for (r, fvs) in enumerate(cons_freevars)
    append!(adj[r], fvs)
end

# ── Running Hopcroft-Karp matching ───────────────────────────────────────
println("--- Hopcroft-Karp matching ---")
t0 = time()
match_size, pairL, pairR = hopcroft_karp(NC, N, adj)
println("  Max matching: $(match_size) / $(min(NC,N))")
println("  Matching time: $(round(time()-t0, digits=2))s")

unmatched_cons = findall(iszero, pairL)
unmatched_vars = findall(iszero, pairR)

println("  Unmatched constraints: $(length(unmatched_cons))")
println("  Unmatched free vars:   $(length(unmatched_vars))")

# ── Classify unmatched constraints ───────────────────────────────────────
if !isempty(unmatched_cons)
    println("\n--- Unmatched constraint analysis ---")
    # Count by constraint type
    type_counts = Dict{String,Int}()
    for ci in unmatched_cons
        typ = con_names[ci]
        type_counts[typ] = get(type_counts, typ, 0) + 1
    end
    for (typ, cnt) in sort(collect(type_counts), by=x->-x[2])
        println("  $cnt × $typ")
    end

    # Show first 10 unmatched constraints in detail
    println("\n  First 10 unmatched constraints:")
    for ci in unmatched_cons[1:min(10,length(unmatched_cons))]
        fvs = cons_freevars[ci]
        vr_names = [JuMP.name(freeref[fv]) for fv in fvs]
        println("    #$ci $(con_names[ci]): $(length(fvs)) free vars — $(join(vr_names[1:min(3,end)], ", "))$(length(fvs) > 3 ? "…" : "")")
    end
end

# ── Examine unmatched variables ──────────────────────────────────────────
if !isempty(unmatched_vars)
    println("\n--- Unmatched free variable analysis ---")
    # Group by variable name prefix (before first '[')
    name_counts = Dict{String,Int}()
    for vi in unmatched_vars
        nm = JuMP.name(freeref[vi])
        base = split(nm, '[')[1]
        name_counts[base] = get(name_counts, base, 0) + 1
    end
    for (base, cnt) in sort(collect(name_counts), by=x->-x[2])
        println("  $cnt × $base")
    end
end

# ── Column-norm analysis ─────────────────────────────────────────────────
# The system is structurally full-rank; the 637 rank deficiency must be
# numerical (tiny column norms from small share parameters).
println("\n--- Column norm analysis ---")
# Build the Jacobian at benchmark
freecols = [pos[JuMP.index(v).value] for v in freeref]
allv = JuMP.all_variables(m)
allvidx = [JuMP.index(v).value for v in allv]
Jval = zeros(length(adj))  # placeholder — we need the actual Jacobian
# Instead, compute column norms of the Jacobian sparsity pattern:
# For each column, count how many constraints reference it
col_counts = zeros(Int, N)
for fvs in cons_freevars
    for fv in fvs; col_counts[fv] += 1; end
end
println("  Column reference counts (col count → freq):")
bins = Dict{Int,Int}()
for c in col_counts; bins[c] = get(bins, c, 0) + 1; end
for (cnt, freq) in sort(collect(bins))
    println("    $cnt constraints ref this var → $freq variables")
end

# Variables referenced by very few constraints
rare = findall(<=(3), col_counts)
println("\n  Variables referenced by ≤3 constraints: $(length(rare))")
for vi in rare[1:min(10, length(rare))]
    nm = JuMP.name(freeref[vi])
    println("    #$vi: $nm  (ref'd by $(col_counts[vi]) constraints)")
end

println("\n" * "=" ^ 64)
