"""
Why does `numeraire = :exrate` come out one free variable short of square?

`verify_coalprice_reference.jl` died with

    Newton needs a square system but got 84135 equations vs 84136 free variables

and `apply_swaps!` cannot be the cause: it validates equal lengths and
all-or-nothing exogeneity on each side, then frees exactly one variable per
variable it pins. Swaps are count-neutral by construction. The equation count is
also closure-independent — no constraint in `build_model_full!` looks at what is
fixed. So the culprit has to be the numeraire branch at the end of
`initialize_model!`, and the cheapest way to name it is to build the model under
each numeraire and diff the free-variable sets by family.

Two builds at ~29s each. Prints the count, the difference, and the name of every
variable whose fixed/free status differs between the two.

Run: julia --project=IndotermJulia IndotermJulia/test/diag_numeraire_square.jl
"""

include(joinpath(@__DIR__, "pipeline_cache.jl"))
using Printf, JuMP

agg6, params = cached_pipeline(6)
bmk = benchmark_levels(params)

"""
Free-variable names (family + linear index) under one numeraire choice, PLUS the
orphan/dead-row census `solve_newton!` performs before it counts. That census is
the part a naive free-vs-equations tally misses: `solve_newton!` pins every free
variable that appears in no constraint and deletes every all-fixed row, so two
closures with identical raw counts can still differ in what it finally solves.
"""
function freeset(numeraire::Symbol)
    m, vars = build_model_full!(agg6, params)
    initialize_model!(m, vars; bmk_levels=bmk, numeraire=numeraire)
    free = Set{String}()
    fam = Dict{String,Int}()
    # name lookup by JuMP index, so orphans can be reported by family not by ref
    byidx = Dict{Int64,String}()
    for (nm, v) in vars, (k, vr) in enumerate(v isa AbstractArray ? vec(collect(v)) : [v])
        if !JuMP.is_fixed(vr)
            push!(free, "$nm[$k]")
            fam[nm] = get(fam, nm, 0) + 1
            byidx[vr.index.value] = "$nm[$k]"
        end
    end
    ncon = sum(length(all_constraints(m, F, S))
               for (F, S) in list_of_constraint_types(m) if !(F <: VariableRef))

    # replicate solve_newton!'s scan (it dedups refs by JuMP index; `vars` can
    # alias the same ref under two names, which is why `free` may over-count)
    freeid = Dict{Int64,Int}(); freeref = JuMP.VariableRef[]
    for (_, v) in vars, vr in (v isa AbstractArray ? v : (v,))
        JuMP.is_fixed(vr) && continue
        iv = vr.index.value
        haskey(freeid, iv) && continue
        push!(freeref, vr); freeid[iv] = length(freeref)
    end
    seen = falses(length(freeref)); ndead = 0
    for (F, S) in list_of_constraint_types(m)
        F <: JuMP.VariableRef && continue
        for con in all_constraints(m, F, S)
            buf = IndotermJulia._collect_free!(Int[], JuMP.constraint_object(con).func, freeid)
            isempty(buf) ? (ndead += 1) : (for id in buf; seen[id] = true; end)
        end
    end
    orphans = [get(byidx, freeref[i].index.value, "?") for i in eachindex(freeref) if !seen[i]]
    @printf("\n:%s — raw free %d (deduped %d), equations %d, orphans %d, dead rows %d\n",
            numeraire, length(free), length(freeref), ncon, length(orphans), ndead)
    @printf("   → solve_newton! will see %d equations vs %d unknowns\n",
            ncon - ndead, length(freeref) - length(orphans))
    isempty(orphans) || println("   orphans: ", join(sort(orphans), ", "))
    return free, fam, ncon
end

fg, famg, ncg = freeset(:gdppi)
fe, fame, nce = freeset(:exrate)

@printf("\n%-10s %12s %12s\n", "numeraire", "free vars", "equations")
@printf("%-10s %12d %12d  (diff %+d)\n", ":gdppi", length(fg), ncg, length(fg) - ncg)
@printf("%-10s %12d %12d  (diff %+d)\n", ":exrate", length(fe), nce, length(fe) - nce)

println("\nfree under :exrate but NOT under :gdppi:")
for s in sort(collect(setdiff(fe, fg))); println("  ", s); end
println("free under :gdppi but NOT under :exrate:")
for s in sort(collect(setdiff(fg, fe))); println("  ", s); end

println("\nper-family counts that differ:")
for nm in sort(collect(union(keys(famg), keys(fame))))
    a, b = get(famg, nm, 0), get(fame, nm, 0)
    a == b || @printf("  %-14s :gdppi %6d   :exrate %6d\n", nm, a, b)
end
println("\nDone.")
