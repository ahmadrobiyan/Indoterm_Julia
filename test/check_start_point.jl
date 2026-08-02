push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP

println("--- Running data pipeline (Steps 0-4) ---")
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
params = prepare_parameters!(agg_r.agg)

println("--- build_model_full! ---")
t0 = time()
m, vars = build_model_full!(agg_r.agg, params)
println("Build: $(round(time()-t0, digits=1))s, vars=$(num_variables(m)), cons=$(num_constraints(m; count_variable_in_set_constraints=false))")

println("--- initialize_model! (benchmark, all shocks = 0.0) ---")
initialize_model!(m, vars)

# What Ipopt actually sees at the initial point: the fixed value for closure
# variables (fix() wins over whatever start_value was set before fixing), the
# start value for everything else.
point(v::VariableRef) = JuMP.is_fixed(v) ? JuMP.fix_value(v) : JuMP.start_value(v)

# Collapse "xlab[45,3]" -> "xlab[.,.]" so the ~124k bad constraints (all built
# via anonymous [c=1:na, d=1:nr]-style containers, so JuMP.name(con) is always
# "") group down into a handful of distinct equation signatures we can map back
# to a build_equations.jl line, instead of drowning in per-index duplicates.
canonicalize(s::AbstractString) = replace(s, r"\[[0-9, ]+\]" => "[.]")

println("--- Evaluating every nonlinear/affine constraint at the initial point ---")
n_checked = 0
patterns = Dict{String,Int}()
examples = Dict{String,String}()
for (F, S) in list_of_constraint_types(m)
    F <: VariableRef && continue  # skip plain variable bound/fix constraints
    for con in all_constraints(m, F, S)
        global n_checked += 1
        local val
        try
            val = JuMP.value(point, con)
        catch e
            val = e
        end
        fval = val isa Number ? val : NaN
        if !(val isa Number) || !isfinite(fval)
            raw = string(JuMP.constraint_object(con).func)
            key = canonicalize(raw)
            patterns[key] = get(patterns, key, 0) + 1
            if !haskey(examples, key)
                examples[key] = raw
            end
        end
    end
end
n_bad = sum(values(patterns); init=0)
println("Checked $n_checked constraints; $n_bad non-finite (or erroring) at the initial point")
println("Distinct equation signatures among the bad ones: $(length(patterns))")
for (key, cnt) in sort(collect(patterns); by=last, rev=true)
    println("--- count=$cnt ---")
    println("  pattern:  $key")
    println("  example:  $(examples[key])")
end
