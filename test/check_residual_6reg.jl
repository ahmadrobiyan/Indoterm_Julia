push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP

println("--- pipeline @ 25×6 ---")
nat = read_national_data(); regsup = read_regsupp_data(); distgone = read_distgone_data()
reg0_r = build_reg0!(nat, regsup)
reg1_r = build_reg1!(reg0_r.reg0, regsup, distgone)
reg2_r = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
ras_r  = ras_balance!(reg2_r.reg2)
pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone)
prem_r = build_premod!(pstra_r.pstras, regsup, reg0_r.elast)
agg6 = aggregate_regions!(aggregate_model!(prem_r.premod).agg, REG_MAP_34_to_6, 6)
params = prepare_parameters!(agg6)

m, vars = build_model_full!(agg6, params)
println("vars=$(num_variables(m)), cons=$(num_constraints(m; count_variable_in_set_constraints=false))")
initialize_model!(m, vars; bmk_levels=benchmark_levels(params))

# Value Ipopt actually sees at the start point: fixed value if fixed, else start.
point(v::VariableRef) = JuMP.is_fixed(v) ? JuMP.fix_value(v) : JuMP.start_value(v)
canon(s::AbstractString) = replace(s, r"\[[0-9, ]+\]" => "[.]")

# For an equality constraint  f(x) == 0  (or f(x) == rhs), the residual is
# value(point, con) - rhs. JuMP normalizes to func-in-set; use the set constant.
const MOI = JuMP.MOI

println("--- residuals at benchmark start point ---")
n = 0
worst = 0.0; worst_str = ""
bysig = Dict{String,Tuple{Int,Float64,String}}()  # sig => (count_over_tol, max_abs, example)
TOL = 1e-6
for (F, S) in list_of_constraint_types(m)
    F <: VariableRef && continue
    for con in all_constraints(m, F, S)
        global n += 1
        co = JuMP.constraint_object(con)
        rhs = 0.0
        s = co.set
        if s isa MOI.EqualTo; rhs = s.value
        elseif s isa MOI.LessThan; rhs = s.upper
        elseif s isa MOI.GreaterThan; rhs = s.lower
        end
        local lhs
        try; lhs = JuMP.value(point, con); catch; lhs = NaN; end
        r = abs(lhs - rhs)
        if !isfinite(r) || r > TOL
            sig = canon(string(co.func))
            cnt, mx, ex = get(bysig, sig, (0, 0.0, ""))
            newmx = (isfinite(r) && r > mx) ? r : mx
            newex = (cnt == 0) ? string(co.func) : ex
            bysig[sig] = (cnt + 1, newmx, newex)
        end
        if isfinite(r) && r > worst
            global worst, worst_str = r, string(co.func)
        end
    end
end

nbad = sum(t -> t[1], values(bysig); init=0)
println("checked $n constraints; $nbad violate |resid|>$TOL at the benchmark point")
println("worst finite residual: $worst")
println("distinct violating equation signatures: $(length(bysig))")
println()
for (sig, (cnt, mx, ex)) in sort(collect(bysig); by = kv -> kv[2][2], rev=true)
    println("── count=$cnt   max|resid|=$(round(mx, sigdigits=4))")
    println("   sig: $sig")
    println("   ex : $(first(ex, 240))")
end
