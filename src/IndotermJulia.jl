module IndotermJulia

using NamedArrays
using JuMP, Ipopt

include("prepare_sets.jl")
include("read_data.jl")
include("build_reg0!.jl")
include("build_reg1!.jl")
include("build_reg2!.jl")
include("ras_balance!.jl")
include("build_pstras!.jl")
include("build_premod!.jl")
include("aggregation_data.jl")
include("aggregate_model!.jl")
include("aggregate_regions!.jl")
include("ces_helper.jl")
include("prepare_parameters.jl")
include("build_equations.jl")
include("build_model!.jl")
include("build_dynamics!.jl")
include("build_macros!.jl")
include("benchmark_levels.jl")
include("initialize_model!.jl")
include("solve_model!.jl")
include("solve_newton!.jl")
include("schur_linsolve.jl")
include("closures.jl")
include("continuation.jl")
include("arclength.jl")
include("calculate_gdp.jl")
include("run_model!.jl")
# after run_model!.jl (needs `pct`) and closures.jl (needs `_side_name`) — the
# scenario constants call both at definition time.
include("scenarios.jl")
include("regional_confidence.jl")

export COM, IND, SRC, OCC, MAR, REG, DST, ORG, PRD, HOU
export read_national_data, read_regsupp_data, read_distgone_data
export build_reg0!, build_reg1!, build_reg2!, ras_balance!, build_pstras!, build_premod!
export Reg0Result, Reg1Result, Reg2Result, RasResult, PstrasResult, PremodResult
export aggregate_model!, AGGCOM, SEC_MAP_185_to_25
export aggregate_regions!, REG6, REG_MAP_34_to_6
export AggResult
export ces, ces_calibrate
export prepare_parameters!
export build_model!, build_model_full!
export build_dynamics!, update_dynamics!
export build_macros!, MAINMACROS, SELMACROS
export benchmark_levels
export initialize_model!, BASE_CLOSURE_SCALARS, BASE_CLOSURE_ARRAYS, DYNAMIC_ONLY_CLOSURE_NAMES, ZERO_START_NAMES
export solve_model!
export solve_newton!, NewtonResult, diagnose_jacobian, euler_continuation!
export continuation_solve!, ContinuationResult
export arclength_solve!, ArclengthResult
export apply_swaps!, TERM_CMF_SWAPS, TERM_SR_SWAPS, COALPRICE_SWAPS, TERM_LR_SWAPS_MATCHED
export calculate_gdp, gdp_both_sides_check, gdp_change_consistency, print_gdp_report
export GDPReport, GDPINCCAT, GDPEXPCAT
export run_model!, ScenarioResult, pct, logpct
export Scenario, describe, TERM_CMF_REFERENCE, TERM_CMF_NO_DELUNITY, COALPRICE_REFERENCE
export HILIRISASI_BAN, HILIRISASI_SMELTER_JAVA, HILIRISASI_SMELTER_EAST, HILIRISASI_FULL
export regional_confidence_report, print_regional_confidence_report, REGIONAL_FRAGILITY_THRESHOLD

end # module IndotermJulia
