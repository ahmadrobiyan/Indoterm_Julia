module IndotermJulia

using NamedArrays

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

export COM, IND, SRC, OCC, MAR, REG, DST, ORG, PRD, HOU
export read_national_data, read_regsupp_data, read_distgone_data
export build_reg0!, build_reg1!, build_reg2!, ras_balance!, build_pstras!, build_premod!
export Reg0Result, Reg1Result, Reg2Result, RasResult, PstrasResult, PremodResult
export aggregate_model!, AGGCOM, SEC_MAP_185_to_25
export AggResult

end # module IndotermJulia
