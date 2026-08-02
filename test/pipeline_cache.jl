"""
Pipeline cache for 25×N diagnostics.

The Steps 0-4 data pipeline costs ~250 s and produces a *deterministic* result — it is pure
array transformation, no solve. Paying it once per session instead of once per experiment is
the difference between a 7-minute and a sub-minute debug loop (see PLAN.md "Session findings
(2026-07-29)" §A6: 429 s elapse before `solve_newton!` is even entered).

Usage from any diagnostic script:

    include(joinpath(@__DIR__, "pipeline_cache.jl"))
    agg6, params = cached_pipeline(6)

Pass `refresh=true` to force a rebuild after changing anything in `build_reg*`/`ras_balance!`/
`build_pstras!`/`build_premod!`/`aggregate_*`/`prepare_parameters!`. The cache does NOT cover
`build_model_full!` (a JuMP model holds pointers into a live MOI backend and does not
round-trip through `Serialization`); that remains ~178 s per run.
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using Serialization

const CACHE_DIR = joinpath(@__DIR__, "..", "data")

cache_path(nr::Int) = joinpath(CACHE_DIR, "cache_$(nr)reg.jls")

"""
    cached_pipeline(nr::Int=6; refresh::Bool=false, rmap=nothing) -> (agg, params)

Run Steps 0-4 down to the `nr`-region aggregated database and its derived parameters, caching
the result to `data/cache_<nr>reg.jls`.

`nr == 6` (`REG_MAP_34_to_6`) and `nr == 34` (no region aggregation) resolve their map
automatically. Any other `nr` REQUIRES an explicit `rmap` — before 2026-08-01 this function
passed `REG_MAP_34_to_6` regardless of `nr`, so `cached_pipeline(12)` allocated 12 region
slots, wrote into only the first 6, and serialized a database whose groups 7-12 were
identically zero as if it were valid. `aggregate_regions!` now rejects that partition
outright; this signature makes the map and the region count travel together so they cannot
disagree in the first place.
"""
function cached_pipeline(nr::Int=6; refresh::Bool=false,
                         rmap::Union{Nothing,Vector{Int}}=nothing)
    if rmap === nothing && !(nr in (6, 34))
        error("cached_pipeline($nr): no default region map for nr=$nr — pass `rmap=` " *
              "explicitly (see REG_MAP_34_to_12 in test/verify_aggregation_consistency.jl)")
    end
    path = cache_path(nr)
    if !refresh && isfile(path)
        t = @elapsed obj = open(deserialize, path)
        println("pipeline cache HIT  ($(round(t, digits=1))s)  $path")
        return obj.agg, obj.params
    end

    println("pipeline cache MISS — running Steps 0-4 at 25×$nr (expect ~250s)")
    t0 = time()
    _t(msg) = println("  T+$(round(time() - t0, digits=1))s  $msg")

    _t("read data")
    nat = read_national_data(); regsup = read_regsupp_data(); distgone = read_distgone_data()
    _t("build_reg0!");   reg0_r  = build_reg0!(nat, regsup)
    _t("build_reg1!");   reg1_r  = build_reg1!(reg0_r.reg0, regsup, distgone)
    _t("build_reg2!");   reg2_r  = build_reg2!(reg1_r.reg1, regsup; reg0_diag=reg0_r.diag)
    _t("ras_balance!");  ras_r   = ras_balance!(reg2_r.reg2)
    _t("build_pstras!"); pstra_r = build_pstras!(ras_r.ras, reg1_r.reg1, distgone)
    _t("build_premod!"); prem_r  = build_premod!(pstra_r.pstras, regsup, reg0_r.elast)
    _t("aggregate")
    agg34 = aggregate_model!(prem_r.premod).agg
    _map = rmap === nothing ? REG_MAP_34_to_6 : rmap
    agg = nr == 34 ? agg34 : aggregate_regions!(agg34, _map, nr)
    _t("prepare_parameters!")
    params = prepare_parameters!(agg)

    mkpath(CACHE_DIR)
    open(f -> serialize(f, (; agg, params)), path, "w")
    _t("cached → $path")
    return agg, params
end
