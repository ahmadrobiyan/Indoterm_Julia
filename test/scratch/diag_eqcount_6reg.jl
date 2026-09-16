# Diagnostic: localise the ~6,900-equation drop (90,460 → 83,541) at 25×6.
#
# The zero-flow guards added in gate 5c-ii have the shape
#     if FLOW > 1e-10;  @constraint(...)  else  fix(var, ...; force=true)  end
# so every extra cell falling into the `else` branch removes exactly ONE
# equation and ONE unknown. A symmetric drop is therefore the fingerprint of
# guards firing on more cells than expected.
#
# This script counts variables already FIXED immediately after
# `build_model_full!` — i.e. before `initialize_model!` applies the closure —
# which isolates guard-pins from closure-exogenous fixes. It also caches the
# pipeline output so repeat runs skip the ~2-3 min data stage.

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IndotermJulia
using JuMP, Serialization

const CACHE = joinpath(@__DIR__, "..", "..", "agg6_cache.jls")

function get_inputs()
    if isfile(CACHE)
        println("--- loading cached pipeline output ---")
        return deserialize(CACHE)
    end
    println("--- running data pipeline (Steps 0-4) ---")
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
    agg6 = aggregate_regions!(agg_r.agg, REG_MAP_34_to_6, 6)
    params = prepare_parameters!(agg6)
    serialize(CACHE, (agg6, params))
    return (agg6, params)
end

function main()
    agg6, params = get_inputs()

    println("--- build_model_full! @ 25×6 ---")
    t0 = time()
    m, vars = build_model_full!(agg6, params)
    nv = num_variables(m)
    nc = num_constraints(m; count_variable_in_set_constraints=false)
    println("Build: $(round(time()-t0, digits=1))s, vars=$nv, cons=$nc")

    # ── variables already fixed at BUILD time = guard pins ──────────────────
    pinned = Dict{String,Int}()
    total_cells = Dict{String,Int}()
    npin = 0
    for (nm, v) in vars
        for vr in (v isa AbstractArray ? vec(v) : [v])
            total_cells[nm] = get(total_cells, nm, 0) + 1
            if JuMP.is_fixed(vr)
                pinned[nm] = get(pinned, nm, 0) + 1
                npin += 1
            end
        end
    end

    println()
    println("Variables FIXED at build time (guard pins, pre-closure): $npin")
    println(rpad("family", 20), rpad("pinned", 10), rpad("of total", 10), "pct")
    for (nm, cnt) in sort(collect(pinned), by = x -> -x[2])
        tot = total_cells[nm]
        println(rpad(nm, 20), rpad(cnt, 10), rpad(tot, 10),
                round(100 * cnt / tot, digits=1), "%")
    end

    println()
    println("free (unfixed) vars at build time: ", nv - npin)
    println("constraints:                       ", nc)
    println("gap (free - cons):                 ", nv - npin - nc)
    println()
    println("Reference counts from PLAN.md at 25×6:")
    println("  gate 5c-i  (core)          88,599")
    println("  gate 6     (+dynamics)     89,911")
    println("  gate 6c    (+macros)       90,460")
    println("  today                      ", nc)
    println("  delta vs 90,460            ", nc - 90_460)
end

main()
