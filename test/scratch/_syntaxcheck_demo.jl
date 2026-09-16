src = read(joinpath(@__DIR__, "demo_scenario_development.jl"), String)
Meta.parseall(src)
println("syntax OK — ", length(split(src, '\n')), " lines")
