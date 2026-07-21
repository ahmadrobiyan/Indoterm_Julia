push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "IndotermJulia.jl"))
println("Module loaded OK")
println("COM size = $(length(COM))")
println("REG size = $(length(REG))")
