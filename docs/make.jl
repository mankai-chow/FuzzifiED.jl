push!(LOAD_PATH,"../src/")

include("../src/FuzzifiED.jl")
using Documenter
using LinearAlgebra, ITensors
using .FuzzifiED
using .FuzzifiED.ITensorSupport

makedocs(sitename = "FuzzifiED.jl") 