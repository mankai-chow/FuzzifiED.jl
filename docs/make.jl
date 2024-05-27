# julia --color=yes --project make.jl
push!(LOAD_PATH,"../src/")

include("../src/FuzzifiED.jl")
using Documenter
using LinearAlgebra
using ITensors
using WignerSymbols
using SparseArrays
using .FuzzifiED
using .FuzzifiED.ITensorsSupport
using .FuzzifiED.Models

makedocs(sitename = "FuzzifiED.jl", 
    pages = ["Home" => "index.md", 
        "Example" => "example.md",
        "Reference" => "reference.md",
        "Releases" => "releases.md"]) 