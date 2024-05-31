# julia --color=yes --project make.jl
push!(LOAD_PATH,"../src/")

include("../src/FuzzifiED.jl")
using Documenter
using ITensors
using WignerSymbols
using SparseArrays
using .FuzzifiED

makedocs(sitename = "FuzzifiED.jl", 
    pages = ["Home" => "index.md", 
        "Example" => "example.md",
        "Core functions" => "core.md",
        "ITensors support" => "itensors.md",
        "Models" => "models.md",
        "Releases" => "releases.md"]) 