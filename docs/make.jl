# julia --color=yes --project make.jl && rm -r publish/*/* ; mv -f build/* publish
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
        "Built-in models" => "models.md",
        "Archieved interfaces" => "archieve.md",
        "Releases" => "releases.md"],
    format = Documenter.HTML(assets = ["assets/themes/serif.css"])) 