# julia --color=yes --project make.jl && rm -r publish/*/* ; mv -f build/* publish
push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"../ext/")

using Documenter
using WignerSymbols
using SparseArrays
using ITensors
using ITensorMPS
using FuzzifiED

makedocs(sitename = "FuzzifiED.jl", 
    pages = ["Home" => "index.md", 
        "Introduction" => "intro.md",
        "Example" => "example.md",
        "Core functions" => "core.md",
        "ITensors support" => "itensors.md",
        "Built-in models" => "models.md",
        "Archieved interfaces" => "archieve.md",
        "Releases" => "releases.md"],
    format = Documenter.HTML(assets = ["assets/themes/serif.css"], repolink = "https://github.com/mankai-chow/FuzzifiED.jl")
) 