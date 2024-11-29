# julia --color=yes --project make.jl && rm -r publish/*/* ; mv -f build/* publish
push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"../ext/")

using Documenter
using WignerSymbols
using SparseArrays
using ITensors
using ITensorMPS
using CUDA
using LinearAlgebra
using HDF5
using FuzzifiED
using FuzzifiED.Fuzzifino

makedocs(sitename = "FuzzifiED.jl", 
    pages = ["Home" => "index.md", 
        "Introduction" => "intro.md",
        "Examples" => "example.md",
        "Core functions" => "core.md",
        "Built-in models" => "models.md",
        "ITensor extension" => "itensors.md",
        "Other extensions" => "extension.md", 
        "Fuzzifino" => "fuzzifino.md",
        "Archieved interfaces" => "archieve.md",
        "Releases" => "releases.md"],
    format = Documenter.HTML(assets = ["assets/themes/serif.css"], repolink = "https://github.com/mankai-chow/FuzzifiED.jl")
) 