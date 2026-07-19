#=
julia -O3 --color=yes make.jl
=#
push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"../ext/")

using Documenter
using FuzzifiED
using FuzzifiED.Fuzzifino
using FuzzifiEDFullRotation

makedocs(sitename = "FuzzifiEDFullRotation.jl", 
    pages = ["Home" => "index.md", 
        "Principle" => "principle.md",
        "Tutorial" => "tutorial.md",
        "Interface" => "interface.md",
        "Releases" => "releases.md"],
    format = Documenter.HTML(
        assets = ["assets/serif.css", "assets/favicon.ico"], 
        repolink = "https://github.com/FuzzifiED/FuzzifiEDFullRotation.jl",
        footer = "Powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) and the [Julia Programming Language](https://julialang.org/). Copyright (c) 2026 Zheng Zhou (周正) and contributors."
    )
)
