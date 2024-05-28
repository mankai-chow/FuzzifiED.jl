# FuzzifiED.jl

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculation on the fuzzy sphere, and also facilitates the DMRG calculations by ITensors. It can also be used for generic fermion models. 

This package is developped by Zheng Zhou (周正) at Perimeter Institute. If you use this package, please mention in the acknowledgement. If you have any questions, please contact at [zzhou@pitp.ca](mailto:zzhou@pitp.ca).

## Installation

Install the package with the commands
```julia
julia> using Pkg; Pkg.add(url="https://github.com/mankai-chow/FuzzifiED.jl.git")
```
Include at the start of your Julia script
```julia
using FuzzifiED
```

## Outline 

```@contents
Pages = [
    "example.md",
    "reference.md"
]
Depth = 2
```

## Index 

```@index
Pages = ["reference.md"]
```