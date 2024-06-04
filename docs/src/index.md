# FuzzifiED.jl v

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculation on the fuzzy sphere, and also facilitates the DMRG calculations by ITensors. It can also be used for generic fermion models. 

This package is currently only supported on 64-bit environments. 

This package is developped by Zheng Zhou (周正) at Perimeter Institute and collaborators. If this package is helpful in your research, we would appreciate it if you mention in the acknowledgement. If you have any questions, please contact at [zzhou@pitp.ca](mailto:zzhou@pitp.ca).

## Install

To install the package, please first enter Julia by entering in the command line `julia`, and then enter the commands
```julia
julia> using Pkg; 
julia> Pkg.add(url="https://github.com/mankai-chow/FuzzifiED_jll.jl.git")
julia> Pkg.add(url="https://github.com/mankai-chow/FuzzifiED.jl.git")
```
Include at the start of your Julia script
```julia
using FuzzifiED
```


## Outline 

```@contents
Pages = [
    "example.md",
    "core.md",
    "itensors.md",
    "models.md"
]
Depth = 2
```

## Index 

```@index
Pages = [
    "core.md",
    "itensors.md",
    "models.md"
]
```