# FuzzifiED.jl

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculation on the fuzzy sphere. 

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

### ITensor support

This package also supports importing the `Site` and `OpSum` objects from `ITensors` library, to use that, include also 
```julia
using ITensors 
using FuzzifiED.ITensorsSupport
```

### Built-in models 

To use the built-in models including the Ising model and ``\mathrm{Sp}(N)`` model, include also 
```julia
using FuzzifiED.Models
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