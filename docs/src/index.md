# FuzzifiED.jl

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculation on the fuzzy sphere. Install the package with the commands
```julia
julia> using Pkg; Pkg.add(url="https://github.com/mankai-chow/FuzzifiED.jl.git")
```
Include at the start of your Julia script
```julia
using LinearAlgebra
using FuzzifiED
```

#### ITensor support

This package also supports importing the `Site` and `OpSum` objects from `ITensors` library, to use that, include also 
```julia
using ITensors 
using FuzzifiED.ITensorsSupport
```
