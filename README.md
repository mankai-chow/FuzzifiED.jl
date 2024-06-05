# FuzzifiED.jl Version 0.5.10

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculation on the fuzzy sphere, and also facilitates the DMRG calculations by ITensors. It can also be used for generic fermion models. 

Documentations can be found at [http://docs.fuzzified.world/](http://docs.fuzzified.world/). (For Perimeter researchers, for now please use a VPN to access the documentation.)

To Apple Silicon (Mac M1, M2, etc.) users: you may want to use [Julia](https://julialang.org/downloads/) for x86 (Intel or Rosetta) instead of macOS (Apple Silicon). 

To install the package, please enter the following commands in Julia :
```julia
julia> using Pkg; 
julia> Pkg.add(url="https://github.com/mankai-chow/FuzzifiED_jll.jl.git")
julia> Pkg.add(url="https://github.com/mankai-chow/FuzzifiED.jl.git")
```
Include at the start of your Julia script
```julia
using FuzzifiED
```

If this package is helpful in your research, we would appreciate it if you mention in the acknowledgement. If you have any questions, please contact Zheng Zhou (周正) at [zzhou@pitp.ca](mailto:zzhou@pitp.ca).
