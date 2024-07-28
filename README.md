# FuzzifiED.jl Version 0.8.2

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculation on the fuzzy sphere, and also facilitates the DMRG calculations by ITensors. It can also be used for generic fermion models. Using this package, you can reproduce almost all the ED results in fuzzy sphere works ; for detail, see [« The fuzzified world »](https://www.fuzzified.world/fuzzified-world).

Documentations can be found at [http://docs.fuzzified.world/](http://docs.fuzzified.world/).

Download Julia at [this link](https://julialang.org/downloads/). 

To Apple Silicon (Mac M1, M2, etc.) users : The support for AArch64 architecture is not guarenteed. You may want to use Julia for x86 (Intel or Rosetta) instead of macOS (Apple Silicon).

To install the package, please enter the following commands in Julia :
```julia
using Pkg ;
Pkg.add(url="https://github.com/mankai-chow/FuzzifiED_jll.jl.git") ;
Pkg.add(url="https://github.com/mankai-chow/FuzzifiED.jl.git") ;
```
Include at the start of your Julia script
```julia
using FuzzifiED
```

If this package is helpful in your research, we would appreciate it if you mention in the acknowledgement. If you have any questions, please contact Zheng Zhou (周正) at [zzhou@pitp.ca](mailto:zzhou@pitp.ca).
