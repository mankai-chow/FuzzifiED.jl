# FuzzifiED.jl Version 0.10.4

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculations on the fuzzy sphere, and also facilitates the DMRG calculations by ITensors. It can also be used for generic fermionic and bosonic models. Using this package, you can reproduce almost all the ED results in fuzzy sphere works.

Documentations can be found at [http://docs.fuzzified.world/](http://docs.fuzzified.world/).

Download Julia at [this link](https://julialang.org/downloads/). 

To install the package, please enter the following commands in Julia :
```julia
using Pkg
Pkg.add(url="https://github.com/mankai-chow/FuzzifiED_jll.jl.git")
Pkg.add(url="https://github.com/mankai-chow/FuzzifiED.jl.git")
```
Include at the start of your Julia script
```julia
using FuzzifiED
```

Apple Silicon (Mac M1, M2, _etc._) users may experience trouble using `GetEigensystem` function. Another function `GetEigensystemKrylov` is provided instead. See documentation for detail.

If this package is helpful in your research, we would appreciate it if you mention in the acknowledgement. If you have any questions, please contact Zheng Zhou (周正) at [physics@zhengzhou.page](mailto:physics@zhengzhou.page).
