# FuzzifiED.jl

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculations on the fuzzy sphere, and also facilitates the DMRG calculations by ITensors. It can also be used for generic fermion models. Using this package, you can reproduce almost all the ED results in fuzzy sphere works.

If this package is helpful in your research, we would appreciate it if you mention in the acknowledgement. If you have any questions, please contact Zheng Zhou (周正) at [physics@zhengzhou.page](mailto:physics@zhengzhou.page).

## Install

To install the package, please first enter Julia by entering in the command line `julia`, and then enter the commands
```julia
using Pkg
Pkg.add(url="https://github.com/mankai-chow/FuzzifiED_jll.jl.git")
Pkg.add(url="https://github.com/mankai-chow/FuzzifiED.jl.git")
```
Include at the start of your Julia script
```julia
using FuzzifiED
```

## Tips 

- Download Julia at [this link](https://julialang.org/downloads/). 
- Apple Silicon (Mac M1, M2, _etc._) users may experience trouble using `GetEigensystem` function. Another function [`GetEigensystemKrylov`](@ref) is provided instead. 
- Jupyter Notebook is highly recommended as it allows you to run Julia (and Python) just like running a Mathematica notebook. _N.b._, you may need to install the package `IJulia` by hand to use Jupyter notebook ; in Jupyter notebooks, you may need to define how many threads OpenMP uses by hand in `FuzzifiED.NumThreads`.
- The package is under active development, so certain interfaces may get changed, superceded or obsolete. We are sorry for any possible inconvenience. 
- For the DMRG calculation, due to the change of interface in update of ITensors v0.7, now both package `ITensors` and `ITensorMPS` need to be installed. 

## Outline 

```@contents
Pages = [
    "intro.md",
    "example.md",
    "core.md",
    "models.md",
    "itensors.md",
    "extension.md",
    "fuzzifino.md",
    "archieve.md",
    "releases.md"
]
Depth = 2
```

## Index 

```@index
Pages = [
    "core.md",
    "models.md",
    "itensors.md",
    "extension.md",
    "fuzzifino.md",
    "archieve.md"
]
```