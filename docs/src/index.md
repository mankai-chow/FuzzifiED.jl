# FuzzifiED.jl

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculation on the fuzzy sphere, and also facilitates the DMRG calculations by ITensors. It can also be used for generic fermion models. 

If this package is helpful in your research, we would appreciate it if you mention in the acknowledgement. If you have any questions, please contact Zheng Zhou (周正) at [zzhou@pitp.ca](mailto:zzhou@pitp.ca).

## Install

Download Julia at [this link](https://julialang.org/downloads/). 

To Apple Silicon (Mac M1, M2, etc.) users : The support for AArch64 architecture is not guarenteed. You may want to use Julia for x86 (Intel or Rosetta) instead of macOS (Apple Silicon).

To install the package, please first enter Julia by entering in the command line `julia`, and then enter the commands
```julia
using Pkg ;
Pkg.add(url="https://github.com/mankai-chow/FuzzifiED_jll.jl.git") ;
Pkg.add(url="https://github.com/mankai-chow/FuzzifiED.jl.git") ;
```
Include at the start of your Julia script
```julia
using FuzzifiED
```
Jupyter Notebook is highly recommended as it allows you to run Julia (and Python) just like running a Mathematica notebook.

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
