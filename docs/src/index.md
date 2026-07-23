# FuzzifiEDFullRotation.jl 

FuzzifiED Full Rotation is an extension of the package [FuzzifiED](https://docs.fuzzified.world). It enables exact diagonalization (ED) on the fuzzy sphere with the implementation of full $\mathrm{SO}(3)$ spatial rotation symmetry and some non-abelian flavour symmetry through partitioning the Hilbert space into segments and compose these segments. It allows the calculation 

## Installation

If you have the permission to the package, you first need to authenticate GitHub CLI by running the following command in the terminal
```bash
gh auth login
```
and then run the following command in the Julia REPL (read-eval-print loop) (To enter Julia REPL, simply type `julia` in the command line) 
```julia
using Pkg
Pkg.add(url="https://github.com/FuzzifiED/FuzzifiEDFullRotation.jl")
```
To use the package, include at the start of the Julia script
```julia
using FuzzifiED
using FuzzifiEDFullRotation
```
To suppress the log into `stderr`, load the `Logging` module.
```julia
using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn))
```

## Outline 

```@contents
Pages = [
    "formalism.md",
    "tutorial.md",
    "interface.md"
]
Depth = 2
```