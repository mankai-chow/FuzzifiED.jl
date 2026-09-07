# SO(3)lver 

SO(3)lver is an extension of the package [FuzzifiED](https://docs.fuzzified.world). It enables exact diagonalization (ED) on the fuzzy sphere with the implementation of full $\mathrm{SO}(3)$ spatial rotation symmetry and some non-Abelian flavour symmetry through partitioning the Hilbert space into segments and compose these segments. 

The implementation of $\mathrm{SO}(3)$ symmetry greatly reduces the demand for the memory and enables calculation for much larger size. _E. g._, with the $\mathrm{SO}(3)$ resolution, the Hilbert space of the ground-state sector for $N_m=22$ has a dimension $2.07\times 10^7$, and the Hamiltonian requires a memory of 76.5 Go ; as a comparison, without the $\mathrm{SO}(3)$ resolution, the Hilbert space for $N_m=18$ has a similar dimension $2.84\times 10^7$, and the Hamiltonian requires a memory of 63.9 Go.

SO(3)lver can apply to a lot of current fuzzy-sphere research. Its flexible design also makes it straightforward to adapt to new models. We provide [a collection of examples](@ref Examples-Using-SO(3)lver).

## Installation

If you have the permission to the package, you first need to authenticate GitHub CLI by running the following command in the terminal
```bash
gh auth login
```
and then run the following command in the Julia REPL (read-eval-print loop) (To enter Julia REPL, simply type `julia` in the command line) 
```julia
using Pkg
Pkg.add(url="https://github.com/FuzzifiED/SO3lver.jl")
```
To use the package, include at the start of the Julia script
```julia
using FuzzifiED
using FuzzifiED.SO3lver
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

## Index 

```@index
Pages = [
    "interface.md"
]
```
