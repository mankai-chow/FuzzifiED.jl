# FuzzifiED.jl

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculations on the fuzzy sphere, and also facilitates the DMRG calculations by ITensor. It can also be used for generic fermionic and bosonic models. Using this package, you can reproduce almost all the ED results in fuzzy sphere works.

If this package is helpful in your research, we would appreciate it if you mention in the acknowledgement. If you have any questions, please contact Zheng Zhou (周正) at [fuzzified@zhengzhou.page](mailto:fuzzified@zhengzhou.page).

## Install

To install the package, please first enter Julia by entering in the command line `julia`, and then enter the commands
```julia
using Pkg ; Pkg.add("FuzzifiED")
```
Include at the start of your Julia script
```julia
using FuzzifiED
```

## Useful information

- Download Julia at [this link](https://julialang.org/downloads/). 
- We have registered FuzzifiED at Julia General Registry ! The package regisitry may have some delay, to bring up to date, use `Pkg.Registry.update()`, or install from the github repositories 
```Julia
using Pkg
Pkg.add(url="https://github.com/FuzzifiED/FuzzifiED_jll.jl")
Pkg.add(url="https://github.com/FuzzifiED/FuzzifiED.jl")
```
- We are migrating the GitHub repositories from the personal account to the organisation account `FuzzifiED`. The old repositories will still be accessible. 
- Jupyter Notebook is highly recommended as it allows you to run Julia (and Python) just like running a Mathematica notebook. _N.b._, you may need to install the package `IJulia` by hand to use Jupyter notebook ; in Jupyter notebooks, you may need to define how many threads OpenMP uses by hand in `FuzzifiED.NumThreads`.
- The package is under active development, so certain interfaces may get changed, superceded or obsolete. We are sorry for any possible inconvenience. 
- For the DMRG calculation, due to the change of interface in update of ITensors v0.7, now both package `ITensors` and `ITensorMPS` need to be installed. 
- The supporting Fortran code for ED is stored at the GitHub Repo at [FuzzifiED/FuzzifiED_Fortran](https://github.com/FuzzifiED/FuzzifiED_Fortran).

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
    "releases.md"
]
Depth = 2
```

## References

* __[Zhu 2022]__ Uncovering conformal symmetry in the 3d Ising transition : state-operator correspondence from a quantum fuzzy sphere regularisation, Wei Zhu, Chao Han, Emilie Huffman, Johannes S. Hofmann, and Yin-Chen He, [arXiv:2210.13482](https://arxiv.org/abs/2210.13482), [Phys. Rev. X __13__, 021009 (2023)](https://doi.org/10.1103/PhysRevX.13.021009).
* __[Hu 2023Mar]__ Operator product expansion coefficients of the 3d Ising criticality via quantum fuzzy sphere, Liangdong Hu, Yin-Chen He, and Wei Zhu, [arXiv:2303.08844](https://arxiv.org/abs/2303.08844), [Phys. Rev. Lett __131__, 031601 (2023)](https://doi.org/10.1103/PhysRevLett.131.031601).
* __[Han 2023Jun]__ Conformal four-point correlators of the 3d Ising transition via the quantum fuzzy sphere, Chao Han, Liangdong Hu, Wei Zhu, and Yin-Chen He, [arXiv:2306.04681](https://arxiv.org/abs/2306.04681), [Phys. Rev. B __108__, 235123 (2023)](https://doi.org/10.1103/PhysRevB.108.235123).
* __[Zhou 2023]__ The ``\mathrm{SO}(5)`` deconfined phase transition under the fuzzy sphere microscope: approximate conformal symmetry, pseudo-criticality, and operator spectrum, Zheng Zhou, Liangdong Hu, Wei Zhu, and Yin-Chen He, [arXiv:2306.16435](https://arxiv.org/abs/2306.16435), [Phys. Rev. X __14__, 021044 (2024)](https://doi.org/10.1103/PhysRevX.14.021044).
* __[Lao 2023]__ 3d Ising CFT and exact diagonalisation on icosahedron : the power of conformal perturbation theory, Bing-Xin Lao, and Slava Rychkov [arXiv:2307.02540](https://arxiv.org/abs/2307.02540), [SciPost Phys. __15__, 243 (2023)](https://doi.org/10.21468/SciPostPhys.15.6.243).
* __[Hu 2023Aug]__ Solving conformal defects in 3d conformal field theory using fuzzy sphere regularisation, Liangdong Hu, Yin-Chen He, and Wei Zhu, [arXiv:2308.01903](https://arxiv.org/abs/2308.01903), [Nat. Commun. __15__, 3659 (2024)](https://doi.org/10.1038/s41467-024-47978-y).
* __[Hofmann 2024]__ Quantum Monte Carlo simulation of the 3d Ising transition on the fuzzy sphere, Johannes S. Hofmann, Florian Goth, Wei Zhu, Yin-Chen He, and Emilie Huffman, [arXiv:2310.19880](https://arxiv.org/abs/2310.19880), [SciPost Phys. Core __7__, 028 (2024)](https://doi.org/10.21468/SciPostPhysCore.7.2.028).
* __[Han 2023Dec]__ Conformal operator content of the Wilson-Fisher transition on fuzzy sphere bilayers, Chao Han, Liangdong Hu, and Wei Zhu, [arXiv:2312.04047](https://arxiv.org/abs/2312.04047), [Phys. Rev. B __110__, 115113 (2024)](https://doi.org/10.1103/PhysRevB.110.115113).
* __[Zhou 2024Jan]__ The ``g``-function and defect changing operators from wavefunction overlap on a fuzzy sphere, Zheng Zhou, Davide Gaiotto, Yin-Chen He, Yijian Zou, [arXiv:2401.00039](https://arxiv.org/abs/2401.00039), [SciPost Phys. __17__, 021 (2024)](https://doi.org/10.21468/SciPostPhys.17.1.021).
* __[Hu 2024]__ Entropic ``F``-function of 3d Ising conformal field theory via the fuzzy sphere regularisation, Liangdong Hu, Wei Zhu, and Yin-Chen He, [arXiv:2401.17362](https://arxiv.org/abs/2401.17362).
* __[Cuomo 2024]__ Impurities with a cusp : general theory and 3d Ising, Gabriel Cuomo, Yin-Chen He, Zohar Komargodski, [arXiv:2406.10186](https://arxiv.org/abs/2406.10186). 
* __[Zhou 2024Jul]__ Studying the 3d Ising surface CFTs on the fuzzy sphere, Zheng Zhou, and Yijian Zou, [arXiv:2407.15914](https://arxiv.org/abs/2407.15914).
* __[Dedushenko 2024]__ Ising BCFTs from the fuzzy hemisphere, Mykola Dedushenko, [arXiv:2407.15948](https://arxiv.org/abs/2407.15948).
* __[Fardelli 2024]__ Constructing the infrared conformal generators on the fuzzy sphere, Giulia Fardelli, A. Liam Fitzpatrick, and Emanuel Katz, [arXiv:2409.02998](https://arxiv.org/abs/2409.02998).
* __[Fan 2024]__ Note on explicit construction of conformal generators on the fuzzy sphere, Ruihua Fan, [arXiv:2409.08257](https://arxiv.org/abs/2409.08257).
* __[Zhou 2024Oct]__ A new series of 3d CFTs with ``\mathrm{Sp}(N)`` global symmetry on fuzzy sphere, Zheng Zhou, and Yin-Chen He, [arXiv:2410.00087](https://arxiv.org/abs/2410.00087).
* __[Voinea 2024]__ Regularising 3d conformal field theories via anyons on the fuzzy sphere, [arXiv:2411.15299](https://arxiv.org/abs/2411.15299).

## Index 

```@index
Pages = [
    "core.md",
    "models.md",
    "itensors.md",
    "extension.md",
    "fuzzifino.md"
]
```