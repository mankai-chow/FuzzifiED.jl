# FuzzifiED.jl

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculation on the fuzzy sphere, and also facilitates the DMRG calculations by ITensors. It can also be used for generic fermion models. Using this package, you can reproduce almost all the ED results in fuzzy sphere works ; for detail, see [« The fuzzified world »](https://www.fuzzified.world/fuzzified-world).

If this package is helpful in your research, we would appreciate it if you mention in the acknowledgement. If you have any questions, please contact Zheng Zhou (周正) at [physics@zhengzhou.page](mailto:physics@zhengzhou.page).

## Install

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

## Tips 

- Download Julia at [this link](https://julialang.org/downloads/). 
- To Apple Silicon (Mac M1, M2, etc.) users : The support for AArch64 architecture is not guarenteed. You may want to use Julia for x86 (Intel or Rosetta) instead of macOS (Apple Silicon).
- Jupyter Notebook is highly recommended as it allows you to run Julia (and Python) just like running a Mathematica notebook. _N.b._, you may need to install the package `IJulia` by hand to use Jupyter notebook ; in Jupyter notebooks, you may need to define how many threads OpenMP uses by hand in `FuzzifiED.NumThreads` 
- We are currently trying to train our own GPT [« Getting FuzzifiED »](https://chatgpt.com/g/g-WvSuxsXus-getting-fuzzified), but the performance is not ideal. Please help us out. 
- The package is under active development, so certain interfaces may get changed, superceded or obsolete. We are sorry for any possible inconvenience. 

## Outline 

```@contents
Pages = [
    "example.md",
    "core.md",
    "itensors.md",
    "models.md",
    "archieve.md"
]
Depth = 2
```

## What is new (0.9)

- Add angular modes observables.
- Allow input of initial vectors for diagonalisation. (We acknowledge Andrew Fitzpatrick for the suggestion.)
- Add the example of Ising generators and Sp(3) CFT.

## References

For a more detailed summary of the background, please visit [« The fuzzified world »](https://www.fuzzified.world/fuzzified-world)

* Uncovering Conformal Symmetry in the 3D Ising Transition: State-Operator Correspondence from a Quantum Fuzzy Sphere Regularization, Wei Zhu, Chao Han, Emilie Huffman, Johannes S. Hofmann, and Yin-Chen He, [Phys. Rev. X 13, 021009 (2023)](https://doi.org/10.1103/PhysRevX.13.021009).
* Operator Product Expansion Coefficients of the 3D Ising Criticality via Quantum Fuzzy Sphere, Liangdong Hu, Yin-Chen He, and Wei Zhu, [Phys. Rev. Lett 131, 031601 (2023)](https://doi.org/10.1103/PhysRevLett.131.031601).
* Conformal four-point correlators of the 3D Ising transition via the quantum fuzzy sphere, Chao Han, Liangdong Hu, Wei Zhu, and Yin-Chen He, [Phys. Rev. B 108, 235123 (2023)](https://doi.org/10.1103/PhysRevB.108.235123).
* The ``\mathrm{SO}(5)`` Deconfined Phase Transition under the Fuzzy Sphere Microscope: Approximate Conformal Symmetry, Pseudo-Criticality, and Operator Spectrum, Zheng Zhou, Liangdong Hu, Wei Zhu, and Yin-Chen He, [Phys. Rev. X 14, 021044 (2024)](https://doi.org/10.1103/PhysRevX.14.021044).
* Solving Conformal Defects in 3D Conformal Field Theory using Fuzzy Sphere Regularization, Liangdong Hu, Yin-Chen He, and Wei Zhu, [Nat. Commun. 15, 3659 (2024)](https://doi.org/10.1038/s41467-024-47978-y).
* Quantum Monte Carlo Simulation of the 3D Ising Transition on the Fuzzy Sphere, Johannes S. Hofmann, Florian Goth, Wei Zhu, Yin-Chen He, and Emilie Huffman, [SciPost Phys. Core 7, 028 (2024)](https://doi.org/10.21468/SciPostPhysCore.7.2.028).
* Conformal Operator Content of the Wilson-Fisher Transition on Fuzzy Sphere Bilayers, Chao Han, Liangdong Hu, and Wei Zhu, [Phys. Rev. B 110, 115113 (2024)](https://doi.org/10.1103/PhysRevB.110.115113).
* The ``g``-function and Defect Changing Operators from Wavefunction Overlap on a Fuzzy Sphere, Zheng Zhou, Davide Gaiotto, Yin-Chen He, Yijian Zou, [SciPost Phys. 17, 021 (2024)](https://doi.org/10.21468/SciPostPhys.17.1.021).
* Entropic ``F``-function of 3D Ising conformal field theory via the fuzzy sphere regularization, Liangdong Hu, Wei Zhu, and Yin-Chen He, [arXiv : 2401.17362](https://arxiv.org/abs/2401.17362).
* Impurities with a cusp: general theory and 3d Ising, Gabriel Cuomo, Yin-Chen He, Zohar Komargodski, [arXiv : 2406.10186](https://arxiv.org/abs/2406.10186). 
* Studying the 3d Ising surface CFTs on the fuzzy sphere, Zheng Zhou, and Yijian Zou, [arXiv : 2407.15914](https://arxiv.org/abs/2407.15914).
* Ising BCFTs from the fuzzy hemisphere, Mykola Dedushenko, [arXiv : 2407.15948](https://arxiv.org/abs/2407.15948).
* Constructing the Infrared Conformal Generators on the Fuzzy Sphere, Giulia Fardelli, A. Liam Fitzpatrick, and Emanuel Katz, [arXiv : 2409.02998](https://arxiv.org/abs/2409.02998).
* Note on explicit construction of conformal generators on the fuzzy sphere, Ruihua Fan, [arXiv : 2409.08257](https://arxiv.org/abs/2409.08257).
* A new series of 3D CFTs with ``\mathrm{Sp}(N)`` global symmetry on fuzzy sphere, Zheng Zhou, and Yin-Chen He, [arXiv : 2410.00087](https://arxiv.org/abs/2410.00087).

## Index 

```@index
Pages = [
    "core.md",
    "itensors.md",
    "models.md",
    "archieve.md"
]
```