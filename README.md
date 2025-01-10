# FuzzifiED.jl Version 1.0.0

Since its proposal, the fuzzy sphere regularisation has made significant contribution to the study of 3d CFTs. This Julia package FuzzifiED is aimed at simplifying the numerical calculations on the fuzzy sphere. It facilitates the exact diagonalisation (ED) calculations as well as the density matrix renormalisation group (DMRG) with the help of ITensor. It can also be used for generic fermionic and bosonic models. This package features the following characteristics : 

* Versatality : FuzzifiED can help reproduce almost all the ED and DMRG results in fuzzy sphere works, and it is easy and flexible for the adaption to new models. 
* Usability : Julia interfaces make the code intuitive and short. To help the users get started, we have also provided a collection of examples.
* Efficiency : FuzzifiED can produce results on reasonable system sizes within minutes.
* Open source : The code for FuzzifiED is fully open source. 

Documentations can be found at [docs.fuzzified.world](https://docs.fuzzified.world/). A PDF version is provided at [this link](https://docs.fuzzified.world/assets/FuzzifiED_Documentation.pdf).

If you have any questions, please contact Zheng Zhou (周正) at [physics@zhengzhou.page](mailto:physics@zhengzhou.page).

Download Julia at [this link](https://julialang.org/downloads/). 

To install the package, run the following command in the Julia REPL (read-eval-print loop) (To enter Julia REPL, simply type `julia` in the command line) 
```julia
using Pkg ; Pkg.add("FuzzifiED")
```
To use the package, include at the start of the Julia script
```julia
using FuzzifiED
```
To obtain the documentation for an interface, type `?` followed by the keyword in the Julia REPL, _e.g._, `?Confs`.

If this package is helpful in your research, we would appreciate it if you mention in the acknowledgement. We have also provided a BibTeX file that includes all the works on the fuzzy sphere works at [this link](https://docs.fuzzified.world/assets/bib_fuzzy.bib).