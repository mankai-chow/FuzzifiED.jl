# ITensors support

This package also supports convertion with the `Site` and `OpSum` objects from `ITensors` library and management of DMRG sweeping process. To use these functions, include
```julia
using ITensors
```
at the heading of the Julia script.

## Format conversion

The `Sites` objects in `ITensors` can be converted to a `Confs` object with the conserved quantities extracted. 
```@docs
ConfsFromSites(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, qn_s :: QN)
ConfsFromSites(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, cf_ref :: Vector{Int64})
```
Conversely, the  `Sites` objects in `ITensors` can also be generated from a set of diagonal quantum numbers 
```@docs
SitesFromQnu(; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)])
```
If the number of quantum numbers are too many, it can be truncated by 
```@docs
TruncateQnu(; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], trunc_lth :: Int64 = 3, trunc_wt :: Vector{Int64} = [ 2 ^ (i - trunc_lth) for i = trunc_lth : length(qnu_o)]) 
```

The `OpSum` objects in `ITensors` can be converted with the collection of `Term`'s
```@docs
TermsFromOpSum(opsum :: Sum{Scaled{ComplexF64, Prod{Op}}})
OpSumFromTerms(tms :: Vector{Term})
```

## Easy sweep

This tool kit facilitates the management of DMRG process. It automatically records the intermediate results and recover these results if a job is stopped and run again on HPC. It also manages the gradual increase of maximal bond dimensions and the determination of convergence by the criteria of energy. These functions require the package [`ITensorMPOConstruction`](https://github.com/ITensor/ITensorMPOConstruction.jl), which can be installed by 
```julia
julia> using Pkg; Pkg.add(url="https://github.com/ITensor/ITensorMPOConstruction.jl.git")
```

```@docs
EasySweep(id :: String, hmt :: MPO, st00 :: MPS ; path :: String = "./", dim_list :: Vector{Int64} = [1000,2000,3000,4000,5000,6000], proj :: Vector{String} = String[], e_tol1 :: Float64 = 1E-6, e_tol :: Float64 = 1E-7, cutoff :: Vector{Float64} = [1E-9], maxdim0 :: Vector{Int64} = [10,20,50,100,200,500], noise0 :: Vector{Float64} = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], noise :: Vector{Float64} = [1E-6,2E-7,5E-8,1E-8,0], nsweeps :: Int64 = 10, weight :: Float64 = 100.0, observer :: AbstractObserver = EasySweepObserver(e_tol1))
SweepOne(id :: String, hmt :: MPO, st0 :: MPS, dim1 :: Int64 ; path :: String = "./", cutoff :: Vector{Float64} = [1E-9], maxdim :: Vector{Int64} = [dim1], nsweeps :: Int64 = 10, noise :: Vector{Float64} = [1E-6,1E-7,0], proj :: Vector{String} = String[], e_tol :: Float64 = 1E-6, weight :: Float64 = 100.0, observer :: AbstractObserver = EasySweepObserver(e_tol))
GetMpoSites(id :: String, tms :: Union{Vector{Term}, Sum{Scaled{ComplexF64, Prod{Op}}}} ; path :: String = "./", qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], mpo_method :: Function = MPO)
GetMpo(id :: String, tms :: Union{Vector{Term}, Sum{Scaled{ComplexF64, Prod{Op}}}}, sites :: Vector{<:Index} ; path :: String = "./", mpo_method :: Function = MPO)
```