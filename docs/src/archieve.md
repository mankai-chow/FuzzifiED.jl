# Archieved interfaces

The following functions are used in the historical versions and might be obsolete in the future versions.

## Core 

```@docs
Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2))
Basis(cfs :: Confs, qnz_s :: Vector{ComplexF64} ; cyc :: Vector{Int64}, perm_o :: Vector{Any}, ph_o :: Vector{Any}, fac_o :: Vector{Any})
GetEntSpec(st :: Vector{<:Number}, bs0 :: Basis, qnu_s_lst :: Vector{Any}, qnz_s_lst :: Vector{Any} ; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], cyc :: Vector{Int64}, perm_o :: Vector{Any}, ph_o :: Vector{Any}, fac_o :: Vector{Any}, amp_oa :: Vector{<:Number}, amp_ob :: Vector{<:Number} = sqrt.(1 .- abs.(amp_oa .^ 2)))
```

## ITensor support 

```@docs
TruncateQnu(; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], trunc_lth :: Int64 = 3, trunc_wt :: Vector{Int64} = [ 2 ^ (i - trunc_lth) for i = trunc_lth : length(qnu_o)]) 
SitesFromQnu(; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)])
GetMPOSites(id :: String, tms :: Union{Vector{Term}, Sum{Scaled{ComplexF64, Prod{Op}}}} ; path :: String = "./", qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], mpo_method :: Function = MPO)
```

## Built-in models 

```@docs
GetLzQnu(nm :: Int64, nf :: Int64)
GetLzConfs(nm :: Int64, nf :: Int64, ne :: Int64 ; lz :: Float64 = 0.0)
GetLzZnQnu(nm :: Int64, nf :: Int64)
GetLzZnConfs(nm :: Int64, nf :: Int64, ne :: Int64 ; lz :: Float64 = 0.0, zn :: Int64 = 0)
GetIsingBasis(cfs :: Confs ; qn_p :: Int64 = 0, qn_z :: Int64 = 0, qn_r :: Int64 = 0)
GetSnBasis(cfs :: Confs, nf :: Int64 ; qn_r :: Int64 = 0, perm :: Vector = [], qn_z :: Vector{<:Number} = Number[]) 
GetIsingIntTerms(nm :: Int64 ; ps_pot :: Vector = [1.])
GetXPolTerms(nm :: Int64)
GetZPolTerms(nm :: Int64)
```