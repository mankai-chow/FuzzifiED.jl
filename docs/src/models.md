# Built-in models

## General models on fuzzy sphere

```@docs 
GetIntMatrix(nm :: Int64, ps_pot :: Vector{Number})
GetL2Terms(nm :: Int64, nf :: Int64)
GetDenIntTerms(nm :: Int64, nf :: Int64 ; ps_pot :: Vector{<:Number} = [1.0], mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a'))
GetPairIntTerms(nm :: Int64, nf :: Int64 ; ps_pot :: Vector{<:Number} = [1.0], mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a'))
GetPolTerms(nm :: Int64, nf :: Int64 ; mat :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf))
```

## Spherical observables

FuzzifiED supports local observables on sphere that can be decomposed into angular components ``\Phi(\Omega)=\sum_{lm}\Phi_{lm}Y^{(s)}_{lm}``
```@docs
SphereObs
```
It can be initialised with the following methods 
```@docs
SphereObs(s2 :: Int64, l2m :: Int64, get_comp :: Function)
SphereObs(s2 :: Int64, l2m :: Int64, cmps :: Dict{Tuple{Int64, Int64}, Vector{Term}})
```
The following methods explicitly calculates and stores each component
```@docs
StoreComps!(obs :: SphereObs)
StoreComps(obs :: SphereObs)
```
The multiplication, addition, conjugate and Laplacian operation of an observable is supported 
```@docs
*(fac :: Number, obs :: SphereObs) 
+(obs1 :: SphereObs, obs2 :: SphereObs) 
adjoint(obs :: SphereObs)
*(obs1 :: SphereObs, obs2 :: SphereObs)
Laplacian(obs :: SphereObs)
```
The observables can be evaluated either at an angular component or at a real-space point.
```@docs
GetComponent(obs :: SphereObs, l :: Number, m :: Number)
GetPointValue(obs :: SphereObs, θ :: Float64, ϕ :: Float64)
```
Two types of operators, _viz._ electrons and density operators are built-in.
```@docs
Electron(nf :: Int64, nm :: Int64, f :: Int64)
Density(nf :: Int64, nm :: Int64 ; mat :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf))
```

## Ising model and Potts model

The following methods is especially helpful for Ising model and Potts model.

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

## ``\mathrm{Sp}(N)`` model

The following methods is especially helpful for ``\mathrm{Sp}(N)`` model and Potts model.

```@docs
GetSpnConfs(nm :: Int64, nf :: Int64, ne :: Int64 ; lz :: Float64 = 0.0, sz :: Vector{Int64} = fill(0, div(nf, 2)))
GetSpnBasis(cfs :: Confs, nf :: Int64 ; qn_p :: Int64 = 0, qn_r :: Int64 = 0, qn_z :: Vector{Int64} = fill(0, div(nf, 2)), qn_x :: Vector{Int64} = fill(0, div(nf, 4)))
GetSpnPairIntTerms(nm :: Int64, nf :: Int64 ; ps_pot :: Vector{<:Number} = [1.])
GetSpnC2Terms(nm :: Int64, nf :: Int64) 
```