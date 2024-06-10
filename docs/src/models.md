# Built-in models

## Diagonal quantum numbers on fuzzy sphere

The following diagonal quantum numbers (symmetry charges) on fuzzy sphere are built in in FuzzifiED.

```@docs
GetNeQNDiag(no :: Int64)
GetLz2QNDiag(nm :: Int64, nf :: Int64)
GetFlavQNDiag(nm :: Int64, nf :: Int64, qf :: Union{Dict{Int64, Int64}, Vector{Int64}}, id :: Int64 = 1, modul :: Int64 = 1)
GetZnfChargeQNDiag(nm :: Int64, nf :: Int64)
GetPinOrbQNDiag(no :: Int64, pin_o :: Vector{Int64}, id :: Int64 = 1) 
```

## Off-diagonal quantum numbers on fuzzy sphere

The following off-diagonal quantum numbers (transformations) on fuzzy sphere are built in in FuzzifiED.

```@docs
GetParityQNOffd(nm :: Int64, nf :: Int64, permf :: Union{Dict{Int64, Int64}, Vector{Vector{Int64}}, Vector{Int64}} = Dict{Int64, Int64}(), fac :: Union{Dict{Int64, <: Number}, Vector{<: Number}} = Dict{Int64, ComplexF64}()) 
GetFlavPermQNOffd(nm :: Int64, nf :: Int64, permf :: Union{Dict{Int64, Int64}, Vector{Vector{Int64}}, Vector{Int64}}, fac :: Union{Dict{Int64, <: Number}, Vector{<: Number}} = Dict{Int64, ComplexF64}())
GetRotyQNOffd(nm :: Int64, nf :: Int64)
```

## Operators on fuzzy sphere

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
Electron(nm :: Int64, nf :: Int64, f :: Int64)
Density(nm :: Int64, nf :: Int64 ; mat :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf))
```