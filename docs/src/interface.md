# Interface 

## Segment Spaces

The segment Hilber space is stored in a [`SegSpace`](@ref) (when fully fermionic) or [`SSegSpace`](@ref) (when bosonic or mixed), both as subtypes of [`AbstractSegSpace`](@ref). 

```@docs
SegSpace
SSegSpace
AbstractSegSpace
```

They can be built by
```@docs
BuildSegSpace
```

## Composed Spaces

```@docs
CompSpace
BuildCompSpace
```

Several methods help comparing, composing sectors, and finding coupling channels of angular momenta.
```@docs
EquivSec
ComposeSec
FindCouplingChannels
```

## Coupling Decompositions

```@docs
CoupleDecomp
CoupleDecomps
Base.:+(cpd1 :: CoupleDecomps, cpd2 :: CoupleDecomps)
```

Several methods are used to generate contact couplings and operators that act only on one segment. 
```@docs
ContactCouple
SingleSegCouple
InsertSegment
```

Several methods help recoupling pseudo-potentials, converting them to coupling channels, and preparing an assembled operator.
```@docs
RecouplePsPot
ConvPsPot
```

## Segmented Operators

```@docs
SegOperator
BuildSegOperator
BuildSegOperators
```

## Composed Operators

```@docs
CompOperator
BuildCompOperator
Base.:*(cpop :: CompOperator{T}, std :: Vector{T}) where T <: Union{Float64, ComplexF64}
GetEigensystem(cpop :: CompOperator{ComplexF64}, nst :: Int64)
```
