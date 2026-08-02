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

Spaces with different $C_2$ can be built simultaneously.
```@docs
BuildSegSpaces
```

## Composite Spaces

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

Several methods help re-coupling angular momenta and converting pseudo-potentials to coupling channel.
```@docs
RecoupleAngMom
ConvPsPot
```

## Segmented Operators

```@docs
SegOperator
BuildSegOperator
BuildSegOperators
```

## Composite Operators

```@docs
CompOperator
BuildCompOperator
Base.:*(cpop :: CompOperator{T}, std :: Vector{T}) where T <: Union{Float64, ComplexF64}
Base.Matrix(cpop :: CompOperator{T}) where T <: Union{Float64, ComplexF64}
GetEigensystem(cpop :: CompOperator{ComplexF64}, nst :: Int64)
```
