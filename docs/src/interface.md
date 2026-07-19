# Interface 

## Segment Spaces

```@docs
SegSpace
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
Base.:+(cpd1 :: Union{CoupleDecomp, Vector{CoupleDecomp}}, cpd2 :: Union{CoupleDecomp, Vector{CoupleDecomp}})
```

Several methods are used to generate contact couplings and operators that act only on one segment. 
```@docs
ContactCouple
SingleSegCouple
```

Several methods help recoupling pseudo-potantials, and converting pseudo-potentials to coupling channels
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
