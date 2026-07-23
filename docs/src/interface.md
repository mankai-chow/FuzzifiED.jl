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

## The Boson-Fermion Mixture

For bosonic degrees of freedom — and for mixed fermion–boson systems such as gauge theories with monopole/flux modes — FuzzifiEDFullRotation provides a parallel set of types and functions built on the [Fuzzifino](https://docs.fuzzified.world/fuzzifino/) module of FuzzifiED. Their names carry an `S` prefix : `BuildSSegSpace`, `SCoupleDecomp`, `ContactSCouple`, `SingleSegSCouple`, `BuildSSegOperators`, `BuildSCompSpace` and `BuildSCompOperator`. The workflow is identical to the fermionic one described above ; a segment may then be a bosonic mode, and fermionic and bosonic segments are coupled by the same $\mathrm{SO}(3)$ recoupling. See the [`majorana_spectrum.jl`](@ref List-of-Examples), [`u1_2_higgs_spectrum.jl`](@ref List-of-Examples) and [`qed_cs52_nf1_spectrum.jl`](@ref List-of-Examples) examples.
