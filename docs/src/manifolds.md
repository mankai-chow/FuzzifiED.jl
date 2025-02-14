# Fuzzy manifolds

FuzzyManifolds is a module to support calculations on other manifolds regularised by lowest Landau level besides sphere. The supported manifolds include torus $T^2$. 

## Torus 

```@docs
GetTorusLz2QNDiag(nm :: Int64, nf :: Int64) 
GetTorusIntMatrix(nm :: Int64, lx :: Number, ps_pot :: Vector{<:Number})
GetTorusDenIntTerms(nm :: Int64, nf :: Int64, lx :: Number, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf), mat_b :: Matrix{<:Number} = Matrix(mat_a'))
GetTorusPairIntTerms(nm :: Int64, nf :: Int64, lx :: Number, ps_pot :: Vector{<:Number}, mat_a :: Matrix{<:Number}, mat_b :: Matrix{<:Number} = Matrix(mat_a'))
```

## Related examples 

* [`circle_ising.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/circle_ising.jl) calculates the spectrum of 2d Ising CFT on a fuzzy thin torus. This example reproduces Figure 4 and Tables I--III in [Han 2025](@ref Han2025)
