# Reference

## Basics 

### Configurations
```@docs
Confs
Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2))
```

### Basis
```@docs
Basis
Basis(cfs :: Confs, qnz_s :: Vector{ComplexF64}, cyc :: Vector{Int64}, perm_o :: Vector{Any}, ph_o :: Vector{Any}, fac_o :: Vector{Any})
Basis(cfs :: Confs)
```
### Term

```@docs
Term
*(fac :: Number, tms :: Vector{Term})
+(tms1 :: Vector{Term}, tms2 :: Vector{Term})
*(tms1 :: Vector{Term}, tms2 :: Vector{Term})
adjoint(tms :: Vector{Term})
```

### Operator

```@docs
Operator
Operator(bsd :: Basis, bsf :: Basis, terms :: Vector{Term} ; red_q :: Int64 = 0, sym_q :: Int64 = 0)
*(op :: Operator, st_d :: Vector{ComplexF64})
*(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, op :: Operator, st_d :: Vector{ComplexF64})
```

### Sparse matrix

```@docs
OpMat
OpMat(op :: Operator)
*(mat :: OpMat, st_d :: Vector{ComplexF64})
*(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat, st_d :: Vector{ComplexF64})
GetEigensystem(mat :: OpMat, nst :: Int64 ; tol :: Float64 = 1E-8)
```

## ITensors support

```@docs
ConfsFromSites(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, qn_s :: QN)
ConfsFromSites(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, cf_ref :: Vector{Int64})
OperatorFromOpSum(bsd :: Basis, bsf :: Basis, opsum :: Sum{Scaled{ComplexF64, Prod{Op}}} ; red_q :: Int64 = 0, sym_q :: Int64 = 0)
```

## Built-in models

```@docs
GetInteractionMatrix(nm :: Int64, ps_pot :: Vector{Number})
GenerateL2Terms(nm :: Int64, nf :: Int64)
```

### Ising model

```@docs
GenerateIsingConfs(nm :: Int64 ; ne :: Int64 = 0, lz :: Float64 = 0.0)
GenerateIsingBasis(cfs :: Confs ; PH :: Int64 = 0, Z2 :: Int64 = 0, Ry :: Int64 = 0)
GenerateIsingHamiltonianTerms(nm :: Int64 ; ps_pot :: Vector{Number} = [4.75, 1], fld_h :: Number = 3.16)
```