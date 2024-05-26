## Reference

### Configurations
```@docs
Confs
Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2))
ConfsFromSites(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, qn_s :: QN)
ConfsFromSites(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, cf_ref :: Vector{Int64})
```

### Basis
```@docs
Basis
Basis(cfs :: Confs, qnz_s :: Vector{ComplexF64}, cyc :: Vector{Int64}, perm_o :: Vector{Any}, ph_o :: Vector{Any}, fac_o :: Vector{Any})
Basis(cfs :: Confs)
```

### Operator

```@docs
Term
Operator
Operator(bsd :: Basis, bsf :: Basis, terms :: Vector{Term} ; red_q :: Int64 = 0, sym_q :: Int64 = 0)
OperatorFromOpSum(bsd :: Basis, bsf :: Basis, opsum :: Sum{Scaled{ComplexF64, Prod{Op}}} ; red_q :: Int64 = 0, sym_q :: Int64 = 0)
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