"""
    mutable struct Basis

# Fields
* `cfs :: Confs` is the basis with only conserved quantities generated in the last step ;
* `dim :: Int64` is the dimension of the basis ;
* `szz :: Int64` records the maximum size ``\\max m_g`` of groups;
* `cfgr :: Vector{Int64}` is a vector of length `cfs.ncf` and records which group ``|I\\rangle`` each configuration ``|i\\rangle`` belong to ;
* `cffac :: Vector{ComplexF64}` is a vector of length `cfs.ncf` and records the coefficients ``\\lambda_i`` ;
* `grel :: Matrix{Int64}` is a `szz```\\times```dim` matrix that records the configurations in each group ``|i_{I1}\\rangle,\\dots,|i_{Im_I}\\rangle``
* `grsz :: Vector{Int64}` is a vector of length `dim` that records the size ``m_I`` of each group.
"""
mutable struct Basis
    cfs :: Confs
    dim :: Int64
    szz :: Int64
    cfgr :: Vector{Int64}
    cffac :: Vector{ComplexF64}
    grel :: Array{Int64, 2}
    grsz :: Vector{Int64}
end 

"""
    function Basis(cfs :: Confs, qnz_s :: Vector{ComplexF64}, cyc :: Vector{Int64}, perm_o :: Vector{Vector{Int64}}, ph_o :: Vector{Vector{Int64}}, fac_o :: Vector{Vector{ComplexF64}}) :: Basis

# Arguments 

* `cfs :: Confs` is the configuration set with only conserved quantities generated in the last step ;
* `qnz_s :: Vector{ComplexF64}` is a vector of length the same as the number of discrete symmetries ``N_Z`` that records the eigenvalue of each transformation ;
* `cyc :: Vector{Int64}` records the cycle of each transformation. For ``\\mathbb{Z}_n`` symmetry, record ``n`` ;
* `perm_o :: Vector{Vector{Int64}}` records the permutation ``\\pi_o``. It has ``N_Z`` elements and each of its elements is a vector of length ``N_o``. 
* `ph_o :: Vector{Vector{Int64}}` records ``p_o`` to determine whether or not to perform a particle-hole transformation. It has ``N_Z`` elements and each of its elements is a vector of length ``N_o``. 
* `fac_o :: Vector{Vector{ComplexF64}}` records the factor ``p_o`` is determine whether or not to perform a particle-hole transformation. It has ``N_Z`` elements and each of its elements is a vector of length ``N_o``. 
"""
function Basis(cfs :: Confs, qnz_s :: Vector{ComplexF64}, cyc :: Vector{Int64}, perm_o :: Vector{Any}, ph_o :: Vector{Any}, fac_o :: Vector{Any})
    nqnz = length(qnz_s)
    perm_o_mat = reduce(hcat, perm_o)
    ph_o_mat = reduce(hcat, ph_o)
    fac_o_mat = reduce(hcat, fac_o)
    dim_ref = Ref{Int64}(0)
    cfgr = Array{Int64, 1}(undef, cfs.ncf)
    cffac = Array{ComplexF64, 1}(undef, cfs.ncf)
    szz = prod([ abs(qnz_s[i]) < 1E-8 ? 1 : cyc[i] for i = 1 : nqnz ])
    @ccall LibpathFuzzifiED.bs_mp_generate_bs_cfgr_(cfs.no :: Ref{Int64}, cfs.nor :: Ref{Int64}, cfs.ncf :: Ref{Int64}, cfs.lid :: Ref{Int64}, cfs.rid :: Ref{Int64}, cfs.conf :: Ref{Int64}, nqnz :: Ref{Int64}, qnz_s :: Ref{ComplexF64}, cyc :: Ref{Int64}, perm_o_mat :: Ref{Int64}, ph_o_mat :: Ref{Int64}, fac_o_mat :: Ref{ComplexF64}, szz :: Ref{Int64}, dim_ref :: Ref{Int64}, cfgr :: Ref{Int64}, cffac :: Ref{ComplexF64}) :: Nothing
    dim = dim_ref[]
    grel = Array{Int64, 2}(undef, szz, dim)
    grsz = Array{Int64, 1}(undef, dim)
    @ccall LibpathFuzzifiED.bs_mp_generate_bs_grel_(cfs.ncf :: Ref{Int64}, szz :: Ref{Int64}, dim :: Ref{Int64}, cfgr :: Ref{Int64}, grel :: Ref{Int64}, grsz :: Ref{Int64}) :: Nothing
    return Basis(cfs, dim, szz, cfgr, cffac, grel, grsz)
end 

"""
    function Basis(cfs :: Confs) :: Basis

Generate a basis from the configurations without applying the ``\\mathbb{Z}_2`` symmetries
"""
function Basis(cfs :: Confs)
    nqnz = 1
    dim = cfs.ncf
    szz = 1 
    cfgr = collect(1 : dim)
    cffac = fill(ComplexF64(1), dim)
    grel = reshape(collect(1 : dim), 1, :)
    grsz = fill(1, dim)
    return Basis(cfs, dim, szz, cfgr, cffac, grel, grsz)
end 


"""
    GetConfWeight(bs :: Basis, st :: Vector{ComplexF64}, cf :: Int64) :: ComplexF64
    GetConfWeight(bs :: Basis, st :: Vector{Float}, cf :: Int64) :: ComplexF64

looks up a the weight of a configuration in a state. 

# Arguments 

* ``bs :: Basis`` is the basis of the state ; 
* ``st :: Vector{ComplexF64}`` or ``st :: Vector{Float}`` is a vector of length `bs.dim` that stores the state ; 
* ``cf :: Int64`` stores the configuration to be looked-up expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th orbital in the `i`-th configuration is occupied ; if the bit is 0, then the orbital is empty. 

# Output

* The weight of the configuration in the state 

"""
function GetConfWeight(bs :: Basis, st :: Vector, cf :: Int64)
    id = GetConfId(bs.cfs, cf)
    return st[bs.cfgr[id]] * bs.cffac[id]
end