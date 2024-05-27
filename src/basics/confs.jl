"""
    mutable struct Confs

This type stores all the configurations that respects the conserved quantities and also a table to inversely look up the index from the configuration. 

# Fields

* `no :: Int64` is the number of orbitals
* `ncf :: Int64` is the number of configurations 
* `conf :: Vector{Int64}` is an array of length `ncf` containing all the configurations. Each configuration is expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th orbital in the `i`-th configuration is occupied ; if the bit is 0, then the orbital is empty. 
* `nor :: Int64`, `lid :: Vector{Int64}` and `rid :: Vector{Int64}` contain the information of Lin table that is used to inversely look up the index `i` from the configuration. 
"""
mutable struct Confs
    no :: Int64 
    nor :: Int64
    ncf :: Int64
    conf :: Array{Int64,1} 
    lid :: Array{Int64,1}
    rid :: Array{Int64,1}
end 


"""
    function Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Vector{Int64}} ; nor :: Int64) :: Confs

generates the configurations that has the quantum numbers given by `qnu_s` of certain conserved quantities specified by `qnu_o :: Vector{Vector{Int64}}`

```math
Q_i=\\sum_{o=1}^{N_o}q_{io}n_o
```

where ``i=1,\\dots,N_U`` is the index of conserved quantities, ``o`` is the index of orbital, ``n_o=c^\\dagger_oc_o``, and ``q_o`` is a set of coefficients that must be non negative integer valued. 

# Arguments

* `no :: Int64` is the number of orbitals ``N_o`` ;
* `qnu_s :: Vector{Int64}` is the set of ``Q_i`` for the selected configurations ;
* `qnu_o :: Vector{Vector{Int64}}` is the set of ``q_{io}`` for each quantum number and for each orbital. It should contain ``N_U`` elements and each element should be a vector of length ``N_o``. 
* (`nor :: Int64` is the number of less significant bits used to generate the Lin table.)

# Output

* `cfs :: Confs` is a [`Confs`](@ref) object
"""
function Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2))
    # qnu_o :: Vector{Vector{Int64}}
    nqnu = length(qnu_s)
    lid = Array{Int64, 1}(undef, 2 ^ (no - nor) + 1)
    ref_ncf = Ref{Int64}(0)
    qnu_o_mat = reduce(hcat, qnu_o)
    @ccall LibpathFuzzifiED.cfs_mp_count_cfs_(no :: Ref{Int64}, nor :: Ref{Int64}, nqnu :: Ref{Int64}, qnu_s :: Ref{Int64}, qnu_o_mat :: Ref{Int64}, ref_ncf :: Ref{Int64}, lid :: Ref{Int64}) :: Nothing
    ncf = ref_ncf[]
    rid = Array{Int64, 1}(undef, 2 ^ nor + 1)
    conf = Array{Int64, 1}(undef, ncf)
    @ccall LibpathFuzzifiED.cfs_mp_generate_cfs_(no :: Ref{Int64}, nor :: Ref{Int64}, nqnu :: Ref{Int64}, qnu_s :: Ref{Int64}, qnu_o_mat :: Ref{Int64}, ncf :: Ref{Int64}, lid :: Ref{Int64}, rid :: Ref{Int64}, conf :: Ref{Int64}) :: Nothing
    return Confs(no, nor, ncf, conf, lid, rid)
end 

"""
    GetConfId(cfs :: Confs, cf :: Int64) :: Int64

inversely look up the index from the configuration

# Arguments

* `cfs :: Confs` stores the configurations.
* `cf :: Int64` stores the configuration to be looked-up expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th orbital in the `i`-th configuration is occupied ; if the bit is 0, then the orbital is empty. 

# Output
* `id :: Int64` is the id of the configuration such that `cfs.conf[id] == cf`

"""
function GetConfId(cfs :: Confs, cf :: Int64)
    cf_l = cf >> cfs.nor
    cf_r = cf & (1 << cfs.nor - 1)
    return cfs.lid[cf_l + 1] + cfs.rid[cf_r + 1]
end 