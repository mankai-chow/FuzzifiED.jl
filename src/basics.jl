"""
    Confs

stores all the configurations that respects the U(1) conserved quantities. 

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
    Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Vector{Int64}} ; nor :: Int64) :: Confs

# Arguments

* `no :: Int64` is the number of orbitals ``N_o`` ;
* `qnu_s :: Vector{Int64}` is the set of ``Q_i`` for the selected configurations ;
* `qnu_o :: Vector{Vector{Int64}}` is the set of ``q_{io}`` for each quantum number and for each orbital. It should contain ``N_U`` elements and each element should be a vector of length ``N_o``. 
* (`nor :: Int64` is the number of less significant bits used to generate the Lin table.)
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
    Basis

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
    Basis(cfs :: Confs, qnz_s :: Vector{ComplexF64}, cyc :: Vector{Int64}, perm_o :: Vector{Vector{Int64}}, ph_o :: Vector{Vector{Int64}}, fac_o :: Vector{Vector{ComplexF64}}) :: Basis

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
    Basis(cfs :: Confs) :: Basis

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
    Operator
"""
mutable struct Operator 
    bsd :: Basis 
    bsf :: Basis
    red_q :: Int64
    sym_q :: Int64
    ntm :: Int64
    nc :: Int64
    cstr :: Array{Int64, 2}
    fac :: Vector{ComplexF64}
end 

"""
    Operator(bsd :: Basis, bsf :: Basis, cstr_vec :: Vector{Vector{Int64}}, fac :: Vector{ComplexF64} ; red_q :: Int64, sym_q :: Int64) :: Operator

# Arguments
* `bsd :: Basis` is the basis of the initial state ;
* `bsf :: Basis` is the basis of the final state ;
* `cstr_vec :: Vector{Vector{Integer}}` records the ``c`` and ``c^\\dagger`` string of each term. A term with ``l`` operators ``c^{(p_{t1})}_{o_{t1}}c^{(p_{t2})}_{o_{t2}}\\dots c^{(p_{tl})}_{o_{tl}}`` correspond to a length-``2l`` vector ``(p_{t1},o_{t1},p_{t2},o_{t2},\\dots p_{tl},o_{tl})``. Note that each element can have different length. For constant term, input `[-1, -1]`; 
* `fac :: Vector{ComplexF64}` corresponds to the factor ``U_t`` in each term.
* `red_q :: Int64` is a flag that records whether or not the conversion to a sparse martrix can be simplified : if `bsd` and `bsf` have exactly the same quantum numbers, and the operator fully respects the symmetries, and all the elements in `bsd.cffac` and `bsf.cffac` has the same absolute value, then `red_q = 1` ; otherwise `red_q = 0` ; 
* `sym_q :: Int64` records the symmetry of the operator : if the matrix is Hermitian, then `sym_q = 1` ; if it is symmetric, then `sym_q = 2` ; otherwise `sym_q = 0`. 
"""
function Operator(bsd :: Basis, bsf :: Basis, cstr_vec :: Vector{Any}, fac :: Vector{ComplexF64} ; red_q :: Int64 = 0, sym_q :: Int64 = 0)
    ntm = length(cstr_vec)
    nc = div(maximum(length.(cstr_vec)), 2)
    cstr_vec_eq = [ [tm ; fill(-1, 2 * nc - length(tm))] for tm in cstr_vec]
    cstr = reduce(hcat, cstr_vec_eq)
    return Operator(bsd, bsf, red_q, sym_q, ntm, nc, cstr, fac)
end

"""
    OpMat

# Fields
* `dimd :: Int64` and `dimf :: Int64` are the number of columns and rows of the matrix ;
* `symq :: Int64` records whether or not the matrix is Hermitian or symmetric ;
* `nel :: Int64` records the number of elements ;
* `colptr :: Vector{Int64}`, `rowid :: Vector{Int64}` and `elval :: Vector{ComplexF64}` records the elements of the sparse matrix as in the `SparseMatrixCSC` elements of Julia. 
"""
mutable struct OpMat
    dimd :: Int64
    dimf :: Int64
    sym_q :: Int64
    nel :: Int64 
    colptr :: Vector{Int64}
    rowid :: Vector{Int64}
    elval :: Vector{ComplexF64}
end 

"""
    OpMat(op :: Operator) :: OpMat
"""
function OpMat(op :: Operator)
    colptr = Array{Int64, 1}(undef, op.bsd.dim + 1)
    nel_ref = Ref{Int64}(0)
    @ccall LibpathFuzzifiED.op_mp_count_op_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
        op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
        op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
        op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstr :: Ref{Int64}, op.fac :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}) :: Nothing
    nel = nel_ref[]
    rowid = Array{Int64, 1}(undef, nel)
    elval = Array{ComplexF64, 1}(undef, nel)
    @ccall LibpathFuzzifiED.op_mp_generate_op_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
        op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
        op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
        op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstr :: Ref{Int64}, op.fac :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, rowid :: Ref{Int64}, elval :: Ref{ComplexF64}) :: Nothing
    return OpMat(op.bsd.dim, op.bsf.dim, op.sym_q, nel, colptr, rowid, elval)
end

"""
    GetEigensystem(mat :: OpMat, nst :: Int64 ; tol :: Float64 = 1E-8)

# Output
* A length-`nst` array recording the eigenvalues, and 
* A `dimd```\\times```nst` matrix where every column records an eigenstate. 
"""
function GetEigensystem(mat :: OpMat, nst :: Int64 ; tol :: Float64 = 1E-8)
    eigval = Array{ComplexF64, 1}(undef, nst + 1)
    eigvec = Array{ComplexF64, 2}(undef, mat.dimd, nst)
    @ccall LibpathFuzzifiED.diag_mp_diagonalisation_(mat.dimd :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, nst :: Ref{Int64}, tol :: Ref{Float64}, eigval :: Ref{ComplexF64}, eigvec :: Ref{ComplexF64}) :: Nothing
    return eigval[1 : end - 1], eigvec
end 

"""
    *(op :: Operator, st_d :: Vector{ComplexF64})

Measure the action of an operator on a state. `st_d` must be of length `op.bsd.dim`. Returns a vector of length `op.bsf.dim` that represents the final state. 

Note that sometimes it is needed to transform a state from one basis to another. This can be done by constructing an identity operator. 
```julia
stf = Operator(bsd, bsf, [[-1, -1]], [ComplexF64(1)]) * std
```
"""
function *(op :: Operator, st_d :: Vector{ComplexF64})
    st_f = Array{ComplexF64, 1}(undef, op.bsf.dim)
    @ccall LibpathFuzzifiED.op_mp_action_op_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
    op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
    op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstr :: Ref{Int64}, op.fac :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, st_d :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}) :: Nothing
    return st_f
end

"""
    *(mat :: OpMat, st_d :: Vector{ComplexF64})

Measure the action of a sparse matrix on a state. `st_d` must be of length `mat.dimd`. Returns a vector of length `mat.dimf` that represents the final state. 
"""
function *(mat :: OpMat, st_d :: Vector{ComplexF64})
    st_f = Array{ComplexF64, 1}(undef, mat.dimf)
    @ccall LibpathFuzzifiED.diag_mp_vec_prod_(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, st_d :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}) :: Nothing
    return st_f
end

"""
    *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, op :: Operator, st_d :: Vector{ComplexF64}) :: ComplexF64

Measuring the inner product between two states and an operator. `st_d` must be of length `op.bsd.dim` and `st_fp` must be of length `op.bsf.dim`, and `st_fp` must be an adjoint. 
"""
function *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, op :: Operator, st_d :: Vector{ComplexF64})
    ovl_ref = Ref{ComplexF64}(0)
    @ccall LibpathFuzzifiED.op_mp_overlap_op_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
    op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
    op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstr :: Ref{Int64}, op.fac :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, st_d :: Ref{ComplexF64}, st_fp' :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}) :: Nothing 
    return ovl_ref[]
end

"""
    *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat, st_d :: Vector{ComplexF64}) :: ComplexF64

Measuring the inner product between two states and a sparse matrix. `st_d` must be of length `mat.dimd` and `st_fp` must be of length `mat.dimf`, and `st_fp` must be an adjoint. 
"""
function *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat, st_d :: Vector{ComplexF64})
    ovl_ref = Ref{ComplexF64}(0)
    @ccall LibpathFuzzifiED.diag_mp_scal_prod_(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, st_d :: Ref{ComplexF64}, st_fp' :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}) :: Nothing
    return ovl_ref[]
end