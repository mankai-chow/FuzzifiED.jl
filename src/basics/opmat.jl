"""
    mutable struct OpMat

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
    function OpMat(op :: Operator) :: OpMat
"""
function OpMat(op :: Operator)
    colptr = Array{Int64, 1}(undef, op.bsd.dim + 1)
    nel_ref = Ref{Int64}(0)
    @ccall LibpathFuzzifiED.op_mp_count_op_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
        op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
        op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
        op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}) :: Nothing
    nel = nel_ref[]
    rowid = Array{Int64, 1}(undef, nel)
    elval = Array{type, 1}(undef, nel)
    if (type == ComplexF64) 
        @ccall LibpathFuzzifiED.op_mp_generate_op_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
            op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
            op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
            op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, rowid :: Ref{Int64}, elval :: Ref{ComplexF64}) :: Nothing
    elseif (type == Float64) 
        @ccall LibpathFuzzifiED.op_mp_generate_op_re_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
        op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
        op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
        op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, rowid :: Ref{Int64}, elval :: Ref{Float64}) :: Nothing
    end
    return OpMat{type}(op.bsd.dim, op.bsf.dim, op.sym_q, nel, colptr, rowid, elval)
end


"""
    function GetEigensystem(mat :: OpMat, nst :: Int64 ; tol :: Float64 = 1E-8)

# Output
* A length-`nst` array recording the eigenvalues, and 
* A `dimd```\\times```nst` matrix where every column records an eigenstate. 
"""
function GetEigensystem(mat :: OpMat, nst :: Int64 ; tol :: Float64 = 1E-8)
    eigval = Array{ComplexF64, 1}(undef, nst + 1)
    eigvec = Array{ComplexF64, 2}(undef, mat.dimd, nst)
    @ccall LibpathFuzzifiED.diag_mp_diagonalisation_(mat.dimd :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, nst :: Ref{Int64}, tol :: Ref{Float64}, ncv :: Ref{Int64}, eigval :: Ref{ComplexF64}, eigvec :: Ref{ComplexF64}) :: Nothing
    return eigval[1 : end - 1], eigvec
end 
function GetEigensystem(mat :: OpMat{Float64}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = max(2 * nst, nst + 10))
    eigval = Array{Float64, 1}(undef, nst + 1)
    eigvec = Array{Float64, 2}(undef, mat.dimd, nst)
    @ccall LibpathFuzzifiED.diag_re_mp_diagonalisation_re_(mat.dimd :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{Float64}, nst :: Ref{Int64}, tol :: Ref{Float64}, ncv :: Ref{Int64}, eigval :: Ref{Float64}, eigvec :: Ref{Float64}) :: Nothing
    return eigval[1 : end - 1], eigvec
end 

"""
    function *(mat :: OpMat, st_d :: Vector{ComplexF64})

Measure the action of a sparse matrix on a state. `st_d` must be of length `mat.dimd`. Returns a vector of length `mat.dimf` that represents the final state. 
"""
function *(mat :: OpMat, st_d :: Vector{ComplexF64})
    st_f = Array{ComplexF64, 1}(undef, mat.dimf)
    @ccall LibpathFuzzifiED.diag_mp_vec_prod_(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, st_d :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}) :: Nothing
    return st_f
end
function *(mat :: OpMat{Float64}, st_d :: Vector{Float64})
    st_f = Array{Float64, 1}(undef, mat.dimf)
    @ccall LibpathFuzzifiED.diag_re_mp_vec_prod_re_(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{Float64}, real.(st_d) :: Ref{Float64}, st_f :: Ref{Float64}) :: Nothing
    return st_f
end


"""
    function *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat, st_d :: Vector{ComplexF64}) :: ComplexF64

Measuring the inner product between two states and a sparse matrix. `st_d` must be of length `mat.dimd` and `st_fp` must be of length `mat.dimf`, and `st_fp` must be an adjoint. 
"""
function *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat, st_d :: Vector{ComplexF64})
    ovl_ref = Ref{ComplexF64}(0)
    @ccall LibpathFuzzifiED.diag_mp_scal_prod_(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, st_d :: Ref{ComplexF64}, st_fp' :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}) :: Nothing
    return ovl_ref[]
end
function *(st_fp :: LinearAlgebra.Adjoint{Float64, Vector{Float64}}, mat :: OpMat{Float64}, st_d :: Vector{Float64})
    ovl_ref = Ref{Float64}(0)
    @ccall LibpathFuzzifiED.diag_re_mp_scal_prod_re_(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{Float64}, st_d :: Ref{Float64}, st_fp' :: Ref{Float64}, ovl_ref :: Ref{Float64}) :: Nothing
    return ovl_ref[]
end