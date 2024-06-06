"""
    mutable struct OpMat{Complex64}
    mutable struct OpMat{Float64}

This type stores a sparse matrix in the same form as `SparseMatrixCSC` in `SparseArrays`. If the matrix is Hermitian or symmetric, only the lower triangle is stored. 

# Fields
* `dimd :: Int64` and `dimf :: Int64` are the number of columns and rows of the matrix ;
* `symq :: Int64` records whether or not the matrix is Hermitian or symmetric ;
* `nel :: Int64` records the number of elements ;
* `colptr :: Vector{Int64}`, `rowid :: Vector{Int64}` and `elval :: Vector{ComplexF64}` records the elements of the sparse matrix as in the `SparseMatrixCSC` elements of Julia. 
"""
mutable struct OpMat{T}
    dimd :: Int64
    dimf :: Int64
    sym_q :: Int64
    nel :: Int64 
    colptr :: Vector{Int64}
    rowid :: Vector{Int64}
    elval :: Vector{T}
end 

"""
    function OpMat(op :: Operator ; type :: DataType = ComplexF64) :: OpMat{type}

Generates the sparse matrix from the operator

# Arguments 

* `op :: Operator` is the operator ;
* `type :: DataType` specifies the type of the matrix. It can either be `ComplexF64` or `Float64`. 
"""
function OpMat(op :: Operator ; type :: DataType = ComplexF64)
    colptr = Array{Int64, 1}(undef, op.bsd.dim + 1)
    nel_ref = Ref{Int64}(0)
    @ctrlstd @ccall FuzzifiED_jll.LibpathFuzzifiED.__op_MOD_count_op(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
        op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
        op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
        op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}) :: Nothing
    nel = nel_ref[]
    rowid = Array{Int64, 1}(undef, nel)
    elval = Array{type, 1}(undef, nel)
    if (type == ComplexF64) 
        @ctrlstd @ccall FuzzifiED_jll.LibpathFuzzifiED.__op_MOD_generate_op(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
            op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
            op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
            op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, rowid :: Ref{Int64}, elval :: Ref{ComplexF64}) :: Nothing
    elseif (type == Float64) 
        @ctrlstd @ccall FuzzifiED_jll.LibpathFuzzifiED.__op_MOD_generate_op_re(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
        op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
        op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
        op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, rowid :: Ref{Int64}, elval :: Ref{Float64}) :: Nothing
    end
    return OpMat{type}(op.bsd.dim, op.bsf.dim, op.sym_q, nel, colptr, rowid, elval)
end


"""
    function GetEigensystem(mat :: OpMat{ComplexF64}, nst :: Int64 ; tol :: Float64, ncv :: Int64) :: Tuple{Vector{ComplexF64}, Matrix{ComplexF64}}
    function GetEigensystem(mat :: OpMat{Float64}, nst :: Int64 ; tol :: Float64, ncv :: Int64) :: Tuple{Vector{Float64}, Matrix{Float64}}

calls the Arpack package to calculate the lowest eigenstates of sparse matrix. 

# Arguments 

* `mat :: OpMat{ComplexF64}` is the matrix ;
* `nst :: Int64` is the number of eigenstates to be calculated ;
* `tol :: Float64` is the tolerence for the Arpack process. The default value is `1E-8` ;
* `ncv :: Int64` is an euxiliary parameter needed in the Arpack process. The default value is `max(2 * nst, nst + 10)`

# Output

* A length-`nst` array that has the same type as `mat` recording the eigenvalues, and 
* A `dimd`\\*`nst` matrix that has the same type as `mat` where every column records an eigenstate. 
"""
function GetEigensystem(mat :: OpMat{ComplexF64}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = max(2 * nst, nst + 10))
    eigval = Array{ComplexF64, 1}(undef, nst + 1)
    eigvec = Array{ComplexF64, 2}(undef, mat.dimd, nst)
    @ctrlstd @ccall FuzzifiED_jll.LibpathFuzzifiED.__diag_MOD_diagonalisation(mat.dimd :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, nst :: Ref{Int64}, tol :: Ref{Float64}, ncv :: Ref{Int64}, eigval :: Ref{ComplexF64}, eigvec :: Ref{ComplexF64}) :: Nothing
    return eigval[1 : end - 1], eigvec
end 
function GetEigensystem(mat :: OpMat{Float64}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = max(2 * nst, nst + 10))
    eigval = Array{Float64, 1}(undef, nst + 1)
    eigvec = Array{Float64, 2}(undef, mat.dimd, nst)
    @ctrlstd @ccall FuzzifiED_jll.LibpathFuzzifiED.__diag_re_MOD_diagonalisation_re(mat.dimd :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{Float64}, nst :: Ref{Int64}, tol :: Ref{Float64}, ncv :: Ref{Int64}, eigval :: Ref{Float64}, eigvec :: Ref{Float64}) :: Nothing
    return eigval[1 : end - 1], eigvec
end 

"""
    function *(mat :: OpMat{ComplexF64}, st_d :: Vector{ComplexF64}) :: Vector{ComplexF64}
    function *(mat :: OpMat{Float64}, st_d :: Vector{Float64}) :: Vector{Float64}

Measure the action of a sparse matrix on a state. `st_d` must be of length `mat.dimd`. Returns a vector of length `mat.dimf` that represents the final state. 
"""
function *(mat :: OpMat{ComplexF64}, st_d :: Vector{ComplexF64})
    st_f = Array{ComplexF64, 1}(undef, mat.dimf)
    @ctrlstd @ccall FuzzifiED_jll.LibpathFuzzifiED.__diag_MOD_vec_prod(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, st_d :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}) :: Nothing
    return st_f
end
function *(mat :: OpMat{Float64}, st_d :: Vector{Float64})
    st_f = Array{Float64, 1}(undef, mat.dimf)
    @ctrlstd @ccall FuzzifiED_jll.LibpathFuzzifiED.__diag_re_MOD_vec_prod_re(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{Float64}, real.(st_d) :: Ref{Float64}, st_f :: Ref{Float64}) :: Nothing
    return st_f
end


"""
    function *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat{ComplexF64}, st_d :: Vector{ComplexF64}) :: ComplexF64
    function *(st_fp :: LinearAlgebra.Adjoint{Float64, Vector{Float64}}, mat :: OpMat{Float64}, st_d :: Vector{Float64}) :: Float64

Measuring the inner product between two states and a sparse matrix. `st_d` must be of length `mat.dimd` and `st_fp` must be of length `mat.dimf`, and `st_fp` must be an adjoint. 
"""
function *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat{ComplexF64}, st_d :: Vector{ComplexF64})
    ovl_ref = Ref{ComplexF64}(0)
    @ctrlstd @ccall FuzzifiED_jll.LibpathFuzzifiED.__diag_MOD_scal_prod(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, st_d :: Ref{ComplexF64}, st_fp' :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}) :: Nothing
    return ovl_ref[]
end
function *(st_fp :: LinearAlgebra.Adjoint{Float64, Vector{Float64}}, mat :: OpMat{Float64}, st_d :: Vector{Float64})
    ovl_ref = Ref{Float64}(0)
    @ctrlstd @ccall FuzzifiED_jll.LibpathFuzzifiED.__diag_re_MOD_scal_prod_re(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{Float64}, st_d :: Ref{Float64}, st_fp' :: Ref{Float64}, ovl_ref :: Ref{Float64}) :: Nothing
    return ovl_ref[]
end


"""
    SparseMatrixCSCFromOpMat(mat :: OpMat{ComplexF64}) :: SparseMatrixCSC{Int64,ComplexF64}
    SparseMatrixCSCFromOpMat(mat :: OpMat{Float64}) :: SparseMatrixCSC{Int64,Float64}

converts the `OpMat` objects to a `SparseMatrixCSC` object in the `SparseArrays` package.
"""
function SparseMatrixCSCFromOpMat(mat :: OpMat)
    matcsc1 = SparseMatrixCSC(mat.dimf, mat.dimd, mat.colptr, mat.rowid, mat.elval)
    if (mat.sym_q == 0) 
        return matcsc1
    elseif (mat.sym_q == 1)
        return matcsc1 + adjoint(matcsc1) - spdiagm(diag(matcsc1))
    elseif (mat.sym_q == 2)
        return matcsc1 + transpose(matcsc1) - spdiagm(diag(matcsc1))
    end
end


"""
    OpMat(matcsc :: SparseMatrixCSC{Int64,ComplexF64}) :: OpMat{ComplexF64}
    OpMat(matcsc :: SparseMatrixCSC{Int64,Float64}) :: OpMat{Float64}

converts the `SparseMatrixCSC` object in the `SparseArrays` package to an `OpMat` objects.
"""
function OpMat(matcsc :: SparseMatrixCSC)
    return OpMat{typeof(matcsc.nzval[1])}(matcsc.n, matcsc.m, 0, length(matcsc.rowval), matcsc.colptr, matcsc.rowval, matcsc.nzval)
end


"""
    MatrixFromOpMat(mat :: OpMat{ComplexF64}) :: SparseMatrixCSC{Int64,ComplexF64}
    MatrixFromOpMat(mat :: OpMat{Float64}) :: SparseMatrixCSC{Int64,Float64}

converts the `OpMat` objects to a full matrix.
"""
function MatrixFromOpMat(mat :: OpMat)
    return Matrix(SparseMatrixCSCFromOpMat(mat))
end