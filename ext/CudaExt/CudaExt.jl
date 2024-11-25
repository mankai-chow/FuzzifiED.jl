module CudaExt

using CUDA
using CUDA.CUSPARSE
using SparseArrays
using FuzzifiED
using KrylovKit

import FuzzifiED:GetEigensystemCuda

"""
    CUSPARSE.CuSparseMatrixCSC(mat :: OpMat{ComplexF64})
    CUSPARSE.CuSparseMatrixCSC(mat :: OpMat{Float64})

converts the `OpMat` objects to a `CuSparseMatrixCSC` object in the `CUDA.CUSPARSE` package.
"""
function CUDA.CUSPARSE.CuSparseMatrixCSC(mat :: OpMat{T}) where T <: Union{Float64, ComplexF64}
    if (mat.sym_q ≠ 0) println("It works better for CUDA to take `sym_q = 0`.") end
    mat_csc = SparseMatrixCSC(mat)
    return CuSparseMatrixCSC(mat_csc)
end

"""
    function GetEigensystemCuda(mat :: OpMat{ComplexF64}, nst :: Int64 ; initvec :: Vector{ComplexF64}, num_th :: Int64, disp_std :: Bool, kwargs...) :: Tuple{Vector{ComplexF64}, Matrix{ComplexF64}}
    function GetEigensystemCuda(mat :: OpMat{Float64}, nst :: Int64 ; tol :: Float64, ncv :: Int64, initvec :: Vector{Float64}, num_th :: Int64, disp_std :: Bool, kwargs...) :: Tuple{Vector{Float64}, Matrix{Float64}}

This method uses Julia `KrylovKit` package to calculate the lowest eigenstates of sparse matrix. The sparse matrix multiplication is realised by `CUDA.CUSPARSE`..

# Arguments 

* `mat :: OpMat{ComplexF64}` or `mat :: OpMat{Float64}` is the matrix ;
* `nst :: Int64` is the number of eigenstates to be calculated ;
* `initvec` is the initial vector. Facultative, a random initialisation `CUDA.rand(T, mat.dimd)` by default ;
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 
* `kwargs...` is the options that will directly sent into `eigsolve`, this includes `tol :: Float64` and `krylovdim`, see its [documentation](https://jutho.github.io/KrylovKit.jl/stable/man/eig/#KrylovKit.eigsolve) for detail.

# Output

* A length-`nst` array that has the same type as `mat` recording the eigenvalues, and 
* A `dimd`\\*`nst` matrix that has the same type as `mat` where every column records an eigenstate. 
"""
function GetEigensystemCuda(mat :: OpMat{T}, nst :: Int64 ; initvec = CUDA.rand(T, mat.dimd), disp_std = !SilentStd, kwargs...) where T <: Union{Float64, ComplexF64}
    mat_cu = CUSPARSE.CuSparseMatrixCSC(mat) 
    enrg, eigs, info = eigsolve(x -> mat_cu * x, initvec, nst, :SR ; kwargs...)
    if (disp_std) print(info) end
    return enrg, eigs
end

end 