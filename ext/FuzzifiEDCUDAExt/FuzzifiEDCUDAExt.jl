module FuzzifiEDCUDAExt

using CUDA
using CUDA.CUSPARSE
using SparseArrays
using FuzzifiED
using KrylovKit

import FuzzifiED.GetEigensystemCuda
import FuzzifiED.SilentStd

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
    GetEigensystemCuda(mat :: OpMat{ComplexF64}, nst :: Int64 ; initvec :: Vector{ComplexF64}, num_th :: Int64, disp_std :: Bool, kwargs...) :: Tuple{Vector{ComplexF64}, CuArray{ComplexF64, 2, CUDA.DeviceMemory}}
    GetEigensystemCuda(mat :: OpMat{Float64}, nst :: Int64 ; initvec :: Vector{Float64}, num_th :: Int64, disp_std :: Bool, kwargs...) :: Tuple{Vector{Float64}, CuArray{Float64, 2, CUDA.DeviceMemory}}

This method uses Julia `KrylovKit` package to calculate the lowest eigenstates of sparse matrix. The sparse matrix multiplication is realised by `CUDA.CUSPARSE`. For an example, refer to [`ising_spectrum_cuda.jl`](https://github.com/FuzzifiED/FuzzifiED.jl/blob/main/examples/ising_spectrum_cuda.jl).

# Arguments 

* `mat :: OpMat{ComplexF64}` or `mat :: OpMat{Float64}` is the matrix ;
* `nst :: Int64` is the number of eigenstates to be calculated ;
* `tol :: Float64` is the tolerence for the KrylovKit process. The default value is `1E-8` ;
* `ncv :: Int64` is the maximum dimension of the Krylov subspace. The default value is `max(2 * nst, nst + 10)`. If `krylovdim` is also given, `ncv` will not be used.
* `initvec` is the initial vector. Facultative, a random initialisation `CUDA.rand(T, mat.dimd)` by default ;
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 
* `kwargs...` is the options that will directly sent into `eigsolve`, see its [documentation](https://jutho.github.io/KrylovKit.jl/stable/man/eig/#KrylovKit.eigsolve) for detail.

# Output

* A length-`nst` array that has the same type as `mat` recording the eigenvalues, and 
* A `dimd`\\*`nst` matrix that has the same type as `mat` where every column records an eigenstate. 
"""
function GetEigensystemCuda(mat :: OpMat{T}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = max(2 * nst, nst + 10), initvec = CUDA.rand(T, mat.dimd), disp_std = !SilentStd, kwargs...) where T <: Union{Float64, ComplexF64}
    kwargs1 = haskey(kwargs, :krylovdim) ? kwargs : (kwargs..., krylovdim = ncv)
    mat_cu = CUSPARSE.CuSparseMatrixCSC(mat) 
    eigval, eigvec, info = eigsolve(x -> mat_cu * x, initvec, nst, :SR ; tol, kwargs1...)
    if (disp_std) print(info) end
    return Vector{T}(eigval), CuArray{T, 2, CUDA.DeviceMemory}(hcat(eigvec...))
end

end 