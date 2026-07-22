export SCompOperator, BuildSCompOperator
import FuzzifiED: GetEigensystem


"""
    SCompOperator{Float64}
    SCompOperator{ComplexF64}

The mutable type `SCompOperator` represents a composite operator — such as the Hamiltonian — acting on a [SCompSpace](@ref SCompSpace) of definite total angular momentum. It combines the reduced matrix elements of the per-part [SSegOperators](@ref SSegOperator) with the ``9j`` recoupling coefficients that relate the coupled basis of the initial and final composite spaces. The operator is stored in a block-structured form. It can be multiplied formally to a state.

# Fields

* `cpspd :: SCompSpace` and `cpspf :: SCompSpace` are the initial and final composite spaces.
* `nd :: Int64` is the number of decomposition channels of the operator.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}`` carried by the operator (``0`` for a scalar such as the Hamiltonian).
* `coeff :: Vector{T}` is the coefficient of each decomposition channel.
* `sgop :: Matrix{SSegOperator}` is the ``N_p×N_d`` matrix of segment operators produced by [BuildSSegOperators](@ref BuildSSegOperators).
* `colptr :: Matrix{Int64}` and `rowid :: Vector{Vector{Int64}}` store the allowed blocks of sectors for each channel. `colptr[:, d]` and `rowid[d]` bear the format of a CSC sparse matrix.
* `idel :: Vector{Matrix{Int64}}` locates the matrix element block in the SSegOperators for each segment. It takes three indices `idel[d][p, e]` where `d` is the channel index, `e` is the element index, and `p` is the part index.
* `mat9j :: Vector{Vector{Matrix{Float64}}}` stores the precomputed ``9j`` recoupling coefficients, including the fermionic reordering sign. It takes four indices `mat9j[d][e][ich, jch]` where `ich`, `jch` is the channel indices within the sectors specified by `e`.
"""
mutable struct SCompOperator{T <: Union{Float64, ComplexF64}}
    cpspd :: SCompSpace
    cpspf :: SCompSpace
    nd :: Int64
    ltot :: Int64
    coeff :: Vector{T}
    sgop :: Matrix{SSegOperator}
    colptr :: Matrix{Int64}
    rowid :: Vector{Vector{Int64}}
    idel :: Vector{Matrix{Int64}}
    mat9j :: Vector{Vector{Matrix{Float64}}}
end


"""
    BuildSCompOperator(cpspd :: SCompSpace{T}[, cpspf :: SCompSpace{T}], cpd :: Vector{SCoupleDecomp}, sgop :: Matrix{SSegOperator}, ltot :: Int64) :: SCompOperator

constructs a [SCompOperator](@ref SCompOperator) from the composite spaces, the coupling decompositions `cpd` and the segment operators `sgop`. It computes and stores the ``9j`` recoupling coefficient between every pair of initial and final coupling channels and every decomposition channel together with the sign arising from fermion parity.

# Arguments

* `cpspd :: SCompSpace{T}` is the initial composite space.
* `cpspf :: SCompSpace{T}` is the final composite space. Facultative, the same as `cpspd` by default.
* `cpd :: Vector{SCoupleDecomp}` is the list of coupling decompositions ; it must be the same one used to build `sgop`.
* `sgop :: Matrix{SSegOperator}` is the matrix of segment operators from [BuildSSegOperators](@ref BuildSSegOperators).
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}`` carried by the operator. Facultative, ``0`` (a scalar) by default.

# Output

* `cpop :: SCompOperator` is the resulting composite operator.
"""
function BuildSCompOperator(cpspd :: SCompSpace{T}, cpspf :: SCompSpace{T}, cpd :: Vector{SCoupleDecomp}, sgop :: Matrix{SSegOperator}, ltot :: Int64 = 0) where T <: Union{Float64, ComplexF64}
    id_tc = vcat([ fill(i, length(cpd[i].ch)) for i in eachindex(cpd)]...)
    ch = vcat([ cpd[i].ch for i in eachindex(cpd) ]...)
    coeff = vcat([ cpd[i].coeff for i in eachindex(cpd) ]...)
    nd = length(id_tc)
    np = cpspd.np

    colptr = zeros(Int64, size(cpspd.idsec, 2) + 1, nd)
    colptr[1, :] .= 1
    rowid = [ Int64[] for d = 1 : nd ]
    idel = [ Matrix{Int64}(undef, np, 0) for d = 1 : nd]
    Threads.@threads :greedy for d = 1 : nd
        index = 0
        for jsec in axes(cpspd.idsec, 2)
            idsecj = cpspd.idsec[:, jsec]
            idel_rng = [ sgop[p, d].colptr[idsecj[p]] : sgop[p, d].colptr[idsecj[p] + 1] - 1 for p = 1 : np ]
            for ideli in Iterators.product(idel_rng...)
                idseci = [ sgop[p, d].rowid[ideli[p]] for p = 1 : np]
                isec_rng = searchsorted(eachcol(cpspf.idsec), idseci)
                isempty(isec_rng) && continue
                isec = isec_rng[1]

                index += 1
                push!(rowid[d], isec)
                idel[d] = hcat(idel[d], collect(ideli))
            end
            colptr[jsec + 1, d] = index + 1
        end
    end

    mat9j = [ Vector{Matrix{Float64}}(undef, colptr[end, d] - 1) for d = 1 : nd ]
    Threads.@threads :greedy for (jsec, d) in collect(Iterators.product(axes(cpspd.idsec, 2), 1 : nd))
        idsecj = cpspd.idsec[:, jsec]
        for e = colptr[jsec, d] : colptr[jsec + 1, d] - 1
            isec = rowid[d][e]
            pfh = [ mod(cpd[id_tc[d]].sec[1, p], 2) for p = 1 : np ]
            pfj = [ mod(cpspd.sgsp[p].sec[1, idsecj[p]], 2) for p = 1 : np]
            pftot = sum([pfj[p] * sum(pfh[p + 1 : end]) for p = 1 : np]) % 2
            chh = ch[d]
            mat = Matrix{Float64}(undef, length(cpspf.chs[isec]), length(cpspd.chs[jsec]))
            for j in eachindex(cpspd.chs[jsec])
                chj = cpspd.chs[jsec][j]
                for i in eachindex(cpspf.chs[isec])
                    chi = cpspf.chs[isec][i]
                    fac = 1
                    for p = 2 : np
                        fac = fac * float(d9j(chj[2,p-1], chj[1,p], chj[2,p],  chh[2,p-1], chh[1,p], chh[2,p],  chi[2,p-1], chi[1,p], chi[2,p])) * √((chj[2,p]+1) * (chh[2,p]+1) * (chi[2,p]+1))
                    end
                    (pftot == 1) && (fac = -fac)
                    mat[i, j] = fac
                end
            end
            mat9j[d][e] = mat
        end
    end
    @info "FINISH GENERATING COMP OPERATOR"

    return SCompOperator{T}(cpspd, cpspf, nd, ltot, coeff, sgop, colptr, rowid, idel, mat9j)
end
BuildSCompOperator(cpspd :: SCompSpace{T}, cpd :: Vector{SCoupleDecomp}, sgop :: Matrix{SSegOperator}, ltot :: Int64 = 0) where T <: Union{Float64, ComplexF64} = BuildSCompOperator(cpspd, cpspd, cpd, sgop, ltot)


"""
    *(cpop :: SCompOperator{T}, std :: Vector{T}) :: Vector{T}

applies the composite operator `cpop` to a state `std` of the initial composite space and returns the resulting state of the final composite space. The action is evaluated block by block : for every decomposition channel and every pair of coupling channels it takes the Kronecker product of the corresponding per-part reduced matrix element blocks, weighted by the channel coefficient and the ``9j`` recoupling factor.
"""
function Base.:*(cpop :: SCompOperator{T}, std :: Vector{T}) where T <: Union{Float64, ComplexF64}
    th_lock = ReentrantLock()
    stf = zeros(T, cpop.cpspf.dim)
    np = cpop.cpspd.np
    Threads.@threads :greedy for (jsec, d) in collect(Iterators.product(axes(cpop.cpspd.idsec, 2), 1 : cpop.nd))
        stf1 = zeros(T, cpop.cpspf.dim)
        idsecj = cpop.cpspd.idsec[:, jsec]
        coeff = cpop.coeff[d]

        jrng_sg = Vector{UnitRange{Int64}}(undef, np)
        irng_sg = Vector{UnitRange{Int64}}(undef, np)
        for e = cpop.colptr[jsec, d] : cpop.colptr[jsec + 1, d] - 1
            idel_sg = cpop.idel[d][:, e]
            isec = cpop.rowid[d][e]
            idseci = cpop.cpspf.idsec[:, isec]
            for jch in eachindex(cpop.cpspd.chs[jsec])
                jrng = cpop.cpspd.ptr_st[jsec][jch] : cpop.cpspd.ptr_st[jsec][jch + 1] - 1
                for p = 1 : np
                    lj = cpop.cpspd.chs[jsec][jch][1, p]
                    idlj = cpop.cpspd.sgsp[p].l_lookup[idsecj[p]][lj]
                    jrng_sg[p] = (cpop.cpspd.sgsp[p].ptr_st[idsecj[p]][idlj] + 1 : cpop.cpspd.sgsp[p].ptr_st[idsecj[p]][idlj + 1]) .- cpop.cpspd.sgsp[p].ptr_st[idsecj[p]][1]
                end
                for ich in eachindex(cpop.cpspf.chs[isec])
                    fac9j = cpop.mat9j[d][e][ich, jch]
                    abs(fac9j) < √eps(Float64) && continue
                    irng = cpop.cpspf.ptr_st[isec][ich] : cpop.cpspf.ptr_st[isec][ich + 1] - 1
                    for p = 1 : np
                        li = cpop.cpspf.chs[isec][ich][1, p]
                        idli = cpop.cpspf.sgsp[p].l_lookup[idseci[p]][li]
                        irng_sg[p] = (cpop.cpspf.sgsp[p].ptr_st[idseci[p]][idli] + 1 : cpop.cpspf.sgsp[p].ptr_st[idseci[p]][idli + 1]) .- cpop.cpspf.sgsp[p].ptr_st[idseci[p]][1]
                    end
                    @views stf1[irng] .+= (coeff * fac9j) .* (⊗([cpop.sgop[p, d].elmat[idel_sg[p]][irng_sg[p], jrng_sg[p]] for p = 1 : np]...) * std[jrng])
                end
            end
        end
        lock(th_lock) do
            stf .+= stf1
        end
    end
    return stf
end


"""
    GetEigensystem(cpop :: SCompOperator{T}, nst :: Int64 ; tol :: Float64, ncv :: Int64, initvec :: Vector{T}, kwargs...) :: Tuple{Vector{T}, Matrix{T}}

computes the lowest `nst` eigenvalues and eigenvectors of the composite operator `cpop` through `KrylovKit.eigsolve`. This yields the spectrum resolved by angular momentum and flavour symmetries within segments.

# Arguments

* `cpop :: SCompOperator{T}` is the composite operator (usually the Hamiltonian).
* `nst :: Int64` is the number of eigenpairs to compute.
* `tol :: Float64` is the tolerance of the eigensolver. Facultative, `1E-8` by default.
* `ncv :: Int64` is the dimension of the Krylov subspace. Facultative, `max(2 * nst, nst + 10)` by default.
* `initvec :: Vector{T}` is the initial vector. Facultative, a random vector by default.
* `kwargs...` are further keyword arguments forwarded to `eigsolve`.

# Output

* `eigval :: Vector{T}` is the vector of the `nst` lowest eigenvalues.
* `eigvec :: Matrix{T}` is the matrix whose columns are the corresponding eigenvectors.
"""
function FuzzifiED.GetEigensystem(cpop :: SCompOperator{T}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = max(2 * nst, nst + 10), initvec = rand(T, cpop.cpspd.dim), kwargs...) where T <: Union{ComplexF64,Float64}
    eigval, eigvec, info = eigsolve(x -> cpop * x, initvec, nst, :SR ; tol, krylovdim = ncv, verbosity = 2, kwargs...)
    return Vector{T}(eigval), Matrix{T}(hcat(eigvec...))
end
