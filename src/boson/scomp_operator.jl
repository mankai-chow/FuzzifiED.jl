export SCompOperator, BuildSCompOperator
import FuzzifiED: GetEigensystem

"""
    SCompOperator{Float64}
    SCompOperator{ComplexF64}

The mutable type `SCompOperator` represents a composite operator — such as the Hamiltonian — acting on a [SCompSpace](@ref SCompSpace) of definite total angular momentum. It combines the reduced matrix elements of the per-part [SSegOperators](@ref SSegOperator) with the ``9j`` recoupling coefficients that relate the coupled basis of the initial and final composite spaces. The operator is never materialised as a dense matrix ; instead `*` applies it to a state on the fly, which is used by [GetEigensystem](@ref) to obtain the low-lying spectrum through a Krylov method.

# Fields

* `cpspd :: SCompSpace` and `cpspf :: SCompSpace` are the initial and final composite spaces.
* `nd :: Int64` is the number of decomposition channels of the operator.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}`` carried by the operator (``0`` for a scalar such as the Hamiltonian).
* `coeff :: Vector{T}` is the coefficient of each decomposition channel.
* `sgop :: Matrix{SSegOperator}` is the ``N_p×N_d`` matrix of segment operators produced by [BuildSSegOperators](@ref BuildSSegOperators).
* `mat9j :: Array{Float64, 3}` stores the precomputed ``9j`` recoupling coefficients, indexed by `[final channel, initial channel, decomposition]`, including the fermionic reordering sign.
"""
mutable struct SCompOperator{T <: Union{Float64, ComplexF64}}
    cpspd :: SCompSpace
    cpspf :: SCompSpace
    nd :: Int64
    ltot :: Int64
    coeff :: Vector{T}
    sgop :: Matrix{SSegOperator}
    mat9j :: Array{Float64, 3}
end

"""
    BuildSCompOperator(cpspd :: SCompSpace{T}, cpspf :: SCompSpace{T}, cpd :: Vector{SCoupleDecomp}, sgop :: Matrix{SSegOperator}, ltot :: Int64) :: SCompOperator
    BuildSCompOperator(cpspd :: SCompSpace{T}, cpd :: Vector{SCoupleDecomp}, sgop :: Matrix{SSegOperator}, ltot :: Int64) :: SCompOperator

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

    mat9j = Array{Float64}(undef, cpspf.nch, cpspd.nch, nd)
    Threads.@threads :greedy for (isec, jsec, d) in collect(Iterators.product(axes(cpspf.idsec, 2), axes(cpspd.idsec, 2), 1 : nd))
        chh = ch[d]
        pfh = mod.(cpd[id_tc[d]].sec[1, :], 2)
        idsecj = cpspd.idsec[:, jsec]
        pfj = [ mod(cpspd.sgsp[p].sec[1, idsecj[p]], 2) for p = 1 : np]
        pftot = sum([pfj[p] * sum(pfh[p + 1 : end]) for p = 1 : np]) % 2

        jst = cpspd.ptr_ch[jsec] - 1
        ist = cpspf.ptr_ch[isec] - 1
        for j in eachindex(cpspd.chs[jsec])
            chj = cpspd.chs[jsec][j]
            for i in eachindex(cpspf.chs[isec])
                chi = cpspf.chs[isec][i]
                fac = 1
                for p = 2 : np
                    fac = fac * float(d9j(chj[2,p-1], chj[1,p], chj[2,p],  chh[2,p-1], chh[1,p], chh[2,p],  chi[2,p-1], chi[1,p], chi[2,p])) * √((chj[2,p]+1) * (chh[2,p]+1) * (chi[2,p]+1))
                end
                (pftot == 1) && (fac = -fac)
                mat9j[i + ist, j + jst, d] = fac
            end
        end
    end

    return SCompOperator{T}(cpspd, cpspf, nd, ltot, coeff, sgop, mat9j)
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
        idel_rng = [ cpop.sgop[p, d].colptr[idsecj[p]] : cpop.sgop[p, d].colptr[idsecj[p] + 1] - 1 for p = 1 : np ]

        jrng_sg = Vector{UnitRange{Int64}}(undef, np)
        irng_sg = Vector{UnitRange{Int64}}(undef, np)
        for idel in Iterators.product(idel_rng...)
            idseci = [ cpop.sgop[p, d].rowid[idel[p]] for p = 1 : np]
            isec = searchsortedfirst(axes(cpop.cpspf.idsec, 2), idseci, lt = (k, t) -> isless(@view(cpop.cpspf.idsec[:, k]), t))
            (isec > size(cpop.cpspf.idsec, 2) || @view(cpop.cpspf.idsec[:, isec]) != idseci) && continue
            for jch in eachindex(cpop.cpspd.chs[jsec])
                jrng = cpop.cpspd.ptr_st[jsec][jch] : cpop.cpspd.ptr_st[jsec][jch + 1] - 1
                for p = 1 : np
                    lj = cpop.cpspd.chs[jsec][jch][1, p]
                    idlj = cpop.cpspd.sgsp[p].l_lookup[idsecj[p]][lj]
                    jrng_sg[p] = (cpop.cpspd.sgsp[p].ptr_st[idsecj[p]][idlj] + 1 : cpop.cpspd.sgsp[p].ptr_st[idsecj[p]][idlj + 1]) .- cpop.cpspd.sgsp[p].ptr_st[idsecj[p]][1]
                end
                for ich in eachindex(cpop.cpspf.chs[isec])
                    fac9j = cpop.mat9j[cpop.cpspf.ptr_ch[isec] - 1 + ich, cpop.cpspd.ptr_ch[jsec] - 1 + jch, d]
                    abs(fac9j) < √eps(Float64) && continue
                    irng = cpop.cpspf.ptr_st[isec][ich] : cpop.cpspf.ptr_st[isec][ich + 1] - 1
                    for p = 1 : np
                        li = cpop.cpspf.chs[isec][ich][1, p]
                        idli = cpop.cpspf.sgsp[p].l_lookup[idseci[p]][li]
                        irng_sg[p] = (cpop.cpspf.sgsp[p].ptr_st[idseci[p]][idli] + 1 : cpop.cpspf.sgsp[p].ptr_st[idseci[p]][idli + 1]) .- cpop.cpspf.sgsp[p].ptr_st[idseci[p]][1]
                    end
                    @views stf1[irng] .+= (coeff * fac9j) .* (⊗([cpop.sgop[p, d].elmat[idel[p]][irng_sg[p], jrng_sg[p]] for p = 1 : np]...) * std[jrng])
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
    GetEigensystem(cpop :: SCompOperator{T}, nst :: Int64 ; tol :: Float64, ncv :: Int64, initvec :: Vector{T}, disp_std :: Bool, kwargs...) :: Tuple{Vector{T}, Matrix{T}}

computes the lowest `nst` eigenvalues and eigenvectors of the composite operator `cpop` by feeding its matrix–vector product `x -> cpop * x` to the Lanczos/Arnoldi method of `KrylovKit.eigsolve`. Since `cpop` already acts within a definite total angular momentum sector, this directly yields the low-lying spectrum resolved by ``l``.

# Arguments

* `cpop :: SCompOperator{T}` is the composite operator (usually the Hamiltonian).
* `nst :: Int64` is the number of eigenpairs to compute.
* `tol :: Float64` is the tolerance of the eigensolver. Facultative, `1E-8` by default.
* `ncv :: Int64` is the dimension of the Krylov subspace. Facultative, `max(2 * nst, nst + 10)` by default.
* `initvec :: Vector{T}` is the initial vector. Facultative, a random vector by default.
* `disp_std :: Bool` controls whether the log is displayed. Facultative, `!FuzzifiED.SilentStd` by default.
* `kwargs...` are further keyword arguments forwarded to `eigsolve`.

# Output

* `eigval :: Vector{T}` is the vector of the `nst` lowest eigenvalues.
* `eigvec :: Matrix{T}` is the matrix whose columns are the corresponding eigenvectors.
"""
function FuzzifiED.GetEigensystem(cpop :: SCompOperator{T}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = max(2 * nst, nst + 10), initvec = rand(T, cpop.cpspd.dim), disp_std = !FuzzifiED.SilentStd, kwargs...) where T <: Union{ComplexF64,Float64}
    kwargs1 = haskey(kwargs, :krylovdim) ? kwargs : (kwargs..., krylovdim = ncv)
    eigval, eigvec, info = eigsolve(x -> cpop * x, initvec, nst, :SR ; tol, kwargs1...)
    #print(info)
    return Vector{T}(eigval), Matrix{T}(hcat(eigvec...))
end
