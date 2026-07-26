export SegOperator, BuildSegOperator, BuildSegOperators


"""
    SegOperator{Float64}
    SegOperator{ComplexF64}

The mutable type `SegOperator` stores the action of a spherical-symmetric operator from an initial [SegSpace](@ref SegSpace) to a final SegSpace. By virtue of the Wigner—Eckart theorem, the full ``m``-dependence is factored out and only the reduced matrix elements ``⟨l_2‖[Φ]_l‖l_1⟩`` — independent of ``l^z`` — need to be kept, and 
```math
    ⟨l_2m_2|[Φ]_{lm}|l_1m_1⟩=(-1)^{l_2-m_2}\\begin{pmatrix}l_2&l&l_1\\\\-m_2&m&m_1\\end{pmatrix}⟨l_2‖[Φ]_l‖l_1⟩
```
The matrix ``⟨\\{Q\\}_2l_2α_2‖[Φ]_l‖\\{Q\\}_1l_1α_1⟩`` is stored in blocks of QNDiag — for given sectors ``\\{Q\\}_{\\{12\\}}``, the matrix elements ``M_{l_2α_2,l_1α_1}`` are stored. 

# Fields

* `sgspd :: SegSpace` and `sgspf :: SegSpace` are the initial and final segment spaces.
* `colptr :: Vector{Int64}` and `rowid :: Vector{Int64}` store the allowed blocks of sectors ``\\{Q\\}_{12}`` in the format of a sparse matrix.
* `elmat :: Vector{Matrix{T}}` stores, for each block, the reduced matrix elements. 
"""
mutable struct SegOperator{T <: Union{Float64, ComplexF64}}
    colptr :: Vector{Int64}
    rowid :: Vector{Int64}
    elmat :: Vector{Matrix{Matrix{T}}}
end


_SegBasis(:: SegSpace, cfs) = Basis(cfs)
_SegBasis(:: SSegSpace, cfs) = SBasis(cfs)
_SegOperator(:: SegSpace, bsd, bsf, tms) = Operator(bsd, bsf, tms)
_SegOperator(:: SSegSpace, bsd, bsf, tms) = SOperator(bsd, bsf, tms)
_SegIdentity(:: SegSpace) = one(AngModes)
_SegIdentity(:: SSegSpace) = one(SAngModes)


"""
    BuildSegOperator(sgspd :: SegSpace[, sgspf :: SegSpace], amd :: AngModes, ll :: Int64, secop :: Vector{Int64} ; eltype :: Type, num_th :: Int64) :: SegOperator
    BuildSegOperator(sgspd :: SSegSpace[, sgspf :: SSegSpace], amd :: SAngModes, ll :: Int64, secop :: Vector{Int64} ; eltype :: Type, num_th :: Int64) :: SegOperator

constructs a [SegOperator](@ref SegOperator) for the spherical spherical-symmetric operator `amd` of rank `ll` acting on a single segment. For every pair of sectors related by the quantum number shift `secop`, and every pair of ``l``-multiplets allowed by the triangle rule, it computes the reduced matrix element from the full matrix element via the Wigner—Eckart theorem by dividing out the phase and the ``3j``-symbol. When the ``3j``-symbol vanishes (for ``m_1=m_2=0`` and ``l>0``) the reduced matrix element is instead recovered from the ``m=1`` components, using the ``L^z=1`` states ``L^+|l,0⟩=\\sqrt{l(l+1)}|l,1⟩`` stored in the segment space.

# Arguments

* `sgspd :: SegSpace` or `sgspd :: SSegSpace` is the initial segment space.
* `sgspf :: SegSpace` or `sgspf :: SSegSpace` is the final segment space. Facultative, the same as `sgspd` by default.
* `amd :: AngModes` or `amd :: SAngModes` is the spherical spherical-symmetric operator.
* `ll :: Int64` is twice the rank ``2l`` of the spherical-symmetric operator.
* `secop :: Vector{Int64}` is the change of quantum numbers induced by the operator ; a final sector matches an initial sector when `secd .+ secop` is equivalent to it.
* `eltype :: Type` is the type of the matrix elements, either `Float64` or `ComplexF64`. Facultative, `FuzzifiED.ElementType` by default.
* `num_th :: Int64` is the number of threads. Facultative, ``1`` by default.

# Output

* `sgop :: SegOperator` is the resulting segment operator.
"""
function BuildSegOperator(sgspd :: AbstractSegSpace, sgspf :: AbstractSegSpace, amd :: Union{AngModes, SAngModes, Symbol}, ll :: Int64, secop :: Vector{Int64} ; eltype = FuzzifiED.ElementType, num_th = FuzzifiED.NumThreads)
    (amd === :Identity) && (amd = _SegIdentity(sgspd))
    index = 0
    colptr = zeros(Int64, size(sgspd.sec, 2) + 1)
    colptr[1] = 1
    rowid = Int64[]
    elmat = Matrix{Matrix{eltype}}[]
    modul = sgspd.sec_modul

    for j in axes(sgspd.sec, 2)
        secd = sgspd.sec[:, j]
        bsd = _SegBasis(sgspd, sgspd.cfs[j])
        std = sgspd.sts[j]
        md = secd[2]
        if (md == 0) 
            bsd1 = _SegBasis(sgspd, sgspd.cfs1[j])
            std1 = sgspd.sts1[j]
        end
        
        for i in axes(sgspf.sec, 2)
            secf = sgspf.sec[:, i]
            EquivSec(secd .+ secop, secf, modul) || continue
            index += 1
            push!(rowid, i)
            mf = secf[2]
            mm = mf - md

            bsf = _SegBasis(sgspf, sgspf.cfs[i])
            stf = sgspf.sts[i]
            tms = GetComponent(amd, ll/2, mm/2)
            op = _SegOperator(sgspd, bsd, bsf, tms)
            op_mat = Matrix(OpMat(op ; num_th))
            hmt_block = stf' * op_mat * std

            if (md == 0 && mf == 0 && ll > 0) # when 3j could vanish
                tms1 = GetComponent(amd, ll/2, mm/2 - 1)
                op1 = _SegOperator(sgspd, bsd1, bsf, tms1)
                op_mat1 = Matrix(OpMat(op1 ; num_th))
                hmt_block1 = stf' * op_mat1 * std1
            end

            hmt_mat = Matrix{Matrix{eltype}}(undef, length(sgspf.l_rng[i]), length(sgspf.l_rng[j]))

            for jl in eachindex(sgspd.l_rng[j])
                ld = sgspd.l_rng[j][jl]
                rngj = (sgspd.ptr_st[j][jl] + 1 : sgspd.ptr_st[j][jl + 1]) .- sgspd.ptr_st[j][1]
                for il in eachindex(sgspf.l_rng[i])
                    lf = sgspf.l_rng[i][il]
                    (ll < abs(ld - lf) || ll > ld + lf) && continue 
                    rngi = (sgspf.ptr_st[i][il] + 1 : sgspf.ptr_st[i][il + 1]) .- sgspf.ptr_st[i][1]
                    fac3j = wigner3j(lf/2, ll/2, ld/2, -mf/2, mm/2, md/2)
                    ((lf - mf) % 4 == 2) && (fac3j = -fac3j)
                    if (fac3j ≠ 0) 
                        hmt_mat[il, jl] = hmt_block[rngi, rngj] / fac3j
                    else
                        fac3j1 = wigner3j(lf/2, ll/2, ld/2, -mf/2, mm/2 - 1, md/2 + 1)
                        ((lf - mf) % 4 == 2) && (fac3j1 = -fac3j1)
                        hmt_mat[il, jl] = hmt_block1[rngi, rngj] / √(ld/2 * (ld/2 + 1)) / fac3j1
                    end
                end
            end
            push!(elmat, hmt_mat)
        end
        colptr[j + 1] = index + 1
    end
    return SegOperator{eltype}(colptr, rowid, elmat)
end
BuildSegOperator(sgspd :: AbstractSegSpace, amd :: Union{AngModes, SAngModes, Symbol}, ll :: Int64, secop :: Vector{Int64} ; eltype = FuzzifiED.ElementType, num_th = 1) = BuildSegOperator(sgspd, sgspd, amd, ll, secop ; eltype, num_th)


"""
    BuildSegOperators(sgspd :: Vector{SegSpace{T}}[, sgspf :: Vector{SegSpace{T}}], cpd :: CoupleDecomps :: Matrix{SegOperator}

constructs, in parallel, all the [SegOperators](@ref SegOperator) required to assemble the composite operators described by the coupling decompositions `cpd`. 

# Arguments

* `sgspd :: Vector{SegSpace{T}}` is the list of the initial segment spaces.
* `sgspf :: Vector{SegSpace{T}}` is the list of the final segment spaces. Facultative, the same as `sgspd` by default.
* `cpd :: CoupleDecomps` is the list of coupling decompositions, _e. g._, an assembled Hamiltonian.
* `p_rng :: Vector{Int64}`. When specified, only the SegOperators of the specified parts will be generated. It must be of the same length as `sgspd`. An array ``1:N_p`` by default.

# Output

* `sgop :: Matrix{SegOperator}` is a matrix of segment operators of size ``N_p×N_d``, where ``N_p`` is the number of parts and ``N_d`` the total number of channels of `cpd`. It is passed together with the same `cpd` to [BuildCompOperator](@ref BuildCompOperator).
"""
function BuildSegOperators(sgspd :: Vector{<:AbstractSegSpace{T}}, sgspf :: Vector{<:AbstractSegSpace{T}}, cpd :: CoupleDecomps ; p_rng :: Vector{Int64} = collect(eachindex(sgspd))) where T <: Union{Float64, ComplexF64}
    nd = length(cpd)
    np = length(p_rng)
    sgop = Matrix{SegOperator}(undef, np, nd)
    # Threads.@threads :greedy for (ip, d) in collect(Iterators.product(1 : np, 1 : nd))
    for d = 1 : nd, ip = 1 : np
        p = p_rng[ip]
        amd = cpd[d].amd[p]
        secop = cpd[d].sec[:, p]
        ll = cpd[d].ch[1, p]
        sgop[ip, d] = BuildSegOperator(sgspd[ip], sgspf[ip], amd, ll, secop)
    end
    @info "FINISH BUILDING $np * $nd SEG OPERATORS"
    return sgop
end
BuildSegOperators(sgspd :: Vector{<:AbstractSegSpace{T}}, cpd :: CoupleDecomps ; p_rng :: Vector{Int64} = collect(eachindex(sgspd))) where T <: Union{Float64, ComplexF64} = BuildSegOperators(sgspd, sgspd, cpd ; p_rng)
