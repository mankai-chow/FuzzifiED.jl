export SSegOperator, BuildSSegOperator, BuildSSegOperators

"""
    SSegOperator{Float64}
    SSegOperator{ComplexF64}

The mutable type `SSegOperator` stores the action of a spherical-symmetric operator from an initial [SSegSpace](@ref SSegSpace) to a final SSegSpace. By virtue of the Wigner—Eckart theorem, the full ``m``-dependence is factored out and only the reduced matrix elements ``⟨l_2‖[Φ]_l‖l_1⟩`` — independent of ``l^z`` — need to be kept, and
```math
    ⟨l_2m_2|[Φ]_{lm}|l_1m_1⟩=(-1)^{l_2-m_2}\\begin{pmatrix}l_2&l&l_1\\\\-m_2&m&m_1\\end{pmatrix}⟨l_2‖[Φ]_l‖l_1⟩
```
The matrix ``⟨\\{Q\\}_2l_2α_2‖[Φ]_l‖\\{Q\\}_1l_1α_1⟩`` is stored in blocks of SQNDiag — for given sectors ``\\{Q\\}_{12}``, the matrix elements ``M_{l_2α_2,l_1α_1}`` are stored.

# Fields

* `sgspd :: SSegSpace` and `sgspf :: SSegSpace` are the initial and final segment spaces.
* `colptr :: Vector{Int64}` and `rowid :: Vector{Int64}` store the allowed blocks of sectors ``\\{Q\\}_{12}`` in the format of a sparse matrix.
* `elmat :: Vector{Matrix{T}}` stores, for each block, the reduced matrix elements.
"""
mutable struct SSegOperator{T <: Union{Float64, ComplexF64}}
    sgspd :: SSegSpace
    sgspf :: SSegSpace
    colptr :: Vector{Int64}
    rowid :: Vector{Int64}
    elmat :: Vector{Matrix{T}}
end

"""
    BuildSSegOperator(sgspd :: SSegSpace, sgspf :: SSegSpace, amd :: SAngModes, ll :: Int64, secop :: Vector{Int64} ; eltype :: Type, num_th :: Int64) :: SSegOperator

constructs a [SSegOperator](@ref SSegOperator) for the spherical spherical-symmetric operator `amd` of rank `ll` acting on a single segment. For every pair of sectors related by the quantum number shift `secop`, and every pair of ``l``-multiplets allowed by the triangle rule, it computes the reduced matrix element from the full matrix element via the Wigner—Eckart theorem by dividing out the phase and the ``3j``-symbol. When the ``3j``-symbol vanishes (for ``m_1=m_2=0`` and ``l>0``) the reduced matrix element is instead recovered from the ``m=1`` components, using the ``L^z=1`` states ``L^+|L,0⟩=\\sqrt{l(l+1)}|l,1⟩`` stored in the segment space.

# Arguments

* `sgspd :: SSegSpace` and `sgspf :: SSegSpace` are the initial and final segment spaces.
* `amd :: SAngModes` is the spherical spherical-symmetric operator.
* `ll :: Int64` is twice the rank ``2l`` of the spherical-symmetric operator.
* `secop :: Vector{Int64}` is the change of quantum numbers induced by the operator ; a final sector matches an initial sector when `secd .+ secop` is equivalent to it.
* `eltype :: Type` is the type of the matrix elements, either `Float64` or `ComplexF64`. Facultative, `FuzzifiED.ElementType` by default.
* `num_th :: Int64` is the number of threads. Facultative, ``1`` by default.

# Output

* `sgop :: SSegOperator` is the resulting segment operator.
"""
function BuildSSegOperator(sgspd :: SSegSpace, sgspf :: SSegSpace, amd :: SAngModes, ll :: Int64, secop :: Vector{Int64} ; eltype = FuzzifiED.ElementType, num_th = 1)
    index = 0
    colptr = zeros(Int64, size(sgspd.sec, 2) + 1)
    colptr[1] = 1
    rowid = Int64[]
    elmat = Matrix{eltype}[]
    modul = sgspd.sec_modul

    for j in axes(sgspd.sec, 2)
        secd = sgspd.sec[:, j]
        bsd = SBasis(sgspd.cfs[j])
        std = sgspd.sts[j]
        md = secd[2]
        if (md == 0)
            bsd1 = SBasis(sgspd.cfs1[j])
            std1 = sgspd.sts1[j]
        end

        for i in axes(sgspf.sec, 2)
            secf = sgspf.sec[:, i]
            EquivSec(secd .+ secop, secf, modul) || continue
            index += 1
            push!(rowid, i)
            mf = secf[2]
            mm = mf - md

            bsf = SBasis(sgspf.cfs[i])
            stf = sgspf.sts[i]
            tms = GetComponent(amd, ll/2, mm/2)
            op = SOperator(bsd, bsf, tms)
            op_mat = Matrix(OpMat(op ; num_th))
            hmt_block = stf' * op_mat * std

            if (md == 0 && mf == 0 && ll > 0) # when 3j could vanish
                tms1 = GetComponent(amd, ll/2, mm/2 - 1)
                op1 = SOperator(bsd1, bsf, tms1)
                op_mat1 = Matrix(OpMat(op1 ; num_th))
                hmt_block1 = stf' * op_mat1 * std1
            end

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
                        hmt_block[rngi, rngj] /= fac3j
                    else
                        fac3j1 = wigner3j(lf/2, ll/2, ld/2, -mf/2, mm/2 - 1, md/2 + 1)
                        ((lf - mf) % 4 == 2) && (fac3j1 = -fac3j1)
                        hmt_block[rngi, rngj] = hmt_block1[rngi, rngj] / √(ld/2 * (ld/2 + 1)) / fac3j1
                    end
                end
            end
            push!(elmat, hmt_block)
        end
        colptr[j + 1] = index + 1
    end
    return SSegOperator{eltype}(sgspd, sgspf, colptr, rowid, elmat)
end

"""
    BuildSSegOperators(sgspd :: Vector{SSegSpace{T}}[, sgspf :: Vector{SSegSpace{T}}], cpd :: Vector{SCoupleDecomp}) :: Matrix{SSegOperator}

constructs, in parallel, all the [SSegOperators](@ref SSegOperator) required to assemble the composite operators described by the coupling decompositions `cpd`.

# Arguments

* `sgspd :: Vector{SSegSpace{T}}` is the list of the initial segment spaces.
* `sgspf :: Vector{SSegSpace{T}}` is the list of the final segment spaces. Facultative, the same as `sgspd` by default.
* `cpd :: Vector{SCoupleDecomp}` is the list of coupling decompositions, _e. g._, an assembled Hamiltonian.

# Output

* `sgop :: Matrix{SSegOperator}` is a matrix of segment operators of size ``N_p×N_d``, where ``N_p`` is the number of parts and ``N_d`` the total number of channels of `cpd`. It is passed together with the same `cpd` to [BuildSCompOperator](@ref BuildSCompOperator).
"""
function BuildSSegOperators(sgspd :: Vector{SSegSpace{T}}, sgspf :: Vector{SSegSpace{T}}, cpd :: Vector{SCoupleDecomp}) where T <: Union{Float64, ComplexF64}
    id_tc = vcat([ fill(i, length(cpd[i].ch)) for i in eachindex(cpd)]...)
    ch = vcat([ cpd[i].ch for i in eachindex(cpd) ]...)
    coeff = vcat([ cpd[i].coeff for i in eachindex(cpd) ]...)
    nd = length(id_tc)
    np = length(sgspd)
    sgop = Matrix{SSegOperator}(undef, np, nd)
    Threads.@threads :greedy for (p, d) in collect(Iterators.product(1 : np, 1 : nd))
        amd = cpd[id_tc[d]].amd[p]
        secop = cpd[id_tc[d]].sec[:, p]
        ll = ch[d][1, p]
        sgop[p, d] = BuildSSegOperator(sgspd[p], sgspf[p], amd, ll, secop)
    end
    return sgop
end
BuildSSegOperators(sgspd :: Vector{SSegSpace{T}}, cpd :: Vector{SCoupleDecomp}) where T <: Union{Float64, ComplexF64} = BuildSSegOperators(sgspd, sgspd, cpd)
