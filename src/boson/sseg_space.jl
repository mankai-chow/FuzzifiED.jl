export SSegSpace, BuildSSegSpace

"""
    SSegSpace{Float64}
    SSegSpace{ComplexF64}

The mutable type `SSegSpace` stores the Hilbert space of a single segment (part) of a boson–fermion mixed system, block-diagonalised into multiplets of definite total angular momentum ``l`` and, facultatively, definite flavour Casimir ``C_2``.
```math
    |\\{Q\\}C_2,lm,α⟩
```
where α is the multiplicity of the sector. For each multiplet only one state with representative``m`` is  is stored ; the remaining partners of the multiplet are recovered on demand through the ladder operators. Throughout this type the angular momenta are stored as twice their value (_i. e._ ``2l``, ``2m``) so that they remain integers.

# Fields

* `sec :: Matrix{Int64}` collects the diagonal quantum number (`SQNDiag`) sectors that are diagonalised. It takes two indices `sec[iqn, isec]` where `iqn` is the index of the SQNDiag and `isec` is the index of the sector.
* `l_rng :: Vector{Vector{Int64}}` records, for each sector, the sorted list of the values of ``2l`` that appear. It takes two indices `l_rng[isec][il]`.
* `l_lookup :: Vector{Dict{Int64, Int64}}` gives, for each sector, a dictionary that maps a value of ``2l`` to its index in `l_rng`.
* `ptr_st :: Vector{Vector{Int64}}` records, for each sector, the pointers that delimit the block of states of each ``L`` : the multiplets of angular momentum `l_rng[isec][il]` occupy the columns `ptr_st[isec][il] - ptr_st[isec][1] + 1 : ptr_st[isec][il + 1] - ptr_st[isec][1]` of `sts[isec]`. The pointers are shifted by a cumulative offset so that they are unique across sectors.
* `cfs :: Vector{SConfs}` stores, for each sector, the configurations `SConfs`.
* `cfs1 :: Vector{SConfs}` stores, for each ``m=0`` sector, the configurations of the auxiliary ``m=1`` (_i. e._ ``2m=2``) sector, used when the ``3j``-symbol of the ``m=0`` component vanishes.
* `sts :: Vector{Matrix{T}}` stores, for each sector, the representative multiplet states as its columns, sorted and grouped by ``L`` according to `ptr_st`.
* `sts1 :: Vector{Matrix{T}}` stores, for each ``m=0`` sector, the ``m=1`` components ``L^+`` acting on `sts`, used together with `cfs1` when the ``3j``-symbol vanishes.
"""
mutable struct SSegSpace{T <: Union{Float64, ComplexF64}}
    sec :: Matrix{Int64}
    l_rng :: Vector{Vector{Int64}}
    l_lookup :: Vector{Dict{Int64, Int64}}
    ptr_st :: Vector{Vector{Int64}}
    cfs :: Vector{SConfs}
    cfs1 :: Vector{SConfs}
    sts :: Vector{Matrix{T}}
    sts1 :: Vector{Matrix{T}}
end

"""
    BuildSSegSpace(nof :: Int64, nob :: Int64, nebm :: Vector{Int64}, sec :: Matrix{Int64}, qnd :: Vector{SQNDiag}, tms_lzlp :: Tuple{STerms, STerms}, tms_c2 :: STerms, c2_rng :: Vector{Float64} ; eltype :: Type, num_th :: Int64) :: SSegSpace

constructs a [SSegSpace](@ref SSegSpace) by diagonalising the total angular momentum ``L^2`` — and, facultatively, the flavour Casimir ``C_2`` — within each diagonal quantum number sector, and organising the resulting eigenstates into multiplets.

_N. b._, in including the diagonal quantum numbers, it is required that  the first SQNDiag must contain fermion parity – it must be odd the when state contains odd number of fermions and even when the state constains even number of fermions, and in many cases the total electric charge satisfies this requirement – and the second SQNDiag must be the angular momentum ``2L^z``.

For each sector the operator ``√2L^2+C_2`` is built and diagonalised ; the factor ``√2`` guarantees that the eigenvalues of ``L^2`` and ``C_2`` can be disentangled. Only the multiplets whose ``C_2`` lies within `c2_rng` are retained. For sectors with ``m=0`` the ``m=1`` components ``L^+|l,0⟩=√{l(l+1)}|l,1⟩`` are also computed and stored.

# Arguments

* `nof :: Int64` is the number of fermionic orbitals ``N_{of}`` of the segment.
* `nob :: Int64` is the number of bosonic orbitals ``N_{ob}`` of the segment.
* `nebm :: Vector{Int64}` is the maximal number of bosons allowed in each sector ; it takes one index `nebm[isec]` where `isec` is the index of the sector.
* `sec :: Matrix{Int64}` collects the diagonal quantum number (`SQNDiag`) sectors that are diagonalised. It takes two indices `sec[iqn, isec]` where `iqn` is the index of the SQNDiag and `isec` is the index of the sector.
* `qnd :: Vector{SQNDiag}` is the list of diagonal quantum numbers `SQNDiag`. _N. b._, in the construction, the first SQNDiag must contain fermion parity – it must be odd the when state contains odd number of fermions and even when the state constains even number of fermions, and in many cases the total electric charge satisfies this requirement – and the second SQNDiag must be the ``2L^z`` quantum number, _e. g._, from `GetLz2QNDiag` wrapped in `SQNDiag` (for fermions) or `GetBosonLz2SQNDiag` (for bosons).
* `tms_lzlp :: Tuple{STerms, STerms}` is the pair of terms ``(L^z,L^+)`` from which ``L^2`` is built, _e. g._, from `GetLpLzTerms` converted to `STerms` (for fermions) or `GetBosonLpLzSTerms` (for bosons).
* `tms_c2 :: STerms` is the flavour Casimir ``C_2``. Facultative, no flavour resolution by default.
* `c2_rng :: Vector{Float64}` is the list of allowed eigenvalues of ``C_2`` ; a multiplet is kept when its Casimir is within `1E-4` of one of these values. Facultative, `[0.0]` by default.
* `eltype :: Type` is the type of the matrix elements, either `Float64` or `ComplexF64`. Facultative, `FuzzifiED.ElementType` by default.
* `num_th :: Int64` is the number of threads. Facultative, `FuzzifiED.NumThreads` by default.

# Output

* `sgsp :: SSegSpace` is the resulting [SSegSpace](@ref SSegSpace) object.
"""
function BuildSSegSpace(nof :: Int64, nob :: Int64, nebm :: Vector{Int64}, sec :: Matrix{Int64}, qnd :: Vector{SQNDiag}, tms_lzlp :: Tuple{STerms, STerms}, tms_c2 :: STerms = 0 * one(STerms), c2_rng :: Vector{Float64} = [0.0] ; eltype = FuzzifiED.ElementType, num_th = FuzzifiED.NumThreads)
    nsec = size(sec, 2)
    cfs = Vector{SConfs}(undef, nsec)
    cfs1 = Vector{SConfs}(undef, nsec)
    sts = Vector{Matrix{eltype}}(undef, nsec)
    sts1 = Vector{Matrix{eltype}}(undef, nsec)
    l_rng = [ Int64[] for _ ∈ axes(sec, 2) ]
    ptr_st = [ Int64[] for _ ∈ axes(sec, 2) ]
    l_lookup = [ Dict{Int64, Int64}() for _ ∈ axes(sec, 2) ]
    tms_l2 = GetL2STerms(tms_lzlp)
    BLAS.set_num_threads(num_th)
    for isec ∈ axes(sec, 2)
        seci = sec[:, isec]
        cfs[isec] = SConfs(nof, nob, nebm[isec], seci, qnd ; num_th)
        bs = SBasis(cfs[isec])

        l2c2_mat = Matrix(OpMat(SOperator(bs, √2 * tms_l2 + tms_c2) ; num_th))
        l2c2_val, st = eigen(Hermitian(l2c2_mat))

        l2_mat = OpMat(SOperator(bs, tms_l2) ; num_th)
        l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in axes(st, 2)]
        c2_val = l2c2_val .- √2 .* l2_val
        l_val = round.(Int64, sqrt.(real.(4 * l2_val) .+ 1) .- 1)

        ls = sort(unique(l_val))
        index = 1
        i_rng = Int64[]
        for l in ls
            for i = 1 : bs.dim
                (l_val[i] == l) || continue
                flag = false
                for c2 in c2_rng
                    abs(c2_val[i] - c2) > 1E-4 && continue
                    flag = true
                    break
                end
                flag || continue
                if (isempty(l_rng[isec]) || l_rng[isec][end] ≠ l)
                    push!(l_rng[isec], l)
                    push!(ptr_st[isec], index)
                end
                push!(i_rng, i)
                index += 1
            end
        end
        for il in eachindex(l_rng[isec])
            l_lookup[isec][l_rng[isec][il]] = il
        end
        push!(ptr_st[isec], index)
        sts[isec] = st[:, i_rng]

        if (seci[2] == 0)
            seci1 = deepcopy(seci)
            seci1[2] = 2
            cfs1[isec] = SConfs(nof, nob, nebm[isec], seci1, qnd ; num_th)
            bs1 = SBasis(cfs1[isec])
            lp = SOperator(bs, bs1, tms_lzlp[2])
            lp_mat = Matrix(OpMat(lp))
            sts1[isec] = lp_mat * sts[isec]
        end
        #println("SECTOR $(seci), TOTAL DIMENSION $(bs.dim), SELECTED DIMENSION $(index - 1).")
    end
    BLAS.set_num_threads(1)
    ptr_sec = cumsum([ptr_st[isec][end] - 1 for isec ∈ axes(sec, 2)])
    for isec = 2 : nsec
        ptr_st[isec] .+= ptr_sec[isec - 1]
    end
    return SSegSpace{eltype}(sec, l_rng, l_lookup, ptr_st, cfs, cfs1, sts, sts1)
end
