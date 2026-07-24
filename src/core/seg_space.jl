export AbstractSegSpace, SegSpace, BuildSegSpace, BuildSegSpaces


"""
    AbstractSegSpace{T}

Abstract supertype of the single-segment Hilbert spaces. Its concrete subtypes are [SegSpace](@ref SegSpace) for a fermionic segment (backed by FuzzifiED `Confs`) and [SSegSpace](@ref SSegSpace) for a bosonic / mixed segment (backed by Fuzzifino `SConfs`). Both share the same field layout — only the configuration storage and the operator backend differ — so the composite space, composite operator and coupling machinery are written generically on `AbstractSegSpace` and support an arbitrary mixture of fermionic and bosonic segments.
"""
abstract type AbstractSegSpace{T <: Union{Float64, ComplexF64}} end


"""
    SegSpace{Float64}
    SegSpace{ComplexF64}

The mutable type `SegSpace` stores the Hilbert space of a single segment (part) of the system, diagonalized to have definite total angular momentum ``l`` and, facultatively, definite flavour Casimir ``C_2``. 
```math
    |\\{Q\\}C_2,lm,α⟩
```
where ``α`` is the multiplicity of the sector. For each multiplet only one state with representative ``m`` is stored. Throughout this type the angular momenta are stored as twice their value (_i. e._ ``2l``, ``2m``) so that they remain integers.

# Fields

* `sec :: Matrix{Int64}` collects the diagonal quantum number (`QNDiag`) sectors that are diagonalised. It takes two indices `sec[iqn, isec]` where `iqn` is the index of the QNDiag and `isec` is the index of the sector.
* `sec_modul :: Vector{Int64}` collects the moduli of the QNDiags.
* `l_rng :: Vector{Vector{Int64}}` records, for each sector, the sorted list of the values of ``2l`` that appear. It takes two indices `l_rng[isec][il]`.
* `l_lookup :: Vector{Dict{Int64, Int64}}` gives, for each sector, a dictionary that maps a value of ``2l`` to its index in `l_rng`.
* `ptr_st :: Vector{Vector{Int64}}` records, for each sector, the pointers that delimit the block of states of each ``l`` : the multiplets of angular momentum `l_rng[isec][il]` are numbered `ptr_st[isec][il] : ptr_st[isec][il + 1] - 1` across the space.
* `cfs :: Vector{Confs}` stores, for each sector, the configurations `Confs`.
* `cfs1 :: Vector{Confs}` stores, for each ``m=0`` sector, the configurations of the auxiliary ``m=1`` (_i. e._ ``2m=2``) sector, used when the ``3j``-symbol of the ``m=0`` component vanishes.
* `sts :: Vector{Matrix{T}}` stores, for each sector, the states as its columns. These states are numbered `ptr_st[isec][1] : ptr_st[isec][end] - 1`.
* `sts1 :: Vector{Matrix{T}}` stores, for each ``m=0`` sector, the ``m=1`` components ``L^+|l,0⟩=\\sqrt{l(l+1)}|l,1⟩``, used together with `cfs1` when the ``3j``-symbol vanishes.
"""
mutable struct SegSpace{T <: Union{Float64, ComplexF64}} <: AbstractSegSpace{T}
    sec :: Matrix{Int64}
    sec_modul :: Vector{Int64}
    l_rng :: Vector{Vector{Int64}}
    l_lookup :: Vector{Dict{Int64, Int64}}
    ptr_st :: Vector{Vector{Int64}}
    cfs :: Vector{Confs}
    cfs1 :: Vector{Confs}
    sts :: Vector{Matrix{T}}
    sts1 :: Vector{Matrix{T}}
end


"""
    BuildSegSpace(no :: Int64, sec :: Matrix{Int64}, qnd :: Vector{QNDiag}, tms_lzlp :: Tuple{Terms, Terms}[, tms_c2 :: Terms, c2_rng :: Vector{Float64}][, sec_modul :: Vector{Int64}] ; l2c2_ratio :: Float64, nst_max :: Vector{Int64}, eltype :: Type, num_th :: Int64) :: SegSpace

constructs a [SegSpace](@ref SegSpace) by diagonalising the total angular momentum ``L^2`` — and, facultatively, the flavour Casimir ``C_2`` — within each diagonal quantum number sector, and organising the resulting eigenstates into multiplets.

_N. b._, in including the diagonal quantum numbers, it is required that  the first QNDiag must contain fermion parity — it must be odd the when state contains odd number of fermions and even when the state constains even number of fermions, and in many cases the total electric charge satisfies this requirement — and the second QNDiag must be the angular momentum ``2L^z``.

For each sector the operator ``αL^2+C_2`` is built and diagonalised ; the factor `α`` usually guarantees that the eigenvalues of ``L^2`` and ``C_2`` can be disentangled. Only the multiplets whose ``C_2`` lies within `c2_rng` are retained. For sectors with ``m=0`` the ``m=1`` components ``L^+|l,0⟩=\\sqrt{l(l+1)}|l,1⟩`` are also computed and stored.

# Arguments

* `no :: Int64` is the number of orbitals ``N_o`` of the segment.
* `sec :: Matrix{Int64}` collects the diagonal quantum number (`QNDiag`) sectors that are diagonalised. It takes two indices `sec[iqn, isec]` where `iqn` is the index of the QNDiag and `isec` is the index of the sector.
* `qnd :: Vector{QNDiag}` is the list of diagonal quantum numbers `QNDiag`. _N. b._, in the construction, the first QNDiag must contain fermion parity — it must be odd the when state contains odd number of fermions and even when the state constains even number of fermions, and in many cases the total electric charge satisfies this requirement — and the second QNDiag must be the ``2L^z`` quantum number, _e. g._, from `GetLz2QNDiag`.
* `tms_lzlp :: Tuple{Terms, Terms}` is the pair of terms ``(L^z,L^+)`` from which ``L^2`` is built, _e. g._, from `GetLpLzTerms`.
* `tms_c2 :: Terms` is the flavour Casimir ``C_2``. Facultative, no flavour resolution by default.
* `c2_rng :: Vector{Float64}` is the list of allowed eigenvalues of ``C_2`` ; a multiplet is kept when its Casimir is within `1E-4` of one of these values. Facultative, `[0.0]` by default.
* `sec_modul :: Vector{Int64}` collects the moduli of the QNDiags. Facultative, all 1 by default.
* `l2c2_ratio :: Float64` is the ratio ``α`` that determines ``αL^2+C_2`` to be diagonalised. Facultative, ``\\sqrt{2}`` by default.
* `nst_max :: Vector{Int64}` specifies the maximal number of eigen-states for each sector. For each sector, if the number is `0` or exceeds half the total dimension, then ``αL^2+C_2`` is fully diagonalised ; if the number is non-zero, then the lowest `nst_max[i]` eigen-states of ``αL^2+C_2`` will be generated using Arnoldi. 
* `eltype :: Type` is the type of the matrix elements, either `Float64` or `ComplexF64`. Facultative, `FuzzifiED.ElementType` by default.
* `num_th :: Int64` is the number of threads. Facultative, `FuzzifiED.NumThreads` by default.

# Output

* `sgsp :: SegSpace` is the resulting [SegSpace](@ref SegSpace) object.
"""
function BuildSegSpace(no :: Int64, sec :: Matrix{Int64}, qnd :: Vector{QNDiag}, tms_lzlp :: Tuple{Terms, Terms}, tms_c2 :: Terms = zero(Terms), c2_rng :: Vector{Float64} = [0.0], sec_modul :: Vector{Int64} = ones(Int64, size(sec, 1)) ; l2c2_ratio :: Float64 = √2, nst_max :: Vector{Int64} = zeros(Int64, size(sec, 2)), eltype = FuzzifiED.ElementType, num_th = FuzzifiED.NumThreads)
    nsec = size(sec, 2)
    cfs = Vector{Confs}(undef, nsec)
    cfs1 = Vector{Confs}(undef, nsec)
    sts = Vector{Matrix{eltype}}(undef, nsec)
    sts1 = Vector{Matrix{eltype}}(undef, nsec)
    l_rng = [ Int64[] for _ ∈ axes(sec, 2) ]
    ptr_st = [ Int64[] for _ ∈ axes(sec, 2) ]
    l_lookup = [ Dict{Int64, Int64}() for _ ∈ axes(sec, 2) ]
    tms_l2 = GetL2Terms(tms_lzlp)
    BLAS.set_num_threads(num_th)
    for isec ∈ axes(sec, 2)
        seci = sec[:, isec]
        cfs[isec] = Confs(no, seci, qnd ; num_th)
        bs = Basis(cfs[isec])

        l2c2_mat = OpMat(Operator(bs, l2c2_ratio * tms_l2 + tms_c2) ; num_th)
        nsti = nst_max[isec]
        if (nsti == 0 || nsti > bs.dim / 2)
            l2c2_val, st = eigen(Hermitian(Matrix(l2c2_mat)))
        else 
            l2c2_val, st = GetEigensystem(l2c2_mat, nsti ; num_th)
        end

        l2_mat = OpMat(Operator(bs, tms_l2) ; num_th)
        l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in axes(st, 2)]
        c2_val = l2c2_val .- l2c2_ratio .* l2_val
        l_val = round.(Int64, sqrt.(real.(4 * l2_val) .+ 1) .- 1)

        ls = sort(unique(l_val))
        index = 1
        i_rng = Int64[]
        for l in ls
            for i = 1 : size(st, 2)
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
            cfs1[isec] = Confs(no, seci1, qnd ; num_th)
            bs1 = Basis(cfs1[isec])
            lp = Operator(bs, bs1, tms_lzlp[2])
            lp_mat = Matrix(OpMat(lp))
            sts1[isec] = lp_mat * sts[isec]
        end
    end
    BLAS.set_num_threads(1)
    ptr_sec = cumsum([ptr_st[isec][end] - 1 for isec ∈ axes(sec, 2)])
    for isec = 2 : nsec
        ptr_st[isec] .+= ptr_sec[isec - 1]
    end
    @info "FINISH BUILDING SEG SPACE, TOTAL DIMENSION $(ptr_st[end][end] - 1)"
    return SegSpace{eltype}(sec, sec_modul, l_rng, l_lookup, ptr_st, cfs, cfs1, sts, sts1)
end
BuildSegSpace(no :: Int64, sec :: Matrix{Int64}, qnd :: Vector{QNDiag}, tms_lzlp :: Tuple{Terms, Terms}, modul :: Vector{Int64} ; l2c2_ratio :: Float64 = √2, nst_max :: Vector{Int64} = zeros(Int64, size(sec, 2)), eltype = FuzzifiED.ElementType, num_th = FuzzifiED.NumThreads) = BuildSegSpace(no, sec, qnd, tms_lzlp, zero(Terms), [0.0], modul ; l2c2_ratio, nst_max, eltype, num_th)


"""
    BuildSegSpaces(no :: Int64, sec :: Matrix{Int64}, qnd :: Vector{QNDiag}, tms_lzlp :: Tuple{Terms, Terms}[, tms_c2 :: Terms, c2_rng :: Vector{Float64}][, sec_modul :: Vector{Int64}] ; l2c2_ratio :: Float64, nst_max :: Vector{Int64}, eltype :: Type, num_th :: Int64) :: SegSpace

constructs multiple [SegSpaces](@ref SegSpace) simultaneosly with different list of flavour Casimir ``C_2``.

# Arguments

* `no :: Int64` is the number of orbitals ``N_o`` of the segment.
* `sec :: Matrix{Int64}` collects the diagonal quantum number (`QNDiag`) sectors that are diagonalised. It takes two indices `sec[iqn, isec]` where `iqn` is the index of the QNDiag and `isec` is the index of the sector.
* `qnd :: Vector{QNDiag}` is the list of diagonal quantum numbers `QNDiag`. _N. b._, in the construction, the first QNDiag must contain fermion parity — it must be odd the when state contains odd number of fermions and even when the state constains even number of fermions, and in many cases the total electric charge satisfies this requirement — and the second QNDiag must be the ``2L^z`` quantum number, _e. g._, from `GetLz2QNDiag`.
* `tms_lzlp :: Tuple{Terms, Terms}` is the pair of terms ``(L^z,L^+)`` from which ``L^2`` is built, _e. g._, from `GetLpLzTerms`.
* `tms_c2 :: Terms` is the flavour Casimir ``C_2``.
* `c2_rng :: Vector{Vector{Float64}}` is a collection of lists of allowed eigenvalues of ``C_2`` ; for each list within, a SegSpace is generated. 
* `sec_modul :: Vector{Int64}` collects the moduli of the QNDiags. Facultative, all 1 by default.
* `l2c2_ratio :: Float64` is the ratio ``α`` that determines ``αL^2+C_2`` to be diagonalised. Facultative, ``\\sqrt{2}`` by default.
* `nst_max :: Vector{Int64}` specifies the maximal number of eigen-states for each sector. For each sector, if the number is `0`, then ``αL^2+C_2`` is fully diagonalised ; if the number is non-zero, then the lowest `nst_max[i]` eigen-states of ``αL^2+C_2`` will be generated using Arnoldi. 
* `eltype :: Type` is the type of the matrix elements, either `Float64` or `ComplexF64`. Facultative, `FuzzifiED.ElementType` by default.
* `num_th :: Int64` is the number of threads. Facultative, `FuzzifiED.NumThreads` by default.

# Output

* `sgsp :: SegSpace` is the resulting [SegSpace](@ref SegSpace) object.
"""
function BuildSegSpaces(no :: Int64, sec :: Matrix{Int64}, qnd :: Vector{QNDiag}, tms_lzlp :: Tuple{Terms, Terms}, tms_c2 :: Terms, c2_rng :: Vector{Vector{Float64}}, sec_modul :: Vector{Int64} = ones(Int64, size(sec, 1)) ; l2c2_ratio :: Float64 = √2, nst_max :: Vector{Int64} = zeros(Int64, size(sec, 2)), eltype = FuzzifiED.ElementType, num_th = FuzzifiED.NumThreads)
    nsec = size(sec, 2)
    cfs = Vector{Confs}(undef, nsec)
    cfs1 = Vector{Confs}(undef, nsec)
    sts = [ Vector{Matrix{eltype}}(undef, nsec) for _ in c2_rng ]
    sts1 = [ Vector{Matrix{eltype}}(undef, nsec) for _ in c2_rng ]
    l_rng = [ [ Int64[] for _ ∈ axes(sec, 2) ] for _ in c2_rng ]
    ptr_st = [ [ Int64[] for _ ∈ axes(sec, 2) ] for _ in c2_rng ]
    l_lookup = [ [ Dict{Int64, Int64}() for _ ∈ axes(sec, 2) ] for _ in c2_rng ]
    tms_l2 = GetL2Terms(tms_lzlp)
    BLAS.set_num_threads(num_th)
    for isec ∈ axes(sec, 2)
        seci = sec[:, isec]
        cfs[isec] = Confs(no, seci, qnd ; num_th)
        bs = Basis(cfs[isec])

        l2c2_mat = OpMat(Operator(bs, l2c2_ratio * tms_l2 + tms_c2) ; num_th)
        nsti = nst_max[isec]
        if (nsti == 0 || nsti ≥ bs.dim)
            l2c2_val, st = eigen(Hermitian(Matrix(l2c2_mat)))
        else 
            l2c2_val, st = GetEigensystem(l2c2_mat, nsti ; num_th)
        end

        l2_mat = OpMat(Operator(bs, tms_l2) ; num_th)
        l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in axes(st, 2)]
        c2_val = l2c2_val .- l2c2_ratio .* l2_val
        l_val = round.(Int64, sqrt.(real.(4 * l2_val) .+ 1) .- 1)

        ls = sort(unique(l_val))
        index = ones(Int64, length(c2_rng))
        i_rng = [ Int64[] for _ in c2_rng ]
        for ic2 in eachindex(c2_rng)
            for l in ls
                for i = 1 : size(st, 2)
                    (l_val[i] == l) || continue 
                    flag = false
                    for c2 in c2_rng[ic2]
                        abs(c2_val[i] - c2) > 1E-4 && continue 
                        flag = true 
                        break
                    end
                    flag || continue
                    if (isempty(l_rng[ic2][isec]) || l_rng[ic2][isec][end] ≠ l)
                        push!(l_rng[ic2][isec], l)
                        push!(ptr_st[ic2][isec], index[ic2])
                    end
                    push!(i_rng[ic2], i)
                    index[ic2] += 1
                end
            end
            for il in eachindex(l_rng[ic2][isec])
                l_lookup[ic2][isec][l_rng[ic2][isec][il]] = il
            end
            push!(ptr_st[ic2][isec], index[ic2])
            sts[ic2][isec] = st[:, i_rng[ic2]]
        end

        if (seci[2] == 0)
            seci1 = deepcopy(seci)
            seci1[2] = 2
            cfs1[isec] = Confs(no, seci1, qnd ; num_th)
            bs1 = Basis(cfs1[isec])
            lp = Operator(bs, bs1, tms_lzlp[2])
            lp_mat = Matrix(OpMat(lp))
            for ic2 in eachindex(c2_rng)
                sts1[ic2][isec] = lp_mat * sts[ic2][isec]
            end
        end
    end
    BLAS.set_num_threads(1)
    for ic2 in eachindex(c2_rng)
        ptr_sec = cumsum([ptr_st[ic2][isec][end] - 1 for isec ∈ axes(sec, 2)])
        for isec = 2 : nsec
            ptr_st[ic2][isec] .+= ptr_sec[isec - 1]
        end
    end
    @info "FINISH BUILDING $(length(c2_rng)) SEG SPACES, TOTAL DIMENSIONS $([ptr_st[ic2][end][end] - 1 for ic2 in eachindex(c2_rng)])"
    sgsp = [ SegSpace{eltype}(sec, sec_modul, l_rng[ic2], l_lookup[ic2], ptr_st[ic2], cfs, cfs1, sts[ic2], sts1[ic2]) for ic2 in eachindex(c2_rng)]
    return sgsp
end
