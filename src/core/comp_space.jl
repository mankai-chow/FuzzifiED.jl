export CompSpace, BuildCompSpace
export EquivSec, ComposeSec, FindCouplingChannels


"""
    CompSpace{Float64}
    CompSpace{ComplexF64}

The mutable type `CompSpace` stores the composite Hilbert space obtained by combining the segment spaces of all the parts and projecting onto a definite total angular momentum ``l_{\\text{tot}}``. A basis state of the composite space is specified by a composite sector (which sector each part sits in), a coupling channel (the angular momentum of each part and the intermediate totals along the coupling chain), and the position of the multiplet within each part. _E. g._ for bipartite and tri-partite systems, it may take the form 
```math
\\begin{aligned}
    |Q_{\\{12\\}}C_{2,\\{12\\}},(l_1l_2)lm,α_{\\{12\\}}⟩&=∑_{m_1m_2}|Q_1C_{2,1},l_1m_1,α_1⟩|Q_2C_{2,2},l_2m_2,α_2⟩⟨l_1m_1,l_2m_2|lm⟩\\\\
    |Q_{\\{123\\}}C_{2,\\{123\\}},((l_1l_2)l_{12}l_3)lm,α_{\\{123\\}}⟩&=∑_{m_1m_2m_3}|Q_1C_{2,1},l_1m_1,α_1⟩|Q_2C_{2,2},l_2m_2,α_2⟩|Q_3C_{2,3},l_3m_3,α_3⟩⟨l_1m_1,l_2m_2|l_{12}m_{12}⟩⟨l_{12}m_{12},l_3m_3|lm⟩
\\end{aligned}
```

Angular momenta are stored as twice their value so that they remain integers.

# Fields

* `np :: Int64` is the number of parts.
* `nch :: Int64` is the total number of coupling channels summed over all composite sectors.
* `dim :: Int64` is the total dimension of the composite space.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}``.
* `sgsp :: Vector{<:AbstractSegSpace{T}}` is the list of the [AbstractSegSpaces](@ref AbstractSegSpace) of the parts, which may mix fermionic [SegSpace](@ref SegSpace) and bosonic [SSegSpace](@ref SSegSpace).
* `idsec :: Matrix{Int64}` is the list of composite sector indices. It takes two indices `idsec[p, isec]` where `isec` is the index of the composite sector and `p` is the index of the part. The sector is then given by `sgsp[p].sec[idsec[p, isec]]`.
* `chs :: Vector{Vector{Matrix{Int64}}}` records, for each composite sector, the list of angular momentum coupling channels. Each channel is stored as a ``2×N_p`` matrix ``\\begin{pmatrix}l_1&l_2&⋯&l_p&⋯&l_{N_p}\\\\l_1&l_{12}&⋯&l_{1⋯ p}&⋯&l\\end{pmatrix}``
, where the first row is the angular momentum of each part ``2l_p``, and the second row is the accumulated angular momentum ``2l_{12⋯p}`` of the first ``p`` parts. It takes two indices `chs[isec][ich]` where the `isec` is the index of the composite sector and `ich` is the index of the channel within the sector.
* `ptr_ch :: Vector{Int64}` are the pointers that delimit the channels of each composite sector.
* `ptr_st :: Vector{Vector{Int64}}` records the pointers that delimit the states of each channel.
"""
mutable struct CompSpace{T <: Union{Float64, ComplexF64}}
    np :: Int64
    nch :: Int64
    dim :: Int64
    ltot :: Int64
    sgsp :: Vector{AbstractSegSpace{T}}
    idsec :: Matrix{Int64}
    chs :: Vector{Vector{Matrix{Int64}}}
    ptr_ch :: Vector{Int64}
    ptr_st :: Vector{Vector{Int64}}
end


"""
    BuildCompSpace(sgsp :: Vector{<:AbstractSegSpace{T}}, idsec :: Matrix{Int64}, ltot :: Int64) :: CompSpace
    BuildCompSpace(sgsp :: Vector{<:AbstractSegSpace{T}}, sec_tot :: Vector{Int64}, ltot :: Int64) :: CompSpace

constructs a [CompSpace](@ref CompSpace) from the segment spaces of the parts with total angular momentum `ltot`. For every composite sector it enumerates, through [FindCouplingChannels](@ref FindCouplingChannels), all the ways of coupling the per-part angular momenta into ``l_{\\text{tot}}``, and computes the resulting dimensions and pointers.

# Arguments

* `sgsp :: Vector{<:AbstractSegSpace{T}}` is the list of the [AbstractSegSpaces](@ref AbstractSegSpace) of the parts, which may mix fermionic [SegSpace](@ref SegSpace) and bosonic [SSegSpace](@ref SSegSpace).
* `idsec :: Matrix{Int64}` is the list of composite sector indices. It takes two indices `idsec[p, isec]`.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}``.

# Output

* `cpsp :: CompSpace` is the resulting composite space.
"""
function BuildCompSpace(sgsp :: Vector{<:AbstractSegSpace{T}}, idsec :: Matrix{Int64}, ltot :: Int64 ; disp_std = !FuzzifiED.SilentStd) where T <: Union{Float64, ComplexF64}
    np = length(sgsp)
    chs = [ Matrix{Int64}[] for _ ∈ axes(idsec, 2)]
    ptr_ch = Int64[1]
    ptr_st = [ Int64[] for _ ∈ axes(idsec, 2)]
    index = 1
    for i ∈ axes(idsec, 2)
        idseci = idsec[:, i]
        l_rng = [ sgsp[p].l_rng[idseci[p]] for p = 1 : np]
        for lpt in Iterators.product(l_rng...)
            append!(chs[i], FindCouplingChannels(np, [lpt...], ltot))
        end
        push!(ptr_ch, ptr_ch[end] + length(chs[i]))
        for ich in eachindex(chs[i])
            push!(ptr_st[i], index)
            lpt = chs[i][ich][1, :]
            dimch = 1
            for p = 1 : np
                il = sgsp[p].l_lookup[idseci[p]][lpt[p]]
                dimch *= sgsp[p].ptr_st[idseci[p]][il + 1] - sgsp[p].ptr_st[idseci[p]][il]
            end
            index += dimch
        end
        push!(ptr_st[i], index)
    end
    nch = ptr_ch[end] - 1
    dim = ptr_st[end][end] - 1
    disp_std && @info "FINISH BUILDING COMP SPACE, ANGULAR MOMENTUM $(ltot/2), TOTAL DIMENSION $(dim), # OF CHANNELS $(nch), # OF SECTORS $(size(idsec, 2))"
    return CompSpace{T}(np, nch, dim, ltot, sgsp, idsec, chs, ptr_ch, ptr_st)
end
"""
    BuildCompSpace(sgsp :: Vector{<:AbstractSegSpace{T}}, sec_tot :: Vector{Int64}, ltot :: Int64) :: CompSpace

constructs a [CompSpace](@ref CompSpace) from the segment spaces of the parts with total angular momentum `ltot`. For every composite sector it enumerates, through [FindCouplingChannels](@ref FindCouplingChannels), all the ways of coupling the per-part angular momenta into ``l_{\\text{tot}}``, and computes the resulting dimensions and pointers. The composite sectors are found automatically with [ComposeSec](@ref).

# Arguments

* `sgsp :: Vector{<:AbstractSegSpace{T}}` is the list of the [AbstractSegSpaces](@ref AbstractSegSpace) of the parts, which may mix fermionic [SegSpace](@ref SegSpace) and bosonic [SSegSpace](@ref SSegSpace).
* `sec_tot :: Vector{Int64}` the target total diagonal quantum numbers.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}``.

# Output

* `cpsp :: CompSpace` is the resulting composite space.
"""
function BuildCompSpace(sgsp :: Vector{<:AbstractSegSpace{T}}, sec_tot :: Vector{Int64}, ltot :: Int64 ; disp_std = !FuzzifiED.SilentStd) where T <: Union{Float64, ComplexF64}
    sec_pt = [ sgspi.sec for sgspi in sgsp]
    modul = sgsp[1].sec_modul
    idsec = ComposeSec(sec_tot, sec_pt, modul)
    return BuildCompSpace(sgsp, idsec, ltot ; disp_std)
end


"""
    EquivSec(sec1 :: Vector{Int64}, sec2 :: Vector{Int64}, modul :: Vector{Int64}) :: Bool

tests whether two diagonal quantum number sectors `sec1` and `sec2` are equivalent. The comparison skips the second entry ``L^z`` ; for every other quantum number `i` the entries must agree, either exactly when `modul[i] == 1` or modulo `modul[i]`.
"""
function EquivSec(sec1 :: Vector{Int64}, sec2 :: Vector{Int64}, modul :: Vector{Int64})
    flag = true 
    for i ∈ eachindex(modul)
        (i == 2) && continue
        if (modul[i] == 1) 
            (sec1[i] == sec2[i]) && continue 
        else
            (sec1[i] - sec2[i]) .% modul[i] == 0 && continue 
        end
        flag = false
        break
    end
    return flag
end


"""
    ComposeSec(sec_tot :: Vector{Int64}, sec_pt :: Vector{Matrix{Int64}}[, modul :: Vector{Int64}]) :: Matrix{Int64}

finds every combination of segment sectors into the total sector `sec_tot`.

# Arguments

* `sec_tot :: Vector{Int64}` is the target total diagonal quantum numbers.
* `sec_pt :: Vector{Matrix{Int64}}` lists, for each part, the sectors available in that part. It takes three indices `sec_pt[p][iqn, isec]`.
* `modul :: Vector{Int64}` are the moduli used to match the quantum numbers. Facultative, all ``1`` by default.

# Output

* `idsec_tot :: Matrix{Int64}` is the sorted list of composite sectors, one per column ; each column is a vector of per-part sector indices.
"""
function ComposeSec(sec_tot :: Vector{Int64}, sec_pt :: Vector{Matrix{Int64}}, modul :: Vector{Int64} = fill(1, length(sec_tot)))
    idsec_tot = Vector{Int64}[]
    for isec in Iterators.product(axes.(sec_pt, 2)...)
        seci_tot = sum([ sec_pt[p][:, isec[p]] for p ∈ eachindex(sec_pt) ])
        EquivSec(seci_tot, sec_tot, modul) && push!(idsec_tot, collect(isec))
    end
    idsec_tot = sort(idsec_tot)
    return isempty(idsec_tot) ? Matrix{Int64}(undef, length(sec_pt), 0) : reduce(hcat, idsec_tot)
end


"""
    FindCouplingChannels(np :: Int64, lpt :: Vector{Int64}, ltot :: Int64) :: Vector{Matrix{Int64}}

recursively enumerates every way of coupling the `np` angular momenta `lpt` successively along the chain of parts into the total angular momentum `ltot`. All the angular momenta are given as twice their value.

# Arguments

* `np :: Int64` is the number of parts.
* `lpt :: Vector{Int64}` is the list of the individual angular momenta ``2l_p`` to be coupled.
* `ltot :: Int64` is twice the target total angular momentum ``2l_{\\text{tot}}``.

# Output

* `chs :: Vector{Matrix{Int64}}` is the list of coupling channels.
"""
function FindCouplingChannels(np :: Int64, lpt :: Vector{Int64}, ltot :: Int64)
    chs = Matrix{Int64}[]
    if (np == 1) 
        if (lpt[1] == ltot) 
            push!(chs, [lpt[1] ; lpt[1] ;;])
        end
        return chs
    end
    (sum(lpt) < ltot) && return chs
    for lRem = abs(lpt[end] - ltot) : 2 : lpt[end] + ltot
        chRem = FindCouplingChannels(np - 1, lpt[1 : end - 1], lRem)
        for ch in chRem 
            push!(chs, hcat(ch, [lpt[end], ltot]))
        end
    end
    return chs
end
