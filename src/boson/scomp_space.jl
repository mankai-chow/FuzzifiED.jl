export SCompSpace, BuildSCompSpace
# The helper functions EquivSec, ComposeSec and FindCouplingChannels are type
# agnostic and shared with the fermionic code in `comp_space.jl`.

"""
    SCompSpace{Float64}
    SCompSpace{ComplexF64}

The mutable type `SCompSpace` stores the composite Hilbert space obtained by combining the segment spaces of all the parts and projecting onto a definite total angular momentum ``l_{\\text{tot}}``. A basis state of the composite space is specified by a composite sector (which sector each part sits in), a coupling channel (the angular momentum of each part and the intermediate totals along the coupling chain), and the position of the multiplet within each part. _E. g._ for bipartite and tri-partite systems, it may take the form
```math
\\begin{aligned}
    |\\{Q\\}_{12}C_{2,12},(l_1l_2)lm,α_{12}⟩&=∑_{m_1m_2}|\\{Q\\}_1C_{2,1},l_1m_1,α_1⟩|\\{Q\\}_2C_{2,2},l_2m_2,α_2⟩⟨l_1m_1,l_2m_2|lm⟩\\\\
    |\\{Q\\}_{123}C_{2,123},((l_1l_2)l_{12}l_3)lm,α_{123}⟩&=∑_{m_1m_2m_3}|\\{Q\\}_1C_{2,1},l_1m_1,α_1⟩|\\{Q\\}_2C_{2,2},l_2m_2,α_2⟩|\\{Q\\}_3C_{2,3},l_3m_3,α_3⟩⟨l_1m_1,l_2m_2|l_{12}m_{12}⟩⟨l_{12}m_{12},l_3m_3|lm⟩
\\end{aligned}
```

Angular momenta are stored as twice their value so that they remain integers.

# Fields

* `np :: Int64` is the number of parts.
* `nch :: Int64` is the total number of coupling channels summed over all composite sectors.
* `dim :: Int64` is the total dimension of the composite space.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}``.
* `sgsp :: Vector{SSegSpace{T}}` is the list of the [SSegSpaces](@ref SSegSpace) of the parts.
* `idsec :: Matrix{Int64}` is the list of composite sector indices. It takes two indices `idsec[p, isec]` where `isec` is the index of the composite sector and `p` is the index of the part. The sector is then given by `sgsp[p].sec[idsec[p, isec]]`.
* `chs :: Vector{Vector{Matrix{Int64}}}` records, for each composite sector, the list of angular momentum coupling channels. Each channel is stored as a ``2×N_p`` matrix, where the first row is the angular momentum of each part ``2L_p``, and the second row is the accumulated angular momentum ``2L_{12⋯p}`` of the first ``p`` parts. It takes two indices `chs[isec][ich]` where the `isec` is the index of the composite sector and `ich` is the index of the channel within the sector.
* `ptr_ch :: Vector{Int64}` are the pointers that delimit the channels of each composite sector.
* `ptr_st :: Vector{Vector{Int64}}` records the pointers that delimit the states of each channel.
"""
mutable struct SCompSpace{T <: Union{Float64, ComplexF64}}
    np :: Int64
    nch :: Int64
    dim :: Int64
    ltot :: Int64
    sgsp :: Vector{SSegSpace{T}}
    idsec :: Matrix{Int64}
    chs :: Vector{Vector{Matrix{Int64}}}
    ptr_ch :: Vector{Int64}
    ptr_st :: Vector{Vector{Int64}}
end

"""
    BuildSCompSpace(sgsp :: Vector{SSegSpace{T}}, idsec :: Matrix{Int64}, ltot :: Int64) :: SCompSpace
    BuildSCompSpace(sgsp :: Vector{SSegSpace{T}}, sec_tot :: Vector{Int64}, ltot :: Int64, modul :: Vector{Int64}) :: SCompSpace

constructs a [SCompSpace](@ref SCompSpace) from the segment spaces of the parts with total angular momentum `ltot`. For every composite sector it enumerates, through [FindCouplingChannels](@ref FindCouplingChannels), all the ways of coupling the per-part angular momenta into ``l_{\\text{tot}}``, and computes the resulting dimensions and pointers.

# Arguments

* `sgsp :: Vector{SSegSpace{T}}` is the list of the [SSegSpaces](@ref SSegSpace) of the parts.
* `idsec :: Matrix{Int64}` is the list of composite sector indices. It takes two indices `idsec[p, isec]`.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}``.

In the second form the composite sectors are found automatically with [ComposeSec](@ref ComposeSec) from

* `sec_tot :: Vector{Int64}` the target total diagonal quantum numbers, and
* `modul :: Vector{Int64}` the moduli used to match the quantum numbers. Facultative, all ``1`` by default.

# Output

* `cpsp :: SCompSpace` is the resulting composite space.
"""
function BuildSCompSpace(sgsp :: Vector{SSegSpace{T}}, idsec :: Matrix{Int64}, ltot :: Int64) where T <: Union{Float64, ComplexF64}
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
    return SCompSpace{T}(np, nch, dim, ltot, sgsp, idsec, chs, ptr_ch, ptr_st)
end

function BuildSCompSpace(sgsp :: Vector{SSegSpace{T}}, sec_tot :: Vector{Int64}, ltot :: Int64, modul :: Vector{Int64} = fill(1, length(sec_tot))) where T <: Union{Float64, ComplexF64}
    sec_pt = [ sgspi.sec for sgspi in sgsp]
    idsec = ComposeSec(sec_tot, sec_pt, modul)
    return BuildSCompSpace(sgsp, idsec, ltot)
end
