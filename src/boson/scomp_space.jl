export SCompSpace, BuildSCompSpace
# The helper functions EquivSec, ComposeSec and FindCouplingChannels are type
# agnostic and shared with the fermionic code in `comp_space.jl`.

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
