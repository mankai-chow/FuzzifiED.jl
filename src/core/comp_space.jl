export CompSpace, BuildCompSpace
export EquivSec, ComposeSec, FindCouplingChannels

mutable struct CompSpace{T <: Union{Float64, ComplexF64}}
    np :: Int64
    nch :: Int64
    dim :: Int64
    ltot :: Int64
    sgsp :: Vector{SegSpace{T}}
    idsec :: Vector{Vector{Int64}}
    chs :: Vector{Vector{Matrix{Int64}}}
    ptr_ch :: Vector{Int64}
    ptr_st :: Vector{Vector{Int64}}
end

function BuildCompSpace(sgsp :: Vector{SegSpace{T}}, idsec :: Vector{Vector{Int64}}, ltot :: Int64) where T <: Union{Float64, ComplexF64}
    np = length(sgsp)
    chs = [ Matrix{Int64}[] for _ ∈ idsec]
    ptr_ch = Int64[1]
    ptr_st = [ Int64[] for _ ∈ idsec]
    index = 1
    for i ∈ eachindex(idsec)
        idseci = idsec[i] 
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
    return CompSpace{T}(np, nch, dim, ltot, sgsp, idsec, chs, ptr_ch, ptr_st)
end

function BuildCompSpace(sgsp :: Vector{SegSpace{T}}, sec_tot :: Vector{Int64}, ltot :: Int64, modul :: Vector{Int64} = fill(1, length(sec_tot))) where T <: Union{Float64, ComplexF64}
    sec_pt = [ sgspi.sec for sgspi in sgsp]
    idsec = ComposeSec(sec_tot, sec_pt, modul)
    return BuildCompSpace(sgsp, idsec, ltot)
end

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

function ComposeSec(sec_tot :: Vector{Int64}, sec_pt :: Vector{Vector{Vector{Int64}}}, modul :: Vector{Int64} = fill(1, length(sec_tot)))
    id_sec_tot = Vector{Int64}[]
    for isec in Iterators.product(eachindex.(sec_pt)...)
        seci_tot = sum([ sec_pt[i][isec[i]] for i ∈ eachindex(modul) ])
        EquivSec(seci_tot, sec_tot, modul) && push!(id_sec_tot, collect(isec))
    end
    return sort(id_sec_tot)
end

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