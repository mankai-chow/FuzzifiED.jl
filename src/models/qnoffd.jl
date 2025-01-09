export GetParityQNOffd, GetFlavPermQNOffd, GetRotyQNOffd


function PermDictOrVector(perm_dict :: Dict{Int64, Int64}, nf :: Int64)
    perm = collect(1 : nf)
    for (f, f1) in perm_dict
        perm[f] = f1
    end
    return perm
end
function PermDictOrVector(perm_cyc :: Vector{Vector{Int64}}, nf :: Int64)
    perm = collect(1 : nf)
    for cyc in perm_cyc
        for f = 1 : length(cyc) - 1
            perm[cyc[f]] = cyc[f + 1]
        end
        perm[cyc[end]] = cyc[1]
    end
    return perm
end
PermDictOrVector(perm :: Vector{Int64}, nf :: Int64) = perm

function DictOrVectorPhase(ph_dict :: Dict{Int64, T}, nf :: Int64) where T <: Number
    ph = ones(T, nf)
    for (f, p) in ph_dict 
        ph[f] = p
    end
    return ph
end
DictOrVectorPhase(ph :: Vector{<: Number}, nf :: Int64) = ph


"""
    GetParityQNOffd(nm :: Int64, nf :: Int64[, permf, fac])
    
Return the particle-hole transformation 
```math
    𝒫: c^†_{mf}↦α_fc_{mπ_f}
```
# Arguments 

* `nm :: Int64` and `nf :: Int64` are the number of orbitals and the flavours.
* `permf :: Dict{Int64, Int64}`, `permf :: Vector{Vector{Int64}}` or `Vector{Int64}` gives the flavour permutation ``π_f``. It is either a vector of the cycles, a vector of the target flavours, or a dictionary of the changed elements. _E.g._, a permutation ``1↦4,2↦5,3↦3,4↦1,5↦6,6↦2`` can be expressed as `[4,5,3,1,6,2]`, `[[1,4],[2,5,6]]` or `Dict(1=>4,2=>5,4=>1,5=>6,6=>2)`. Facultative, identity by default. 
* `fac :: Dict{Int64, <: Number}` or `Vector{<: Number}` gives the factor ``α_f``. It is either a vector of all vectors, or a dictionary of all non-unity elements. Facultative, all unity by default. 
"""
function GetParityQNOffd(nm :: Int64, nf :: Int64, permf :: Union{Dict{Int64, Int64}, Vector{Vector{Int64}}, Vector{Int64}} = Dict{Int64, Int64}(), fac :: Union{Dict{Int64, <: Number}, Vector{<: Number}} = Dict{Int64, ComplexF64}()) 
    permf1 = PermDictOrVector(permf, nf)
    fac1 = ComplexF64.(DictOrVectorPhase(fac, nf))
    return QNOffd(
        vcat([permf1[f] + (m - 1) * nf for f = 1 : nf, m = 1 : nm]...),
        true, 
        vcat([fac1[f] for f = 1 : nf, m = 1 : nm]...)
    )
end


"""
    GetFlavPermQNOffd(nm :: Int64, nf :: Int64, permf, fac][, cyc :: Int64])
    
Return the flavour permutaiton transformation 
```math
    𝒵: c^†_{mf}↦α_fc^†_{mπ_f}
```
# Arguments 

* `nm :: Int64` and `nf :: Int64` are the number of orbitals and the flavours.
* `permf :: Dict{Int64, Int64}`, `permf :: Vector{Vector{Int64}}` or `Vector{Int64}` gives the flavour permutation ``π_f``. It is either a vector of the cycles, a vector of the target flavours, or a dictionary of the changed elements. _E.g._, a permutation ``1↦4,2↦5,3↦3,4↦1,5↦6,6↦2`` can be expressed as `[4,5,3,1,6,2]`, `[[1,4],[2,5,6]]` or `Dict(1=>4,2=>5,4=>1,5=>6,6=>2)`. Facultative, identity by default. 
* `fac :: Dict{Int64, <: Number}` or `Vector{<: Number}` gives the factor ``α_f``. It is either a vector of all vectors, or a dictionary of all non-unity elements. Facultative, all unity by default. 
* `cyc :: Int64` is the period of the permutation. 
"""
function GetFlavPermQNOffd(nm :: Int64, nf :: Int64, permf :: Union{Dict{Int64, Int64}, Vector{Vector{Int64}}, Vector{Int64}}, fac :: Union{Dict{Int64, <: Number}, Vector{<: Number}} = Dict{Int64, ComplexF64}(), cyc :: Int64 = 2)
    permf1 = PermDictOrVector(permf, nf)
    fac1 = ComplexF64.(DictOrVectorPhase(fac, nf))
    return QNOffd(
        vcat([permf1[f] + (m - 1) * nf for f = 1 : nf, m = 1 : nm]...),
        vcat([fac1[f] for f = 1 : nf, m = 1 : nm]...),
        cyc
    )
end
GetFlavPermQNOffd(nm :: Int64, nf :: Int64, permf :: Union{Dict{Int64, Int64}, Vector{Vector{Int64}}, Vector{Int64}}, cyc :: Int64) = GetFlavPermQNOffd(nm, nf, permf, Dict{Int64, Int64}(), cyc)


"""
    GetRotyQNOffd(nm :: Int64, nf :: Int64)
    
Return the ``π``-rotation with respect to the ``y``-axis. 
```math
    ℛ_y: c^†_{mf}↦(-)^{m+s}c^†_{-mf}
```
"""
GetRotyQNOffd(nm :: Int64, nf :: Int64) = QNOffd(
    vcat([f + (nm - m) * nf for f = 1 : nf, m = 1 : nm]...), 
    ComplexF64(-1) .^ (collect(0 : nm * nf - 1) .÷ nf)
)