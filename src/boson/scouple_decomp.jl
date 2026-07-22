export SCoupleDecomp
export ContactSCouple, SingleSegSCouple


"""
    SCoupleDecomp

The mutable type `SCoupleDecomp` records the decomposition of couplings into channels of direct-producted spherical-symmetric actions onto each part. _E. g._ for bipartite and tri-partite systems, each coupling channel may take the form
```math
\\begin{aligned}
    [𝒪]_{(l_1l_2)l}&=[𝒪_1]_{l_1}⊗[𝒪_2]_{l_2}&𝒪_{lm}&=[𝒪_1]_{l_1l_1}[𝒪_2]_{l_2m_2}⟨l_1m_1,l_2m_2|lm⟩\\\\
    [𝒪]_{((l_1l_2)l_{12}l_3)l}&=([𝒪_1]_{l_1}⊗[𝒪_2]_{l_2})_{l_{12}}[𝒪_3]_{l_3}&𝒪_{lm}&=[𝒪_1]_{l_1l_1}[𝒪_2]_{l_2m_2}[𝒪_3]_{l_3m_3}⟨l_1m_1,l_2m_2|l_{12}m_{12}⟩⟨l_{12}m_{12},l_3m_3|lm⟩
\\end{aligned}
```

# Fields

* `amd :: Vector{SAngModes}` records, for each part `p`, the spherical-symmetric action on that part, stored as an `SAngModes` object.
* `ch :: Vector{Matrix{Int64}}` is the list of coupling channels.
* `coeff :: Vector{ComplexF64}` is the coefficient of each channel.
* `sec :: Matrix{Int64}` records the change of quantum numbers that the term induces : `sec[iqn, p]` is the shift of the `iqn`-th diagonal quantum number on part `p`.
"""
mutable struct SCoupleDecomp
    amd :: Vector{SAngModes}
    ch :: Vector{Matrix{Int64}}
    coeff :: Vector{ComplexF64}
    sec :: Matrix{Int64} # sec[iqn, p]
end


"""
    SCoupleDecomp(amd :: Vector{SAngModes}, ch :: Matrix{Int64}, sec :: Matrix{Int64}) :: SCoupleDecomp

constructs a single-channel `SCoupleDecomp` with unit coefficient from the operators `amd`, a single coupling channel `ch` and the quantum number shift `sec`.
"""
function SCoupleDecomp(amd :: Vector{SAngModes}, ch :: Matrix{Int64}, sec :: Matrix{Int64})
    return SCoupleDecomp(amd, [ch], [1.0], sec)
end


"""
    cpd1 + cpd2 :: Vector{SCoupleDecomp}
    cpd1 - cpd2 :: Vector{SCoupleDecomp}
    -cpd :: Vector{SCoupleDecomp}
    fac * cpd :: Vector{SCoupleDecomp}

enable the linear combination of coupling decompositions.
"""
function Base.:+(cpd1 :: Union{SCoupleDecomp, Vector{SCoupleDecomp}}, cpd2 :: Union{SCoupleDecomp, Vector{SCoupleDecomp}})
    return [ cpd1 ; cpd2 ]
end
function Base.:*(fac :: Number, cpd :: SCoupleDecomp)
    return SCoupleDecomp(cpd.amd, cpd.ch, fac .* cpd.coeff, cpd.sec)
end
function Base.:*(fac :: Number, cpd :: Vector{SCoupleDecomp})
    return fac .* cpd
end
function Base.:-(cpd :: Union{SCoupleDecomp, Vector{SCoupleDecomp}})
    return (-1) * cpd
end
function Base.:-(cpd1 :: Union{SCoupleDecomp, Vector{SCoupleDecomp}}, cpd2 :: Union{SCoupleDecomp, Vector{SCoupleDecomp}})
    return cpd1 + (-1) * cpd2
end


"""
    CountChannels(cpd :: Vector{SCoupleDecomp}) :: Int64 

return the total number of channels for a list of CoupleDecomps
"""
function CountChannels(cpd :: Vector{SCoupleDecomp})
    return sum([length(cpdi.ch) for cpdi in cpd])
end

function Fuzzifino.SAngModes(obs :: SSphereObs)
    return SAngModes(obs.l2m, obs.get_comp)
end


"""
    ContactSCouple(obs :: Vector{SSphereObs}, sec :: Matrix{Int64}, ltot :: Int64) :: SCoupleDecomp

constructs a [SCoupleDecomp](@ref SCoupleDecomp) for a contact term, _i. e._ the product of one spherical observable per part evaluated at the same point on the sphere
```math
    ∫\\mathrm{d}^2𝐫\\,\\sqrt{4π}Ȳ_{lm}(𝐫)\\,n_1(𝐫)n_2(𝐫)⋯n_{N_p}(𝐫)
```

# Arguments

* `obs :: Vector{SSphereObs}` is the list of spherical observables acting on each segment.
* `sec :: Matrix{Int64}` records the change of quantum numbers, `sec[iqn, p]` for the `iqn`-th quantum number on part `p`.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}`` of the term. Facultative, ``0`` (a scalar) by default.

# Output

* `cpd :: SCoupleDecomp` is the resulting coupling decomposition.
"""
function ContactSCouple(obs :: Vector{SSphereObs}, sec :: Matrix{Int64}, ltot :: Int64 = 0)
    amd = SAngModes.(obs)
    np = length(obs)
    s2 = [ obsi.s2 for obsi in obs ]
    s2_ptsum = cumsum(s2)
    l_rng = [ abs(obsi.s2) : 2 : obsi.l2m for obsi in obs ]
    ch1 = Matrix{Int64}[]
    ch = Matrix{Int64}[]
    coeff = ComplexF64[]
    for li in Iterators.product(l_rng...)
        chi = FindCouplingChannels(np, collect(li), ltot)
        append!(ch1, chi)
    end
    for chi in ch1
        coeffi = 1
        flag = true
        for p = 2 : np
            if (chi[2, p] < abs(s2_ptsum[p]))
                flag = false
                break
            end
            coeffi *= clebschgordan(chi[2, p - 1]/2, -s2_ptsum[p - 1]/2, chi[1, p]/2, -s2[p]/2, chi[2, p]/2, -s2_ptsum[p]/2)
        end
        flag || continue
        coeffi *= √(prod(chi[1, :] .+ 1) / (ltot + 1)) * FuzzifiED.ObsNormRadSq / (4π) ^ (np/2 - 1)
        push!(ch, chi)
        push!(coeff, coeffi)
    end
    return SCoupleDecomp(amd, ch, coeff, sec)
end


"""
    SingleSegSCouple(np :: Int64, p :: Int64, amdp :: SAngModes, l :: Int64, secp :: Vector{Int64}) :: SCoupleDecomp
    SingleSegSCouple(np :: Int64, p :: Int64, tms :: STerms, sec :: Vector{Int64}) :: SCoupleDecomp

constructs a [SCoupleDecomp](@ref SCoupleDecomp) for a term that acts only on a single part ``p`` and as the identity on all the other parts.  

# Arguments

* `np :: Int64` is the number of parts.
* `p :: Int64` is the index of the part on which the operator acts.
* `amdp :: SAngModes` is the spherical tensor operator acting on part ``p`.
* `l :: Int64` is twice the rank ``2l`` of the operator on part `p`.
* `secp :: Vector{Int64}` is the change of quantum numbers on part `p`.

In the second form the operator is a scalar (rank ``0``) given directly as a list of terms `tms :: STerms`, and `sec` is its quantum number shift.

# Output

* `cpd :: SCoupleDecomp` is the resulting coupling decomposition.
"""
function SingleSegSCouple(np :: Int64, p :: Int64, amdp :: SAngModes, l :: Int64, secp :: Vector{Int64})
    amd = [one(SAngModes) for p = 1 : np]
    amd[p] = amdp
    ch = zeros(Int64, 2, np)
    ch[1, p] = l
    ch[2, p : end] .= l
    sec = zeros(Int64, length(secp), np)
    sec[:, p] = secp
    return SCoupleDecomp(amd, [ch], [1], sec)
end
function SingleSegSCouple(np :: Int64, p :: Int64, tms :: STerms, sec :: Vector{Int64})
    amdp = SAngModes(0, Dict((0, 0) => tms))
    return SingleSegSCouple(np, p, amdp, 0, sec)
end


"""
    InsertSegment(np :: Int64, p_rng :: Vector{Int64}, cpd :: SCoupleDecomp) :: CoupleDecomp
    InsertSegment(np :: Int64, p_rng :: Vector{Int64}, cpd :: Vector{CoupleDecomp}) :: Vector{SCoupleDecomp}

Insert segments where the coupling acts as identity to SCoupleDecomp(s)
```math
    H_1⊗H_2⊗⋯↦𝕀⊗⋯⊗H_1⊗𝕀⊗⋯⊗H_2⊗𝕀⊗⋯
```

# Arguments 
* `np :: Int64` is the number of segments after the insertion.
* `p_rng :: Int64` is the positions of the non-trivial segments after the insertion.
* `cpd :: SCoupleDecomp`, `cpd :: Vector{SCoupleDecomp}` is the SCoupleDecomp(s) before the insertion. 
"""
function InsertSegment(np :: Int64, p_rng :: Vector{Int64}, cpd :: SCoupleDecomp)
    amd1 = ones(SAngModes, np)
    amd1[p_rng] = cpd.amd
    ch1 = InsertSegment(np, p_rng, cpd.ch)
    sec1 = zeros(Int64, size(cpd.sec, 1), np)
    sec1[:, p_rng] = cpd.sec
    return SCoupleDecomp(amd1, ch1, cpd.coeff, sec1)
end
InsertSegment(np :: Int64, p_rng :: Vector{Int64}, cpds :: Vector{SCoupleDecomp}) = InsertSegment.(Ref(np), Ref(p_rng), cpds)


"""
    PrepareCouple(cpd :: SCoupleDecomp ; eltype :: Type) :: SCoupleDecomp
    PrepareCouple(cpd :: Vector{SCoupleDecomp} ; eltype :: Type) :: Vector{SCoupleDecomp}

prepares a coupling decomposition for the construction of operators by pre-storing the components of each `SAngModes` and, when `eltype == Float64`, rotating any purely imaginary operator by ``i`` while compensating the phase in the coefficients, so that all the matrix elements can be represented with real numbers. This should be called once on an assembled operator before it is passed to [BuildSSegOperators](@ref BuildSSegOperators).

# Arguments

* `cpd :: SCoupleDecomp` or `Vector{SCoupleDecomp}` is the coupling decomposition to be prepared.
* `eltype :: Type` is the target matrix element type, either `Float64` or `ComplexF64`. Facultative, `FuzzifiED.ElementType` by default.

# Output

* the prepared coupling decomposition, of the same shape as the input.
"""
function PrepareCouple(cpd :: SCoupleDecomp ; eltype = FuzzifiED.ElementType)
    amd1 = SAngModes[]
    ph = 1.0 + 0.0im
    for amdi in cpd.amd
        amdi1 = StoreComps(amdi)
        if (eltype == Float64)
            coeff1 = collect(amdi1.comps)[1][2][1].coeff
            if (abs(coeff1.re / coeff1.im) < 1E-4)
                amdi1 *= 1.0im 
                ph *= -1.0im
            end
        end
        push!(amd1, amdi1)
    end
    return SCoupleDecomp(amd1, cpd.ch, cpd.coeff .* ph, cpd.sec)
end
PrepareCouple(cpd :: Vector{SCoupleDecomp} ; eltype = FuzzifiED.ElementType) = PrepareCouple.(cpd ; eltype)