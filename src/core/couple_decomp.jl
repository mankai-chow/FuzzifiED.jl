export CoupleDecomp
export RecouplePsPot, ConvPsPot, ContactCouple, SingleSegCouple
export PrepareCouple


"""
    CoupleDecomp

The mutable type `CoupleDecomp` records the decomposition of couplings into channels of direct-producted spherical-symmetric actions onto each part. _E. g._ for bipartite and tri-partite systems, each coupling channel may take the form 
```math
\\begin{aligned}
    [𝒪]_{(l_1l_2)l}&=[𝒪_1]_{l_1}⊗[𝒪_2]_{l_2}&𝒪_{(l_1l_2)lm}&=[𝒪_1]_{l_1l_1}[𝒪_2]_{l_2m_2}⟨l_1m_1,l_2m_2|lm⟩\\\\
    [𝒪]_{((l_1l_2)l_{12}l_3)l}&=([𝒪_1]_{l_1}⊗[𝒪_2]_{l_2})_{l_{12}}[𝒪_3]_{l_3}&𝒪_{lm}&=[𝒪_1]_{l_1l_1}[𝒪_2]_{l_2m_2}[𝒪_3]_{l_3m_3}⟨l_1m_1,l_2m_2|l_{12}m_{12}⟩⟨l_{12}m_{12},l_3m_3|lm⟩
\\end{aligned}
```

# Fields

* `amd :: Vector{AngModes}` records, for each part ``p``, the spherical-symmetric action on that part, stored as an `AngModes` object.
* `ch :: Vector{Matrix{Int64}}` is the list of coupling channels. 
* `coeff :: Vector{ComplexF64}` is the coefficient of each channel.
* `sec :: Matrix{Int64}` records the change of quantum numbers that the term induces : `sec[iqn, p]` is the shift of the `iqn`-th diagonal quantum number on part ``p``.
"""
mutable struct CoupleDecomp
    amd :: Vector{AngModes}
    ch :: Vector{Matrix{Int64}}
    coeff :: Vector{ComplexF64}
    sec :: Matrix{Int64} # sec[iqn, p]
end


"""
    CoupleDecomp(amd :: Vector{AngModes}, ch :: Matrix{Int64}, sec :: Matrix{Int64}) :: CoupleDecomp

constructs a single-channel `CoupleDecomp` with unit coefficient from the operators `amd`, a single coupling channel `ch` and the quantum number shift `sec`.
"""
function CoupleDecomp(amd :: Vector{AngModes}, ch :: Matrix{Int64}, sec :: Matrix{Int64})
    return CoupleDecomp(amd, [ch], [1.0], sec)
end

"""
    cpd1 + cpd2 :: Vector{CoupleDecomp}
    cpd1 - cpd2 :: Vector{CoupleDecomp}
    -cpd :: Vector{CoupleDecomp}
    fac * cpd :: Vector{CoupleDecomp}

enable the linear combination of coupling decompositions. 
"""
function Base.:+(cpd1 :: Union{CoupleDecomp, Vector{CoupleDecomp}}, cpd2 :: Union{CoupleDecomp, Vector{CoupleDecomp}})
    return [ cpd1 ; cpd2 ]
end

function Base.:*(fac :: Number, cpd :: CoupleDecomp)
    return CoupleDecomp(cpd.amd, cpd.ch, fac .* cpd.coeff, cpd.sec)
end

function Base.:*(fac :: Number, cpd :: Vector{CoupleDecomp})
    return fac .* cpd
end

function Base.:-(cpd :: Union{CoupleDecomp, Vector{CoupleDecomp}})
    return (-1) * cpd
end

function Base.:-(cpd1 :: Union{CoupleDecomp, Vector{CoupleDecomp}}, cpd2 :: Union{CoupleDecomp, Vector{CoupleDecomp}})
    return cpd1 + (-1) * cpd2
end

"""
    RecouplePsPot(s1 :: Number, s2 :: Number, s3 :: Number, s4 :: Number, ps_pot0 :: Dict) :: Dict
    RecouplePsPot(s :: Number, ps_pot :: Vector{<:Number}) :: Dict

recouples the Haldane pseudopotentials of a two-body interaction ``c^†_1c^†_2c_3c_4`` from the pairing channel ``(12)(34)`` — creations and annihilations each coupled to a pair angular momentum ``j`` — into the density channel ``(14)(23)`` — density modes ``n^{(14)}=c^†_1c_4`` and ``n^{(23)}=c^†_2c_3`` each coupled to a rank ``l``. This is the form consumed by [ConvPsPot](@ref ConvPsPot) and expressed as a product of two density operators. 

# Arguments

* `s1, s2, s3, s4 :: Number` are the single-particle angular momenta of the four operators.
* `ps_pot0 :: Dict` maps a pair angular momentum ``j`` to the pseudopotential ``V_j``.

In the second form the interaction is a single flavour with spin `s`, and `ps_pot :: Vector` lists the pseudopotentials ordered from the largest pair angular momentum ``j=2s`` downwards, _i. e._ `ps_pot[i]` is ``V_{2s+1-i}``.

# Output

* `ps_pot1 :: Dict` maps a density rank ``l`` to the recoupled coefficient ``W_l``, ready to be passed to [ConvPsPot](@ref ConvPsPot).
"""
function RecouplePsPot(s1 :: Number, s2 :: Number, s3 :: Number, s4 :: Number, ps_pot0 :: Dict)
    Lmin = max(abs(s1 - s4), abs(s2 - s3))
    Lmax = min(s1 + s4, s2 + s3)
    ps_pot1 = Dict{Float64, ComplexF64}()
    for L in Lmin : Lmax
        W = 0.0
        for (J, V) in ps_pot0
            (abs(s1 - s2) ≤ J ≤ s1 + s2) || continue
            (abs(s3 - s4) ≤ J ≤ s3 + s4) || continue
            sixj = Float64(wigner6j(s1, s2, J, s3, s4, L))
            sixj == 0.0 && continue
            #iseven(L) || (sixj = -sixj)
            sixj *= (1.0im) ^ Int(-2J)
            W += (2J + 1) * sixj * V
        end
        ps_pot1[L] = W
    end
    return ps_pot1
end
function RecouplePsPot(s :: Number, ps_pot :: Vector{<:Number}) 
    ps_pot0 = Dict([ 2s + 1 - i => ps_pot[i] for i ∈ eachindex(ps_pot)])
    return RecouplePsPot(s, s, s, s, ps_pot0)
end

"""
    ConvPsPot(ps_pot0 :: Dict) :: Tuple{Vector{Matrix{Int64}}, Vector{ComplexF64}}

converts a density-channel pseudopotential — a dictionary mapping a density rank ``l`` to its coefficient ``V_l`` — into the coupling channels and coefficients used by a two-part density–density [CoupleDecomp](@ref CoupleDecomp). For each rank ``l`` it produces the channel ``\\begin{smallmatrix}2l&2l\\\\2l&0\\end{smallmatrix}`` (two rank-``l`` operators coupled to a total scalar) with coefficient ``W_l√{2l+1}(-1)^{l}``.

# Output

* `ch :: Vector{Matrix{Int64}}` is the list of coupling channels.
* `coeff :: Vector{ComplexF64}` is the list of channel coefficients.

These are typically splatted into a [CoupleDecomp](@ref CoupleDecomp) constructor, _e. g._ `CoupleDecomp([n_mod, n_mod], ConvPsPot(ps_pot0)..., sec)`.
"""
function ConvPsPot(ps_pot0 :: Dict)
    ch = Matrix{Int64}[]
    coeff = ComplexF64[]
    for (J, V) in ps_pot0
        J2 = Int64(2 * J)
        push!(ch, [J2  J2 ; J2  0])
        push!(coeff, V * √(J2 + 1) * (1.0im) ^ (J2))
    end
    return ch, coeff
end


function FuzzifiED.AngModes(obs :: SphereObs)
    return AngModes(obs.l2m, obs.get_comp)
end


"""
    ContactCouple(obs :: Vector{SphereObs}, sec :: Matrix{Int64}, ltot :: Int64) :: CoupleDecomp

constructs a [CoupleDecomp](@ref CoupleDecomp) for a contact term, _i. e._ the product of one spherical observable per part evaluated at the same point on the sphere
```math 
    ∫\\mathrm{d}^2𝐫\\,√{4π}Ȳ_{lm}(𝐫)\\,n_1(𝐫)n_2(𝐫)⋯n_{N_p}(𝐫)
```

# Arguments

* `obs :: Vector{SphereObs}` is the list of spherical observables, one per part.
* `sec :: Matrix{Int64}` records the change of quantum numbers, `sec[iqn, p]` for the `iqn`-th quantum number on part ``p``.
* `ltot :: Int64` is twice the total angular momentum ``2l_{\\text{tot}}`` of the term. Facultative, ``0`` (a scalar) by default.

# Output

* `cpd :: CoupleDecomp` is the resulting coupling decomposition.
"""
function ContactCouple(obs :: Vector{SphereObs}, sec :: Matrix{Int64}, ltot :: Int64 = 0)
    amd = AngModes.(obs)
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
        coeffi *= √(prod(chi[1, :] .+ 1) / (ltot + 1)) * FuzzifiED.ObsNormRadSq
        push!(ch, chi)
        push!(coeff, coeffi)
    end
    return CoupleDecomp(amd, ch, coeff, sec)
end

"""
    SingleSegCouple(np :: Int64, p :: Int64, amdp :: AngModes, l :: Int64, secp :: Vector{Int64}) :: CoupleDecomp
    SingleSegCouple(np :: Int64, p :: Int64, tms :: Terms, sec :: Vector{Int64}) :: CoupleDecomp

constructs a [CoupleDecomp](@ref CoupleDecomp) for a term that acts as a nontrivial operator on a single part ``p`` and as the identity on all the other parts. This is used, _e. g._, for a single-part chemical potential or on-site interaction.

# Arguments

* `np :: Int64` is the number of parts.
* `p :: Int64` is the index of the part on which the operator acts.
* `amdp :: AngModes` is the spherical tensor operator acting on part ``p``.
* `l :: Int64` is twice the rank ``2l`` of the operator on part ``p``.
* `secp :: Vector{Int64}` is the change of quantum numbers on part ``p``.

In the second form the operator is a scalar (rank ``0``) given directly as a list of terms `tms :: Terms`, and `sec` is its quantum number shift.

# Output

* `cpd :: CoupleDecomp` is the resulting coupling decomposition.
"""
function SingleSegCouple(np :: Int64, p :: Int64, amdp :: AngModes, l :: Int64, secp :: Vector{Int64})
    amd1 = AngModes(0, Dict((0, 0) => one(Terms)))
    amd = [amd1 for p = 1 : np]
    amd[p] = amdp
    ch = zeros(Int64, 2, np)
    ch[1, p] = l
    ch[2, p : end] .= l 
    sec = zeros(Int64, length(secp), np)
    sec[:, p] = secp
    return CoupleDecomp(amd, [ch], [1], sec)
end
function SingleSegCouple(np :: Int64, p :: Int64, tms :: Terms, sec :: Vector{Int64})
    amdp = AngModes(0, Dict((0, 0) => tms))
    return SingleSegCouple(np, p, amdp, 0, sec)
end

"""
    PrepareCouple(cpd :: CoupleDecomp ; eltype :: Type) :: CoupleDecomp
    PrepareCouple(cpd :: Vector{CoupleDecomp} ; eltype :: Type) :: Vector{CoupleDecomp}

prepares a coupling decomposition for the construction of operators by pre-storing the components of each `AngModes` and, when `eltype == Float64`, rotating any purely imaginary operator by ``i`` while compensating the phase in the coefficients, so that all the matrix elements can be represented with real numbers. This should be called once on an assembled operator before it is passed to [BuildSegOperators](@ref BuildSegOperators).

# Arguments

* `cpd :: CoupleDecomp` or `Vector{CoupleDecomp}` is the coupling decomposition to be prepared.
* `eltype :: Type` is the target matrix element type, either `Float64` or `ComplexF64`. Facultative, `FuzzifiED.ElementType` by default.

# Output

* the prepared coupling decomposition, of the same shape as the input.
"""
function PrepareCouple(cpd :: CoupleDecomp ; eltype = FuzzifiED.ElementType)
    amd1 = AngModes[]
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
    return CoupleDecomp(amd1, cpd.ch, cpd.coeff .* ph, cpd.sec)
end
PrepareCouple(cpd :: Vector{CoupleDecomp} ; eltype = FuzzifiED.ElementType) = PrepareCouple.(cpd ; eltype)