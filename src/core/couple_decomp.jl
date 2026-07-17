export CoupleDecomp
export RecouplePsPot, ConvPsPot, ContactCouple, SingleSegCouple
export PrepareCouple

mutable struct CoupleDecomp
    amd :: Vector{AngModes}
    ch :: Vector{Matrix{Int64}}
    coeff :: Vector{ComplexF64}
    sec :: Matrix{Int64} # sec[iqn, p]
end

function CoupleDecomp(amd :: Vector{AngModes}, ch :: Matrix{Int64}, sec :: Matrix{Int64})
    return CoupleDecomp(amd, [ch], [1.0], sec)
end

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