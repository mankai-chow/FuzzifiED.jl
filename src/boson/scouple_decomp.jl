export SCoupleDecomp
export ContactSCouple, SingleSegSCouple

mutable struct SCoupleDecomp
    amd :: Vector{SAngModes}
    ch :: Vector{Matrix{Int64}}
    coeff :: Vector{ComplexF64}
    sec :: Matrix{Int64} # sec[iqn, p]
end

function SCoupleDecomp(amd :: Vector{SAngModes}, ch :: Matrix{Int64}, sec :: Matrix{Int64})
    return SCoupleDecomp(amd, [ch], [1.0], sec)
end

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

function Fuzzifino.SAngModes(obs :: SSphereObs)
    return SAngModes(obs.l2m, obs.get_comp)
end

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
        coeffi *= √(prod(chi[1, :] .+ 1) / (ltot + 1)) * FuzzifiED.ObsNormRadSq
        push!(ch, chi)
        push!(coeff, coeffi)
    end
    return SCoupleDecomp(amd, ch, coeff, sec)
end

function SingleSegSCouple(np :: Int64, p :: Int64, amdp :: SAngModes, l :: Int64, secp :: Vector{Int64})
    amd1 = SAngModes(0, Dict((0, 0) => one(STerms)))
    amd = [amd1 for p = 1 : np]
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