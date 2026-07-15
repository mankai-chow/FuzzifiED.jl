export RecouplePsPot, ConvPsPot

function RecouplePsPot(s1 :: Number, s2 :: Number, s3 :: Number, s4 :: Number, ps_pot0 :: Dict)
    Lmin = max(abs(s1 - s3), abs(s2 - s4))
    Lmax = min(s1 + s3, s2 + s4)
    ps_pot1 = Dict{Float64, ComplexF64}()
    for L in Lmin : Lmax
        W = 0.0
        for (J, V) in ps_pot0
            (abs(s1 - s2) ≤ J ≤ s1 + s2) || continue
            (abs(s3 - s4) ≤ J ≤ s3 + s4) || continue
            sixj = Float64(wigner6j(s1, s2, J, s4, s3, L))
            sixj == 0.0 && continue
            iseven(L) || (sixj = -sixj)
            #iseven(round(Int, s2 + s3 + J + L)) || (sixj = -sixj)
            W += (2J + 1) * sixj * V
        end
        ps_pot1[L] = W
    end
    return ps_pot1
end
function RecouplePsPot(s :: Number, ps_pot :: Vector{<:Number}) 
    ps_pot0 = Dict([ 2s + 1 - i => ps_pot[i] * (-1) ^ (i - 1) for i ∈ eachindex(ps_pot)])
    return RecouplePsPot(s, s, s, s, ps_pot0)
end

function ConvPsPot(ps_pot0 :: Dict)
    ch = Matrix{Int64}[]
    coeff = ComplexF64[]
    for (J, V) in ps_pot0
        J2 = Int64(2 * J)
        push!(ch, [2J  2J ; 2J  0])
        push!(coeff, V * √(J2 + 1))
    end
    return ch, coeff
end