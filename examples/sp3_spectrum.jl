# This example calculates the spectrum of the Sp(3) CFT. 
# It uses a tri-partition, each with two flavours and a SU(2) symmetry.

using FuzzifiED
using FuzzifiEDFullRotation
using LinearAlgebra
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 5
ne = nm * 2
nf = 6
s = (nm - 1) / 2
qnd_pt = [
    GetNeQNDiag(nm * 2),
    GetLz2QNDiag(nm, 2),
    GetFlavQNDiag(nm, 2, [1, -1])
]
tms_lzlp_pt = GetLpLzTerms(nm, 2)
tms_c2 = GetC2Terms(nm, 2, :SU)
sec_pt = stack([[ne1, ((nm + 1) * ne1) % 2, ne1 % 2] for ne1 = 0 : ne])
sec_tot = [ne, 0, 0]

Δ = GetPairingMod(nm, 2, [0 1;0 0])
n = GetDensityMod(nm, 2, [1 0;0 1])
tms_intra = GetDenIntTerms(nm, 2, [1.0, 0.2643, 0.0652]) - GetPairIntTerms(nm, 2, [0.3798, 0.0, -0.0219], [0 1;0 0])
cpd_nn_inter = 2 * CoupleDecomps([n, n], ConvPsPot(RecouplePsPot(s, [1.0, 0.2643, 0.0652]))..., [0 0;0 0;0 0])
cpd_ΔΔ_inter = CoupleDecomps([Δ', Δ], ConvPsPot(Dict(2s => 0.3798, 2s-2 => -0.0219))..., [2 -2;0 0;0 0]) + 
    CoupleDecomps([Δ, Δ'], ConvPsPot(Dict(2s => 0.3798, 2s-2 => -0.0219))..., [-2 2;0 0;0 0])
cpd_hmt = PrepareCouple(sum([SingleSegCouple(3, p, tms_intra, [0, 0, 0]) for p = 1 : 3]) +
    sum([InsertSegment(3, [p1, p2], cpd_nn_inter) for p1 = 1 : 3 for p2 = p1 + 1 : 3]) -
    sum([InsertSegment(3, [p1, p2], cpd_ΔΔ_inter) for p1 = 1 : 3 for p2 = p1 + 1 : 3])) ;

tms_c2_sg = GetC2Terms(nm, 2, :SU) + (nf - 2) / 4 * GetPolTerms(nm, 2)
cpd_c2_inter = CoupleDecomps([n, n], ConvPsPot(Dict([2s - l => -0.5 for l = 0 : 2s]))..., [0 0;0 0;0 0]) +
    CoupleDecomps([Δ', Δ], ConvPsPot(Dict([2s - l => -1 for l = 0 : 2 : 2s]))..., [2 -2;0 0;0 0]) + 
    CoupleDecomps([Δ, Δ'], ConvPsPot(Dict([2s - l => -1 for l = 0 : 2 : 2s]))..., [-2 2;0 0;0 0])
cpd_c2 = PrepareCouple(sum([SingleSegCouple(3, p, tms_c2_sg, [0, 0, 0]) for p = 1 : 3]) +
    sum([InsertSegment(3, [p1, p2], cpd_c2_inter) for p1 = 1 : 3 for p2 = p1 + 1 : 3]))

s_rng = collect(0 : 2)
c2_rng = [Float64[s * (s + 1)] for s in s_rng]
sgsp = BuildSegSpaces(nm * 2, sec_pt, qnd_pt, tms_lzlp_pt, tms_c2, c2_rng)
result = []
for is in [[1,1,1], [1,1,2], [1,1,3], [1,2,2], [1,2,3], [2,2,2]]
    si = s_rng[is]
    sgop_hmt = BuildSegOperators(sgsp[is], cpd_hmt)
    sgop_c2 = BuildSegOperators(sgsp[is], cpd_c2)
    for l = 0 : 2
        cpsp = BuildCompSpace(sgsp[is], sec_tot, 2l)
        cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
        cpop_c2 = BuildCompOperator(cpsp, cpd_c2, sgop_c2)
        nst = 10 - sum(si) - l
        enrg, st = GetEigensystem(cpop_hmt, nst ; issymmetric = true)
        for i in eachindex(enrg)
            c2_val = st[:, i]' * (cpop_c2 * st[:, i]) / √(2l + 1)
            push!(result, [enrg[i] / √(2l + 1), l, c2_val])
        end
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2 && st[3] ≈ 0, result)[1][1]
spec = unique([ round.([ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ])
display(permutedims(stack(spec)))
