# This example calculates the spectrum of the SO(5) deconfined criticality. 
# It uses a bi-partition into two segments each of two flavours. 
# The implemented symmetry is SU(2)×SU(2) ⊃ SO(5) and the SO(5) representation 
# is resolved by measuring flavour Casimir.

using FuzzifiED
using FuzzifiEDFullRotation
using LinearAlgebra
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 7
ne = nm * 2
nf = 4
s = (nm - 1) / 2
qnd_pt = [
    GetNeQNDiag(nm * 2),
    GetLz2QNDiag(nm, 2),
    GetFlavQNDiag(nm, 2, [1, -1])
]
tms_lzlp_pt = GetLpLzTerms(nm, 2)
tms_c2 = GetC2Terms(nm, 2, :SU)
sec_pt = stack([[ne1, ((nm + 1) * ne1) % 2, 0] for ne1 = 0 : ne])
sec_tot = [ne, 0, 0]

Δ = GetPairingMod(nm, 2, [0 1;0 0])
n = GetDensityMod(nm, 2, [1 0;0 1])
cpd_nn = SingleSegCouple(2, 1, GetDenIntTerms(nm, 2), [0, 0, 0]) + 
    SingleSegCouple(2, 2, GetDenIntTerms(nm, 2), [0, 0, 0]) + 
    2 * CoupleDecomps([n, n], ConvPsPot(RecoupleAngMom(s, [1]))..., [0 0;0 0;0 0])
cpd_ΔΔ = SingleSegCouple(2, 1, GetPairIntTerms(nm, 2, [0 1;0 0]), [0, 0, 0]) + 
    SingleSegCouple(2, 2, GetPairIntTerms(nm, 2, [0 1;0 0]), [0, 0, 0]) + 
    CoupleDecomps([Δ', Δ], ConvPsPot(s, [1])..., [2 -2;0 0;0 0]) + 
    CoupleDecomps([Δ, Δ'], ConvPsPot(s, [1])..., [-2 2;0 0;0 0])
cpd_hmt = cpd_nn - 0.9 * cpd_ΔΔ

tms_c2_sg = GetC2Terms(nm, 2, :SU) + (nf - 2) / 4 * GetPolTerms(nm, 2)
cpd_c2 = SingleSegCouple(2, 1, tms_c2_sg, [0, 0, 0]) + SingleSegCouple(2, 2, tms_c2_sg, [0, 0, 0]) + 
    CoupleDecomps([n, n], ConvPsPot(Dict([2s - l => -0.5 for l = 0 : 2s]))..., [0 0;0 0;0 0]) +
    CoupleDecomps([Δ', Δ], ConvPsPot(Dict([2s - l => -1 for l = 0 : 2 : 2s]))..., [2 -2;0 0;0 0]) + 
    CoupleDecomps([Δ, Δ'], ConvPsPot(Dict([2s - l => -1 for l = 0 : 2 : 2s]))..., [-2 2;0 0;0 0])

s_rng = collect(0 : 2)
c2_rng = [Float64[s * (s + 1)] for s in s_rng]
sgsp = BuildSegSpaces(nm * 2, sec_pt, qnd_pt, tms_lzlp_pt, tms_c2, c2_rng)
result = []
for is = 1 : 3, js = is : 3
    si = s_rng[is]
    sj = s_rng[js]
    sgop_hmt = BuildSegOperators([sgsp[is], sgsp[js]], cpd_hmt)
    sgop_c2 = BuildSegOperators([sgsp[is], sgsp[js]], cpd_c2)
    for l = 0 : 2
        cpsp = BuildCompSpace([sgsp[is], sgsp[js]], sec_tot, 2l)
        cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
        cpop_c2 = BuildCompOperator(cpsp, cpd_c2, sgop_c2)
        enrg, st = GetEigensystem(cpop_hmt, 10 ; issymmetric = true)
        for i in eachindex(enrg)
            c2_val = st[:, i]' * (cpop_c2 * st[:, i]) / √(2l + 1)
            push!(result, [enrg[i] / √(2l + 1), l, c2_val, si, sj])
        end
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2 && st[3] ≈ 0, result)[1][1]
spec = [ round.([ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(stack(spec)))
