# This example calculates the spectrum of the 3D Ising model on the fuzzy sphere 
# at N_m = 12 in the sectors with total angular momentum l = 0, 1, and 2.
# It resolves the Z_2 symmetry as a QNDiag by using XX - Z basis. 

using FuzzifiED
using FuzzifiEDFullRotation
FuzzifiED.ElementType = Float64

nm = 12
ne = nm
s = (nm - 1) / 2

qnd_pt = [ 
    GetNeQNDiag(nm)  GetNeQNDiag(nm) ;
    GetLz2QNDiag(nm, 1)  GetLz2QNDiag(nm, 1) ;
    zero(QNDiag, nm, 2)  GetNeQNDiag(nm, 2)
] 
sec_pt = stack([ [ne1  ne1 ; ((nm + 1) * ne1) % 2  ((nm + 1) * ne1) % 2 ; 0  ne1 % 2 ] for ne1 = 0 : ne])
tms_lzlp_pt = GetLpLzTerms(nm, 1) 

c = GetElectronMod(nm, 1, 1)
cpd_int_0 = CoupleDecomps([ c' * c, c' * c ], ConvPsPot(RecouplePsPot(s, [2.0]))..., zeros(Int64, 3, 2))
cpd_int_1 = CoupleDecomps([ c' * c', c * c], ConvPsPot(Dict(2s - 1 => -1/2))..., [2 -2;0 0;0 0]) + 
    CoupleDecomps([ c * c, c' * c'], ConvPsPot(Dict(2s - 1 => -1/2))..., [-2 2;0 0;0 0]) +
    SingleSegCouple(2, 1, GetDenIntTerms(nm, 1, [0, 1/2]), zeros(Int64, 3)) + 
    SingleSegCouple(2, 2, GetDenIntTerms(nm, 1, [0, 1/2]), zeros(Int64, 3))
cpd_h = SingleSegCouple(2, 1, GetPolTerms(nm, 1), zeros(Int64, 3)) - SingleSegCouple(2, 2, GetPolTerms(nm, 1), zeros(Int64, 3))
cpd_hmt = PrepareCouple(4.75 * cpd_int_0 + cpd_int_1 - 3.16 * cpd_h) ;

sgsp = [ BuildSegSpace(nm, sec_pt[:, p, :], qnd_pt[:, p], tms_lzlp_pt, [1, 1, 2]) for p = 1 : 2]
sgop_hmt = BuildSegOperators(sgsp, cpd_hmt)

result = []
for l in 0 : 2, Z = 0 : 1 
    sec_tot = [ne, 0, Z]
    cpsp = BuildCompSpace(sgsp, sec_tot, 2l)
    cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
    nst = 10 - l
    enrg, st = GetEigensystem(cpop_hmt, nst ; issymmetric = true)
    for i in eachindex(enrg)
        push!(result, [enrg[i] / √(2l + 1), l, Z])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2 && st[3] ≈ 0, result)[1][1]
spec = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(stack(spec)))
