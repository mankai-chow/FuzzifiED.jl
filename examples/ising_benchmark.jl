using FuzzifiED
using FuzzifiEDFullRotation
using BenchmarkTools
FuzzifiED.ElementType = Float64

nm = 16
nf = 2
ne = nm
s = (nm-1)/2

qnd_pt = [ 
    GetNeQNDiag(nm), 
    GetLz2QNDiag(nm, 1)
] 
sec_pt = stack([ [ne1, 0] for ne1 = 0 : 2 : ne])
sec_tot = [ne, 0]

c = GetElectronMod(nm, 1, 1)
cpd_0 = CoupleDecomps([ c' * c, c' * c ], ConvPsPot(RecouplePsPot(s, [4.75 * 2]))..., zeros(Int64, 2, 2))
cpd_1 = CoupleDecomps([ c' * c', c * c], ConvPsPot(Dict(2s - 1 => -1/2))..., [2 -2;0 0]) + 
    CoupleDecomps([ c * c, c' * c'], ConvPsPot(Dict(2s - 1 => -1/2))..., [-2 2;0 0]) +
    SingleSegCouple(2, 1, GetDenIntTerms(nm, 1, [0, 1/2]), zeros(Int64, 2)) + 
    SingleSegCouple(2, 2, GetDenIntTerms(nm, 1, [0, 1/2]), zeros(Int64, 2))
cpd_h = SingleSegCouple(2, 1, GetPolTerms(nm, 1), zeros(Int64, 2)) - SingleSegCouple(2, 2, GetPolTerms(nm, 1), zeros(Int64, 2))
cpd_hmt = PrepareCouple(cpd_0 + cpd_1 - 3.16 * cpd_h) 
tms_lzlp_pt = GetLpLzTerms(nm, 1) 

b = @benchmark begin

sgsp = BuildSegSpace(nm, sec_pt, qnd_pt, tms_lzlp_pt)
sgop_hmt = Matrix{SegOperator}(undef, 2, length(cpd_hmt))
sgop_hmt[1, :] = BuildSegOperators([sgsp], cpd_hmt ; p_rng = [1])
sgop_hmt[2, 1 : nm] = sgop_hmt[1, 1 : nm]
sgop_hmt[2, nm + 1 : 2 : end] = sgop_hmt[1, nm + 2 : 2 : end]
sgop_hmt[2, nm + 2 : 2 : end] = sgop_hmt[1, nm + 1 : 2 : end]

cpsp = BuildCompSpace([sgsp, sgsp], sec_tot, 0)
cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
enrg, st = GetEigensystem(cpop_hmt, 10 ; issymmetric = true)
display(enrg)

end

display(b)

