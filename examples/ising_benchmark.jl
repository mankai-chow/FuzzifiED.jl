# This example provides a benchmark for the package at its best performance. 
# It calculates the 10 lowest Z_2-even states with l = 0.
# Currently, it only supports even nm. 

using FuzzifiED
using FuzzifiEDFullRotation
using BenchmarkTools
FuzzifiED.ElementType = Float64

nm = 16
nf = 2
ne = nm
s = (nm-1)/2
l = 0 
Z = 0

qnd_pt = [ 
    GetNeQNDiag(nm), 
    GetLz2QNDiag(nm, 1)
] 
sec_pt = stack([ [ne1, Z] for ne1 = Z : 2 : ne])
sec_tot = [ne, 0]

c = GetElectronMod(nm, 1, 1)
cpd_int_0 = CoupleDecomps([ c' * c, c' * c ], ConvPsPot(RecoupleAngMom(s, [2.0]))..., zeros(Int64, 2, 2))
cpd_int_1 = CoupleDecomps([ c' * c', c * c], ConvPsPot(s, [0, -1/2])..., [2 -2;0 0]) + 
    CoupleDecomps([ c * c, c' * c'], ConvPsPot(s, [0, -1/2])..., [-2 2;0 0]) +
    SingleSegCouple(2, 1, GetDenIntTerms(nm, 1, [0, 1/2]), zeros(Int64, 2)) + 
    SingleSegCouple(2, 2, GetDenIntTerms(nm, 1, [0, 1/2]), zeros(Int64, 2))
cpd_h = SingleSegCouple(2, 1, GetPolTerms(nm, 1), zeros(Int64, 2)) - SingleSegCouple(2, 2, GetPolTerms(nm, 1), zeros(Int64, 2))
cpd_hmt = 4.75 * cpd_int_0 + cpd_int_1 - 3.16 * cpd_h
tms_lzlp_pt = GetLpLzTerms(nm, 1) 

b = @benchmark begin

sgsp = BuildSegSpace(nm, sec_pt, qnd_pt, tms_lzlp_pt)
sgop_hmt = BuildSegOperators([sgsp, sgsp], cpd_hmt ; ident_seg = [1, 1])
cpsp = BuildCompSpace([sgsp, sgsp], sec_tot, 2l)
cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
enrg, st = GetEigensystem(cpop_hmt, 10 ; issymmetric = true)
display(enrg)

end

display(b)

