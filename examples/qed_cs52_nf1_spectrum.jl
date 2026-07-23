# This example gives the model for QED-Chern-Simons of U(1)_{5/2} with one flavour of Dirac fermion.

using FuzzifiED
using FuzzifiED.Fuzzifino
using FuzzifiEDFullRotation
FuzzifiED.ElementType = Float64

nmc = 5
s = (nmc - 1) / 2
nmf = 3 * nmc - 2
nmb = 2 * nmc - 1
FuzzifiED.ObsNormRadSq = nmc 

qnd_c = [
    SQNDiag(GetNeQNDiag(nmc), 1),
    SQNDiag(GetLz2QNDiag(nmc, 1), 1), 
    zero(SQNDiag, nmc, 1),
    GetBosonNeSQNDiag(nmc, 1)
]
qnd_b = [
    2 * GetBosonNeSQNDiag(1, nmb),
    GetBosonLz2SQNDiag(1, nmb, 1), 
    GetBosonNeSQNDiag(1, nmb),
    SQNDiag(GetNeQNDiag(1), nmb)
]
qnd_f = [
    3 * SQNDiag(GetNeQNDiag(nmf), 1),
    SQNDiag(GetLz2QNDiag(nmf, 1), 1), 
    SQNDiag(GetNeQNDiag(nmf), 1),
    GetBosonNeSQNDiag(nmf, 1)
]

sec_c = stack([ne1, ((nmc + 1) * ne1) % 2, 0, 0] for ne1 = 0 : nmc)
sec_b = stack([2ne1, 0, ne1, 0] for ne1 = 0 : nmc)
sec_f = stack([3ne1, ((nmf + 1) * ne1) % 2, ne1, 0] for ne1 = 0 : nmc) 
nebm_c = [0 for ne1 = 0 : nmc]
nebm_b = [ne1 for ne1 = 0 : nmc]
nebm_f = [0 for ne1 = 0 : nmc]

tms_lzlp_c = STerms.(GetLpLzTerms(nmc, 1))
tms_lzlp_b = GetBosonLpLzSTerms(nmb, 1)
tms_lzlp_f = STerms.(GetLpLzTerms(nmf, 1))
tms_proj_b = GetBosonDenIntSTerms(nmb, 1)
tms_proj_f = STerms(GetDenIntTerms(nmf, 1, [0, 1]))

sgsp_c = BuildSSegSpace(nmc, 1, nebm_c, sec_c, qnd_c, tms_lzlp_c)
sgsp_b = BuildSSegSpace(1, nmb, nebm_b, sec_b, qnd_b, tms_lzlp_b, tms_proj_b)
sgsp_f = BuildSSegSpace(nmf, 1, nebm_f, sec_f, qnd_f, tms_lzlp_f, tms_proj_f)


c = GetFermionSObs(nmc, 1, 1)
f = GetFermionSObs(nmf, 1, 1)
b = GetBosonSObs(nmb, 1, 1)
nc = c' * c 
nf = f' * f 
nb = b' * b

tms_int_b = SimplifyTerms(GetIntegral(nb * nb))
tms_int_1 = SimplifyTerms(GetIntegral(nf * Laplacian(nf)))
tms_pol_f = SimplifyTerms(GetIntegral(nf))

cpd_pol_f = SingleSegSCouple(3, 3, tms_pol_f, zeros(Int64, 4))
cpd_hop = ContactSCouple([c', b', f], [ 1 2 -3 ; 0 0 0 ; 0 1 -1 ; 0 0 0]) -
          ContactSCouple([c, b, f'],  [-1 -2 3 ; 0 0 0 ; 0 -1 1 ; 0 0 0])
cpd_int_e = 4 * ContactSCouple([nc, nb, one(SSphereObs)], zeros(Int64, 4, 3)) +
            6 * ContactSCouple([nc, one(SSphereObs), nf ], zeros(Int64, 4, 3)) +
            12 * ContactSCouple([one(SSphereObs), nb, nf ], zeros(Int64, 4, 3)); 
cpd_hmt = PrepareCouple(cpd_int_e - 0.2 * cpd_hop + 0.2 * cpd_pol_f) # An example of Hamiltonian, not necessarily conformal

result = []
for q = -1 : 1, l = 0 : 2
    sec_tot = [3nmc, 0, nmc + q, 0]
    cpsp = BuildSCompSpace([sgsp_c, sgsp_b, sgsp_f], sec_tot, 2l)
    sgop_hmt = BuildSSegOperators([sgsp_c, sgsp_b, sgsp_f], cpd_hmt)
    cpop_hmt = BuildSCompOperator(cpsp, cpd_hmt, sgop_hmt)
    nst = 10 - abs(q) - l
    enrg, st = GetEigensystem(cpop_hmt, nst ; issymmetric = true)
    for i in eachindex(enrg)
        push!(result, [enrg[i] / √(2l + 1), l, q])
    end
end
 
sort!(result, by = st -> real(st[1])) 
enrg_0 =  filter(st -> st[3] ≈ 0, result)[1][1]
enrg_q = (filter(st -> st[3] ≈ -1, result)[1][1] - filter(st -> st[3] ≈ 1, result)[1][1]) / 2
enrg_1 = (filter(st -> st[2] ≈ 1 && st[3] ≈ 0, result)[1][1] - enrg_0) / 2
spec = [ [ (st[1] - enrg_0 + st[3] * enrg_q) / enrg_1 ; st] for st in result]
sort!(spec, by = st -> real(st[1])) 
display(permutedims(hcat(spec...)))
