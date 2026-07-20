using FuzzifiED
using FuzzifiED.Fuzzifino
using FuzzifiEDFullRotation
FuzzifiED.ElementType = Float64

nmf = 7
nf = 2
nof = nmf * nf
nmb = 2 * nmf - 1
s = (nmf - 1) / 2
FuzzifiED.ObsNormRadSq = nmf

qnd_f = [ 
    SQNDiag(GetNeQNDiag(nof), 1),
    SQNDiag(GetLz2QNDiag(nmf, nf), 1),
    SQNDiag(GetFlavQNDiag(nmf, nf, [1, -1]), 1) , 
    GetBosonNeSQNDiag(nof, 1) 
]
qnd_b = [
    2 * GetBosonNeSQNDiag(1, nmb),
    GetBosonLz2SQNDiag(1, nmb, 1), 
    zero(SQNDiag, 1, nmb),
    SQNDiag(GetNeQNDiag(1), nmb) 
]

sec_f = stack([ [ne1, 0, 0, 0] for ne1 = 0 : 2 : nof ])
sec_b = stack([ [ne1, 0, 0, 0] for ne1 = 0 : 2 : nof ])
sec_tot = [nof, 0, 0, 0]
nebm_f = [0   for ne1 = 0 : 2 : nof ]
nebm_b = [ne1 for ne1 = 0 : 2 : nof ]
tms_lzlp_f = STerms.(GetLpLzTerms(nmf, nf))
tms_lzlp_b = GetBosonLpLzSTerms(nmb, 1)
tms_c2 = STerms(GetC2Terms(nmf, nf, :SU))
tms_proj = GetBosonDenIntSTerms(nmb, 1)

b = GetBosonSObs(nmb, 1, 1)
ff = GetFermionSObs(nmf, 2, 1) * GetFermionSObs(nmf, 2, 2)
nb = GetBosDensitySObs(nmb, 1)
nf = GetFerDensitySObs(nmf, 2)

cpd_hop = ContactSCouple([ff', b], [2 -2; 0 0; 0 0; 0 0]) + ContactSCouple([ff, b'], [-2 2; 0 0; 0 0; 0 0])
cpd_μ = SingleSegSCouple(2, 1, STerms(GetPolTerms(nmf, 2)), [0, 0, 0, 0])
cpd_int_e = SingleSegSCouple(2, 1, GetIntegral(nf * nf), [0, 0, 0, 0]) + 4 * SingleSegSCouple(2, 2, GetIntegral(nb * nb), [0, 0, 0, 0]) + 4 * ContactSCouple([nf, nb], [0 0; 0 0; 0 0; 0 0])

cpd_hmt = cpd_int_e - 0.5 * cpd_hop + 0.312 * cpd_μ 

nst_max = [ zeros(Int64, size(sec_b, 2) - 3) ; 5 .* nmf .^ [2,1,0] ]
sgsp_b = BuildSSegSpace(1, nmb, nebm_b, sec_b, qnd_b, tms_lzlp_b, tms_proj, [0.0] ; l2c2_ratio = 0.1/nmf^2, nst_max)
ss = collect(0 : 2)
c2_rng = [ Float64[s * (s + 1)] for s in ss]
sgsps_f = BuildSSegSpaces(nof, 1, nebm_f, sec_f, qnd_f, tms_lzlp_f, tms_c2, c2_rng)
sgop_hmt = Matrix{SSegOperator}(undef, 2, CountChannels(cpd_hmt))
sgop_hmt[2, :] = BuildSSegOperators([sgsp_b], cpd_hmt ; p_rng = [2])

result = []
for is in eachindex(ss)
    s = ss[is]
    sgsp_f = sgsps_f[is]
    sgop_hmt[1, :] = BuildSSegOperators([sgsp_f], cpd_hmt ; p_rng = [1])
    for l = 0 : 2
        ll = Int64(2l)
        cpsp = BuildSCompSpace([sgsp_f, sgsp_b], sec_tot, ll)
        cpop_hmt = BuildSCompOperator(cpsp, cpd_hmt, sgop_hmt)
        enrg, st = GetEigensystem(cpop_hmt, 10 ; issymmetric = true)
        for i in eachindex(enrg)
            push!(result, [enrg[i] / √(ll + 1), l, s])
        end
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2 && st[3] ≈ 0, result)[1][1]
spec = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(stack(spec)))
