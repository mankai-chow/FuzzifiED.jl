using FuzzifiED
using FuzzifiED.Fuzzifino
using FuzzifiEDFullRotation
FuzzifiED.ElementType = Float64
Base.zero(SQNDiag, nof, nob) = SQNDiag(zeros(Int64, nof), zeros(Int64, nob))

nmf = 10
nmb = nmf - 1
ne = nmf
sf = (nmf - 1) / 2
sb = (nmb - 1) / 2

qnd_f = [ 
    SQNDiag(GetFlavQNDiag(nmf, 1, [1], 1, 2), 1), 
    SQNDiag(GetLz2QNDiag(nmf, 1), 1),
    SQNDiag(GetNeQNDiag(nmf), 1), 
    GetBosonNeSQNDiag(nmf, 1)]
qnd_b = [ 
    zero(SQNDiag, 1, nmb), 
    GetBosonLz2SQNDiag(1, nmb, 1), 
    GetBosonNeSQNDiag(1, nmb),
    SQNDiag(GetNeQNDiag(1), nmb)]
modul = [2, 1, 1, 1]
sec_f = [ [ne1 % 2, ((nmf + 1) * ne1) % 2, ne1, 0] for ne1 = 0 : ne]
sec_b = [ [0, ((nmb + 1) * ne1) % 2, ne1, 0] for ne1 = 0 : ne]
nebm_f = [0 for ne1 = 0 : ne]
nebm_b = [ne1 for ne1 = 0 : ne]
tms_lzlp_f = STerms.(GetLpLzTerms(nmf, 1))
tms_lzlp_b = GetBosonLpLzSTerms(nmb, 1)

amd_ff = GetFermionSMod(nmf, 1, 1) * GetFermionSMod(nmf, 1, 1)
amd_fb = GetFermionSMod(nmf, 1, 1) * GetBosonSMod(nmb, 1, 1)
amd_bb = GetBosonSMod(nmb, 1, 1) * GetBosonSMod(nmb, 1, 1)
amd_nf = GetFerDensitySMod(nmf, 1, [1;;]) 
amd_nb = GetBosDensitySMod(nmb, 1, [1;;]) 
tms_bb  = ContractMod(amd_bb', amd_bb, nmf - 2)

cpd_hop = SCoupleDecomp([amd_ff', amd_bb], ConvPsPot(Dict(nmf - 2 => 1))..., [0 0 ; 0 0 ; 2 -2 ; 0 0]) + SCoupleDecomp([amd_ff, amd_bb'], ConvPsPot(Dict(nmf - 2 => 1))..., [0 0 ; 0 0 ; -2 2 ; 0 0])
cpd_fb = SCoupleDecomp([amd_nf, amd_nb], ConvPsPot(RecouplePsPot(sb, sf, sf, sb, Dict(nmf - 3/2 => 1)))..., [0 0 ; 0 0 ; 0 0 ; 0 0 ])
cpd_bb = SingleSegSCouple(2, 2, tms_bb, [0, 0, 0, 0])
cpd_μ = SingleSegSCouple(2, 1, STerms(GetPolTerms(nmf, 1, [1;;])), [0, 0, 0, 0])
cpd_hmt = PrepareCouple(2.0 * cpd_fb + 1.0 * cpd_bb - 0.3 * cpd_hop)

sgsp_f = BuildSSegSpace(nmf, 1, nebm_f, sec_f, qnd_f, tms_lzlp_f)
sgsp_b = BuildSSegSpace(1, nmb, nebm_b, sec_b, qnd_b, tms_lzlp_b)
sgop_hmt = BuildSSegOperators([sgsp_f, sgsp_b], cpd_hmt, [2, 1, 1, 1])
result = []
for ltot = 0 : 1/2 : 2
    ll = Int64(2ltot)
    sec_tot = [nmf + ll%2, 0, ne, 0]
    cpsp = BuildSCompSpace([sgsp_f, sgsp_b], sec_tot, ll, modul)
    @show ltot, cpsp.dim
    cpop_hmt = BuildSCompOperator(cpsp, cpd_hmt, sgop_hmt)
    enrg, st = GetEigensystem(cpop_hmt, 10)
    for i in eachindex(enrg)
        push!(result, [enrg[i] / √(ll + 1), ltot])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2, result)[1][1]
spec = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(hcat(spec...)))