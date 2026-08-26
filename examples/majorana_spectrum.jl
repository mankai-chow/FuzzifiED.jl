# This example calculates the spectrum of a free Majorana fermion.

using FuzzifiED
using FuzzifiED.Fuzzifino
using FuzzifiED.JackToolkit
using SO3lver
FuzzifiED.ElementType = Float64

function GetBPfJack(nof :: Int64, nob :: Int64, sec :: Vector{Int64}, bs :: SBasis, l2c2_mat :: OpMat, nst :: Int64)
    st0 = GetJackStates(bs, nob, sec[3], 2, 2, sec[2]) 
    l2c2_val, st = OrganizeJackStates(st0, l2c2_mat)
    @info "Sector $sec, # of states $(size(st, 2))"
    return l2c2_val, st
end

nmf = 10
nmb = nmf - 1
ne = nmf
s = (nmf - 1) / 2

qnd_f = [ 
    GetNeQNDiag(nmf, 2),
    GetLz2QNDiag(nmf, 1),
    GetNeQNDiag(nmf)]
qnd_b = [ 
    zero(SQNDiag, 0, nmb), 
    GetBosonLz2SQNDiag(0, nmb, 1), 
    GetBosonNeSQNDiag(0, nmb)]
modul = [2, 1, 1]
sec_f = stack([ [ne1 % 2, ((nmf + 1) * ne1) % 2, ne1] for ne1 = 0 : ne])
sec_b = stack([ [0, ((nmb + 1) * ne1) % 2, ne1] for ne1 = 0 : ne])
nebm_b = [ne1 for ne1 = 0 : ne]

f = GetElectronMod(nmf, 1, 1)
b = GetBosonSMod(nmb, 1, 1)

tms_lzlp_f = GetLzLpTerms(nmf, 1)
tms_lzlp_b = GetBosonLzLpSTerms(nmb, 1)
tms_proj = ContractMod(b' * b' * b', b * b * b, 3(s-1/2))

cpd_hop = CoupleDecomps([f' * f', b * b], ConvPsPot(Dict(2s-1 => 1))..., [0 0 ; 0 0 ; 2 -2]) + CoupleDecomps([f * f, b' * b'], ConvPsPot(Dict(2s-1 => 1))..., [0 0 ; 0 0 ; -2 2])
cpd_fb = CoupleDecomps([f' * f, b' * b], ConvPsPot(RecoupleAngMom(s-1/2, s, s, s-1/2, Dict(2s-1/2 => 1)))..., [0 0 ; 0 0 ; 0 0 ])
cpd_bb = SingleSegCouple(2, 2, ContractMod(b' * b', b * b, 2s-1), [0, 0, 0])
cpd_μ = SingleSegCouple(2, 1, GetPolTerms(nmf, 1, [1;;]), [0, 0, 0])
cpd_hmt = 2.0 * cpd_fb + 1.0 * cpd_bb - 0.3 * cpd_hop

sgsp_f = BuildSegSpace(nmf, sec_f, qnd_f, tms_lzlp_f, modul)

nst_max = [ length(GetJackRoots(nmb, sec_b[3, i], 2, 2 ; fermion = false)) for i in axes(sec_b, 2) ]
sgsp_b = BuildSegSpace(0, nmb, nebm_b, sec_b, qnd_b, tms_lzlp_b, tms_proj, [0.0], modul ; l2c2_ratio = 0.1/nmf^2, diag_method = GetBPfJack, nst_max)
sgop_hmt = BuildSegOperators([sgsp_f, sgsp_b], cpd_hmt)
result = []
for l = 0 : 1/2 : 2
    ll = Int64(2l)
    sec_tot = [nmf + ll%2, 0, ne]
    cpsp = BuildCompSpace([sgsp_f, sgsp_b], sec_tot, ll)
    cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
    enrg, st = GetEigensystem(cpop_hmt, 10 ; issymmetric = true)
    for i in eachindex(enrg)
        push!(result, [enrg[i], l])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2, result)[1][1]
spec = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(stack(spec)))
