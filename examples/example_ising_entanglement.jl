using FuzzifiED
using SpecialFunctions
FuzzifiED.SilentStd = true

nm = 8
cfs = GetLzConfs(nm, 2, nm)
bs = GetIsingBasis(cfs ; qn_p = 1, qn_r = 1, qn_z = 1)
tms_hmt = GetIsingIntTerms(nm ; ps_pot = [4.75, 1.]) - 3.16 * GetXPolTerms(nm)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 5)
@show enrg
st_g = st[:, 1]

amp_oa = vcat([ sqrt(beta_inc(1 + m, nm - m, 0.5)[1]) for f = 1 : 2, m = 0 : nm - 1]...) ;
amp_ob = vcat([ sqrt(beta_inc(1 + m, nm - m, 0.5)[2]) for f = 1 : 2, m = 0 : nm - 1]...) ;
qnz_s_lst = Any[ [[1], [1]], [[-1], [-1]] ]
qnu_s_lst = []
for nea = 0 : nm 
    neb = nm - nea 
    for lza = (nea - min(nea, neb)) * (nm - 1) ÷ 2 : (nea + min(nea, neb)) * (nm - 1) ÷ 2
        lzb = nm * (nm - 1) ÷ 2 - lza 
        push!(qnu_s_lst, [[nea, lza], [neb, lzb]])
    end
end
@time "Calculate entanglement spectrum" ent_spec = GetEntSpec(st_g, bs, qnu_s_lst, qnz_s_lst ; GetLzQnu(nm, 2)..., GetIsingQnz(nm ; qn_z = 1)..., amp_oa, amp_ob)

eig_rho = vcat(values(ent_spec)...)
tr_rho = sum(eig_rho)
@show tr_rho
ent_entropy = -sum(eig_rho .* log.(eig_rho))
@show ent_entropy