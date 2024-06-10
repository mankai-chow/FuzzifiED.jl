using FuzzifiED
using SpecialFunctions
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
≊(x, y) = abs(x - y) < eps(Float32)

nm = 8
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) ]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [1, -1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) ]

tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot = 2 .* [4.75, 1.], mat_a = σ1, mat_b = σ2) - 
    3.16 * GetPolTerms(nm, 2 ; mat = σx) )

cfs = Confs(2 * nm, [nm, 0], qnd)
bs = Basis(cfs, [1, 1, 1], qnf)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 5)
@show enrg
st_g = st[:, 1]

amp_oa = vcat([ sqrt(beta_inc(1 + m, nm - m, 0.5)[1]) for f = 1 : 2, m = 0 : nm - 1]...) ;
amp_ob = vcat([ sqrt(beta_inc(1 + m, nm - m, 0.5)[2]) for f = 1 : 2, m = 0 : nm - 1]...) ;
secf_lst = [ [[1], [1]], [[-1], [-1]] ]
secd_lst = Vector{Vector{Int64}}[]
for nea = 0 : nm 
    neb = nm - nea 
    for lza = -min(nea, neb) * (nm - 1) : 2 : min(nea, neb) * (nm - 1)
        lzb = -lza 
        push!(secd_lst, [[nea, lza], [neb, lzb]])
    end
end
@time "Calculate entanglement spectrum" ent_spec = GetEntSpec(st_g, bs, secd_lst, secf_lst ; qnd_a = qnd, qnf_a = [GetFlavPermQNOffd(nm, 2, [2, 1])], amp_oa, amp_ob)

eig_rho = vcat(values(ent_spec)...)
tr_rho = sum(eig_rho)
@show tr_rho
ent_entropy = -sum(eig_rho .* log.(eig_rho))
@show ent_entropy