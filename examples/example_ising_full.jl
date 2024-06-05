using FuzzifiED
using LinearAlgebra

nm = 8
cfs = GetLzConfs(nm, 2, nm)
bs = GetIsingBasis(cfs ; qn_p = 1, qn_r = 1, qn_z = 1)
tms_hmt = GetIsingIntTerms(nm ; ps_pot = [4.75, 1.]) - 3.16 * GetXPolTerms(nm)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
hmt_mat_full = MatrixFromOpMat(hmt_mat)
enrg, st = eigen(hmt_mat_full)
@show enrg