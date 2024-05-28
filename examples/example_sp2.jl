using FuzzifiED

nm = 7
nf = 4
cfs = GetSpnConfs(nm, nf, 2 * nm)
bs = Basis(cfs)
bs = GetSpnBasis(cfs, nf ; qn_p = 1, qn_r = 1, qn_z = [1, 1], qn_x = [1] )

tms_hmt = GetIdDenIntTerms(nm, nf, [1.]) - GetSpnPairIntTerms(nm, nf, [.9])
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)

tms_l2 = GetL2Terms(nm, nf)
l2 = Operator(bs, bs, tms_l2 ; red_q = 1, sym_q = 1)
l2_mat = OpMat(l2 ; type = Float64)
l2_val = [ st[:, i]' * l2_mat * st[:, i] for i = 1 : length(enrg)]

tms_c2 = GetSpnC2Terms(nm, nf)
c2 = Operator(bs, bs, tms_c2 ; red_q = 1, sym_q = 1)
c2_mat = OpMat(c2 ; type = Float64)
c2_val = [ st[:, i]' * c2_mat * st[:, i] for i = 1 : length(enrg)]

display(reduce(hcat, [["Energy";enrg], ["L^2";l2_val], ["C_2";c2_val]]))