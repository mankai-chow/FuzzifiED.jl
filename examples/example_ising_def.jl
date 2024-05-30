using FuzzifiED

nm = 12
cfs_pp = GetIsingDefConfs(nm, nm ; def_conf = [1, 1])
@show digits(cfs_pp.conf[1], base = 2, pad = 2 * nm)
@show cfs_pp.ncf
bs_pp = GetIsingBasis(cfs_pp ; qn_p = 1, qn_r = 1)
tms_hmt = GetIsingIntTerms(nm, [4.75, 1.]) - 3.16 * GetXPolTerms(nm)
hmt = Operator(bs_pp, bs_pp, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
st_pp = st[:, 1]
@show enrg

cfs_pm = GetIsingDefConfs(nm, nm ; def_conf = [1, -1])
@show digits(cfs_pm.conf[1], base = 2, pad = 2 * nm)
@show cfs_pm.ncf
bs_pm = GetIsingBasis(cfs_pm ; qn_p = 1)
tms_hmt = GetIsingIntTerms(nm, [4.75, 1.]) - 3.16 * GetXPolTerms(nm)
hmt = Operator(bs_pm, bs_pm, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
st_pm = st[:, 1]
@show enrg

cfs_p0 = GetIsingDefConfs(nm, nm ; def_conf = [1, 0])
@show cfs_p0.ncf
bs_p0 = GetIsingBasis(cfs_p0 ; qn_p = 1)
tms_hmt = GetIsingIntTerms(nm, [4.75, 1.]) - 3.16 * GetXPolTerms(nm)
hmt = Operator(bs_p0, bs_p0, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
@show enrg

op_conv = Operator(bs_pp, bs_pm, [Term(1.0, [1, 2 * nm, 0, 2 * nm - 1])])
@show st_pm' * op_conv * st_pp