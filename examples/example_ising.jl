using FuzzifiED

nm = 8
cfs = GetLzConfs(nm, 2, nm)
bs = GetIsingBasis(cfs ; qn_p = 1, qn_r = 1, qn_z = 1)
tms_hmt = GetIsingIntTerms(nm ; ps_pot = [4.75, 1.]) - 3.16 * GetXPolTerms(nm)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
@show enrg

tms_l2 = GetL2Terms(nm, 2)
l2 = Operator(bs, bs, tms_l2 ; red_q = 1, sym_q = 1)
l2_mat = OpMat(l2 ; type = Float64)
l2_val = [ st[:, i]' * l2_mat * st[:, i] for i = 1 : length(enrg)]
@show l2_val

st_T = st[:, 3]
st_L2T = l2_mat * st[:, 3]
@show abs(st_L2T' * st_T) ^ 2 / ((st_T' * st_T) * (st_L2T' * st_L2T))

qnz_s1 = ComplexF64[ 1, -1, 1 ] 
bs1 = GetIsingBasis(cfs ; qn_p = 1, qn_r = 1, qn_z = -1)
hmt = Operator(bs1, bs1, tms_hmt ; red_q = 1, sym_q = 1) 
hmt_mat = OpMat(hmt ; type = Float64)
enrg1, st1 = GetEigensystem(hmt_mat, 10)

st_I = st[:, 1] 
st_e = st[:, 2] 
st_s = st1[:, 1]
tms_nz = SimplifyTerms(GetComponent(Density(nm, 2 ; mat = [ 1 0 ; 0 -1 ]), 0.0, 0.0))
nz = Operator(bs, bs1, tms_nz ; red_q = 1) 
@show abs((st_s' * nz * st_e) / (st_s' * nz * st_I))

tms_ndn = Density(nm, 2 ; mat = [ 1 0 ; 0 -1 ]) * Density(nm, 2 ; mat = [ 1 0 ; 0 -1 ])
@show SimplifyTerms(GetPointValue(tms_ndn, 0.0, 0.0))