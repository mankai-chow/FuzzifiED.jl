using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]

#================================================
IMPLEMENT THE DIAGONAL QNS AND GENERATE THE CONFS
================================================#

nm = 12
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) 
]
cfs = Confs(2 * nm, [nm, 0], qnd)


#======================================================
IMPLEMENT THE OFF-DIAGONAL QNS AND INITIALISE THE BASIS
======================================================#

qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [1, -1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) 
]
bs = Basis(cfs, [1, 1, 1], qnf)


#===========================
RECORD THE HAMILTONIAN TERMS
===========================#

tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot = 2 .* [4.75, 1.], mat_a = σ1, mat_b = σ2) - 
    3.16 * GetPolTerms(nm, 2 ; mat = σx) 
)


#=========================================
GENERATE THE SPARSE MATRIX AND DIAGONALISE
=========================================#

hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
@show enrg


#============================================
MEASURE THE TOTAL ANGULAR MOMENTUM OBSERVABLE
============================================#

tms_l2 = GetL2Terms(nm, 2)
l2 = Operator(bs, bs, tms_l2 ; red_q = 1, sym_q = 1)
l2_mat = OpMat(l2 ; type = Float64)
l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
@show l2_val

st_T = st[:, 3]
st_L2T = l2_mat * st[:, 3]
@show abs(st_L2T' * st_T) ^ 2 / ((st_T' * st_T) * (st_L2T' * st_L2T))


#======================================
MEASURE THE DENSITY OPERATOR OBSERVABLE
======================================#

bs1 = Basis(cfs, [1, -1, 1], qnf)
hmt = Operator(bs1, bs1, tms_hmt ; red_q = 1, sym_q = 1) 
hmt_mat = OpMat(hmt ; type = Float64)
enrg1, st1 = GetEigensystem(hmt_mat, 10)
st_I = st[:, 1] 
st_e = st[:, 2] 
st_s = st1[:, 1]

obs_nz = Density(nm, 2 ; mat = [ 1 0 ; 0 -1 ])
tms_nz = SimplifyTerms(GetComponent(obs_nz, 0.0, 0.0))
nz = Operator(bs, bs1, tms_nz ; red_q = 1) 
@show abs((st_s' * nz * st_e) / (st_s' * nz * st_I))