# This example calculates the g-function of magnetic line defect in 3d Ising model
# using the ovelaps between the bulk, defect ground state and the lowest defect-creation state.
# This example reproduces Figure 6 in arXiv : 2401.00039
# On my table computer, this calculation takes 3.011 s

using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 12
no = nm * 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, 2 .* [4.75, 1.], σ1, σ2) - 
    3.16 * GetPolTerms(nm, 2, σx) )
    
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) 
]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) 
]
cfs_00 = Confs(2 * nm, [nm, 0], qnd)
bs_00 = Basis(cfs_00, [1, 1, 1], qnf)
hmt = Operator(bs_00, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 3)
st_00 = st[:, 1]

qnd_pp = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, 2), 
    GetPinOrbQNDiag(no, [1, no - 1]), 
    GetPinOrbQNDiag(no, [2, no])
]
qnf_pp = [
    GetRotyQNOffd(nm, 2), 
    GetParityQNOffd(nm, 2, [2,1], [-1,1])
]
cfs_pp = Confs(no, [nm, 0, 2, 0], qnd_pp)
bs_pp = Basis(cfs_pp, [1,1], qnf_pp)

hmt = Operator(bs_pp, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 3)
st_pp = st[:, 1]

qnd_p0 = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, 2), 
    GetPinOrbQNDiag(no, [1]), 
    GetPinOrbQNDiag(no, [2])
]
qnf_p0 = [
    GetParityQNOffd(nm, 2, [2,1], [-1,1])
]
cfs_p0 = Confs(no, [nm, 0, 1, 0], qnd_p0)
bs_p0 = Basis(cfs_p0, [1], qnf_p0)

hmt = Operator(bs_p0, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 3)
st_p0 = st[:, 1]

tms_I = [Term(1, [-1, -1])]
@show ovl_p000 = st_p0' * Operator(bs_00, bs_p0, tms_I) * st_00
@show ovl_p0pp = st_p0' * Operator(bs_pp, bs_p0, tms_I) * st_pp 
@show g_fn = (ovl_p000 / ovl_p0pp) ^ 2