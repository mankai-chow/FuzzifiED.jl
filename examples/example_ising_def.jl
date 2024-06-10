using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
≊(x, y) = abs(x - y) < eps(Float32)

nm = 12
no = nm * 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot = 2 .* [4.75, 1.], mat_a = σ1, mat_b = σ2) - 
    3.16 * GetPolTerms(nm, 2 ; mat = σx) )
qnd_pp = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, 2), 
    GetPinOrbQNDiag(no, [1, no - 1]), 
    GetPinOrbQNDiag(no, [2, no])
]
qnf_pp = [
    GetRotyQNOffd(nm, 2), 
    GetParityQNOffd(nm, 2, [2,1],[1,-1])
]
cfs_pp = Confs(no, [nm, 0, 2, 0], qnd_pp)
bs_pp = Basis(cfs_pp, [1,1], qnf_pp)

hmt = Operator(bs_pp, bs_pp, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
st_pp = st[:, 1]
@show enrg

qnd_p0 = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, 2), 
    GetPinOrbQNDiag(no, [1]), 
    GetPinOrbQNDiag(no, [2])
]
qnf_p0 = [
    GetParityQNOffd(nm, 2, [2,1],[1,-1])
]
cfs_p0 = Confs(no, [nm, 0, 1, 0], qnd_p0)
bs_p0 = Basis(cfs_p0, [1], qnf_p0)

hmt = Operator(bs_p0, bs_p0, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
st_p0 = st[:, 1]
@show enrg

qnd_pm = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, 2), 
    GetPinOrbQNDiag(no, [1, no]), 
    GetPinOrbQNDiag(no, [2, no - 1])
]
qnf_pm = [
    GetParityQNOffd(nm, 2, [2,1],[1,-1]),    
    GetRotyQNOffd(nm, 2) * GetFlavPermQNOffd(nm, 2, [2,1])
]
cfs_pm = Confs(no, [nm, 0, 2, 0], qnd_pm)
bs_pm = Basis(cfs_pm, [1, 1], qnf_pm)

hmt = Operator(bs_pm, bs_pm, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
st_pm = st[:, 1]
@show enrg