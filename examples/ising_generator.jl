using FuzzifiED
using FuzzifiEDFullRotation
using LinearAlgebra
FuzzifiED.ElementType = Float64

nm = 12
ne = nm
s = (nm - 1) / 2

qnd_pt = [
    GetNeQNDiag(nm),
    GetLz2QNDiag(nm, 1)
]
sec_pt = stack([ [ne1, ((nm + 1) * ne1) % 2] for ne1 = 0 : ne])
sec_tot = [ne, 0]
tms_lzlp_pt = GetLpLzTerms(nm, 1)


n_mod = GetDensityMod(nm, 1, [1;;]) 
cpd_int = 2 * CoupleDecomps([ n_mod, n_mod ], ConvPsPot(RecouplePsPot((nm-1)/2, [4.75, 1.0]))..., [ 0 0 ; 0 0 ])
c_obs = GetElectronObs(nm, 1, 1) 
cpd_nx = ContactCouple([c_obs', c_obs], [ 1 -1 ; 0 0 ]) - ContactCouple([c_obs, c_obs'], [ -1 1 ; 0 0 ])
cpd_hmt = cpd_int - 3.16 * cpd_nx

sgsp = BuildSegSpace(nm, sec_pt, qnd_pt, tms_lzlp_pt)
sgop_hmt = BuildSegOperators([sgsp, sgsp], cpd_hmt)

cpsp0 = BuildCompSpace([sgsp, sgsp], sec_tot, 0)
cpop_hmt0 = BuildCompOperator(cpsp0, cpd_hmt, sgop_hmt)
enrg0, st0 = GetEigensystem(cpop_hmt0, 10 ; issymmetric = true)

cpsp1 = BuildCompSpace([sgsp, sgsp], sec_tot, 2)
cpop_hmt1 = BuildCompOperator(cpsp1, cpd_hmt, sgop_hmt)
enrg1, st1 = GetEigensystem(cpop_hmt1, 10 ; issymmetric = true)
enrg1 ./= √3

display(enrg0)
display(enrg1)

stI = st0[:, 1]
stσ = st0[:, 2]
stϵ = st0[:, 3]
st∂σ = st1[:, 1]
st∂ϵ = st1[:, 2]

cpd_nx_l1 = ContactCouple([c_obs', c_obs], [ 1 -1 ; 0 0 ], 2) - ContactCouple([c_obs, c_obs'], [ -1 1 ; 0 0 ], 2)
sgop_nx_l1 = BuildSegOperators([sgsp, sgsp], cpd_nx_l1)
cpop_nx_l1 = BuildCompOperator(cpsp0, cpsp1, cpd_nx_l1, sgop_nx_l1)

cpd_int_l1_0 = CoupleDecomps([ n_mod, n_mod ], RecouplePsPot(s, s, s, s, [Int64[4s 4s;4s 2]], [1.0])..., [ 0 0 ; 0 0 ])
sgop_int_l1_0 = BuildSegOperators([sgsp, sgsp], cpd_int_l1_0)
cpop_int_l1_0 = BuildCompOperator(cpsp0, cpsp1, cpd_int_l1_0, sgop_int_l1_0)

cpd_int_l1_1 = CoupleDecomps([ n_mod, n_mod ], RecouplePsPot(s, s, s, s, [Int64[4s-2 4s-2;4s-2 2]], [1.0])..., [ 0 0 ; 0 0 ])
sgop_int_l1_1 = BuildSegOperators([sgsp, sgsp], cpd_int_l1_1)
cpop_int_l1_1 = BuildCompOperator(cpsp0, cpsp1, cpd_int_l1_1, sgop_int_l1_1)

opst = [cpop_nx_l1 * stI, cpop_int_l1_1 * stI, cpop_int_l1_0 * stI]
mat = [sti' * stj for sti in opst, stj in opst]
eigval, eigvec = eigen(mat);

cpd_pk = eigvec[:, 1]' * [cpd_nx_l1, cpd_int_l1_1, cpd_int_l1_0]
sgop_pk = BuildSegOperators([sgsp, sgsp], cpd_pk)
cpop_pk = BuildCompOperator(cpsp0, cpsp1, cpd_pk, sgop_pk)

compare_st(st0, st1) = abs(st0' * st1) ^ 2 / ((st0' * st0) * (st1' * st1)) ;

@show compare_st(cpop_pk * stσ, st∂σ)
@show compare_st(cpop_pk * stϵ, st∂ϵ)