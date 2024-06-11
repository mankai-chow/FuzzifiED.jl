# This example calculates various OPE coefficients at nm = 12
# by taking overlaps between CFT states and density operators and composite.
# On my table computer, this calculation takes 2.434 s

using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]
≊(x, y) = abs(x - y) < eps(Float32)

nm = 12
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
tms_l2 = GetL2Terms(nm, 2)

cfs = Confs(2 * nm, [nm, 0], qnd)
bsp = Basis(cfs, [1, 1, 1], qnf)
bsm = Basis(cfs, [1,-1, 1], qnf)

hmt = Operator(bsp, bsp, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 20)
l2 = Operator(bsp, bsp, tms_l2 ; red_q = 1, sym_q = 1)
l2_mat = OpMat(l2 ; type = Float64)
l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
idl0 = [ i for i in eachindex(enrg) if l2_val[i] ≊ 0 ]
idl2 = [ i for i in eachindex(enrg) if l2_val[i] ≊ 6 ]
stg = st[:, idl0[1]]
stϵ = st[:, idl0[2]]
stT = st[:, idl2[1]]
stϵ1 = st[:, idl0[4]]

hmt = Operator(bsm, bsm, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 20)
l2 = Operator(bsm, bsm, tms_l2 ; red_q = 1, sym_q = 1)
l2_mat = OpMat(l2 ; type = Float64)
l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
idl0 = [ i for i in eachindex(enrg) if l2_val[i] ≊ 0 ]
idl2 = [ i for i in eachindex(enrg) if l2_val[i] ≊ 6 ]
stσ = st[:, idl0[1]]
stσ2 = st[:, idl2[2]]
stσ1 = st[:, idl0[4]]

obs_nz = Density(nm, 2 ; mat = σz)
obs_nx = Density(nm, 2 ; mat = σx)
tms_nz00 = SimplifyTerms(GetComponent(obs_nz, 0.0, 0.0))
tms_nz20 = SimplifyTerms(GetComponent(obs_nz, 2.0, 0.0))
tms_nx00 = SimplifyTerms(GetComponent(obs_nx, 0.0, 0.0))
tms_nx20 = SimplifyTerms(GetComponent(obs_nx, 2.0, 0.0))
tms_Oϵ = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot = 2 .* [4.75, 1.], mat_a = σ1, mat_b = σ2) +
    3.16 * GetPolTerms(nm, 2 ; mat = σx) )

nz00 = Operator(bsp, bsm, tms_nz00 ; red_q = 1) 
nz20 = Operator(bsp, bsm, tms_nz20 ; red_q = 1) 
nx00p = Operator(bsp, bsp, tms_nx00 ; red_q = 1) 
nx00m = Operator(bsm, bsm, tms_nx00 ; red_q = 1)
nx20m = Operator(bsm, bsm, tms_nx20 ; red_q = 1)  
Oϵ00p = Operator(bsp, bsp, tms_Oϵ ; red_q = 1) 
Oϵ00m = Operator(bsm, bsm, tms_Oϵ ; red_q = 1) 

ovl_σz0 = stσ' * nz00 * stg
ovl_ϵx0 = stϵ' * nx00p * stg
ovl_0x0 = stg' * nx00p * stg
ovl_ϵO0 = stϵ' * Oϵ00p * stg
ovl_0O0 = stg' * Oϵ00p * stg
f_σσϵ_1 = abs(stσ' * nz00 * stϵ / ovl_σz0)
f_σσϵ_2 = abs((stσ' * nx00m * stσ - ovl_0x0) / ovl_ϵx0)
f_σσϵ_3 = abs((stσ' * Oϵ00m * stσ - ovl_0O0) / ovl_ϵO0)
f_ϵϵϵ_1 = abs((stϵ' * nx00p * stϵ - ovl_0x0) / ovl_ϵx0)
f_ϵϵϵ_2 = abs((stϵ' * nx00p * stϵ - ovl_0O0) / ovl_ϵO0)
f_σσT = abs(stσ' * nz20 * stT / ovl_σz0) * √(15/8)
f_σσ1T = abs(stσ1' * nz20 * stT / ovl_σz0) * √(15/8)
f_σσ2ϵ1 = abs(stσ2' * nz20 * stϵ1 / ovl_σz0) * √(15/8)
f_TTϵ_1 = abs((stT' * nx00p * stT - ovl_0x0) / ovl_ϵx0)
f_TTϵ_2 = abs((stT' * Oϵ00p * stT - ovl_0O0) / ovl_ϵO0)

@show ovl_σz0
@show ovl_ϵx0
@show ovl_0x0
@show ovl_ϵO0
@show ovl_0O0
@show f_σσϵ_1
@show f_σσϵ_2
@show f_σσϵ_3
@show f_ϵϵϵ_1
@show f_ϵϵϵ_2
@show f_σσT
@show f_σσ1T 
@show f_σσ2ϵ1
@show f_TTϵ_1
@show f_TTϵ_2