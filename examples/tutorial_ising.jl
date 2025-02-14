# This tutorial contains the ED code that uses the built-in models. It
# 1. calculates the lowest eigenstates in the symmetry sector L^z=0 and (𝒫,𝒵,ℛ)=(+,+,+),
# 2. measures their total angular momenta, and 
# 3. calcultes the OPE coefficient f_{σσϵ}=⟨σ|n^z_{00}|ϵ⟩/⟨σ|n^z_{00}|0⟩

using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]

# Implement the diagonal QNs and generate the confs
nm = 12
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) 
]
cfs = Confs(2 * nm, [nm, 0], qnd)

# Implement the off-diagonal QNs and initialise the basis 
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) 
]
bs = Basis(cfs, [1, 1, 1], qnf)

# Record the Hamiltonian terms
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, 2 .* [4.75, 1.], σ1, σ2) - 
    3.16 * GetPolTerms(nm, 2, σx) 
)

# Generate the sparse matrix and diagonalise 
hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 10)
@show enrg

# Measure the total angular momentum
tms_l2 = GetL2Terms(nm, 2)
l2 = Operator(bs, tms_l2)
l2_mat = OpMat(l2)
l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
@show l2_val

# Repeat the calculation for the Z_2-odd sector
bs_m = Basis(cfs, [1, -1, 1], qnf)
hmt_m = Operator(bs_m, bs_m, tms_hmt ; red_q = 1, sym_q = 1) 
hmt_mat_m = OpMat(hmt_m)
enrg_m, st_m = GetEigensystem(hmt_mat_m, 10)
st0 = st[:, 1] 
ste = st[:, 2] 
sts = st_m[:, 1]

# Measure the density operator
obs_nz = GetDensityObs(nm, 2, σz)
tms_nz = SimplifyTerms(GetComponent(obs_nz, 0.0, 0.0))
nz = Operator(bs, bs_m, tms_nz ; red_q = 1) 
f_sse = abs((sts' * nz * ste) / (sts' * nz * st0))
@show f_sse