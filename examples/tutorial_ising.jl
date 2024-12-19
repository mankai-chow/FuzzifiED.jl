# In this tutorial, we use the built-in QNs and operators to calculate the Ising model on fuzzy sphere. 
# We show how to compute the eigenstates and their energies in a sector with given U(1) and Z_2 conserved quantities,
# and how to measure their other quantum numbers like total angular momentum, 
# and verify the eigenstates of Hamiltonian are also their eigenstates. 
# We also show how to construct a density operator and measure the inner products. 

using HDF5
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
    # WARNING : THE -1 FACTOR CAN BE ONLY PUT AT THE SPIN-UP ELECTRONS
    # OTHERWISE AN OVERALL (-1) FACTOR WILL BE PRODUCED AT ODD nm.
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

# Write the sparse matrix into HDF5 file
using HDF5 
f = h5open("data_tmp.h5", "cw")
write(f, "hmt_mat", hmt_mat)
close(f)
f = h5open("data_tmp.h5", "r") 
hmt_mat1 = read(f, "hmt_mat", OpMat{ComplexF64})
close(f)

# Measure the total angular momentum
tms_l2 = GetL2Terms(nm, 2)
l2 = Operator(bs, tms_l2)
l2_mat = OpMat(l2)
l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
@show l2_val

st_T = st[:, 3]
st_L2T = l2_mat * st[:, 3]
@show abs(st_L2T' * st_T) ^ 2 / ((st_T' * st_T) * (st_L2T' * st_L2T))

# Repeat the calculation for the Z_2-odd sector
bs1 = Basis(cfs, [1, -1, 1], qnf)
hmt = Operator(bs1, bs1, tms_hmt ; red_q = 1, sym_q = 1) 
hmt_mat = OpMat(hmt)
enrg1, st1 = GetEigensystem(hmt_mat, 10)
st_I = st[:, 1] 
st_e = st[:, 2] 
st_s = st1[:, 1]

# Measure the density operator
obs_nz = GetDensityObs(nm, 2, σz)
tms_nz = SimplifyTerms(GetComponent(obs_nz, 0.0, 0.0))
nz = Operator(bs, bs1, tms_nz ; red_q = 1) 
@show abs((st_s' * nz * st_e) / (st_s' * nz * st_I))