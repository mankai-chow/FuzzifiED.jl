using FuzzifiED
using WignerSymbols 
using FuzzifiED.Models

nm = 12
# Implement the conserved quantities and generate the confs
@time "Initialise configurations" cfs = GenerateIsingConfs(nm)
@show cfs.ncf
# Implement the discrete symmetries and initialise the basis
@time "Initialise basis" bs = GenerateIsingBasis(cfs ; PH = 1, Ry = 1, Z2 = 1)
@show bs.dim 
# Record the Hamiltonian operator
tms_hmt = GenerateIsingHamiltonianTerms(nm ; ps_pot = Number[4.75, 1.], fld_h = 3.16)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
# Generate the sparse matrix and diagonalise
@time "Initialise the Hamiltonian matrix" hmtmat = OpMat(hmt)
@show hmtmat.nel
@time "Diagonalise Hamiltonian" enrg, st = GetEigensystem(hmtmat, 10)
@show real(enrg)
# Measure the total angular momentum observable
tms_l2 = GenerateL2Terms(nm, 2)
l2 = Operator(bs, bs, tms_l2 ; red_q = 1, sym_q = 1)
@time "Initialise L2" l2_mat = OpMat(l2)
# Calculate the inner product for each eigenstate
@time "Measure L2" l2_val = [ st[:, i]' * l2_mat * st[:, i] for i = 1 : length(enrg)]
@show real(l2_val)

# Verify whether T is an eigenstate of L^2
st_T = st[:, 3]
st_L2T = l2_mat * st[:, 3]
@show abs(st_L2T' * st_T) ^ 2 / ((st_T' * st_T) * (st_L2T' * st_L2T))
# Measure the density operator observable
# Repeat the calculation for the Z_2 odd sector (with subscript 1)
qnz_s1 = ComplexF64[ 1, -1, 1 ] # Change only the discrete quantum numbers and generate the basis
@time "Initialise Basis Z" bs1 = GenerateIsingBasis(cfs ; PH = 1, Ry = 1, Z2 = -1)
@show bs1.dim 
hmt = Operator(bs1, bs1, tms_hmt ; red_q = 1, sym_q = 1) # Generate and diagonalise Hamiltonian in the new basis
@time "Initialise Hamiltonian" hmtmat = OpMat(hmt)
@show hmtmat.nel
@time "Diagonalise Hamiltonian" enrg1, st1 = GetEigensystem(hmtmat, 10)
@show real(enrg1)
# Record the identity, sigma and epsilon states 
st_I = st[:, 1] # ground state
st_e = st[:, 2] # epsilon state
st_s = st1[:, 1]

# Record the density operator n^z
tms_nz = [ Term(isodd(o) ? 1 / nm : -1 / nm, [1, o, 0, o]) for o = 1 : nm * 2]
# The nz operator sends a state in bs (+) to bs1 (-)
nz = Operator(bs, bs1, tms_nz ; red_q = 1)
# Measuring the finite size OPE
@show abs((st_s' * nz * st_e) / (st_s' * nz * st_I))