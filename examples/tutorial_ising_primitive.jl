# This tutorial contains the ED code that uses only the core functions. It
# 1. calculates the lowest eigenstates in the symmetry sector L^z=0 and (𝒫,𝒵,ℛ)=(+,+,+),
# 2. measures their total angular momenta, and 
# 3. calcultes the OPE coefficient f_{σσϵ}=⟨σ|n^z_{00}|ϵ⟩/⟨σ|n^z_{00}|0⟩

using FuzzifiED

let

#================================================
IMPLEMENT THE DIAGONAL QNS AND GENERATE THE CONFS
================================================#

# Inputing the basic setups
nf = 2
nm = 12
no = nf * nm
s = .5 * (nm - 1)
# Record the QNDiag
qnd = [ 
    # Number of electrons Ne
    QNDiag(fill(1, no)), 
    # Twice angular momentum 2Lz
    QNDiag([ (o - 1) ÷ nf * 2 - (nm - 1) for o = 1 : no ])
]
cfs = Confs(no, [nm, 0], qnd)
@show cfs.ncf


#======================================================
IMPLEMENT THE OFF-DIAGONAL QNS AND INITIALISE THE BASIS
======================================================#

# Record a QNOffd by site permutation (and facultatives particle-hole, factor, cycle)
qnf = [ 
    # Parity (Particle-hole)
    QNOffd([ isodd(o) ? o + 1 : o - 1 for o = 1 : no], true, ComplexF64[ isodd(o) ? -1 : 1 for o = 1 : no]),
    # Flavour symmetry
    QNOffd([ isodd(o) ? o + 1 : o - 1 for o = 1 : no]),
    # Y-axis pi-rotation
    QNOffd([ isodd(o) ? no - o : no + 2 - o for o = 1 : no], ComplexF64(-1) .^ (collect(0 : nm * nf - 1) .÷ nf))
]
bs = Basis(cfs, [1, 1, 1], qnf) 
@show bs.dim 


#===========================
RECORD THE HAMILTONIAN TERMS
===========================#

using WignerSymbols
# Input the parameters of the Hamiltonian
ps_pot = [ 4.75, 1. ] * 2.
h = 3.16
tms_hmt = Term[]
# Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m
m = zeros(Int64, 4)
for m[1] = 0 : nm - 1, m[2] = 0 : nm - 1, m[3] = 0 : nm - 1
    m[4] = m[1] + m[2] - m[3]
    (m[4] < 0 || m[4] >= nm) && continue
    f = [0, 1, 1, 0]
    o = m .* nf .+ f .+ 1
    mr = m .- s
    
    # Calculate the matrix element val from pseudopotentials
    val = ComplexF64(0)
    for l in eachindex(ps_pot)
        (abs(mr[1] + mr[2]) > nm - l || abs(mr[3] + mr[4]) > nm - l) && break 
        val += ps_pot[l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, mr[1], mr[2], -mr[1] - mr[2]) * wigner3j(s, s, nm - l, mr[4], mr[3], -mr[3] - mr[4])
    end 
    # Record the interaction term val, "Cdag", o1, "Cdag", o2, "C", o3, "C", o4
    tms_hmt += Terms(val, [1, o[1], 1, o[2], 0, o[3], 0, o[4]])
end 
for m = 0 : nm - 1
    o = m * nf .+ [1, 2]
    # Record the transverse field term
    tms_hmt += Terms(-h, [1, o[1], 0, o[2]])
    tms_hmt += Terms(-h, [1, o[2], 0, o[1]])
end


#=========================================
GENERATE THE SPARSE MATRIX AND DIAGONALISE
=========================================#

hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 10)
@show real(enrg)


#============================================
MEASURE THE TOTAL ANGULAR MOMENTUM OBSERVABLE
============================================#

tms_lz = 
    [ begin m = div(o - 1, nf)
        Term(m - s, [1, o, 0, o])
    end for o = 1 : no ]
tms_lp = 
    [ begin m = div(o - 1, nf)
        Term(sqrt(m * (nm - m)), [1, o, 0, o - nf])
    end for o = nf + 1 : no ]
tms_lm = tms_lp' 
tms_l2 = SimplifyTerms(tms_lz * tms_lz - tms_lz + tms_lp * tms_lm)
# Initialise the L2 operator
l2 = Operator(bs, tms_l2)
@time "Initialise L2" l2_mat = OpMat(l2)
# Calculate the inner product for each eigenstate
@time "Measure L2" l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
@show real(l2_val)


#======================================
MEASURE THE DENSITY OPERATOR OBSERVABLE
======================================#

# Repeat the calculation for the Z_2 odd sector (with subscript 1)
bs_m = Basis(cfs, [ 1, -1, 1 ], qnf) 
hmt_m = Operator(bs_m, bs_m, tms_hmt ; red_q = 1, sym_q = 1) # Generate and diagonalise Hamiltonian in the new basis
hmt_mat_m = OpMat(hmt_m)
enrg_m, st_m = GetEigensystem(hmt_mat_m, 10)
@show real(enrg_m)
# Record the identity, sigma and epsilon states 
st0 = st[:, 1] # ground state
ste = st[:, 2] # epsilon state
sts = st_m[:, 1]

# Record the density operator n^z
tms_nz00 = Term[]
for m = 0 : nm - 1
    o = m * nf .+ [1, 2]
    # Record the transverse field term
    tms_nz00 += Terms( 1 / nm, [1, o[1], 0, o[1]])
    tms_nz00 += Terms(-1 / nm, [1, o[2], 0, o[2]])
end
# The nz operator sends a state in bs (+) to bs_m (-)
nz00 = Operator(bs, bs_m, tms_nz00 ; red_q = 1)
# Measuring the finite size OPE
f_sse = abs((sts' * nz00 * ste) / (sts' * nz00 * st0))
@show f_sse

end