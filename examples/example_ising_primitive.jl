using FuzzifiED

#========================================================
IMPLEMENT THE CONSERVED QUANTITIES AND GENERATE THE CONFS
========================================================#

# Inputing the basic setups
nf = 2
nm = 8
no = nf * nm
s = .5 * (nm - 1)
ne = div(no, 2)
# Initialise the arrays
qnu_s = Vector{Int64}(undef, 0)
qnu_o = []
# Record the number of electrons
push!(qnu_o, fill(1, no)) # qnu_o[1] = [1,1,...,1]
push!(qnu_s, ne) 
# Record the angular momentum
push!(qnu_o, [ div(o - 1, nf) for o = 1 : no ]) # qnu_o[2] = [0,0,1,1,...,7,7] to qnu_o
push!(qnu_s, ne * s) 
# Generate the configurations and print the number
@time "Initialise configurations" cfs = Confs(no, qnu_s, qnu_o)
@show cfs.ncf


#=========================================================
IMPLEMENT THE DISCRETE SYMMETRIES AND INITIALISE THE BASIS
=========================================================#

cyc = [ 2, 2, 2 ] # Input three Z_2 symmetries 
qnz_s = ComplexF64[ 1, 1, 1 ] # Quantum numbers are all positive 
# Initialise the vectors
perm_o = []
ph_o = []
fac_o = []
# Record the parity
push!(perm_o, [ isodd(o) ? o + 1 : o - 1 for o = 1 : no]) # perm_o[1] = [2,1,4,3,...,16,15]
push!(ph_o, fill(1, no)) # ph_o[1] = [1,1,...,1] meaning PH
push!(fac_o, [ isodd(o) ? -1 : 1 for o = 1 : no]) # fac_o[1] = [1,-1,1,-1,...,1,-1]
# Record the flavour symmetry
push!(perm_o, [ isodd(o) ? o + 1 : o - 1 for o = 1 : no]) # perm_o[3] = [2,1,4,3,...,16,15]
push!(ph_o, fill(0, no)) # ph_o[3] = [0,0,...,0] meaning no PH
push!(fac_o, fill(ComplexF64(1), no)) # fac_o[3] = [1,1,...,1]
# Record the pi-rotation
push!(perm_o, [ isodd(o) ? no - o : no + 2 - o for o = 1 : no]) # perm_o[2] = [15,16,13,14,...,1,2]
push!(ph_o, fill(0, no)) # ph_o[2] = [0,0,...,0] meaning no PH
push!(fac_o, fill(ComplexF64(1), no)) # fac_o[2] = [1,1,...,1]
# Generate the basis and print the dimension
@time "Initialise basis" bs = Basis(cfs, qnz_s ; cyc, perm_o, ph_o, fac_o)
@show bs.dim 


#==============================
RECORD THE HAMILTONIAN OPERATOR
===============================#

using WignerSymbols
# Input the parameters of the Hamiltonian
ps_pot = [ 4.75, 1. ] * 2.
h = 3.16
tms_hmt = Vector{Term}(undef, 0)
# Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m3
for m1 = 0 : nm - 1
    f1 = 0
    o1 = m1 * nf + f1 + 1
    m1r = m1 - s
    for m2 = 0 : nm - 1
        f2 = 1
        o2 = m2 * nf + f2 + 1
        m2r = m2 - s
        for m3 = 0 : nm - 1
            f3 = 1
            o3 = m3 * nf + f3 + 1
            m3r = m3 - s
            m4 = m1 + m2 - m3 
            if (m4 < 0 || m4 >= nm) continue end
            f4 = 0
            o4 = m4 * nf + f4 + 1
            m4r = m4 - s
            # Calculate the matrix element val from pseudopotentials
            val = ComplexF64(0)
            for l in eachindex(ps_pot)
                if (abs(m1r + m2r) > nm - l || abs(m3r + m4r) > nm - l) break end 
                val += ps_pot[l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, m1r, m2r, -m1r - m2r) * wigner3j(s, s, nm - l, m4r, m3r, -m3r - m4r)
            end 
            # Record the interaction term val, "Cdag", o1, "Cdag", o2, "C", o3, "C", o4
            push!(tms_hmt, Term(val, [1, o1, 1, o2, 0, o3, 0, o4]))
        end
    end
    o1x = o1 + 1
    # Record the transverse field term
    push!(tms_hmt, Term(-h, [1, o1, 0, o1x]))
    push!(tms_hmt, Term(-h, [1, o1x, 0, o1]))
end
# Generate the Hamiltonian operator
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)


#=========================================
GENERATE THE SPARSE MATRIX AND DIAGONALISE
=========================================#

@time "Initialise the Hamiltonian matrix" hmt_mat = OpMat(hmt)
@show hmt_mat.nel
@time "Diagonalise Hamiltonian" enrg, st = GetEigensystem(hmt_mat, 10)
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
tms_l2 = tms_lz * tms_lz - tms_lz + tms_lp * tms_lm
# Initialise the L2 operator
l2 = Operator(bs, bs, tms_l2 ; red_q = 1, sym_q = 1)
@time "Initialise L2" l2_mat = OpMat(l2)
# Calculate the inner product for each eigenstate
@time "Measure L2" l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
@show real(l2_val)
# Verify whether T is an eigenstate of L^2
st_T = st[:, 3]
st_L2T = l2_mat * st[:, 3]
@show abs(st_L2T' * st_T) ^ 2 / ((st_T' * st_T) * (st_L2T' * st_L2T))


#======================================
MEASURE THE DENSITY OPERATOR OBSERVABLE
======================================#

# Repeat the calculation for the Z_2 odd sector (with subscript 1)
qnz_s1 = ComplexF64[ 1, -1, 1 ] # Change only the discrete quantum numbers and generate the basis
@time "Initialise Basis Z" bs1 = Basis(cfs, qnz_s1 ; cyc, perm_o, ph_o, fac_o) 
@show bs1.dim 
hmt = Operator(bs1, bs1, tms_hmt ; red_q = 1, sym_q = 1) # Generate and diagonalise Hamiltonian in the new basis
@time "Initialise Hamiltonian" hmt_mat = OpMat(hmt)
@show hmt_mat.nel
@time "Diagonalise Hamiltonian" enrg1, st1 = GetEigensystem(hmt_mat, 10)
@show real(enrg1)
# Record the identity, sigma and epsilon states 
st_I = st[:, 1] # ground state
st_e = st[:, 2] # epsilon state
st_s = st1[:, 1]

# Record the density operator n^z
tms_nz = Vector{Term}(undef, 0)
for m1 = 0 : nm - 1
    o1u = 2 * m1 + 1
    o1d = 2 * m1 + 2
    push!(tms_nz, Term( 1 / nm, [1, o1u, 0, o1u]))
    push!(tms_nz, Term(-1 / nm, [1, o1d, 0, o1d]))
end
# The nz operator sends a state in bs (+) to bs1 (-)
nz = Operator(bs, bs1, tms_nz ; red_q = 1)
# Measuring the finite size OPE
@show abs((st_s' * nz * st_e) / (st_s' * nz * st_I))