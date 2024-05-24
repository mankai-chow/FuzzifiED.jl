LibpathFuzzifiED = "./lib_fuzzifi_ed.so"
include("./fuzzifi_ed.jl")


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
push!(qnu_o, vcat(fill(collect(0 : nm - 1), nf)...)) # qnu_o[2] = [0,1,...,7,1,2,...,7] to qnu_o
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
push!(perm_o, [ collect(nm + 1 : 2 * nm) ; collect(1 : nm) ]) # perm_o[1] = [9,10,...,16,1,2,...,8]
push!(ph_o, fill(1, no)) # ph_o[1] = [1,1,...,1] meaning PH
push!(fac_o, [ fill(ComplexF64(1), nm) ; fill(ComplexF64(-1), nm) ]) # fac_o[1] = [1,1,...,1,-1,-1,...,-1]
# Record the flavour symmetry
push!(perm_o, [ collect(nm + 1 : 2 * nm) ; collect(1 : nm) ]) # perm_o[2] = [9,10,...,16,1,2,...,8]
push!(ph_o, fill(0, no)) # ph_o[2] = [0,0,...,0] meaning no PH
push!(fac_o, fill(ComplexF64(1), no)) # fac_o[2] = [1,1,...,1]
# Record the pi-rotation
push!(perm_o, [ collect(nm : -1 : 1) ; collect(2 * nm : -1 : nm + 1) ]) # perm_o[3] = [8,7,...,1,16,15,...,9]
push!(ph_o, fill(0, no)) # ph_o[3] = [0,0,...,0] meaning no PH
push!(fac_o, fill(ComplexF64(1), no)) # fac_o[3] = [1,1,...,1]
# Generate the basis and print the dimension
@time "Initialise basis" bs = Basis(cfs, qnz_s, cyc, perm_o, ph_o, fac_o)
@show bs.dim 


#==============================
RECORD THE HAMILTONIAN OPERATOR
===============================#

using WignerSymbols
# Input the parameters of the Hamiltonian
ps_pot = [ 4.75, 1. ] * 2.
h = 3.16
cstr_hmt = []
fac_hmt = Array{ComplexF64, 1}(undef, 0)
# Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m3
for m1 = 0 : nm - 1
    f1 = 0
    o1 = m1 + f1 * nm + 1
    m1r = m1 - s
    for m2 = 0 : nm - 1
        f2 = 1
        o2 = m2 + f2 * nm + 1
        m2r = m2 - s
        for m3 = 0 : nm - 1
            f3 = 1
            o3 = m3 + f3 * nm + 1
            m3r = m3 - s
            m4 = m1 + m2 - m3 
            if (m4 < 0 || m4 >= nm) continue end
            f4 = 0
            o4 = m4 + f4 * nm + 1
            m4r = m4 - s
            # Calculate the matrix element val from pseudopotentials
            val = ComplexF64(0)
            for l in 1 : length(ps_pot)
                if (abs(m1r + m2r) > nm - l || abs(m3r + m4r) > nm - l) break end 
                val += ps_pot[l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, m1r, m2r, -m1r - m2r) * wigner3j(s, s, nm - l, m4r, m3r, -m3r - m4r)
            end 
            # Record the interaction term val, "Cdag", o1, "Cdag", o2, "C", o3, "C", o4
            push!(cstr_hmt, [1, o1, 1, o2, 0, o3, 0, o4])
            push!(fac_hmt, val)
        end
    end
    o1x = o1 + nm
    # Record the transverse field term
    push!(cstr_hmt, [1, o1, 0, o1x])
    push!(fac_hmt, -h)
    push!(cstr_hmt, [1, o1x, 0, o1])
    push!(fac_hmt, -h)
end
# Generate the Hamiltonian operator
hmt = Operator(bs, bs, 1, 1, cstr_hmt, fac_hmt)


#=========================================
GENERATE THE SPARSE MATRIX AND DIAGONALISE
=========================================#

@time "Initialise the Hamiltonian matrix" hmtmat = OpMat(hmt)
@show hmtmat.nel
@time "Diagonalise Hamiltonian" enrg, st = GetEigensystem(hmtmat, 10)
@show real(enrg)


#============================================
MEASURE THE TOTAL ANGULAR MOMENTUM OBSERVABLE
============================================#

cstr_l2 = []
fac_l2 = Array{ComplexF64, 1}(undef, 0)
for o1 = 1 : no 
    m1 = mod(o1 - 1, nm) 
    # record the -Lz term
    push!(cstr_l2, [1, o1, 0, o1])
    push!(fac_l2, -(m1 - s))
    for o2 = 1 : no 
        m2 = mod(o2 - 1, nm)
        # record the Lz^2 term
        push!(cstr_l2, [1, o2, 0, o2, 1, o1, 0, o1])
        push!(fac_l2, (m1 - s) * (m2 - s))
        if m1 == nm - 1 continue end
        if m2 == 0 continue end 
        # record the L+L- term
        push!(cstr_l2, [1, o1 + 1, 0, o1, 1, o2 - 1, 0, o2])
        push!(fac_l2, sqrt(m2 * (nm - m2) * (m1 + 1) * (nm - m1 - 1)))
    end
end
# Initialise the L2 operator
l2 = Operator(bs, bs, 1, 1, cstr_l2, fac_l2)
@time "Initialise L2" l2_mat = OpMat(l2)
# Calculate the inner product for each eigenstate
@time "Measure L2" l2_val = [ st[:, i]' * l2_mat * st[:, i] for i = 1 : length(enrg)]
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
@time "Initialise Basis Z" bs1 = Basis(cfs, qnz_s1, cyc, perm_o, ph_o, fac_o) 
@show bs1.dim 
hmt = Operator(bs1, bs1, 1, 1, cstr_hmt, fac_hmt) # Generate and diagonalise Hamiltonian in the new basis
@time "Initialise Hamiltonian" hmtmat = OpMat(hmt)
@show hmtmat.nel
@time "Diagonalise Hamiltonian" enrg1, st1 = GetEigensystem(hmtmat, 10)
@show real(enrg1)
# Record the identity, sigma and epsilon states 
st_I = st[:, 1] # ground state
st_e = st[:, 2] # epsilon state
st_s = st1[:, 1]

# Record the density operator n^z
cstr_nz = []
fac_nz = Array{ComplexF64, 1}(undef, 0)
for o1u = 1 : nm
    o1d = o1u + nm
    push!(cstr_nz, [1, o1u, 0, o1u])
    push!(fac_nz, 1 / nm)
    push!(cstr_nz, [1, o1d, 0, o1d])
    push!(fac_nz, -1 / nm)
end
# The nz operator sends a state in bs (+) to bs1 (-)
nz = Operator(bs, bs1, 1, 0, cstr_nz, fac_nz)
# Measuring the finite size OPE
@show abs((st_s' * nz * st_e) / (st_s' * nz * st_I))