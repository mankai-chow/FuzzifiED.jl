using FuzzifiED
using SO3lver
FuzzifiED.ElementType = Float64

# Set-up
nm = 12
ne = nm

# Bipartite into two segments of spin-up and spin-down
# Record the QNDiags, possible sectors for a segment, and the total sector
qnd_pt = [
    GetNeQNDiag(nm),
    GetLz2QNDiag(nm, 1)
]
sec_pt = stack([ [ne1, ((nm + 1) * ne1) % 2] for ne1 = 0 : ne]) # 2L^z consistent to be ℤ or ℤ+1/2
sec_tot = [ne, 0]
tms_lzlp_pt = GetLpLzTerms(nm, 1)

# Build the segment Hilbert space
sgsp = BuildSegSpace(nm, sec_pt, qnd_pt, tms_lzlp_pt)
l = 0 # total angular momentum specified by twice the value 2l
# Assemble the composed Hilbert space
cpsp = BuildCompSpace([sgsp, sgsp], sec_tot, 2l)

# Decompose the Hamiltonian
n_mod = GetDensityMod(nm, 1, [1;;]) 
cpd_int = 2 * CoupleDecomps([ n_mod, n_mod ], ConvPsPot(RecoupleAngMom((nm-1)/2, [4.75, 1.0]))..., [ 0 0 ; 0 0 ]) # Ising interaction from pseudo-potential re-coupling
c_obs = GetElectronObs(nm, 1, 1) 
cpd_nx = ContactCouple([c_obs', c_obs], [ 1 -1 ; 0 0 ]) - ContactCouple([c_obs, c_obs'], [ -1 1 ; 0 0 ]) # Transverse field from contact coupling
cpd_hmt = cpd_int - 3.16 * cpd_nx

# Build the segment operators for each channel in the Hamiltonian
sgop_hmt = BuildSegOperators([sgsp, sgsp], cpd_hmt ; ident_seg = [1, 1])
# Assemble the segment operators into the composed operator for the Hamiltonian
cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
# Diagonalise for the 10 lowest states
enrg, st = GetEigensystem(cpop_hmt, 10 ; issymmetric = true)
enrg /= √(2l + 1) # Convert ⟨l‖H‖l⟩ into energy E = ⟨lm|H|lm⟩
display(enrg)

# This calculation do not resolve ℤ_2, so I, σ, and ϵ are obtained within the same Hilbert space 
stI = st[:, 1] 
stσ = st[:, 2] 
stϵ = st[:, 3]
# Measure f_{σσϵ}=⟨σ|n^z_{00}|ϵ⟩ / ⟨σ|n^z_{00}|𝕀⟩
nz_00 = GetComponent(GetDensityObs(nm, 1), 0, 0)
cpd_nz = SingleSegCouple(2, 1, nz_00, [0, 0]) - SingleSegCouple(2, 2, nz_00, [0, 0])
cpop_nz = BuildCompOperator(cpsp, cpd_nz)
f_σσϵ = abs(stσ' * cpop_nz * stϵ) / abs(stσ' * cpop_nz * stI)
@show f_σσϵ
