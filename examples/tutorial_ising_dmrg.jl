# This tutorial contains the DMRG code that converts the format into ITensor. It
# 1. calculates the lowest eigenstates in the symmetry sector L^z=0 and 𝒵=+,
# 2. measures their total angular momenta (for DMRG, the ground state only), and 
# 3. calcultes the OPE coefficient f_{σσϵ}=⟨σ|n^z_{00}|ϵ⟩/⟨σ|n^z_{00}|0⟩

using FuzzifiED
using ITensors, ITensorMPS
FuzzifiED.ElementType = Float64
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]

nm = 12
nf = 2
no = nm * nf

# Construct sites and MPO with 
sites = GetSites([
    GetNeQNDiag(nm * nf), 
    GetLz2QNDiag(nm, nf),
    GetZnfChargeQNDiag(nm, nf)
])
ps_pot = [4.75, 1.] ./ 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, ps_pot) - 
    GetDenIntTerms(nm, 2, ps_pot, σx) - 
    3.16 * GetPolTerms(nm, nf, σz)
)
os_hmt = OpSum(tms_hmt)
@time hmt = MPO(os_hmt, sites)

# Initial ℤ₂-even state
cfi_p = [ isodd(o) ? 1 : 0 for o = 1 : no ]
sti_p = MPS(sites, string.(cfi_p))
# Initial ℤ₂-odd state
cfi_m = [ (isodd(o) == (o > 2)) ? 1 : 0 for o = 1 : no ]
sti_m = MPS(sites, string.(cfi_m))

# Ground state : lowest ℤ₂-even
E0, st0 = dmrg(hmt, sti_p ; 
    nsweeps = 10, 
    maxdim = [10,20,50,100,200,500], 
    noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], 
    cutoff = [1E-8])
# ϵ-state : second ℤ₂-even
Ee, ste = dmrg(hmt, [st0], sti_p ; 
    nsweeps = 10, 
    maxdim = [10,20,50,100,200,500], 
    noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], 
    cutoff = [1E-8], 
    weight = 100)
# σ-state : lowest ℤ₂-odd
Es, sts = dmrg(hmt, sti_m ; 
    nsweeps = 10, 
    maxdim = [10,20,50,100,200,500], 
    noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], 
    cutoff = [1E-8])

# Measure ground state angular momentum
tms_l2 = GetL2Terms(nm, 2)
l2 = MPO(OpSum(tms_l2), sites)
val_l20 = inner(st0', l2, st0)
@show val_l20

# Measure OPE coefficient f_{σσϵ}
obs_nx = GetDensityObs(nm, 2, σx)
tms_nx00 = SimplifyTerms(GetComponent(obs_nx, 0.0, 0.0))
nx00 = MPO(OpSum(tms_nx00), sites)
f_sse = abs(inner(sts', nx00, ste) / inner(sts', nx00, st0))
@show f_sse