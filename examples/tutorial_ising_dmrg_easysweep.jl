# This tutorial contains the DMRG code that uses the EasySweep extension. It
# 1. calculates the lowest eigenstates in the symmetry sector L^z=0 and 𝒵=+,
# 2. measures their total angular momenta (for DMRG, the ground state only), and 
# 3. calcultes the OPE coefficient f_{σσϵ}=⟨σ|n^z_{00}|ϵ⟩/⟨σ|n^z_{00}|0⟩

# We recommend the unregistered package ITensorMPOConstruction to construct MPO. 
# The following command is used to install the package :
#     using Pkg; Pkg.add(url="https://github.com/ITensor/ITensorMPOConstruction.jl.git")

using FuzzifiED
using ITensors, ITensorMPS, HDF5
using ITensorMPOConstruction
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]

function MyMPO(os, sites)
    operatorNames = [ "I", "C", "Cdag", "N" ]
    opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
    return MPO_new(os, sites ; basis_op_cache_vec = opCacheVec)
end

nm = 12
nf = 2
no = nm * nf

path = "nm_$(nm)_tmp/"
mkpath(path)

# Input the Hamiltonian
ps_pot = [4.75, 1.] ./ 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, ps_pot) - 
    GetDenIntTerms(nm, 2, ps_pot, σx) - 
    3.16 * GetPolTerms(nm, 2, σz)
)
# Input the quantum numbers
qnd = [ 
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf), 
    GetZnfChargeQNDiag(nm, nf) 
]
hmt, sites = GetMPOSites("hmt", tms_hmt, qnd ; path, mpo_method = MyMPO)

# Initial ℤ₂-even state
cfi_p = [ isodd(o) ? 1 : 0 for o = 1 : no ]
sti_p = MPS(sites, string.(cfi_p))
# Initial ℤ₂-odd state
cfi_m = [ (isodd(o) == (o > 2)) ? 1 : 0 for o = 1 : no ]
sti_m = MPS(sites, string.(cfi_m))

# Generate the three state
E0, st0 = EasySweep("0", hmt, sti_p ; path)
Ee, ste = EasySweep("e", hmt, sti_p ; path, proj = ["0"])
Es, sts = EasySweep("s", hmt, sti_m ; path)

# Measure the angular momentum of the ground state
tms_l2 = GetL2Terms(nm, 2)
l2 = GetMPO("l2", tms_l2, sites ; path, mpo_method = MyMPO)
val_l20 = inner(st0', l2, st0)
@show val_l20

# Measure OPE coefficient f_{σσϵ}
obs_nx = GetDensityObs(nm, 2, σx)
tms_nx00 = SimplifyTerms(GetComponent(obs_nx, 0.0, 0.0))
nx00 = GetMPO("nx00", tms_nx00, sites ; path, mpo_method = MyMPO)
f_sse = abs(inner(sts', nx00, ste) / inner(sts', nx00, st0))
@show f_sse