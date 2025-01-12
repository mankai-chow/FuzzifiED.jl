# In this tutorial, we show how to combine FuzzifiED and ITensor to do DMRG calculations on fuzzy sphere. 
# We show how to construct ITensor objects such as Sites and OpSum from FuzzifiED interfaces.

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
cf0 = [ isodd(o) ? 1 : 0 for o = 1 : no ]
st0 = MPS(sites, string.(cf0))
# Initial ℤ₂-odd state
cf1 = [ (isodd(o) == o > 2) ? 1 : 0 for o = 1 : no ]
st1 = MPS(sites, string.(cf1))

# Ground state : lowest ℤ₂-even
EI, stI = dmrg(hmt, st0 ; 
    nsweeps = 10, 
    maxdim = [10,20,50,100,200,500], 
    noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], 
    cutoff = [1E-8])
# ϵ-state : second ℤ₂-even
Ee, ste = dmrg(hmt, [stI], st0 ; 
    nsweeps = 10, 
    maxdim = [10,20,50,100,200,500], 
    noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], 
    cutoff = [1E-8], 
    weight = 100)
# σ-state : lowest ℤ₂-odd
Es, sts = dmrg(hmt, st1 ; 
    nsweeps = 10, 
    maxdim = [10,20,50,100,200,500], 
    noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], 
    cutoff = [1E-8])

# Measure ground state angular momentum
tms_l2 = GetL2Terms(nm, 2)
l2 = MPO(OpSum(tms_l2))
val_l2I = inner(stI', l2, stI)

# Measure OPE coefficient f_{σσϵ}
obs_nx = GetDensityObs(nm, 2, sgx)
tms_nx00 = SimplifyTerms(GetComponent(obs_nx, 0.0, 0.0))
nx00 = MPO(OpSum(tms_nx00))
f_sse = abs(inner(sts', nx00, ste) / inner(sts', nx00, stI))