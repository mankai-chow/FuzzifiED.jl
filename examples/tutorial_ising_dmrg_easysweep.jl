# In this tutorial, we show how to use the tools EasySweep to manage DMRG sweeps and store intermediate results 
# so that the task can be restored if they are interrupted due to the time limit.
# We calculate and store the MPO for Hamiltonian and angular momentum, 
# and calculate two Z_2-even states, ground state and epsilon, and one Z_2 odd state sigma. 
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

path = "nm_$(nm)/"
mkpath(path)

ps_pot = [4.75, 1.] ./ 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, ps_pot) - 
    GetDenIntTerms(nm, 2, ps_pot, σx) - 
    3.16 * GetPolTerms(nm, 2, σz)
)
qnd = [ 
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf), 
    GetZnfChargeQNDiag(nm, nf) 
]
hmt, sites = GetMPOSites("hmt", tms_hmt, qnd ; path, mpo_method = MyMPO)

cf0 = [ isodd(o) ? 1 : 0 for o = 1 : no ]
st0 = MPS(sites, string.(cf0))

Eg, stg = EasySweep("g", hmt, st0 ; path)
Ee, ste = EasySweep("e", hmt, st0 ; path, proj = ["g"])

cf1 = cf0
cf1[1] = 0
cf1[2] = 1
st1 = MPS(sites, string.(cf1))
Es, sts = EasySweep("s", hmt, st1 ; path)

tms_l2 = GetL2Terms(nm, 2)
l2 = GetMPO("l2", tms_l2, sites ; path, mpo_method = MyMPO)
@show inner(stg', l2, stg)