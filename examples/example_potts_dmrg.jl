using FuzzifiED
using ITensors
using LinearAlgebra
using ITensorMPOConstruction
const mat_X = [0 0 1 ; 1 0 0 ; 0 1 0]

function MyMPO(os, sites)
    operatorNames = [ "I", "C", "Cdag", "N" ]
    opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
    return MPO_new(os, sites ; basisOpCacheVec = opCacheVec)
end

nm = 7
nf = 3
no = nm * nf

path = "temp_nm_$(nm)/"
mkpath(path)

ps_pot = [4.75, 1.]
h = 10.0
tms_hmt = GetDenIntTerms(nm, nf ; ps_pot)
            - GetDenIntTerms(nm, nf ; ps_pot, mat_a = mat_X)
            - h * GetPolTerms(nm, nf ; mat = diagm([2, -1, -1]))
hmt, sites = GetMpoSites("hmt", tms_hmt ; path, GetLzZnQnu(nm, nf)..., mpo_method = MyMPO)

cf0 = [ o % 3 == 1 ? 1 : 0 for o = 1 : no ]
st0 = MPS(sites, string.(cf0))
Eg, stg = EasySweep("g", hmt, st0 ; path)