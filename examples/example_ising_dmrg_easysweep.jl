using FuzzifiED
using ITensors
using ITensorMPOConstruction

function MyMPO(os, sites)
    operatorNames = [ "I", "C", "Cdag", "N" ]
    opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
    return mpo = MPO_new(os, sites ; basisOpCacheVec = opCacheVec)
end

nm = 12
nf = 2
no = nm * nf

path = "nm_$(nm)/"
mkpath(path)

sigma_x = [ 0 1 ; 1 0 ]
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot = [4.75, 1.]) - 
    GetDenIntTerms(nm, 2 ; ps_pot = [4.75, 1.], mat_a = sigma_x) - 
    3.16 * GetZPolTerms(nm)
)
hmt, sites = GetMpoSites("hmt", tms_hmt ; path, GetLzZnQnu(nm, 2)..., mpo_method = MyMPO)

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
l2 = GetMpo("l2", tms_l2, sites ; path, mpo_method = MyMPO)
@show inner(stg', l2, stg)