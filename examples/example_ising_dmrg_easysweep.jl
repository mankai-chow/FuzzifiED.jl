using FuzzifiED
using ITensors

nm = 12
nf = 2
no = nm * nf

path = "nm_$(nm)/"
mkpath(path)

tms_hmt = GetIsingXIntTerms(nm, [4.75, 1.]) - 3.16 * GetZPolTerms(nm)
hmt, sites = GetMpoSites("hmt", tms_hmt ; path, GetIsingXQnu(nm)...)

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
l2 = GetMpo("l2", tms_l2, sites ; path, old = true)
@show inner(stg', l2, stg)