using FuzzifiED
using ITensors

nm = 12
nf = 2
no = nm * nf

path = "nm_$(nm)/"
mkpath(path)

sigma_x = [ 0 1 ; 1 0 ]
tms_hmt = SimplifyTerms(GetDenIntTerms(nm, 2 ; ps_pot = [4.75, 1.]) - GetDenIntTerms(nm, 2 ; ps_pot = [4.75, 1.] ; mat_a = sigma_x) - 3.16 * GetZPolTerms(nm))
hmt, sites = GetMpoSites("hmt", tms_hmt ; path, GetLzZnQnu(nm, 2)...)

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
l2 = GetMpo("l2", tms_l2, sites ; path)
@show inner(stg', l2, stg)