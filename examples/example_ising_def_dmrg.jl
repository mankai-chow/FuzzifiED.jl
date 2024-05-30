using FuzzifiED
using ITensors

nm = 12
nf = 2
no = nm * nf

path = "nm_$(nm)/"
mkpath(path)

tms_hmt = GetIsingDefIntTerms(nm, [4.75, 1.]) - 3.16 * GetDefXPolTerms(nm)
qnu = TruncateQnu(; GetIsingDefQnu(nm; def_conf = [1, 1])...)
hmt, sites = GetMpoSites("hmt", tms_hmt ; path, qnu...)

cf0 = [ isodd(o) ? 1 : 0 for o = 1 : no ]
st0 = MPS(sites, string.(cf0))

Epp, stpp = EasySweep("pp", hmt, st0 ; path)

cf0[1] = 0
cf0[2] = 1
st0 = MPS(sites, string.(cf0))

Epm, stpm = EasySweep("pm", hmt, st0 ; path)

os = OpSum()
os += "Cdag", 1, "C", 2 
conv = MPO(os, sites)
@show inner(stpp', conv, stpm)