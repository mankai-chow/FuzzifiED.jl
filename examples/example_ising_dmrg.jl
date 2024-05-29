using FuzzifiED
using ITensors

nm = 12
nf = 2
no = nm * nf

sites = SitesFromQnu(; GetIsingXQnu(nm)...)
tms_hmt = GetIsingXIntTerms(nm, [4.75, 1.]) - 3.16 * GetZPolTerms(nm)
@time mpo_hmt = MPO(OpSumFromTerms(tms_hmt), sites)
cf0 = [ isodd(o) ? 1 : 0 for o = 1 : no ]
st0 = MPS(sites, string.(cf0))

# Calculate the ground state energy by DMRG
Eg, stg = dmrg(mpo_hmt, st0 ; nsweeps = 10, maxdim = [10,20,50,100,200,500], noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], cutoff = [1E-8])
@show Eg

# Use the same objects to do ED and compare the results
cfs = ConfsFromSites(sites, cf0)
bs = Basis(cfs)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
@show enrg