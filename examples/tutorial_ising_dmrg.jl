# In this tutorial, we show how to combine FuzzifiED and ITensor to do DMRG calculations on fuzzy sphere. 
# We show how to construct ITensor objects such as Sites and OpSum from FuzzifiED interfaces, 
# and how to convert ITensor objects back to Confs and Terms to benchmark DMRG results with ED. 

using FuzzifiED
using ITensors
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]

nm = 12
nf = 2
no = nm * nf

# Construct sites and MPO with 
sites = SitesFromQNDiag([
    GetNeQNDiag(nm * nf), 
    GetLz2QNDiag(nm, nf),
    GetZnfChargeQNDiag(nm, nf)
])
ps_pot = [4.75, 1.] ./ 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot) - 
    GetDenIntTerms(nm, 2 ; ps_pot, mat_a = σx) - 
    3.16 * GetPolTerms(nm, nf ; mat = σz)
)
os = OpSumFromTerms(tms_hmt)
@time mpo_hmt = MPO(os, sites)

cf0 = [ isodd(o) ? 1 : 0 for o = 1 : no ]
st0 = MPS(sites, string.(cf0))

# Calculate the ground state energy by DMRG
Eg, stg = dmrg(mpo_hmt, st0 ; nsweeps = 10, maxdim = [10,20,50,100,200,500], noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], cutoff = [1E-8])
@show Eg

# Convert sites and OpSum back to do ED and compare the results
cfs = ConfsFromSites(sites, cf0)
bs = Basis(cfs)
tms_hmt1 = TermsFromOpSum(os)
hmt = Operator(bs, bs, tms_hmt1 ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
@show enrg