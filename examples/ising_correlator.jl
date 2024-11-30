# This example calculates the σσ two-point function on sphere
# and the σσσσ four-point function on sphere, 0 and ∞. 
# This example reproduces Figures 1c and 2a in Phys. Rev. B 108, 235123 (2023)
# On my table computer, this calculation takes 1.650 s

using FuzzifiED
using LegendrePolynomials
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 12
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) ]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [1, -1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) ]

tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, 2 .* [4.75, 1.], σ1, σ2) - 
    3.16 * GetPolTerms(nm, 2, σx) )
tms_l2 = GetL2Terms(nm, 2)

cfs = Confs(2 * nm, [nm, 0], qnd)
bsp  = Basis(cfs, [1, 1, 1], qnf)
bsp1 = Basis(cfs, [1, 1,-1], qnf)
bsm  = Basis(cfs, [1,-1, 1], qnf)
bsm1 = Basis(cfs, [1,-1,-1], qnf)

hmt = Operator(bsp, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 20)
stg = st[:, 1]

hmt = Operator(bsm, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 20)
stσ = st[:, 1]

obs_nz = StoreComps(GetDensityObs(nm, 2, σz))
nzl0p = Dict([ l => 
    Operator(bsp, iseven(l) ? bsm : bsm1, GetComponent(obs_nz, l, 0.0)) 
for l = 0 : nm - 1])
ovl_σz0 = stσ' * nzl0p[0] * stg
cor_gzzg_l = [ begin 
    st_zg = nzl0p[l] * stg
    st_zg' * st_zg / abs(ovl_σz0) ^ 2
end for l = 0 : nm - 1]

nzl0m = Dict([ l => 
    Operator(bsm, iseven(l) ? bsp : bsp1, GetComponent(obs_nz, l, 0.0)) 
for l = 0 : nm - 1])
cor_σzzσ_l = [ begin 
    st_zσ = nzl0m[l] * stσ
    st_zσ' * st_zσ / abs(ovl_σz0) ^ 2
end for l = 0 : nm - 1]

Cor_gzzg(θ :: Float64) = [Pl(cos(θ), l) * (2 * l + 1) for l in 0 : nm - 1]' * cor_gzzg_l
Cor_σzzσ(θ :: Float64) = [Pl(cos(θ), l) * (2 * l + 1) for l in 0 : nm - 1]' * cor_σzzσ_l

display(permutedims(hcat([[θ / π, Cor_gzzg(θ), Cor_σzzσ(θ)] for θ = 0 : π / 10 : 2 * π]...)))
# To plot the correlation functions, uncomment the following lines 
# using Plots 
# θs = 0 : π / 10 : 2 * π
# plot_gzzg = plot(θs, Cor_gzzg.(θs))
# plot_σzzσ = plot(θs, Cor_σzzσ.(θs))