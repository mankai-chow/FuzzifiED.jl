# This example calculates the 1-pt function σ and 2-pt function σϕ of magnetic line defect in 3d Ising model.
# The normalisation of the correlators require bulk data ovl_σzI in `ising_ope.jl`.
# This example reproduces Figure 4 in Nat. Commun. 15, 3659 (2024)
# On my table computer, this calculation takes 0.332 s

using FuzzifiED
using SphericalHarmonics
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

nm = 12
no = nm * 2
qnd = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, 2), 
    GetPinOrbQNDiag(no, [1, no - 1]), 
    GetPinOrbQNDiag(no, [2, no])
]
qnf = [
    GetRotyQNOffd(nm, 2), 
    GetParityQNOffd(nm, 2, [2,1], [-1,1])
]
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, 2 .* [4.75, 1.], σ1, σ2) - 
    3.16 * GetPolTerms(nm, 2, σx) )
    
cfs = Confs(no, [nm, 0, 2, 0], qnd)
bs = Basis(cfs, [1,1], qnf)

hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 3)
stI = st[:, 1]
stϕ = st[:, 2]

obs_nz = StoreComps(GetDensityObs(nm, 2, σz))
nzl0p = Dict([ l => 
    Operator(bs, GetComponent(obs_nz, l, 0.0)) 
for l = 0 : 2 : nm - 1])
ovl_IzI = stI' * nzl0p[0] * stI
cor_IzI_l = [ stI' * nzl0p[l] * stI for l = 0 : 2 : nm - 1]
cor_ϕzϕ_l = [ stI' * nzl0p[l] * stϕ for l = 0 : 2 : nm - 1]
Cor_IzI(θ :: Float64) = [ real(computeYlm(θ, 0, lmax = nm - 1)[(l, 0)]) for l = 0 : 2 : nm - 1]' * cor_IzI_l
Cor_ϕzϕ(θ :: Float64) = [ real(computeYlm(θ, 0, lmax = nm - 1)[(l, 0)]) for l = 0 : 2 : nm - 1]' * cor_ϕzϕ_l

display(permutedims(hcat([[θ / π, Cor_IzI(θ), Cor_ϕzϕ(θ)] for θ = 0 : π / 20 : π]...)))
# To plot the correlation functions, uncomment the following lines 
# using Plots 
# θs = 0 : π / 20 : π
# plot_IzI = plot(θs, Cor_IzI.(θs))
# plot_ϕzϕ = plot(θs, Cor_ϕzϕ.(θs))