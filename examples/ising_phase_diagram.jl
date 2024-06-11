# This example calculates the phase diagram of fuzzy sphere Ising model
# by calculating the order parameter <M^2>. 
# On my table computer, this calculation takes 5.089 s

using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]
≊(x, y) = abs(x - y) < eps(Float32)

nm = 12
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) ]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [1, -1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) ]
cfs = Confs(2 * nm, [nm, 0], qnd)
bs = Basis(cfs, [1, 1, 1], qnf)

hm2_lst = Vector{Float64}[]
for h = 2.0 : 0.2 : 4.0
    tms_hmt = SimplifyTerms(
        GetDenIntTerms(nm, 2 ; ps_pot = 2 .* [4.75, 1.], mat_a = σ1, mat_b = σ2) - 
        h * GetPolTerms(nm, 2 ; mat = σx) )
    hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
    hmt_mat = OpMat(hmt ; type = Float64)
    enrg, st = GetEigensystem(hmt_mat, 5)
    st_g = st[:, 1]

    tms_m = GetPolTerms(nm, 2 ; mat = σz)
    tms_m2 = SimplifyTerms(tms_m * tms_m)
    m2 = Operator(bs, bs, tms_m2 ; red_q = 1, sym_q = 1)

    m2_val = st_g' * m2 * st_g
    push!(hm2_lst, [h, m2_val])
    @show h, m2_val
end

display(permutedims(hcat(hm2_lst...)))