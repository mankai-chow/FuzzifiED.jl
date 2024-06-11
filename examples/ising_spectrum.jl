# This example calculates the spectrum of 3d Ising model on fuzzy sphere at nm = 12.
# For each (P,Z,R) sector, 20 states are calculated.
# On my table computer, this calculation takes 7.549 s

using FuzzifiED
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
≊(x, y) = abs(x - y) < eps(Float32)

nm = 12
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) ]
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) ]

tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot = 2 .* [4.75, 1.], mat_a = σ1, mat_b = σ2) - 
    3.16 * GetPolTerms(nm, 2 ; mat = σx) )
tms_l2 = GetL2Terms(nm, 2)

cfs = Confs(2 * nm, [nm, 0], qnd)

result = []
for P in [1, -1], Z in [1, -1], R in [1, -1]
    bs = Basis(cfs, [P, Z, R], qnf)
    hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
    hmt_mat = OpMat(hmt ; type = Float64)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = Operator(bs, bs, tms_l2 ; red_q = 1, sym_q = 1)
    l2_mat = OpMat(l2 ; type = Float64)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, round.([enrg[i], l2_val[i], P, Z], digits = 6))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≊ 6 && st[3] ≊ 1 && st[4] ≊ 1, result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
for P in (1, -1), Z in (1, -1)
    display(permutedims(hcat(
        filter(st -> st[4] ≊ P && st[5] ≊ Z, result_dim)...
    )))
end