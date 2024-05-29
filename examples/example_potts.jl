using FuzzifiED
using LinearAlgebra
using PrettyTables
const ω = (-1 + √3im) / 2
const mat_Z = diagm([1, ω, ω'])
const mat_h = [ 0 1 1 ; 1 0 1 ; 1 1 0 ]

nm = 7
nf = 3
cfs = GetLzConfs(nm, nf, nm)
ps_pot = [4.75, 1.]
h = 10.0
tms_hmt = GetDenIntTerms(nm, nf ; ps_pot) - GetDenIntTerms(nm, nf ; ps_pot, mat_a = mat_Z) - h * GetPolTerms(nm, nf ; mat = mat_h)
tms_l2 = GetL2Terms(nm, nf)

result = []
Z2(x) = x == 1 ? "+" : x == -1 ? "-" : "0"
Z3(x) = x == 1 ? 0 : 1
for (qn_r, (qn_z2, qn_z3)) in Iterators.product([-1, 1], [(1, 1), (-1, 1), (0, ω)])
    @show qn_r, qn_z2, qn_z3
    bs = GetS3Basis(cfs ; qn_z3, qn_z2, qn_r)

    hmt = Operator(bs, bs, tms_hmt ; red_q = 0, sym_q = 0)
    @time "Generate sparse matrix" hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = Operator(bs, bs, tms_l2 ; red_q = 0, sym_q = 0)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i = 1 : length(enrg)]

    for i in eachindex(enrg) 
        push!(result, [real(enrg[i]), round(real(l2_val[i]), digits = 6), Z2(qn_z2), Z3(qn_z3), Z2(qn_r)])
    end
end

sort!(result)
enrg_0 = result[1][1]
enrg_T = filter(st -> abs(st[2] - 6) < eps(Float32) && st[3] == "+" && st[4] == 0, result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
result_display = permutedims(hcat(result_dim...))

pretty_table(result_display, header = ["Dim", "Energy", "L2", "Z2", "Z3", "Ry"])