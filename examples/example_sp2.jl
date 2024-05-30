using FuzzifiED

nm = 7
nf = 4
cfs = GetSpnConfs(nm, nf, 2 * nm)
tms_hmt = SimplifyTerms(GetDenIntTerms(nm, nf) + .9 * GetSpnPairIntTerms(nm, nf))
tms_l2 = GetL2Terms(nm, nf)
tms_c2 = GetSpnC2Terms(nm, nf)

result = []
ds(x) = x == 1 ? "+" : x == -1 ? "-" : "0"
for (qn_p, qn_r, (qn_z, qn_x)) in Iterators.product([1, -1], [1, -1], [([1, 1], [1]), ([1, 1], [-1]), ([1, -1], [0]), ([-1, -1], [1]), ([-1, -1], [-1])])
    @show qn_p, qn_r, qn_z, qn_x

    bs = GetSpnBasis(cfs, nf ; qn_p, qn_r, qn_z, qn_x)
    hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
    hmt_mat = OpMat(hmt ; type = Float64)
    enrg, st = GetEigensystem(hmt_mat, 10)
    
    l2 = Operator(bs, bs, tms_l2 ; red_q = 1, sym_q = 1)
    l2_mat = OpMat(l2 ; type = Float64)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i = 1 : length(enrg)]
    
    c2 = Operator(bs, bs, tms_c2 ; red_q = 1, sym_q = 1)
    c2_mat = OpMat(c2 ; type = Float64)
    c2_val = [ st[:, i]' * c2_mat * st[:, i] for i = 1 : length(enrg)]

    for i in eachindex(enrg) 
        push!(result, [enrg[i] ; round(l2_val[i], digits = 6) ; round(c2_val[i], digits = 6) ; ds.([qn_p ; qn_r ; qn_z ; qn_x])])
    end
end

sort!(result)
enrg_0 = result[1][1]
enrg_T = filter(st -> abs(st[2] - 6) < eps(Float32) && abs(st[3]) < eps(Float32), result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
result_display = permutedims(hcat(result_dim...))
display(result_display)
# For better output, change the last line to 
# using PrettyTables
# pretty_table(result_display, header = ["Dim", "Energy", "L2", "C2", "PH", "RY", "Z1", "Z2", "X"])