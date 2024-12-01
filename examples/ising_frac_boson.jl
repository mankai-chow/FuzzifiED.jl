# This example calculates the spectrum of 3d Ising model on fuzzy sphere
# for bosons at fractional filling ν = 1/2.
# This example reproduces Figure 12a,b in arXiv:2411.15299.
# On my table computer, this calculation takes 10.186 s.
# We acknowlege Cristian Voinea for his help in reproducing the results. 

using FuzzifiED
using FuzzifiED.Fuzzifino
const σ1 = [ 1 0 ; 0 0 ]
const σ2 = [ 0 0 ; 0 1 ]
const σ0 = [ 1 0 ; 0 1 ]
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < eps(Float32)

let

ne = 7
nm = 2 * ne - 1
nf = 2 
nof = 1 # Fuzzifino can only deal with mixture of bosons and fermions, so we put a single orbital of fermion and keep it empty.
nob = nm * nf 
qnd = [
    SQNDiag([1], fill(1, nob)),
    SQNDiag([1], fill(0, nob)),
    SQNDiag([0], collect(0 : nm * nf - 1) .÷ nf .* 2 .- (nm - 1))
]
qnf = [
    SQNOffd([1], vcat([ [2, 1][f] + (m - 1) * nf for f = 1 : nf, m = 1 : nm]...)),
    SQNOffd([1], vcat([f + (nm - m) * nf for f = 1 : nf, m = 1 : nm]...), 
        ComplexF64[1], ComplexF64(-1) .^ (collect(0 : nm * nf - 1) .÷ nf))
]
cfs = SConfs(nof, nob, ne, [ne, 0, 0], qnd)

ps_pot_proj = [1.0]
ps_pot_int = 2 .* [0.0, 0.53, 0.09]
fld_h = 0.25
int_el_proj = GetIntMatrix(nm, ps_pot_proj)
int_el_int = GetIntMatrix(nm, ps_pot_int)

tms_hmt = STerm[]
for o1 = 1 : nob
    m1 = div(o1 - 1, nf) + 1
    f1 = mod(o1 - 1, nf) + 1
    for o2 = 1 : nob
        m2 = div(o2 - 1, nf) + 1
        f2 = mod(o2 - 1, nf) + 1
        for o3 = 1 : nob
            m3 = div(o3 - 1, nf) + 1
            f3 = mod(o3 - 1, nf) + 1
            m4 = m1 + m2 - m3 
            if (m4 <= 0 || m4 > nm) continue end
            for f4 = 1 : nf 
                o4 = (m4 - 1) * nf + f4
                val  = σ0[f1, f4] * σ0[f2, f3] * int_el_proj[m1, m2, m3]
                val += σ1[f1, f4] * σ2[f2, f3] * int_el_int[m1, m2, m3]
                if (abs(val) < 1E-15) continue end 
                tms_hmt += STerms(val, [1, -o1, 1, -o2, 0, -o3, 0, -o4])
            end
        end
    end
end
for o1 = 1 : nob 
    o2 = (o1 - 1) ⊻ 1 + 1
    tms_hmt += STerms(-fld_h, [1, -o1, 0, -o2])
end 
tms_hmt = SimplifyTerms(tms_hmt)

s = (nm - 1) / 2.0
tms_lz = 
    [ begin m = div(o - 1, nf)
        STerm(m - s, [1, -o, 0, -o])
    end for o = 1 : nob ]
tms_lp = 
    [ begin m = div(o - 1, nf)
        STerm(sqrt(m * (nm - m)), [1, -o, 0, -(o - nf)])
    end for o = nf + 1 : nob ]
tms_lm = tms_lp' 
tms_l2 = SimplifyTerms(tms_lz * tms_lz - tms_lz + tms_lp * tms_lm)

result = []
for Z in [1, -1], R in [1, -1]
    bs = SBasis(cfs, [Z, R], qnf)
    hmt = SOperator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20)

    l2 = SOperator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, round.([enrg[i], l2_val[i], Z], digits = 6))
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 6 && st[3] ≈ 1, result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(hcat(result_dim...)))

end