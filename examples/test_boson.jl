# This example tests the nearest-neighbour tight-binding model0 
# $H=\sum_i(b^\dagger_ib_{i+1}+f^\dagger_if_{i+1}+\mathrm{h.c.})$. 
# The example diagonalises the sector with the number of bosons and fermions both $N_o/2$
# and even under the reflection with respect to a bond center $i\mapsto N_o+1-i$, 
# and measures the total particle number squared $\left[\sum_i(b_i^\dagger b_i+f^\dagger_if_i)\right]^2$

using FuzzifiED
using FuzzifiED.Fuzzifino
FuzzifiED.ElementType = Float64

nof = 8
nob = nof
qnd = [
    SQNDiag(fill(1, nof), fill(1, nob)),
    SQNDiag(fill(1, nof), fill(0, nob))
]
qnf = [
    SQNOffd(collect(nof : -1 : 1), collect(nob : -1 : 1))
]
cfs = SConfs(nof, nob, nof ÷ 2, [nof, nof ÷ 2], qnd)

tms_hmt = SimplifyTerms(sum([
    [STerm(1, [1,-x, 0,-x % nob - 1]), STerm(1, [1,-x % nob - 1, 0,-x]),
     STerm(1, [1, x, 0, x % nob + 1]), STerm(1, [1, x % nob + 1, 0, x])]
    for x = 1 : nob
]))

bs = SBasis(cfs, [1], qnf)
hmt = SOperator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 20)
display(sort(enrg))

tms_n2 = SimplifyTerms(sum(
    [STerm(1, [1, o*s, 0, o*s])]
    for o = 1 : nof, s in [1, -1]
) ^ 2)
n2 = SOperator(bs, tms_n2)
@show st[:, 1]' * n2 * st[:, 1]
