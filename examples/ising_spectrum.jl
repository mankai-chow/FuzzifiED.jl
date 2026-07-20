using FuzzifiED
using FuzzifiEDFullRotation
FuzzifiED.ElementType = Float64

nm = 12
nf = 2
ne = nm

qnd_pt = [ 
    GetNeQNDiag(nm), 
    GetLz2QNDiag(nm, 1)
] 
sec_pt = stack([ [ne1, ((nm + 1) * ne1) % 2] for ne1 = 0 : ne])
sec_tot = [ne, 0]

n_mod = GetDensityMod(nm, 1, [1;;])
c_obs = GetElectronObs(nm, 1, 1)
cpd_int = 2 * CoupleDecomp([ n_mod, n_mod ], ConvPsPot(RecouplePsPot((nm-1)/2, [4.75, 1.0]))..., [ 0 0 ; 0 0 ])
cpd_nx = ContactCouple([c_obs', c_obs], [ 1 -1 ; 0 0 ]) - ContactCouple([c_obs, c_obs'], [ -1 1 ; 0 0 ])

cpd_hmt = PrepareCouple(cpd_int - 3.16 * cpd_nx)
tms_lzlp_pt = GetLpLzTerms(nm, 1) 

sgsp = BuildSegSpace(nm, sec_pt, qnd_pt, tms_lzlp_pt)
sgop_hmt = BuildSegOperators([sgsp, sgsp], cpd_hmt)
result = []
for l in 0 : 2
    ll = Int64(2l)
    cpsp = BuildCompSpace([sgsp, sgsp], sec_tot, ll)
    cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
    enrg, st = GetEigensystem(cpop_hmt, 10 ; issymmetric = true)
    for i in eachindex(enrg)
        push!(result, [enrg[i] / √(2l + 1), l])
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2, result)[2][1]
spec = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(stack(spec)))
