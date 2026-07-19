using FuzzifiED
using FuzzifiEDFullRotation
using LinearAlgebra
FuzzifiED.ElementType = Float64
Base.zero(QNDiag, no) = QNDiag(zeros(Int64, no))

nm = 6
nf = 3
noc = nm * nf
nmf = nm * 3 - 2
FuzzifiED.ObsNormRadSq = nm

qnd_pt = [
    GetNeQNDiag(noc)  nf * GetNeQNDiag(nmf) ;
    GetLz2QNDiag(nm, nf)  GetLz2QNDiag(nmf, 1) ;
    GetFlavQNDiag(nm, nf, [1, -1, 0])  zero(QNDiag, nmf) ;
    GetFlavQNDiag(nm, nf, [1, 1, -2])  zero(QNDiag, nmf)
]
tms_lzlp = [GetLpLzTerms(nm, nf), GetLpLzTerms(nmf, 1)]
tms_c2 = GetC2Terms(nm, nf, :SU)
tms_proj = GetDenIntTerms(nmf, 1, [0, 1])

c_obs = [ GetElectronObs(nm, nf, f) for f = 1 : nf ]
f_obs = GetElectronObs(nmf, 1, 1)
ccc_obs = prod(c_obs)
nf_obs = GetDensityObs(nmf, 1)
nc_obs = GetDensityObs(nm, nf)

cpd_μ = SingleSegCouple(2, 1, GetPolTerms(nm, nf, Matrix(I, nf, nf)), [0, 0, 0, 0])
cpd_U1 = SingleSegCouple(2, 2, SimplifyTerms(GetIntegral(nf_obs * Laplacian(nf_obs))), [0, 0, 0, 0])
cpd_t = -ContactCouple([ccc_obs, f_obs'], [-3 3 ; 0 0 ; 0 0 ; 0 0]) +
    ContactCouple([ccc_obs', f_obs], [3 -3 ; 0 0 ; 0 0 ; 0 0])
cpd_U0 = 6 * ContactCouple([nc_obs, nf_obs], [0 0 ; 0 0 ; 0 0 ; 0 0]) + 
    SingleSegCouple(2, 1, GetIntegral(nc_obs * nc_obs), [0, 0, 0, 0]) + 
    9 * SingleSegCouple(2, 2, GetIntegral(nf_obs * nf_obs), [0, 0, 0, 0])

cpd_hmt = 0.5 * cpd_U0 + cpd_U1 + 0.5 * cpd_t + 0.085 * cpd_μ

sec_f = stack([ [ne, ((nmf + 1) * ne) % 2, 0, 0] for ne = 0 : nf : noc])
sgsp_f = BuildSegSpace(nmf, sec_f, qnd_pt[:, 2], tms_lzlp[2], tms_proj, [0.0])
result = []
for (c2, sz) in [(0,0), (3,1), (6,1), (8,2), (12,2)]
    sec_c = stack([ [ne, ((nmf + 1) * ne) % 2, 2 * sz, 0] for ne = 0 : nf : noc])
    sec_tot = [ noc, 0, 2 * sz, 0 ]
    sgsp_c = BuildSegSpace(noc, sec_c, qnd_pt[:, 1], tms_lzlp[1], tms_c2, [Float64(c2)])
    for l = 0 : 2
        cpsp = BuildCompSpace([sgsp_c, sgsp_f], sec_tot, 2l)
        @show c2, l, cpsp.dim 
        sgop_hmt = BuildSegOperators([sgsp_c, sgsp_f], cpd_hmt)
        cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
        enrg, st = GetEigensystem(cpop_hmt, 10)
        for i in eachindex(enrg)
            push!(result, [enrg[i] / √(2l + 1), l, c2])
        end
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2 && st[3] ≈ 0, result)[1][1]
spec = [ round.([ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(stack(spec)))
