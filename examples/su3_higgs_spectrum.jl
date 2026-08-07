# This example calculates the spectrum of the SU(3)_1-Higgs theory 
# realized as a transition between a ν=1/3 Laughlin state and a ν=3 fIQH state.

using FuzzifiED
using FuzzifiED.JackToolkit
using SO3lver
using LinearAlgebra
FuzzifiED.ElementType = Float64

function GetLaughlin13Jack(no :: Int64, sec :: Vector{Int64}, bs :: Basis, l2c2_mat :: OpMat, nst :: Int64)
    st0 = GetJackStates(bs, no, sec[1] ÷ 3, 1, 3, sec[2]) 
    l2c2_val, st = OrganizeJackStates(st0, l2c2_mat)
    @info "Sector $sec, # of states $(size(st, 2))"
    return l2c2_val, st
end

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
tms_lzlp = [GetLzLpTerms(nm, nf), GetLzLpTerms(nmf, 1)]
tms_c2 = GetC2Terms(nm, nf, :SU)
tms_proj = GetDenIntTerms(nmf, 1, [0, 1])

c_obs = [ GetElectronObs(nm, nf, f) for f = 1 : nf ]
f_obs = GetElectronObs(nmf, 1, 1)
ccc_obs = prod(c_obs)
nf_obs = GetDensityObs(nmf, 1)
nc_obs = GetDensityObs(nm, nf)

cpd_μ = SingleSegCouple(2, 1, GetPolTerms(nm, nf, Matrix(I, nf, nf)), [0, 0, 0, 0])
cpd_t = -ContactCouple([ccc_obs, f_obs'], [-3 3 ; 0 0 ; 0 0 ; 0 0]) +
    ContactCouple([ccc_obs', f_obs], [3 -3 ; 0 0 ; 0 0 ; 0 0])
cpd_U0 = 6 * ContactCouple([nc_obs, nf_obs], [0 0 ; 0 0 ; 0 0 ; 0 0]) + 
    SingleSegCouple(2, 1, GetIntegral(nc_obs * nc_obs), [0, 0, 0, 0])
cpd_hmt = cpd_U0 + cpd_t - 0.11 * cpd_μ

sec_f = stack([ [ne, ((nmf + 1) * ne) % 2, 0, 0] for ne = 0 : nf : noc])
nst_max = [ length(GetJackRoots(nmf, sec_f[1, i] ÷ 3, 1, 3, sec_f[2, i])) for i in axes(sec_f, 2) ]
sgsp_f = BuildSegSpace(nmf, sec_f, qnd_pt[:, 2], tms_lzlp[2], tms_proj, [0.0] ; l2c2_ratio = 0.1/nm^2, diag_method = GetLaughlin13Jack, nst_max)
sgop_hmt = Matrix{SegOperator}(undef, 2, length(cpd_hmt))
sgop_hmt[2, :] = BuildSegOperators([sgsp_f], cpd_hmt ; p_rng = [2])

result = []
szc2s = [(0, [0]), (1, [3, 6]), (2, [8, 12])]
for (sz, c2s) in szc2s
    sec_c = stack([ [ne, ((nmf + 1) * ne) % 2, 2 * sz, 0] for ne = 0 : nf : noc])
    sec_tot = [ noc, 0, 2 * sz, 0 ]
    c2_rng = [ Float64[c2] for c2 in c2s]
    sgsps_c = BuildSegSpaces(noc, sec_c, qnd_pt[:, 1], tms_lzlp[1], tms_c2, c2_rng)
    for ic2 in eachindex(c2s)
        sgsp_c = sgsps_c[ic2]
        c2 = c2s[ic2]
        sgop_hmt[1, :] = BuildSegOperators([sgsp_c], cpd_hmt ; p_rng = [1])
        for l = 0 : 2
            cpsp = BuildCompSpace([sgsp_c, sgsp_f], sec_tot, 2l)
            cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
            enrg, st = GetEigensystem(cpop_hmt, 10 ; issymmetric = true)
            for i in eachindex(enrg)
                push!(result, [enrg[i] / √(2l + 1), l, c2])
            end
        end
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2 && st[3] ≈ 0, result)[1][1]
spec = [ round.([ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(stack(spec)))
