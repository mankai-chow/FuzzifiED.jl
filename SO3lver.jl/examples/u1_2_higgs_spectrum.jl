# This example calculates the spectrum of the U(1)_2-Higgs theory realized as a transition between a ν=1/2 bosonic Laughlin state and a ν=2 fIQH state.

using FuzzifiED
using FuzzifiED.Fuzzifino
using FuzzifiED.JackToolkit
using SO3lver
FuzzifiED.ElementType = Float64

function GetLaughlin12Jack(nof :: Int64, nob :: Int64, sec :: Vector{Int64}, bs :: SBasis, l2c2_mat :: OpMat, nst :: Int64)
    st0 = GetJackStates(bs, nob, sec[1] ÷ 2, 1, 2, sec[2]) 
    l2c2_val, st = OrganizeJackStates(st0, l2c2_mat)
    @info "Sector $sec, # of states $(size(st, 2))"
    return l2c2_val, st
end

nmf = 7
nf = 2
nof = nmf * nf
nmb = 2 * nmf - 1
s = (nmf - 1) / 2
FuzzifiED.ObsNormRadSq = nmf

qnd_f = [ 
    GetNeQNDiag(nof),
    GetLz2QNDiag(nmf, nf),
    GetFlavQNDiag(nmf, nf, [1, -1])
]
qnd_b = [
    2 * GetBosonNeSQNDiag(0, nmb),
    GetBosonLz2SQNDiag(0, nmb, 1), 
    zero(SQNDiag, 0, nmb)
]

sec_f = stack([ [ne1, 0, 0] for ne1 = 0 : 2 : nof ])
sec_b = stack([ [ne1, 0, 0] for ne1 = 0 : 2 : nof ])
sec_tot = [nof, 0, 0]
nebm_b = [ne1 for ne1 = 0 : 2 : nof ]
tms_lzlp_f = GetLzLpTerms(nmf, nf)
tms_lzlp_b = GetBosonLzLpSTerms(nmb, 1)
tms_c2 = GetC2Terms(nmf, nf, :SU)
tms_proj = GetBosonDenIntSTerms(nmb, 1)

b = GetBosonSObs(nmb, 1, 1)
ff = GetElectronObs(nmf, 2, 1) * GetElectronObs(nmf, 2, 2)
nb = GetBosDensitySObs(nmb, 1)
nf = GetDensityObs(nmf, 2)

cpd_hop = ContactCouple([ff', b], [2 -2; 0 0; 0 0]) + ContactCouple([ff, b'], [-2 2; 0 0; 0 0])
cpd_μ = SingleSegCouple(2, 1, GetPolTerms(nmf, 2), [0, 0, 0])
cpd_int_e = SingleSegCouple(2, 1, GetIntegral(nf * nf), [0, 0, 0]) + 4 * SingleSegCouple(2, 2, GetIntegral(nb * nb), [0, 0, 0]) + 4 * ContactCouple([nf, nb], [0 0; 0 0; 0 0])

cpd_hmt = cpd_int_e - 0.5 * cpd_hop + 0.312 * cpd_μ 

nst_max = [ length(GetJackRoots(nmb, sec_b[1, i] ÷ 2, 1, 2, sec_b[2, i] ; fermion = false)) for i in axes(sec_b, 2) ]
sgsp_b = BuildSegSpace(0, nmb, nebm_b, sec_b, qnd_b, tms_lzlp_b, tms_proj, [0.0] ; l2c2_ratio = 0.1/nmf^2, diag_method = GetLaughlin12Jack, nst_max)
ss = collect(0 : 2)
c2_rng = [ Float64[s * (s + 1)] for s in ss]
sgsps_f = BuildSegSpaces(nof, sec_f, qnd_f, tms_lzlp_f, tms_c2, c2_rng)
sgop_hmt = Matrix{SegOperator}(undef, 2, length(cpd_hmt))
sgop_hmt[2, :] = BuildSegOperators([sgsp_b], cpd_hmt ; p_rng = [2])

result = []
for is in eachindex(ss)
    s = ss[is]
    sgsp_f = sgsps_f[is]
    sgop_hmt[1, :] = BuildSegOperators([sgsp_f], cpd_hmt ; p_rng = [1])
    for l = 0 : 2
        cpsp = BuildCompSpace([sgsp_f, sgsp_b], sec_tot, 2l)
        cpop_hmt = BuildCompOperator(cpsp, cpd_hmt, sgop_hmt)
        enrg, st = GetEigensystem(cpop_hmt, 10 ; issymmetric = true)
        for i in eachindex(enrg)
            push!(result, [enrg[i], l, s])
        end
    end
end

sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≈ 2 && st[3] ≈ 0, result)[1][1]
spec = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
display(permutedims(stack(spec)))
