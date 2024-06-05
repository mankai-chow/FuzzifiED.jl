using FuzzifiED
using LinearAlgebra

# Configuration setup 
nm = 7
no = 4 * nm
ne = 2 * nm
s = .5 * (nm - 1)

# Diagonal quantum numbers 
qnu_s = Int64[ ne, Int(ne * s), ne, ne]
qnu_o = []
push!(qnu_o, fill(1, no)) 
push!(qnu_o, [ div(o - 1, 4) for o = 1 : no ]) 
push!(qnu_o, vcat([ mod(f1, 4) == 0 ? 2 : (mod(f1, 4) == 2 ? 0 : 1) for f1 = 0 : 3, m1 = 1 : nm]...))
push!(qnu_o, vcat([ mod(f1, 4) == 1 ? 2 : (mod(f1, 4) == 3 ? 0 : 1) for f1 = 0 : 3, m1 = 1 : nm]...))
cfs = Confs(no, qnu_s, qnu_o)

# Off-diagonal quantum numbers
cyc = fill(2, 4)
perm_o = []
ph_o = []
fac_o = []
push!(perm_o, [begin 
        f1 = mod(o - 1, 4)
        o + (f1 < 2 ? 2 : -2)
    end for o = 1 : no ])
push!(ph_o, fill(1, no))
push!(fac_o, [ ComplexF64(1) * (mod(o - 1, 4) < 2 ? 1 : -1) for o = 1 : no ])

push!(perm_o, [begin 
        f1 = mod(o - 1, 4) ; m1 = div(o - 1, 4)
        1 + f1 + 4 * (nm - 1 - m1)
    end for o = 1 : no])
push!(ph_o, fill(0, no)) 
push!(fac_o, fill(ComplexF64(1), no)) 

push!(perm_o, [begin 
        f1 = mod(o - 1, 4)
        o + 2 * (f1 == 0 ? 1 : (f1 == 2 ? -1 : 0))
    end for o = 1 : no])
push!(ph_o, fill(0, no))
push!(fac_o, [ ComplexF64(1) * (mod(o - 1, 4) == 0 ? -1 : 1) for o = 1 : no ])

push!(perm_o, [begin 
    f1 = mod(o - 1, 4)
    o + 2 * (f1 == 1 ? 1 : (f1 == 3 ? -1 : 0))
end for o = 1 : no])
push!(ph_o, fill(0, no))
push!(fac_o, [ ComplexF64(1) * (mod(o - 1, 4) == 1 ? -1 : 1) for o = 1 : no ])

# Hamiltonian
ps_pot_u = [ 1.0]
ps_pot_v = [-0.9]
Ω = [ i - j == 2 ? 1 : 0 for i = 1 : 4, j = 1 : 4 ]
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 4 ; ps_pot = ps_pot_u, mat_a = Matrix{Float64}(I, 4, 4)) + GetPairIntTerms(nm, 4 ; ps_pot = ps_pot_v, mat_a = Ω)
)
# Total Angular momentum
tms_l2 = GetL2Terms(nm, 4)
# Flavour quadratic Casimir
# Do not ask me why C_2 looks like this. It just works.
tms_c2 = SimplifyTerms(
    GetDenIntTerms(nm, 4 ; ps_pot = [isodd(m) ? -.5 : 0 for m = 1 : nm], mat_a = Matrix{Float64}(I, 4, 4)) + 
    GetPairIntTerms(nm, 4 ; ps_pot = [isodd(m) ? -1.0 : 0 for m = 1 : nm], mat_a = Ω)
)
push!(tms_c2, Term(ne + .25 * ne ^ 2, [-1, -1]))
tms_c2 = SimplifyTerms(tms_c2)

# Go through all the Z_2 sectors
result = []
ds(x) = x == 1 ? "+" : x == -1 ? "-" : "0"
for qn_p in [1, -1], qn_r in [1, -1], (qn_z1, qn_z2) in [(1, 1), (1, -1), (-1, 1)]
    @show qn_p, qn_r, qn_z1, qn_z2
    qnz_s = ComplexF64[ qn_p, qn_r, qn_z1, qn_z2 ]

    bs = Basis(cfs, qnz_s, cyc, perm_o, ph_o, fac_o)
    hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
    hmt_mat = OpMat(hmt ; type = Float64)
    enrg, st = GetEigensystem(hmt_mat, 10)
    
    l2 = Operator(bs, bs, tms_l2 ; red_q = 1, sym_q = 1)
    l2_mat = OpMat(l2 ; type = Float64)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
    
    c2 = Operator(bs, bs, tms_c2 ; red_q = 1, sym_q = 1)
    c2_mat = OpMat(c2 ; type = Float64)
    c2_val = [ st[:, i]' * c2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg) 
        push!(result, [enrg[i] ; round(l2_val[i], digits = 6) ; round(c2_val[i], digits = 6) ; ds.(qnz_s)])
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