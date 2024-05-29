
function GetLzQnu(nm :: Int64, nf :: Int64)
    no = nf * nm
    qnu_o = []
    push!(qnu_o, fill(1, no)) 
    push!(qnu_o, [ div(o - 1, nf) for o = 1 : no ]) 
    qnu_name = ["N_e", "L_z"]
    modul = [1, 1]
    return (qnu_o = qnu_o, qnu_name = qnu_name, modul = modul)
end

function GetLzConfs(nm :: Int64, nf :: Int64, ne :: Int64 ; lz :: Float64 = 0.0)
    no = nf * nm
    s = .5 * (nm - 1)
    qnu_s = Int64[ne, ne * s + lz]
    qnu = GetLzQnu(nm, nf)
    # Generate the configurations and print the number
    return Confs(no, qnu_s, qnu.qnu_o)
end

function GetS3Basis(cfs :: Confs ; qn_z3 :: Number = 0, qn_z2 :: Int64 = 0, qn_r = 0)
    no = cfs.no
    nf = 3
    nm = no ÷ nf
    qn_r1 = qn_r
    if (nm % 4 >= 2) qn_r1 = -qn_r end
    cyc = Vector{Int64}(undef, 0)
    qnz_s = Vector{ComplexF64}(undef, 0)
    perm_o = []
    ph_o = []
    fac_o = []
    if qn_z3 != 0
        push!(perm_o, [begin 
                f1 = mod(o - 1, nf)
                o + (f1 == 0 ? 2 : -1)
            end for o = 1 : no])
        push!(ph_o, fill(0, no))
        push!(fac_o, fill(ComplexF64(1), no))
        push!(qnz_s, qn_z3)
        push!(cyc, 3)
    end
    if qn_z2 != 0
        push!(perm_o, [begin 
                f1 = mod(o - 1, nf)
                o + (f1 == 0 ? 1 : (f1 == 1 ? -1 : 0))
            end for o = 1 : no])
        push!(ph_o, fill(0, no))
        push!(fac_o, fill(ComplexF64(1), no))
        push!(qnz_s, qn_z2)
        push!(cyc, 2)
    end
    if qn_r != 0
        push!(perm_o, [begin 
                f1 = mod(o - 1, nf) ; m1 = div(o - 1, nf)
                1 + f1 + nf * (nm - 1 - m1)
            end for o = 1 : no])
        push!(ph_o, fill(0, no)) 
        push!(fac_o, fill(ComplexF64(1), no)) 
        push!(qnz_s, qn_r1)
        push!(cyc, 2)
    end
    # Generate the basis and print the dimension
    return Basis(cfs, qnz_s, cyc, perm_o, ph_o, fac_o)
end

function GetDenIntTerms(nm :: Int64, nf :: Int64 ; ps_pot :: Vector = [1.0], mat_a :: Matrix = Matrix{Float64}(I, nf, nf), mat_b :: Matrix = Matrix(mat_a'))
    no = nm * nf
    int_el = GetIntMatrix(nm, ps_pot)
    tms = Vector{Term}(undef, 0)
    # Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m3
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for o2 = 1 : no 
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf) + 1
            if (o1 == o2) continue end
            # if (f1 < f2) continue end # f1 >= f2
            # if (f1 == f2 && m1 <= m2) continue end 
            for o3 = 1 : no
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf) + 1
                if (abs(mat_b[f2, f3]) < 1E-13) continue end 
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                for f4 = 1 : nf 
                    if (abs(mat_a[f1, f4]) < 1E-13) continue end
                    o4 = (m4 - 1) * nf + f4
                    if (o3 == o4) continue end
                    val = mat_a[f1, f4] * mat_b[f2, f3] * int_el[m1, m2, m3]
                    if (abs(val) < 1E-15) continue end 
                    push!(tms, Term(val, [1, o1, 1, o2, 0, o3, 0, o4]))
                end
            end
        end
    end
    return tms
end

function GetPolTerms(nm :: Int64, nf :: Int64 ; mat :: Matrix = Matrix{Float64}(I, nf, nf))
    no = nm * nf
    tms = Vector{Term}(undef, 0)
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for f2 = 1 : nf 
            if abs(mat[f1, f2]) < 1E-13 continue end
            o2 = (m1 - 1) * nf + f2
            push!(tms, Term(mat[f1, f2], [1, o1, 0, o2]))
        end
    end
    return tms
end