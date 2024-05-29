""" 
    function GetIsingQnu(nm :: Int64) :: @NamedTuple{qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}}

returns the diagonal quantum numbers, _i.e._, particle number ``N_e`` and angular momentum ``L_z``, of the fuzzy sphere Ising model. 

# Arguments 

- `nm :: Int64` is the number of orbitals ; 

# Output

A named tuple with three elements that can be directly fed into [`SitesFromQnu`](@ref)

- `qnu_o :: Vector{Vector{Int64}}` stores the charge of each orbital under each conserved quantity. See [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)))) for detail.
- `qnu_name :: Vector{String}` stores the name of each quantum number.
- `modul :: Vector{Int64}` stores the modulus of each quantum number, 1 if no modulus. 

"""
function GetIsingDefQnu(nm :: Int64 ; def_conf :: Vector{Int64} = [1, 1])
    nf = 2
    no = nf * nm
    qnu_o = []
    push!(qnu_o, fill(1, no)) 
    push!(qnu_o, [ div(o - 1, nf) for o = 1 : no ]) 
    qnu_name = ["N_e", "L_z"]
    if (def_conf[1] != 0)
        push!(qnu_o, [ 1 ; 0 ; fill(0, no - 4) ; 0 ; 0 ]) 
        push!(qnu_o, [ 0 ; 1 ; fill(0, no - 4) ; 0 ; 0 ]) 
        append!(qnu_name, ["N_def_NUp", "N_def_NDn"])
    end
    if (def_conf[2] != 0)
        push!(qnu_o, [ 0 ; 0 ; fill(0, no - 4) ; 1 ; 0 ]) 
        push!(qnu_o, [ 0 ; 0 ; fill(0, no - 4) ; 0 ; 1 ]) 
        append!(qnu_name, ["N_def_SUp", "N_def_SDn"])
    end
    modul = [1, 1, 1, 1, 1, 1]
    return (qnu_o = qnu_o, qnu_name = qnu_name, modul = modul)
end

"""
    function GetIsingConfs(nm :: Int64 ; ne :: Int64 = 0, lz :: Float64) :: Confs

Return the configurations with conserved particle number ``N_e`` and angular momentum ``L_z``.
        
# Arguments

- `nm :: Int64` is the number of orbitals ``2s+1``.
- `ne :: Int64` is the number of electrons.
- `lz :: Float64` is the angular momentum. Facultive, 0 by default. 
"""
function GetIsingDefConfs(nm :: Int64, ne :: Int64 ; lz :: Float64 = 0.0, def_conf :: Vector{Int64} = [1, 1])
    nf = 2
    no = nf * nm
    s = .5 * (nm - 1)
    qnu_s = Int64[ne ; ne * s + lz]
    if def_conf[1] != 0
        append!(qnu_s, def_conf[1] == 1 ? [1, 0] : [0, 1])
    end
    if def_conf[2] != 0
        append!(qnu_s, def_conf[2] == 1 ? [1, 0] : [0, 1])
    end
    qnu = GetIsingDefQnu(nm ; def_conf)
    # Generate the configurations and print the number
    return Confs(no, qnu_s, qnu.qnu_o)
end

"""
    function GetIsingIntTerms(nm :: Int64, ps_pot :: Vector) :: Vector{Term}

Returns the terms for the ising interaction 

```math
\\sum_{m_1m_2m_3m_4}2U_{m_1m_2m_3m_4}c^{\\dagger}_{m_1\\uparrow}c^{\\dagger}_{m_2\\downarrow}c_{m_3\\downarrow}c_{m_4\\uparrow}
```

from the pseudopotentials. 

# Arguments 

- `nm :: Int64` is the number of orbitals ``2s+1``.
- `ps_pot :: Vector{Number}` is the pseudopotential of Ising interaction.
"""
function GetIsingDefIntTerms(nm :: Int64, ps_pot :: Vector ; def_conf :: Vector{Int64} = [1, 1] )
    nf = 2
    no = nm * 2
    int_el = GetIntMatrix(nm, ps_pot)
    tms_ising = Vector{Term}(undef, 0)
    def_orb = []
    if def_conf[1] != 0 push!(def_orb, 1) end
    if def_conf[2] != 0 push!(def_orb, nm) end
    # Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m3
    for m1 = 1 : nm
        f1 = 0
        o1 = (m1 - 1) * nf + f1 + 1
        for m2 = 1 : nm
            f2 = 1
            o2 = (m2 - 1) * nf + f2 + 1
            for m3 = 1 : nm
                f3 = 1
                o3 = (m3 - 1) * nf + f3 + 1
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                if (m1 in def_orb || m4 in def_orb) 
                    if (m1 != m4) continue end
                end
                if (m2 in def_orb || m3 in def_orb) 
                    if (m2 != m3) continue end
                end
                f4 = 0
                o4 = (m4 - 1) * nf + f4 + 1
                push!(tms_ising, Term(int_el[m1, m2, m3] * 2., [1, o1, 1, o2, 0, o3, 0, o4]))
            end
        end
    end
    return tms_ising
end


function GetDefXPolTerms(nm :: Int64 ; def_conf :: Vector{Int64} = [1, 1])
    orb_d = 1 
    orb_f = 2 * nm
    if (def_conf[1] != 0) orb_d += 2 end
    if (def_conf[2] != 0) orb_f -= 2 end
    return [ Term(1, [1, isodd(i) ? i + 1 : i - 1, 0, i]) for i = orb_d : orb_f]
end