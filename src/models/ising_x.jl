""" 
    function GetIsingQnu(nm :: Int64) :: @NamedTuple{qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}}

returns the diagonal quantum numbers, _i.e._, particle number ``N_e``, angular momentum ``L_z`` and spin-parity ``S_z``, of the fuzzy sphere Ising model. The Ising model is written in the basis of 
```math
c_+^\\dagger=\\frac{c_\\uparrow^\\dagger+c_\\downarrow^\\dagger}{\\sqrt{2}},\\quad c_-^\\dagger=\\frac{c_\\uparrow^\\dagger-c_\\downarrow^\\dagger}{\\sqrt{2}}
```
and the parity of the ``\\mathbb{Z}_2``-odd particle numbers is conserved. 

# Arguments 

- `nm :: Int64` is the number of orbitals ; 

# Output

A named tuple with three elements that can be directly fed into [`SitesFromQnu`](@ref)

- `qnu_o :: Vector{Vector{Int64}}` stores the charge of each orbital under each conserved quantity. See [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)))) for detail.
- `qnu_name :: Vector{String}` stores the name of each quantum number.
- `modul :: Vector{Int64}` stores the modulus of each quantum number, 1 if no modulus. 

"""
function GetIsingXQnu(nm :: Int64)
    nf = 2
    no = nf * nm
    qnu_o = []
    push!(qnu_o, fill(1, no)) 
    push!(qnu_o, [ div(o - 1, nf) for o = 1 : no ]) 
    push!(qnu_o, [ isodd(o) ? 1 : 0 for o = 1 : no ]) 
    qnu_name = ["N_e", "L_z", "S_z"]
    modul = [1, 1, 2]
    return (qnu_o = qnu_o, qnu_name = qnu_name, modul = modul)
end

"""
    function GetIsingConfs(nm :: Int64 ; ne :: Int64 = 0, lz :: Float64) :: Confs

Return the configurations with conserved particle number ``N_e``, angular momentum ``L_z`` and spin-parity ``S_z`` of the fuzzy sphere Ising model. 
        
# Arguments

- `nm :: Int64` is the number of orbitals ``2s+1``.
- `ne :: Int64` is the number of electrons.
- `lz :: Float64` is the angular momentum. Facultive, 0 by default. 
"""
function GetIsingXConfs(nm :: Int64, ne :: Int64 ; lz :: Float64 = 0.0, sz :: Int64 = 0)
    nf = 2
    no = nf * nm
    s = .5 * (nm - 1)
    qnu_s = Int64[ne, ne * s + lz, sz]
    qnu = GetIsingQnuX(nm)
    # Generate the configurations and print the number
    return Confs(no, qnu_s, qnu.qnu_o)
end

"""
    function GetIsingIntTerms(nm :: Int64, ps_pot :: Vector) :: Vector{Term}

Returns the terms for the ising interaction.

```math
\\sum_{m_1m_2m_3m_4}2U_{m_1m_2m_3m_4}c^{\\dagger}_{m_1\\uparrow}c^{\\dagger}_{m_2\\downarrow}c_{m_3\\downarrow}c_{m_4\\uparrow}
```
written in the basis of ``c_\\pm`` from the pseudopotentials. 

# Arguments 

- `nm :: Int64` is the number of orbitals ``2s+1``.
- `ps_pot :: Vector{Number}` is the pseudopotential of Ising interaction.
"""
function GetIsingXIntTerms(nm :: Int64, ps_pot :: Vector)
    nf = 2
    no = nm * 2
    int_el = GetIntMatrix(nm, ps_pot)
    tms_ising = Vector{Term}(undef, 0)
    # Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m3
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf)
        for o2 = o1 + 1 : no # consider only o1 < o2, o3 > o4
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf)
            for o3 = 1 : no 
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf)
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                f4 = mod(f1 + f2 + f3, 2)
                o4 = (m4 - 1) * nf + f4 + 1
                if (o3 <= o4) continue end
                V1 = int_el[m1, m2, m3] * (f1 == f4 ? 1 : -1) -int_el[m2, m1, m3] * (f2 == f4 ? 1 : -1)
                push!(tms_ising, Term(V1, [1, o1, 1, o2, 0, o3, 0, o4]))
            end
        end
    end
    return tms_ising
end