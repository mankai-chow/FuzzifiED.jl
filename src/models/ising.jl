"""
    function GetLzQnu(nm :: Int64, nf :: Int64) :: @NamedTuple{qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}}

returns the diagonal quantum numbers, _i.e._, particle number ``N_e`` and angular momentum ``L_z+sN_e``
```math
\\begin{aligned}
    N_e&=∑_o n_o\\\\
    L_z+sN_e&=∑_{mf}(m+s)n_o
\\end{aligned}
```

# Arguments 

- `nm :: Int64` is the number of orbitals ; 
- `nf :: Int64` is the number of flavours ; 

# Output

A named tuple with three elements that can be directly fed into [`SitesFromQnu`](@ref)

- `qnu_o :: Vector{Vector{Int64}}` stores the charge of each orbital under each conserved quantity. See [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)))) for detail.
- `qnu_name :: Vector{String}` stores the name of each quantum number.
- `modul :: Vector{Int64}` stores the modulus of each quantum number, 1 if no modulus. 

"""
function GetLzQnu(nm :: Int64, nf :: Int64)
    no = nf * nm
    qnu_o = []
    push!(qnu_o, fill(1, no)) 
    push!(qnu_o, [ div(o - 1, nf) for o = 1 : no ]) 
    qnu_name = ["N_e", "L_z"]
    modul = [1, 1]
    return (qnu_o = qnu_o, qnu_name = qnu_name, modul = modul)
end

"""
    function GetLzZnQnu(nm :: Int64, nf :: Int64) :: @NamedTuple{qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}}

returns the diagonal quantum numbers, _i.e._, particle number ``N_e``, angular momentum ``L_z+sN_e`` and flavour charge ``Z_{N_f}``
```math
\\begin{aligned}
    N_e&=∑_o n_o\\\\
    L_z+sN_e&=∑_{mf}(m+s)n_{mf}\\\\
    Z_{N_f}&=∑_{m,f=0}^{N_f-1}fn_{mf}\\ \\mathrm{mod}\\ N_f
\\end{aligned}
```

# Arguments 

- `nm :: Int64` is the number of orbitals ; 
- `nf :: Int64` is the number of flavours ; 

# Output

A named tuple with three elements that can be directly fed into [`SitesFromQnu`](@ref)

- `qnu_o :: Vector{Vector{Int64}}` stores the charge of each orbital under each conserved quantity. See [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)))) for detail.
- `qnu_name :: Vector{String}` stores the name of each quantum number.
- `modul :: Vector{Int64}` stores the modulus of each quantum number, 1 if no modulus. 

"""
function GetLzZnQnu(nm :: Int64, nf :: Int64)
    no = nf * nm
    qnu_o = []
    push!(qnu_o, fill(1, no)) 
    push!(qnu_o, [ (o - 1) ÷ nf for o = 1 : no ]) 
    push!(qnu_o, [ (o - 1) % nf for o = 1 : no ]) 
    qnu_name = ["N_e", "L_z", "Z_n"]
    modul = [1, 1, nf]
    return (qnu_o = qnu_o, qnu_name = qnu_name, modul = modul)
end


"""
    function GetLzConfs(nm :: Int64, nf :: Int64, ne :: Int64 ; lz :: Float64) :: Confs

Return the configurations with conserved particle number ``N_e`` and angular momentum ``L_z``.
        
# Arguments

- `nm :: Int64` is the number of orbitals ``2s+1``.
- `nf :: Int64` is the number of flavours ; 
- `ne :: Int64` is the number of electrons.
- `lz :: Float64` is the angular momentum. Facultative, 0 by default. 
"""
function GetLzConfs(nm :: Int64, nf :: Int64, ne :: Int64 ; lz :: Float64 = 0.0)
    no = nf * nm
    s = .5 * (nm - 1)
    qnu_s = Int64[ne, ne * s + lz]
    qnu = GetLzQnu(nm, nf)
    # Generate the configurations and print the number
    return Confs(no, qnu_s, qnu.qnu_o)
end


"""
    function GetLzZnConfs(nm :: Int64, nf :: Int64, ne :: Int64 ; lz :: Float64, zn :: Int64 = 0) :: Confs

Return the configurations with conserved particle number ``N_e``, angular momentum ``L_z`` and flavour charge ``Z_{N_f}``.
        
# Arguments

- `nm :: Int64` is the number of orbitals ``2s+1``.
- `nf :: Int64` is the number of flavours ; 
- `ne :: Int64` is the number of electrons.
- `lz :: Float64` is the angular momentum. Facultative, 0 by default. 
- `zn :: Float64` is the flavour charge. Facultative, 0 by default. 
"""
function GetLzZnConfs(nm :: Int64, nf :: Int64, ne :: Int64 ; lz :: Float64 = 0.0, zn :: Int64 = 0)
    no = nf * nm
    s = .5 * (nm - 1)
    qnu_s = Int64[ne, ne * s + lz, zn]
    qnu = GetLzZnQnu(nm, nf)
    # Generate the configurations and print the number
    return Confs(no, qnu_s, qnu.qnu_o ; modul = qnu.modul)
end

"""
    function GetIsingQnz(nm :: Int64 ; qn_p :: Int64 = 0, qn_z :: Int64 = 0, qn_r :: Int64 = 0) :: @NamedTuple{cyc, perm_o, ph_o, fac_o}

Return the off-diagonal quantum numbers of the parity ``\\mathscr{P}``, flavour symmetry ``𝒵`` and ``π``-rotation along ``y``-axis ``\\mathscr{R}`` from the configurations already generated. Quantum numbers set to zero signify that they are not conserved. 

# Arguments

- `nm :: Int64` is the number of orbitals.
- `qn_p :: Int64` is quantum number for parity transformation. Facultative, 0 by default.
- `qn_z :: Int64` is the particle quantum number for ``ℤ_2``-flavour transformation. Facultative, 0 by default.
- `qn_r :: Int64` is the quantum number for  ``π`` rotation along ``y``-axis compared with the ground state. Facultative, 0 by default.
"""
function GetIsingQnz(nm :: Int64 ; qn_p :: Int64 = 0, qn_z :: Int64 = 0, qn_r :: Int64 = 0)
    no = nm * 2
    cyc = Vector{Int64}(undef, 0)
    perm_o = []
    ph_o = []
    fac_o = []
    if qn_p != 0
        push!(perm_o, [ isodd(o) ? o + 1 : o - 1 for o = 1 : no]) 
        push!(ph_o, fill(1, no))
        push!(fac_o, [ isodd(o) ? -1 : 1 for o = 1 : no])
        push!(cyc, 2)
    end
    if qn_z != 0
        push!(perm_o, [ isodd(o) ? o + 1 : o - 1 for o = 1 : no])
        push!(ph_o, fill(0, no))
        push!(fac_o, fill(ComplexF64(1), no))
        push!(cyc, 2)
    end
    if qn_r != 0
        push!(perm_o, [ isodd(o) ? no - o : no + 2 - o for o = 1 : no])
        push!(ph_o, fill(0, no)) 
        push!(fac_o, fill(ComplexF64(1), no)) 
        push!(cyc, 2)
    end
    # Generate the basis and print the dimension
    return (cyc = cyc, perm_o = perm_o, ph_o = ph_o, fac_o = fac_o)
end

"""
    function GetIsingBasis(cfs :: Confs ; qn_p :: Int64 = 0, qn_z :: Int64 = 0, qn_r :: Int64 = 0) :: Basis

Return the basis with conserved parity ``\\mathscr{P}``, flavour symmetry ``𝒵`` and ``π``-rotation along ``y``-axis ``\\mathscr{R}`` from the configurations already generated. Quantum numbers set to zero signify that they are not conserved. 

# Arguments

- `cfs :: Confs` is the configurations generated by [`GetLzConfs`](@ref) or [`GetLzZnConfs`](@ref).
- `qn_p :: Int64` is quantum number for parity transformation. Facultative, 0 by default.
- `qn_z :: Int64` is the particle quantum number for ``ℤ_2``-flavour transformation. Facultative, 0 by default.
- `qn_r :: Int64` is the quantum number for  ``π`` rotation along ``y``-axis compared with the ground state. Facultative, 0 by default.
"""
function GetIsingBasis(cfs :: Confs ; qn_p :: Int64 = 0, qn_z :: Int64 = 0, qn_r :: Int64 = 0)
    nm = cfs.no ÷ 2
    qn_r1 = qn_r
    if (mod(nm, 4) >= 2) qn_r1 = -qn_r end
    qnz_s = Vector{ComplexF64}(undef, 0)
    if qn_p != 0 push!(qnz_s, qn_p) end
    if qn_z != 0 push!(qnz_s, qn_z) end
    if qn_r != 0 push!(qnz_s, qn_r1) end
    # Generate the basis and print the dimension
    return Basis(cfs, qnz_s ; GetIsingQnz(nm ; qn_p, qn_z, qn_r)...)
end

"""
    function GetIsingIntTerms(nm :: Int64 ; ps_pot :: Vector) :: Vector{Term}

Returns the terms for the ising interaction 

```math
∑_{m_1m_2m_3m_4}2U_{m_1m_2m_3m_4}c^{†}_{m_1\\uparrow}c^{†}_{m_2\\downarrow}c_{m_3\\downarrow}c_{m_4\\uparrow}
```

from the pseudopotentials. 

# Arguments 

- `nm :: Int64` is the number of orbitals ``2s+1``.
- `ps_pot :: Vector{Number}` is the pseudopotential of Ising interaction.
"""
function GetIsingIntTerms(nm :: Int64 ; ps_pot :: Vector = [1.])
    return GetDenIntTerms(nm, 2 ; ps_pot = 2 .* ps_pot, mat_a = diagm([1, 0]), mat_b = diagm([0, 1]))
end

"""
    function GetXPolTerms(nm :: Int64)

Returns the terms for the density operator ``n^x_{l=0,m=0}``

# Arguments 

- `nm :: Int64` is the number of orbitals.
"""
function GetXPolTerms(nm :: Int64)
    return GetPolTerms(nm, 2 ; mat = [ 0 1 ; 1 0 ])
end

"""
    function GetZPolTerms(nm :: Int64)

Returns the terms for the density operator ``n^z_{l=0,m=0}``

# Arguments 

- `nm :: Int64` is the number of orbitals.
"""
function GetZPolTerms(nm :: Int64)
    return GetPolTerms(nm, 2 ; mat = [ 1 0 ; 0 -1 ])
end