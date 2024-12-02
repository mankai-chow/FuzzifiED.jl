"""
    SQNDiag

The mutable type `SQNDiag` records the information of a diagonal ``\\mathrm{U}(1)`` or ``ℤ_p`` quantum number in the form of a symmetry charge
```math
Q=∑_{o=1}^{N_{of}}q_{f,o}n_{f,o}+∑_{o=1}^{N_{ob}}q_{b,o}n_{b,o}
```
or
```math
Q=∑_{o=1}^{N_{of}}q_{f,o}n_{f,o}+∑_{o=1}^{N_{ob}}q_{b,o}n_{b,o}\\ \\mathrm{mod}\\ p
```
where ``i=1,…,N_U`` is the index of quantum number, ``o`` is the index of site, ``N_{of}`` and ``N_{ob}`` are the number of fermionic and bosonic sites, ``n_{f,o}=f^†_of_o``, ``n_{b,o}=b^†_ob_o``, and ``q_{f,o},q_{b,o}`` are a set of symmetry charges that must be integer valued.

# Fields 

- `name :: String` is the name of the diagonal quantum number 
- `chargef :: Vector{Int64}` is the symmetry charge ``q_{f,o}`` of each site
- `chargeb :: Vector{Int64}` is the symmetry charge ``q_{b,o}`` of each site
- `modul :: Vector{Int64}` is the modulus ``p``, set to 1 for ``\\mathrm{U}(1)`` SQNDiags. 

# Initialisation 

It can be initialised by the following method
```julia
SQNDiag([name :: String, ]chargef :: Vector{Int64}, chargeb :: Vector{Int64}[, modul :: Int64]) :: SQNDiag
```
The arguments `name` and `modul` are facultative. By default `name` is set to `\"QN\"` and `modul` is set to 1. 
"""
mutable struct SQNDiag
    name :: String
    chargef :: Vector{Int64}
    chargeb :: Vector{Int64}
    modul :: Int64
    SQNDiag(name :: String, chargef :: Vector{Int64}, chargeb :: Vector{Int64}, modul :: Int64 = 1) = new(name, chargef, chargeb, modul)
    SQNDiag(chargef :: Vector{Int64}, chargeb :: Vector{Int64}, modul :: Int64 = 1) = new("QN", chargef, chargeb, modul)
end

"""
    SQNOffd

The mutable type `SQNOffd` records the information of an off-diagonal ``ℤ_p`` quantum number in the form of a discrete transformation
```math
𝒵:\\ f_o↦ α_{f,o}^* f^{(p_{f,o})}_{π_{f,o}},  f_o^†↦α_{f,o} c^{(1-p_{f,o})}_{π_{f,o}},  b_o^†↦α_{b,o} b^†_{π_{b,o}}
```
where we use a notation ``c^{(1)}=c^†`` and ``c^{0}=c`` for convenience, ``π_{f,o},π_{b,o}`` are permutations of ``1,…,N_{of}`` or ``N_{ob}``, ``α_{f,o},α_{b,o}`` are coefficients, and ``p_{f,o}`` specified whether or not particle-hole transformation is performed for the fermionic site. Note that one must guarentee that all these transformations commute with each other and also commute with the diagonal QNs. 

# Arguments 

- `permf :: Vector{Int64}` is a length-``N_{of}`` vector that records the fermion permutation ``π_{f,o}``.
- `permb :: Vector{Int64}` is a length-``N_{ob}`` vector that records the boson permutation ``π_{b,o}``.
- `phf :: Vector{Int64}` is a length-``N_{of}`` vector that records ``p_{f,o}`` to determine whether or not to perform a particle-hole transformation
- `facf :: Vector{ComplexF64}` is a length-``N_{of}`` vector that records the factor ``α_{f,o}`` in the transformation.
- `facb :: Vector{ComplexF64}` is a length-``N_{ob}`` vector that records the factor ``α_{b,o}`` in the transformation.
- `cyc :: Int64` is the cycle ``p``. 

# Initialisation 

It can be initialised by the following method
```julia
SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}[, phf :: Vector{Int64}][, facf :: Vector{ComplexF64}, facb :: Vector{ComplexF64}][, cyc :: Int64]) :: SQNOffd
SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}, phf_q :: Bool[, fac :: Vector{ComplexF64}, facb :: Vector{ComplexF64}]) :: SQNOffd
```
The arguments `phf`, `facf`, `facb` and `cyc` are facultative. By default `ph` is set all 0, `facf`, `facb` is set to all 1 and `cyc` is set to 2. If `phf_q` is a bool and true, then `ph` is set to all 1. 

"""
mutable struct SQNOffd
    permf :: Vector{Int64}
    permb :: Vector{Int64}
    phf :: Vector{Int64}
    facf :: Vector{ComplexF64}
    facb :: Vector{ComplexF64}
    cyc :: Int64
    SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}, phf :: Vector{Int64} = [0 for o in permf], facf :: Vector{ComplexF64} = [ComplexF64(1) for o in permf], facb :: Vector{ComplexF64} = [ComplexF64(1) for o in permb], cyc :: Int64 = 2) = new(permf, permb, phf, facf, facb, cyc)
    SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}, facf :: Vector{ComplexF64}, facb :: Vector{ComplexF64}, cyc :: Int64 = 2) = new(permf, permb, [0 for o in permf], facf, facb, cyc)
    SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}, cyc :: Int64) = new(permf, permb, [0 for o in permf], [ComplexF64(1) for o in permf], [ComplexF64(1) for o in permb], cyc)
    SQNOffd(permf :: Vector{Int64}, permb :: Vector{Int64}, phf_q :: Bool, facf :: Vector{ComplexF64} = [ComplexF64(1) for o in permf], facb :: Vector{ComplexF64} = [ComplexF64(1) for o in permb]) = new(permf, permb, [phf_q ? 1 : 0 for o in permf], facf, facb, 2)
end