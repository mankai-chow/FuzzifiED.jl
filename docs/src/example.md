# FuzzifiED explained in an example

In this example, we will illustrate how to use `FuzzifiED` to calculate the spectrum of Ising model on fuzzy sphere and how to calculate the OPE coefficient ``\lambda_{\sigma\sigma\epsilon}`` by measuring the expectation value of the density operator ``n^z``. 

The examples can be found in the directory [`examples`](https://github.com/mankai-chow/FuzzifiED.jl/tree/main/examples). Three versions of this example is provided. The first does not use the built-in example and is stored in [`ising_primitive.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_primitive.jl) ; the second uses the interfaces with ITensors and is stored in [`example_ising_itensors.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_itensors.jl) ; the third uese the built-in example for the Ising model and is stored in [`example_ising.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising.jl). The explanations following applied mainly to the primitive version. 

In addition, an example of how `FuzzifiED` can facilitate DMRG calculation is given. Two versions of the DMRG example is provided. The first uses `MPO` and `dmrg` functions of the ITensors package and is stored in [`ising_dmrg.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_dmrg.jl). The second uses the [`EasySweep`](@ref) function in the package which further wraps the `dmrg` function to facilitate the management of sweeps and is stored in [`ising_dmrg_easysweep.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_dmrg_easysweep.jl). 

We also append in the end a list of 

## Exact diagonalisation with FuzzifiED

### Implement the diagonal quantum numbers and generate the configurations

`FuzzifiED` supports diagonal quantum numbers (QNU) in the form of 

```math
Q_i=\sum_{o=1}^{N_o}q_{io}n_o \quad \mathrm{$U(1)$ symmetry} \quad\mathrm{or}\quad Q_i=\sum_{o=1}^{N_o}q_{io}n_o\ \mathrm{mod}\ P_i \quad \mathrm{Discrete $Z_p$ symmetry}
```

where ``i=1,\dots,N_U`` is the index of diagonal quantum numbers, ``o`` is the index of orbital, ``n_o=c^\dagger_oc_o``, and ``q_o`` is a set of coefficients that must be non negative integer valued. (A list of ``q_o`` with both positive and negative entries can be adapted by shifting every elements by a same value)

The function used to implement the diagonal quantum numbers and generate all the configurations (_i.e._, direct product states) is [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2))). There are two diagonal quantum numbers in the Ising model, _viz._ the particle number and the angular momentum

```math
\begin{aligned}
Q_1&=N_e,& q_{1,m\sigma}&=1\\
Q_2&=L_z+N_es,&q_{2,m\sigma}&=m+s
\end{aligned}
```

where the orbital index ``o`` contains both ``m`` and ``\sigma``. In the code, we store the orbitals with the same $m$ together, _viz._ we store the spin-up orbitals in odd ``o=1,3,\dots,2N_m-1`` and the spin-down orbitals in even ``o=2,4,\dots,2N_m``. Thus, if we want to look at the ``L_z=0`` half-filled sector

```julia
# Inputing the basic setups
nf = 2
nm = 8
no = nf * nm
s = .5 * (nm - 1)
ne = div(no, 2)
# Initialise the arrays
qnu_s = Vector{Int64}(undef, 0)
qnu_o = []
# Record the number of electrons
push!(qnu_o, fill(1, no)) # qnu_o[1] = [1,1,...,1]
push!(qnu_s, ne) 
# Record the angular momentum
push!(qnu_o, [ div(o - 1, nf) for o = 1 : no ]) # qnu_o[2] = [0,0,1,1,...,7,7] to qnu_o
push!(qnu_s, ne * s) 
# Generate the configurations and print the number
@time "Initialise configurations" cfs = Confs(no, qnu_s, qnu_o)
@show cfs.ncf
```

#### ITensor support

The quantum numbers can also be imported from the `Sites` objects in `ITensors`. This can be done using the [`ConfsFromSites`](@ref) function.

```julia
# Overload the ITensors type "Fermion"
function ITensors.space( :: SiteType"Fermion" ; m1 :: Int = 0)
    return [
        QN(("Nf", 0, -1), ("Lz",  0)) => 1
        QN(("Nf", 1, -1), ("Lz", m1)) => 1
    ]
end
# Initialise the sites
sites = [ siteind("Fermion", m1 = div(o - 1, nf)) for o :: Int = 1 : no]
qn_s = QN(("Nf", ne), ("Lz", Int(ne * s)))
@time "Initialise configurations" cfs = ConfsFromSites(sites, qn_s)
# Alternatively, one can initialise the configuration quantum number
# cf_ref = [o <= ne ? 1 : 0 for o = 1 : no]
# @time "Initialise configurations" cfs = ConfsFromSites(sites, cf_ref)
@show cfs.ncf
```

#### Built-in model 
Using The built-in Ising model, the process above can be done in one line with the method [`GetLzConfs`](@ref).
```julia
cfs = GetLzConfs(nm, 2)
```

### Implement the discrete symmetries and initialise the basis

`FuzzifiED` supports discrete ``\mathbb{Z}_n`` symmetries in the form of 

```math
\mathscr{Z}:\ c_o\to \alpha_o^* c^{(p_o)}_{\pi_o},\quad c_o^\dagger\to \alpha_o c^{(1-p_o)}_{\pi_o}
```

where we use a notation ``c^{(1)}=c^\dagger`` and ``c^{0}=c`` for convenience, where ``\pi_o`` is a permutation of ``1,\dots N_o``, ``\alpha_o`` is a coefficient, and ``p_o`` specified whether or not particle-hole transformation is performed for the orbital. Note that one must guarentee that all these transformations commute with each other and also commute with the diagonal quantum numbers. 

After implementing these symmetries, a state in the new basis should look like 

```math
|I\rangle=\lambda_{i_{I1}}|i_{I1}\rangle+\lambda_{i_{I2}}|i_{I2}\rangle+\cdots+\lambda_{i_{Im_I}}|i_{Im_I}\rangle
```

where the ``|i\rangle``'s are configurations in the `Confs` generated in the last section. The ``|I\rangle`` is a linear combination, and can be regarded as a grouping of ``m_I`` configurations.

The function used to implement the discrete symmetries is [`Basis`](@ref). There are three ``\mathbb{Z}_2`` transformations in the Ising model, _viz._ the particle-hole transformation ``\mathscr{P}``, the ``\pi``-rotation along the ``y``-axis ``\mathscr{R}_y`` and the flavour (Ising) symmetry ``\mathscr{Z}``

```math
\begin{aligned}
    \mathscr{P}:c^\dagger_{\sigma m}&\to\sigma c_{-\sigma,m}\\
    \mathscr{Z}:c^\dagger_{\sigma m}&\to c^\dagger_{-\sigma,m}\\
    \mathscr{R}_y:c^\dagger_{\sigma m}&\to c^\dagger_{\sigma,-m}\\
\end{aligned}
```

Thus, if we want to look at the all-positive sector
```julia
cyc = [ 2, 2, 2 ] # Input three Z_2 symmetries 
qnz_s = ComplexF64[ 1, 1, 1 ] # Quantum numbers are all positive 
# Initialise the vectors
perm_o = []
ph_o = []
fac_o = []
# Record the parity
push!(perm_o, [ isodd(o) ? o + 1 : o - 1 for o = 1 : no]) # perm_o[1] = [2,1,4,3,...,16,15]
push!(ph_o, fill(1, no)) # ph_o[1] = [1,1,...,1] meaning PH
push!(fac_o, [ isodd(o) ? -1 : 1 for o = 1 : no]) # fac_o[1] = [1,-1,1,-1,...,1,-1]
# Record the flavour symmetry
push!(perm_o, [ isodd(o) ? o + 1 : o - 1 for o = 1 : no]) # perm_o[3] = [2,1,4,3,...,16,15]
push!(ph_o, fill(0, no)) # ph_o[3] = [0,0,...,0] meaning no PH
push!(fac_o, fill(ComplexF64(1), no)) # fac_o[3] = [1,1,...,1]
# Record the pi-rotation
push!(perm_o, [ isodd(o) ? no - o : no + 2 - o for o = 1 : no]) # perm_o[2] = [15,16,13,14,...,1,2]
push!(ph_o, fill(0, no)) # ph_o[2] = [0,0,...,0] meaning no PH
push!(fac_o, fill(ComplexF64(1), no)) # fac_o[2] = [1,1,...,1]
# Generate the basis and print the dimension
@time "Initialise basis" bs = Basis(cfs, qnz_s, cyc, perm_o, ph_o, fac_o)
@show bs.dim 
```

Note that if no discrete symmetry is needed, one can simply put instead `bs = Basis(conf)`

#### Built-in model 
Using The built-in Ising model, the process above can be done in one line with the method [`GetIsingBasis`](@ref).
```julia
bs = GetIsingBasis(cfs ; qn_p = 1, qn_r = 1, qn_z = 1)
```

### Record the Hamiltonian operator
The operator here refers to the sum of product of ``c`` and ``c^\dagger``'s in the form 

```math
\Phi=\sum_{t=1}^{N_t}U_tc^{(p_{t1})}_{o_{t1}}c^{(p_{t2})}_{o_{t2}}\dots c^{(p_{tl})}_{o_{tl}}
```

where ``c^{(0)}=c`` and ``c^{(1)}=c^\dagger``. Here the operator string sum is recorded together with the basis of the initial state and the basis of the final state. 

This can be generated by the [`Operator`](@ref) function. The Hamiltonian for the fuzzy sphere Ising model
```math
H=\sum_{m_1m_2m_3m_4}U_{m_1m_2m_3m_4}\delta_{m_1+m_2,m_3+m_4}c^\dagger_{m_1\uparrow}c^\dagger_{m_2\downarrow}c_{m_3\downarrow}c_{m_4\uparrow}-h\sum_m(c^\dagger_{m\uparrow}c_{m\downarrow}+\mathrm{h.c.})
```
can be recorded with 

```julia
using WignerSymbols
# Input the parameters of the Hamiltonian
ps_pot = [ 4.75, 1. ] * 2.
h = 3.16
tms_hmt = Vector{Term}(undef, 0)
# Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m3
for m1 = 0 : nm - 1
    f1 = 0
    o1 = m1 * nf + f1 + 1
    m1r = m1 - s
    for m2 = 0 : nm - 1
        f2 = 1
        o2 = m2 * nf + f2 + 1
        m2r = m2 - s
        for m3 = 0 : nm - 1
            f3 = 1
            o3 = m3 * nf + f3 + 1
            m3r = m3 - s
            m4 = m1 + m2 - m3 
            if (m4 < 0 || m4 >= nm) continue end
            f4 = 0
            o4 = m4 * nf + f4 + 1
            m4r = m4 - s
            # Calculate the matrix element val from pseudopotentials
            val = ComplexF64(0)
            for l in 1 : length(ps_pot)
                if (abs(m1r + m2r) > nm - l || abs(m3r + m4r) > nm - l) break end 
                val += ps_pot[l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, m1r, m2r, -m1r - m2r) * wigner3j(s, s, nm - l, m4r, m3r, -m3r - m4r)
            end 
            # Record the interaction term val, "Cdag", o1, "Cdag", o2, "C", o3, "C", o4
            push!(tms_hmt, Term(val, [1, o1, 1, o2, 0, o3, 0, o4]))
        end
    end
    o1x = o1 + 1
    # Record the transverse field term
    push!(tms_hmt, Term(-h, [1, o1, 0, o1x]))
    push!(tms_hmt, Term(-h, [1, o1x, 0, o1]))
end
# Generate the Hamiltonian operator
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
```

#### ITensor support

Alternatively, one can generate the operator using an `OpSum` object instead of `cstr_vec` and `fac` using the function [`TermsFromOpSum`](@ref).

For the Hamiltonian of Ising model,
```julia

using WignerSymbols
# Input the parameters of the Hamiltonian
ps_pot = [ 4.75, 1. ] * 2.
h = 3.16
global ops_hmt = OpSum()
# Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m3
for m1 = 0 : nm - 1
    f1 = 0
    o1 = m1 * nf + f1 + 1
    m1r = m1 - s
    for m2 = 0 : nm - 1
        f2 = 1
        o2 = m2 * nf + f2 + 1
        m2r = m2 - s
        for m3 = 0 : nm - 1
            f3 = 1
            o3 = m3 * nf + f3 + 1
            m3r = m3 - s
            m4 = m1 + m2 - m3 
            if (m4 < 0 || m4 >= nm) continue end
            f4 = 0
            o4 = m4 * nf + f4 + 1
            m4r = m4 - s
            # Calculate the matrix element val from pseudopotentials
            val = ComplexF64(0)
            for l in 1 : length(ps_pot)
                if (abs(m1r + m2r) > nm - l || abs(m3r + m4r) > nm - l) break end 
                val += ps_pot[l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, m1r, m2r, -m1r - m2r) * wigner3j(s, s, nm - l, m4r, m3r, -m3r - m4r)
            end 
            # Record the interaction term
            global ops_hmt += val, "Cdag", o1, "Cdag", o2, "C", o3, "C", o4
        end
    end
    o1x = o1 + 1
    # Record the transverse field term
    global ops_hmt += -h, "Cdag", o1, "C", o1x
    global ops_hmt += -h, "Cdag", o1x, "C", o1
end
# Generate the Hamiltonian operator
tms_hmt = TermsFromOpSum(ops_hmt)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
```

#### Built-in model 
Using The built-in Ising model, the process above can be done in one line with the method [`GetIsingIntTerms`](@ref) and [`GetXPolTerms`](@ref)
```julia
tms_hmt = GetIsingIntTerms(nm ; ps_pot = [4.75, 1.]) - 3.16 * GetXPolTerms(nm)
```

### Generate the sparse matrix and diagonalise

After specifying the Hamiltonian, we then use the [`OpMat`](@ref) to generate a sparse matrix from the operator. To get the 10 lowest eigenstates and their energies
```julia
@time "Initialise the Hamiltonian matrix" hmt_mat = OpMat(hmt)
@show hmt_mat.nel
@time "Diagonalise Hamiltonian" enrg, st = GetEigensystem(hmt_mat, 10)
@show real(enrg)
```

We also note that matrices with real elements can be generated with the option `type = Float64` in the `OpMat` function. 

### Measuring the angular momentum

We can measure the inner product of a final state, an operator or its matrix and an initial state or the action of an operator or its matrix on a state by directly using the [`*`](@ref) operator. To measure the total angular momentum ``L^2`` by definition

```math
\begin{aligned}
L^2&=L^+L^-+(L^z)^2-L^z\\
L^z&=\sum_{\sigma m}mc^\dagger_mc_m\\
L^\pm&=\sum_{\sigma m}\sqrt{(s\mp m)(s\pm m+1)}c^\dagger_{m\pm 1}c_m
\end{aligned}
```

The construction of the operator can be simplified by the [addition](@ref +(tms1 :: Vector{Term}, tms2 :: Vector{Term})), [multiplication](@ref *(tms1 :: Vector{Term}, tms2 :: Vector{Term})), [Hermitian conjugate](@ref adjoint(tms :: Vector{Term})) and [simplification](@ref SimplifyTerms) of terms. The following code measures the angular momentum of each eigenstate and verify whether ``|T\rangle`` is an eigenstate of ``L^2`` by measuring 

```math
|L^2T\rangle=L^2|T\rangle,\quad\frac{|\langle T|L^2T\rangle|^2}{\langle T|T\rangle\langle L^2T|L^2T\rangle}
```

```julia
tms_lz = 
    [ begin m = div(o - 1, nf)
        Term(m - s, [1, o, 0, o])
    end for o = 1 : no ]
tms_lp = 
    [ begin m = div(o - 1, nf)
        Term(sqrt(m * (nm - m)), [1, o, 0, o - nf])
    end for o = nf + 1 : no ]
tms_lm = tms_lp' 
tms_l2 = SimplifyTerms(tms_lz * tms_lz - tms_lz + tms_lp * tms_lm)
# Initialise the L2 operator
l2 = Operator(bs, bs, tms_l2 ; red_q = 1, sym_q = 1)
@time "Initialise L2" l2_mat = OpMat(l2)
# Calculate the inner product for each eigenstate
@time "Measure L2" l2_val = [ st[:, i]' * l2_mat * st[:, i] for i = 1 : length(enrg)]
@show real(l2_val)
# Verify whether T is an eigenstate of L^2
st_T = st[:, 3]
st_L2T = l2_mat * st[:, 3]
@show abs(st_L2T' * st_T) ^ 2 / ((st_T' * st_T) * (st_L2T' * st_L2T))
```

#### Built-in model 
Using The built-in Ising model, the generation of the terms in total angular momentum can be done in one line with the method [`GetL2Terms`](@ref)
```julia
tms_l2 = GetL2Terms(nm, 2)
```
This method also applies to other models on fuzzy sphere. 

### Measuring the density operator

Similar process can be used to calculate the OPE coefficient by measuring the density operator, by definition 

```math
\lambda_{\sigma\sigma\epsilon}=\frac{\langle\sigma|n^z_{00}|\epsilon\rangle}{\langle\sigma|n^z_{00}|\mathbb{I}\rangle},\quad n^z_{00}=\frac{1}{N_m}\sum_{\sigma m}\sigma c^\dagger_{\sigma m}c_{\sigma m}
```

To do that, we need to first repeat the calculation in the ``\mathbb{Z}_2``-odd sector
```julia
# Repeat the calculation for the Z_2 odd sector (with subscript 1)
qnz_s1 = ComplexF64[ 1, -1, 1 ] # Change only the discrete quantum numbers and generate the basis
@time "Initialise Basis Z" bs1 = Basis(cfs, qnz_s1, cyc, perm_o, ph_o, fac_o) 
@show bs1.dim 
hmt = Operator(bs1, bs1, tms_hmt ; red_q = 1, sym_q = 1) # Generate and diagonalise Hamiltonian in the new basis
@time "Initialise Hamiltonian" hmt_mat = OpMat(hmt)
@show hmt_mat.nel
@time "Diagonalise Hamiltonian" enrg1, st1 = GetEigensystem(hmt_mat, 10)
@show real(enrg1)
# Record the identity, sigma and epsilon states 
st_I = st[:, 1] # ground state
st_e = st[:, 2] # epsilon state
st_s = st1[:, 1]
```
and then generate and measure the density operator
```julia
# Record the density operator n^z
tms_nz = [ Term(isodd(o) ? 1 / nm : -1 / nm, [1, o, 0, o]) for o = 1 : no]
# The nz operator sends a state in bs (+) to bs1 (-)
nz = Operator(bs, bs1, tms_nz ; red_q = 1)
# Measuring the finite size OPE
@show abs((st_s' * nz * st_e) / (st_s' * nz * st_I))
```

#### Built-in model 
Using The built-in Ising model, the generation of the terms in density operator can be done in one line with the method [`GetZPolTerms`](@ref)
```julia
tms_nz = GetZPolTerms(nm) 
```
This method also applies to other models on fuzzy sphere. 

## DMRG calculations with FuzzifiED

In this example, we calculate the ground state MPS of the fuzzy sphere Ising model by DMRG in ITensors, and use the same objects to do ED calculation and compare the result. 

We first set-up the calculation by
```julia
using FuzzifiED
using ITensors
nm = 12
nf = 2
no = nm * nf
```
The first object we need is the sites. FuzzifiED overloads the fermion type and supports direct generation of the sites from the diagonal quantum numbers by the function [`SitesFromQnu`](@ref). The diagonal quantum numbers of the Ising model is built in and can be generated by the function [`GetLzZnQnu`](@ref). The Ising model is expressed in the basis of ``XX-Z`` so that the flavour symmetry is diagonal. 
```julia
sites = SitesFromQnu(; GetLzZnQnu(nm, 2)...)
```
We then generate the terms of the Hamiltonian using the built-in functions, convert it to `OpSum` type in by the [`OpSumFromTerms`](@ref), and convert it to MPO using the `MPO` function in ITensors
```julia
sigma_x = [ 0 1 ; 1 0 ]
tms_hmt = SimplifyTerms(GetDenIntTerms(nm ; ps_pot = [4.75, 1.]) - GetDenIntTerms(nm ; ps_pot = [4.75, 1.], mat_a = sigma_x) - 3.16 * GetZPolTerms(nm))
@time mpo_hmt = MPO(OpSumFromTerms(tms_hmt), sites)
```
We then use the all-up state as the initial state. In FuzzifiED, the occupied and empty sites are expressed by 0 and 1, while they are expressed by `"0"` and `"1"` in ITensors, so a conversion to string is needed. 
```julia
cf0 = [ isodd(o) ? 1 : 0 for o = 1 : no ]
st0 = MPS(sites, string.(cf0))
```
After that, the Hamiltonian MPO and the initial state MPS can be used for input for DMRG calculation. 
```julia
Eg, stg = dmrg(mpo_hmt, st0 ; nsweeps = 10, maxdim = [10,20,50,100,200,500], noise = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], cutoff = [1E-8])
@show Eg
```
We then convert these objects for the ED calculate. The configurations can be generated from Sites and a reference configuration by the function [`ConfsFromSites`](@ref). 
```julia
cfs = ConfsFromSites(sites, cf0)
bs = Basis(cfs)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1)
hmt_mat = OpMat(hmt ; type = Float64)
enrg, st = GetEigensystem(hmt_mat, 10)
@show enrg
```

### Use EasySweep to manage DMRG sweeps

EasySweep facilitates the management of DMRG process by automatically recording the intermediate results and recovering these results if a job is stopped and run again on HPC. It also manages the gradual increase of maximal bond dimensions and the determination of convergence by the criteria of energy. Before the calculation, we need to define a method to generate MPO from OpSum and Sites. We suggest using `MPO_new` from package `ITensorMPOConstruction`, which can be installed through 
```julia
julia> Pkg.add(url="https://github.com/ITensor/ITensorMPOConstruction.jl.git"); 
```
```julia
using ITensorMPOConstruction
function MyMPO(os, sites)
    operatorNames = [ "I", "C", "Cdag", "N" ]
    opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
    mpo = MPO_new(os, sites ; basisOpCacheVec = opCacheVec)
end
```
We also need to specify a path where the results are stored. 
```julia
path = "nm_$(nm)/"
mkpath(path)
```
The Hamiltonian MPO and the sites can be either generated or read from file by the function [`GetMpoSites`](@ref). The input is the quantum numbers and the terms or OpSum ; the output is the MPO and the sites
```julia
sigma_x = [ 0 1 ; 1 0 ]
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2 ; ps_pot = [4.75, 1.]) - 
    GetDenIntTerms(nm, 2 ; ps_pot = [4.75, 1.], mat_a = sigma_x) - 
    3.16 * GetZPolTerms(nm)
)
hmt, sites = GetMpoSites("hmt", tms_hmt ; path, GetLzZnQnu(nm, 2)..., mpo_method = MyMPO)
```
After generating the initial state MPS, the DMRG calculation of the states ``\mathbb{I}`` and ``\epsilon`` can be done by the [`EasySweep`](@ref) function.  
```julia
cf0 = [ isodd(o) ? 1 : 0 for o = 1 : no ]
st0 = MPS(sites, string.(cf0))
Eg, stg = EasySweep("g", hmt, st0 ; path)
Ee, ste = EasySweep("e", hmt, st0 ; path, proj = ["g"])
```
The total angular momentum can be measured by generating the MPO of ``L^2`` and measure the inner product 
```julia
tms_l2 = GetL2Terms(nm, 2)
l2 = GetMpo("l2", tms_l2, sites ; path)
@show inner(stg', l2, stg)
```
The ``\mathbb{Z}_2``-odd ``\sigma`` state can be calculated similarly.
```julia
cf1 = cf0
cf1[1] = 0
cf1[2] = 1
st1 = MPS(sites, string.(cf1))
Es, sts = EasySweep("s", hmt, st1 ; path)
```

## List of examples

The following examples of FuzzifiED can be found in the repository [`examples`](https://github.com/mankai-chow/FuzzifiED.jl/tree/main/examples).

* [`example_ising.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising.jl) does the ED calculation of Ising model through the built-in models. 
* [`example_ising_primitive.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_primitive.jl) does the ED calculation of Ising model through the primitive functions.
* [`example_ising_itensors.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_itensors.jl) does the ED calculation of Ising model by the Sites and OpSum objects in ITensors.
* [`example_ising_dmrg.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_dmrg.jl) does the DMRG calculation of Ising model through the `dmrg` function in ITensors.
* [`example_ising_dmrg_easysweep.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_dmrg_easysweep.jl) does the DMRG calculation of Ising model through the `EasySweep` function which wraps ITensors.
* [`example_ising_def.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_def.jl) does the ED calculation of Ising model with magnetic line defect or defect creation or changing operators.
* [`example_ising_def_dmrg.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_def_dmrg.jl) does the DMRG calculation of Ising model with magnetic line defect or defect changing operators. 
* [`example_sp2.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_sp2.jl) does the ED calculation of the ``\mathrm{SO}(5)`` DQCP model.
* [`example_potts.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_potts.jl) does the ED calculation of the three-states Potts model.
* [`example_potts_dmrg.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_potts_dmrg.jl) does the DMRG calculation of the three-states Potts model in the X basis.
