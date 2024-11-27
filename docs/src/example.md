# FuzzifiED explained in an example

In this example, we will illustrate how to use `FuzzifiED` to calculate the spectrum of Ising model on fuzzy sphere and how to calculate the OPE coefficient ``\lambda_{\sigma\sigma\epsilon}`` by measuring the expectation value of the density operator ``n^z``. 

The examples can be found in the directory [`examples`](https://github.com/mankai-chow/FuzzifiED.jl/tree/main/examples). Two versions of this example is provided. The first uses the built-in functions for quantum numbers and operators to calculate the observables and is stored in [`tutorial_ising.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/tutorial_ising.jl). The second does not use the built-in example and is stored in [`tutorial_ising_primitive.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/tutorial_ising_primitive.jl). 

In addition, an example of how `FuzzifiED` can facilitate DMRG calculation is given. The example contains two versions. The first uses `MPO` and `dmrg` functions of the ITensors package and is stored in [`tutorial_ising_dmrg.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/tutorial_ising_dmrg.jl). The second uses the [`EasySweep`](@ref) function in the package which further wraps the `dmrg` function to facilitate the management of sweeps and is stored in [`tutorial_ising_dmrg_easysweep.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/tutorial_ising_dmrg_easysweep.jl). 

We also append in the end [a list of given examples](#List-of-examples) at the end of the page.

## Exact diagonalisation with FuzzifiED

### Implement the diagonal quantum numbers and generate the configurations

`FuzzifiED` supports ``\mathrm{U}(1)`` diagonal quantum numbers (QNDiag) 
```math
Q_i=\sum_{o=1}^{N_o}q_{io}n_o
```
or ``\mathbb{Z}_n`` diagonal quantum numbers with period ``P_i``
```math
Q_i=\sum_{o=1}^{N_o}q_{io}n_o\ \mathrm{mod}\ P_i
```
where ``i=1,\dots,N_U`` is the index of diagonal quantum numbers, ``o`` is the index of orbital, ``n_o=c^\dagger_oc_o``, and ``q_o`` is a set of coefficients that must be integer-valued. 

There are two diagonal quantum numbers in the Ising model, _viz._ the particle number and the angular momentum
```math
\begin{aligned}
Q_1&=N_e,& q_{1,m\sigma}&=1\\
Q_2&=2L_z,&q_{2,m\sigma}&=2m
\end{aligned}
```
where the orbital index ``o`` contains both ``m`` and ``\sigma``. In the code, we store the orbitals with the same ``m`` together, _viz._ we store the spin-up orbitals in odd ``o=1,3,\dots,2N_m-1`` and the spin-down orbitals in even ``o=2,4,\dots,2N_m``. We look at the ``L_z=0`` half-filled sector. 

These two diagonal quantum numbers are built in with functions [`GetNeQNDiag`](@ref) and [`GetLz2QNDiag`](@ref) and the configurations can be generated with function [`Confs`](@ref 
Confs(no :: Int64, secd :: Vector{Int64}, qnd :: Vector{QNDiag} ; nor :: Int64 = div(no, 2), num_th :: Int64 = NumThreads, disp_std :: Bool = !SilentStd)). 

```julia
# Inputing the basic setups
nm = 12
qnd = [ 
    GetNeQNDiag(2 * nm), 
    GetLz2QNDiag(nm, 2) 
]
cfs = Confs(2 * nm, [nm, 0], qnd)
```

#### Primitive version 

Alternatively, you can put in the diagonal quantum numbers [`QNDiag`](@ref) by hand by specifying the symmetry charges ``q_o`` of the orbitals, and facultatively the name and the modulus
```julia
nf = 2
nm = 8
no = nf * nm
s = .5 * (nm - 1)
ne = div(no, 2)
# Record the QNDiag
qnd = [ 
    # Number of electrons Ne
    QNDiag(fill(1, no)), 
    # Twice angular momentum 2Lz
    QNDiag([ (o - 1) ÷ nf * 2 - (nm - 1) for o = 1 : no ])
]
cfs = Confs(no, [ne, 0], qnd)
```

### Implement the off-diagonal ``\mathbb{Z}_p`` symmetries and initialise the basis

`FuzzifiED` supports off diagonal discrete ``\mathbb{Z}_p`` symmetries in the form of 

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
These three symmetries are built in with functions [`GetParityQNOffd`](@ref), [`GetFlavPermQNOffd`](@ref) and [`GetRotyQNOffd`](@ref). Thus, if we want to look at the all-positive sector
```julia
qnf = [ 
    GetParityQNOffd(nm, 2, [2, 1], [-1, 1]), 
    GetFlavPermQNOffd(nm, 2, [2, 1]), 
    GetRotyQNOffd(nm, 2) 
]
bs = Basis(cfs, [1, 1, 1], qnf)
# The second argument gives the eigenvalues under the transformations, for Z_2 put 1,-1 ; for Z_n put exp(2im*pi*q/p)
```
If no discrete symmetry is needed, one can simply put instead `bs = Basis(conf)`

#### Primitive version 

Alternatively, you can put in the diagonal quantum numbers [`QNOffd`](@ref) by hand by specifying the permutations ``\pi_o`` of the orbitals, and facultatively the particle-hole transformation ``p_o``, the factor ``\alpha_o`` and the cycle. 
```julia
# Record a QNOffd by orbital permutation (and facultatives particle-hole, factor, cycle)
qnf = [ 
    # Parity (Particle-hole)
    QNOffd([ isodd(o) ? o + 1 : o - 1 for o = 1 : no], true, ComplexF64[ isodd(o) ? -1 : 1 for o = 1 : no]),
    # Flavour symmetry
    QNOffd([ isodd(o) ? o + 1 : o - 1 for o = 1 : no]),
    # Y-axis pi-rotation
    QNOffd([ isodd(o) ? no - o : no + 2 - o for o = 1 : no], ComplexF64(-1) .^ (collect(0 : nm * nf - 1) .÷ nf))
]
bs = Basis(cfs, [1, 1, 1], qnf) 
```
Please refer to the [documenation](@ref QNDiag) if you want to implement ``\mathbb{Z}_{n>2}`` symmetries. 

### Record the Hamiltonian terms

The operators in the form of the sum of product of ``c`` and ``c^\dagger``'s are supported 

```math
\Phi=\sum_{t=1}^{N_t}U_tc^{(p_{t1})}_{o_{t1}}c^{(p_{t2})}_{o_{t2}}\dots c^{(p_{tl})}_{o_{tl}}
```

where ``c^{(0)}=c`` and ``c^{(1)}=c^\dagger``. Here the operator string sum is recorded together with the basis of the initial state and the basis of the final state. 

The Hamiltonian for the fuzzy sphere Ising model
```math
H=\sum_{m_1m_2m_3m_4}U_{m_1m_2m_3m_4}\delta_{m_1+m_2,m_3+m_4}c^\dagger_{m_1\uparrow}c^\dagger_{m_2\downarrow}c_{m_3\downarrow}c_{m_4\uparrow}-h\sum_m(c^\dagger_{m\uparrow}c_{m\downarrow}+\mathrm{h.c.})
```
can be recorded as the sum of a density-density interaction and a polarisation term, which are built-in with the functions [`GetDenIntTerms`](@ref) and [`GetPolTerms`](@ref). We then use [`Operator`](@ref) function to record it together with its basis. 
```julia
const σ1 = [  1  0 ;  0  0 ]
const σ2 = [  0  0 ;  0  1 ]
const σx = [  0  1 ;  1  0 ]
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, 2 .* [4.75, 1.], σ1, σ2) - 
    3.16 * GetPolTerms(nm, 2, σx) )
hmt = Operator(bs, tms_hmt)
```

#### Primitive version

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
            for l in eachindex(ps_pot)
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
hmt = Operator(bs, tms_hmt)
```

### Generate the sparse matrix and diagonalise

After specifying the Hamiltonian, we then use [Operator](@ref Operator(bsd :: Basis, bsf :: Basis, terms :: Vector{Term} ; red_q :: Int64 = 0, sym_q :: Int64 = 0)) to record also the basis and use [`OpMat`](@ref) to generate a sparse matrix from the operator. To get the 10 lowest eigenstates and their energies
```julia
hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 10)
```
Matrices with real elements can be generated by specifying `OpMat{Float64}(hmt)` or setting `FuzzifiED.ElementType = Float64`. 

### Write the sparse matrix into HDF5 file 

It is sometimes needed to write the sparse matrix into file to avoid extra effort to generate it again in another calculation. 

```julia
using HDF5 
f = h5open("data.h5", "cw")
write(f, "hmt_mat", hmt_mat)
close(f)
```

To read from the file

```julia
f = h5open("data.h5", "r") 
hmt_mat = read(f, "hmt_mat", OpMat{ComplexF64})
close(f)
```

Apart from `OpMat`, the supported types for writing include `Confs`, `Basis`, `Vector{Term}`, `Operator`, `OpMat{ComplexF64}` and `OpMat{Float64}`. 

### Measuring the angular momentum

We can measure the inner product of a final state, an operator or its matrix and an initial state or the action of an operator or its matrix on a state by directly using the [`*`](@ref) operator. The total angular momentum ``L^2`` is built-in with function [`GetL2Terms`](@ref). The following code measures the angular momentum of each eigenstate and verify whether ``|T\rangle`` is an eigenstate of ``L^2`` by measuring 

```math
|L^2T\rangle=L^2|T\rangle,\quad\frac{|\langle T|L^2T\rangle|^2}{\langle T|T\rangle\langle L^2T|L^2T\rangle}
```

```julia
tms_l2 = GetL2Terms(nm, 2)
l2 = Operator(bs, tms_l2)
l2_mat = OpMat(l2)
l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]
@show l2_val

st_T = st[:, 3]
st_L2T = l2_mat * st[:, 3]
@show abs(st_L2T' * st_T) ^ 2 / ((st_T' * st_T) * (st_L2T' * st_L2T))
```

#### Primitive version

By definition
```math
\begin{aligned}
L^2&=L^+L^-+(L^z)^2-L^z\\
L^z&=\sum_{\sigma m}mc^\dagger_mc_m\\
L^\pm&=\sum_{\sigma m}\sqrt{(s\mp m)(s\pm m+1)}c^\dagger_{m\pm 1}c_m
\end{aligned}
```

The construction of the operator can be simplified by the [addition](@ref +(tms1 :: Vector{Term}, tms2 :: Vector{Term})), [multiplication](@ref *(tms1 :: Vector{Term}, tms2 :: Vector{Term})), [Hermitian conjugate](@ref adjoint(tms :: Vector{Term})) and [simplification](@ref SimplifyTerms) of terms. 

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
```

### Go through all the sectors

One can repeat the calculation for all the sectors and records the results
```julia
result = []
for P in [1, -1], Z in [1, -1], R in [1, -1]
    bs = Basis(cfs, [P, Z, R], qnf)
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 10)

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    for i in eachindex(enrg)
        push!(result, round.([enrg[i], l2_val[i], P, Z], digits = 6))
    end
end
```
We then sort the eigenstates, find the energy of ground state and stress tensor, and calibrate the scaling dimensions. 
```julia
sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_T = filter(st -> st[2] ≊ 6 && st[3] ≊ 1 && st[4] ≊ 1, result)[1][1]
result_dim = [ [ 3 * (st[1] - enrg_0) / (enrg_T - enrg_0) ; st] for st in result ]
for P in (1, -1), Z in (1, -1)
    display(permutedims(hcat(
        filter(st -> st[4] ≊ P && st[5] ≊ Z, result_dim)...
    )))
end
```

### Measuring the density operator

Similar process can be used to calculate the OPE coefficient by measuring the density operator, by definition 

```math
\lambda_{\sigma\sigma\epsilon}=\frac{\langle\sigma|n^z_{00}|\epsilon\rangle}{\langle\sigma|n^z_{00}|\mathbb{I}\rangle},\quad n^z_{00}=\frac{1}{N_m}\sum_{\sigma m}\sigma c^\dagger_{\sigma m}c_{\sigma m}
```

To do that, we need to first repeat the calculation in the ``\mathbb{Z}_2``-odd sector
```julia
bs1 = Basis(cfs, [1, -1, 1], qnf)
hmt = Operator(bs1, bs1, tms_hmt ; red_q = 1, sym_q = 1) 
hmt_mat = OpMat(hmt)
enrg1, st1 = GetEigensystem(hmt_mat, 10)
st_I = st[:, 1] 
st_e = st[:, 2] 
st_s = st1[:, 1]
```

The [`SphereObs`](@ref) type stores the information of a local observable on the sphere. In particular, the [electron](@ref Electron) and [density operators](@ref Density) are built-in. The addition and multiplication of observables are enabled. It can be evaluated at a certain point [GetPointValue](@ref) or angular component with [`GetComponent`](@ref). 
```julia
obs_nz = Density(nm, 2, [ 1 0 ; 0 -1 ])
tms_nz = SimplifyTerms(GetComponent(obs_nz, 0.0, 0.0))
nz = Operator(bs, bs1, tms_nz ; red_q = 1) 
@show abs((st_s' * nz * st_e) / (st_s' * nz * st_I))
```

#### Primitive version

Alternatively, one can write out the terms of density operator explicitly.
```julia
# Record the density operator n^z
tms_nz = [ Term(isodd(o) ? 1 / nm : -1 / nm, [1, o, 0, o]) for o = 1 : no]
# The nz operator sends a state in bs (+) to bs1 (-)
nz = Operator(bs, bs1, tms_nz ; red_q = 1)
# Measuring the finite size OPE
@show abs((st_s' * nz * st_e) / (st_s' * nz * st_I))
```

## DMRG calculations with FuzzifiED

In this example, we calculate the ground state MPS of the fuzzy sphere Ising model by DMRG in ITensors, and use the same objects to do ED calculation and compare the result. 

We first set-up the calculation by
```julia
using FuzzifiED
using ITensors
const σx = [  0  1 ;  1  0 ]
const σz = [  1  0 ;  0 -1 ]

nm = 12
nf = 2
no = nm * nf
```
The first object we need is the sites. FuzzifiED overloads the fermion type and supports direct generation of the sites from the diagonal quantum numbers by the function [`SitesFromQNDiag`](@ref). The Ising model is expressed in the basis of ``XX-Z`` so that the flavour symmetry is diagonal. 
```julia
sites = SitesFromQNDiag([
    GetNeQNDiag(nm * nf), 
    GetLz2QNDiag(nm, nf),
    GetZnfChargeQNDiag(nm, nf)
])
```
We then generate the terms of the Hamiltonian using the built-in functions, convert it to `OpSum` type in by the [`OpSumFromTerms`](@ref), and convert it to MPO using the `MPO` function in ITensors
```julia
ps_pot = [4.75, 1.] ./ 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, ps_pot) - 
    GetDenIntTerms(nm, 2, ps_pot, σx) - 
    3.16 * GetPolTerms(nm, nf, σz)
)
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
hmt = Operator(bs, tms_hmt)
hmt_mat = OpMat(hmt)
enrg, st = GetEigensystem(hmt_mat, 10)
@show enrg
```

### Use EasySweep to manage DMRG sweeps

EasySweep facilitates the management of DMRG process by automatically recording the intermediate results and recovering these results if a job is stopped and run again on HPC. It also manages the gradual increase of maximal bond dimensions and the determination of convergence by the criteria of energy. Before the calculation, we need to define a method to generate MPO from OpSum and Sites. We suggest using `MPO_new` from package `ITensorMPOConstruction`, which can be installed through 
```julia
using Pkg ; Pkg.add(url="https://github.com/ITensor/ITensorMPOConstruction.jl.git"); 
```
```julia
using ITensorMPOConstruction
function MyMPO(os, sites)
    operatorNames = [ "I", "C", "Cdag", "N" ]
    opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
    return MPO_new(os, sites ; basisOpCacheVec = opCacheVec)
end
```
We also need to specify a path where the results are stored. 
```julia
path = "nm_$(nm)/"
mkpath(path)
```
The Hamiltonian MPO and the sites can be either generated or read from file by the function [`GetMPOSites`](@ref). The input is the quantum numbers and the terms or OpSum ; the output is the MPO and the sites
```julia
ps_pot = [4.75, 1.] ./ 2
tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, 2, ps_pot) - 
    GetDenIntTerms(nm, 2, ps_pot, σx) - 
    3.16 * GetPolTerms(nm, 2, σz)
)
qnd = [ 
    GetNeQNDiag(no), 
    GetLz2QNDiag(nm, nf), 
    GetZnfChargeQNDiag(nm, nf) 
]
hmt, sites = GetMPOSites("hmt", tms_hmt, qnd ; path, mpo_method = MyMPO)
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
l2 = GetMPO("l2", tms_l2, sites ; path)
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

The examples of FuzzifiED can be found in the repository [`examples`](https://github.com/mankai-chow/FuzzifiED.jl/tree/main/examples). Apart from the tutorials that we have introduced above 

* [`tutorial_ising.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/tutorial_ising.jl) does the ED calculation of Ising model through the built-in models. 
* [`tutorial_ising_primitive.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/tutorial_ising_primitive.jl) does the ED calculation of Ising model through the primitive functions.
* [`tutorial_ising_dmrg.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/tutorial_ising_dmrg.jl) does the DMRG calculation of Ising model through the `dmrg` function in ITensors.
* [`tutorial_ising_dmrg_easysweep.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/tutorial_ising_dmrg_easysweep.jl) does the DMRG calculation of Ising model through the `EasySweep` function which wraps ITensors.

We offer a series of other examples that reproduces various achievements of fuzzy sphere. For a more detailed summary of the background, see the [Review of existing work](@ref). 

* [`ising_spectrum.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_spectrum.jl) calculates the spectrum of 3d Ising model on fuzzy sphere at ``N_m = 12``. For each ``(P,Z,R)`` sector, 20 states are calculated. This example reproduces Table I and Figure 4 in [Zhu 2022](@ref References).
* [`ising_phase_diagram.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_phase_diagram.jl) calculates the phase diagram of fuzzy sphere Ising modelby calculating the order parameter ``\langle M^2\rangle``. This example reproduces Figure 3 in [Zhu 2022](@ref References).
* [`ising_ope.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_ope.jl) calculates various OPE coefficients at ``N_m = 12`` by taking overlaps between CFT states and density operators and composite. This example reproduces Figure 2 and Table I in [Hu 2023Mar](@ref References).
* [`ising_correlator.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_correlator.jl) calculates the ``σσ`` two-point function on sphere and the ``σσσσ`` four-point function on sphere, 0 and ``∞``. This example reproduces Figures 1c and 2a in [Han 2023Jun](@ref References).
* [`ising_optimisation.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_optimisation.jl) defines a cost function as the square sum of the deviations of descendants and stress tensor to evaluate the conformal symmetry for Ising model and minimises this cost function to find the best parameter.
* [`ising_full_spectrum.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_full_spectrum.jl) calculates the full spectrum of 3d Ising model on fuzzy sphere at ``N_m = 10`` for sector ``(P,Z,R) = (1,1,1)``.
* [`ising_space_entangle.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_space_entangle.jl) calculates the entanglement entropy of the Ising ground state along the real space cut of ``θ = 0.500π`` and ``0.499π`` respectively, and use these two data to extract finite size ``F``-function without sustracting the IQHE contribution. This example reproduces Figures 3 in [Hu 2024](@ref References).
* [`ising_orbital_entangle.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_space_entangle.jl) calculates the entanglement entropy of the Ising ground state along the orbital space cut at ``m = 0``, and also the entanglement spectrum in the half-filled ``l_z = 0, 1`` and  both ``\mathbb{Z}_2`` sectors.
* [`ising_generator.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_generator.jl) examines the quality of conformal symmetry at ``N_m = 12`` by examining the matrix elements of conformal generators ``P^z + K^z`` and compare the states ``(P^z + K^z)|Φ⟩`` with the CFT expectations. This example reproduces Figure 7 in [Fardelli 2024](@ref References).
* [`defect_spectrum.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/defect_spectrum.jl) calculates the spectrum of magnetic line defect in 3d Ising model in ``l_z = 0, P = ±1`` and ``l_z = 1`` sectors, calibrated by bulk ``T``. This example reproduces Table I in [Hu 2023Aug](@ref References).
* [`defect_correlator.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/defect_correlator.jl) calculates the 1-pt function ``σ`` and 2-pt function ``σ\hat{ϕ}`` of magnetic line defect in 3d Ising model. The normalisation of the correlators require extra bulk data. This example reproduces Figure 4 in [Hu 2023Aug](@ref References).
* [`defect_changing.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/defect_changing.jl) calculates the spectrum of the defect creation and changing operators of the magnetic line defect in 3d Ising model. This example reproduces Table 2 and Figure 5 in [Zhou 2024Jan](@ref References).
* [`defect_overlap.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/defect_overlap.jl) calculates the ``g``-function of magnetic line defect in 3d Ising model using the ovelaps between the bulk, defect ground state and the lowest defect-creation state. This example reproduces Figure 6 in [Zhou 2024Jan](@ref References).
* [`cusp_dim.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/cusp_dim.jl) calculates the scaling dimension of the cusp of the magnetic line defect in 3d Ising model as a function of the angle ``θ``. This example reproduces Table 2, upper panel in [Cuomo 2024](@ref References).
* [`surface_ordinary_spectrum.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/surface_ordinary_spectrum.jl) calculates the spectrum of ordinary surface CFT in 3d Ising model calibrated by surface displacement operator ``D`` in the orbital boundary scheme. This example reproduces Figures 3 and 4 in [Zhou 2024Jul](@ref References).
* [`surface_normal_spectrum.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/surface_normal_spectrum.jl) calculates the spectrum of normal surface CFT in 3d Ising model calibrated by surface displacement operator ``D`` in the orbital boundary scheme. This example reproduces Figure 5 in [Zhou 2024Jul](@ref References).
* [`o3_wf_spectrum.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/o3_wf_spectrum.jl) calculates the spectrum of ``\mathrm{O}(3)`` Wilson-Fisher CFT using the bilayer Heisenberg model. This example reproduces Table I and Figure 2 in [Hu 2023Dec](@ref References).
* [`so5_spectrum.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/so5_spectrum.jl) calculates the spectrum of SO(5) DQCP on fuzzy sphere. This example reproduces Table II in [Zhou 2023](@ref References).
* [`sp3_spectrum.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/sp3_spectrum.jl) calculates the spectrum of Sp(3) CFT on fuzzy sphere. This example reproduces Table I in [Zhou 2024Oct](@ref References).
* [`ising_spectrum_krylov.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_spectrum_krylov.jl) calculates the spectrum of 3d Ising model on fuzzy sphere by calling the eigsolve function in KrylovKit.jl instead of Arpack.
* [`ising_spectrum_cuda.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/ising_spectrum_cuda.jl) calculates the spectrum of 3d Ising model on fuzzy sphere for one sector by performing the sparse matrix multiplication on CUDA.