# FuzzifiED explained in an example

In this example, we will illustrate how to use  `FuzzifiED` to calculate the spectrum of Ising model on fuzzy sphere and how to calculate the OPE coefficient ``\lambda_{\sigma\sigma\epsilon}`` by measuring the expectation value of the density operator ``n^z``. 

## Implement the conserved quantities and generate the configurations

`FuzzifiED` supports conserved quantities in the form of 

```math
Q_i=\sum_{o=1}^{N_o}q_{io}n_o
```

where ``i=1,\dots,N_U`` is the index of conserved quantities, ``o`` is the index of orbital, ``n_o=c^\dagger_oc_o``, and ``q_o`` is a set of coefficients that must be non negative integer valued. 

The function used to implement the conserved quantities and generate all the configurations (_i.e._, direct product states) is [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2))). There are two conserved quantities in the Ising model, _viz._ the particle number and the angular momentum

```math
\begin{aligned}
Q_1&=N_e,& q_{1,m\sigma}&=1\\
Q_2&=L_z+N_es,&q_{2,m\sigma}&=m+s
\end{aligned}
```

where the orbital index ``o`` contains both ``m`` and ``\sigma``. In the code, we store the spin-up orbitals in ``o=1,\dots,N_m`` and the spin-down orbitals in ``o=N_m+1,\dots,2N_m``. Thus, if we want to look at the ``L_z=0`` half-filled sector

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

### ITensor support

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

## Implement the discrete symmetries and initialise the basis

`FuzzifiED` supports discrete ``\mathbb{Z}_n`` symmetries in the form of 

```math
\mathscr{Z}:\ c_o\to \alpha_o^* c^{(p_o)}_{\pi_o},\quad c_o^\dagger\to \alpha_o c^{(1-p_o)}_{\pi_o}
```

where we use a notation ``c^{(1)}=c^\dagger`` and ``c^{0}=c`` for convenience, where ``\pi_o`` is a permutation of ``1,\dots N_o``, ``\alpha_o`` is a coefficient, and ``p_o`` specified whether or not particle-hole transformation is performed for the orbital. Note that one must guarentee that all these transformations commute with each other and also commute with the conserved quantities. 

After implementing these symmetries, a state in the new basis should look like 

```math
|I\rangle=\lambda_{i_{I1}}|i_{I1}\rangle+\lambda_{i_{I2}}|i_{I2}\rangle+\cdots+\lambda_{i_{Im_I}}|i_{Im_I}\rangle
```

where the ``|i\rangle``'s are configurations in the `Confs` generated in the last section. The ``|I\rangle`` is a linear combination, and can be regarded as a grouping of ``m_I`` configurations.

The function used to implement the discrete symmetries is [`Basis`](@ref). There are three ``\mathbb{Z}_2`` transformations in the Ising model, _viz._ the particle-hole transformation ``\mathscr{P}``, the ``\pi``-rotation along the ``y``-axis ``\mathscr{R}_y`` and the flavour symmetry ``\mathscr{Z}``

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

## Record the Hamiltonian operator
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

### ITensor support

Alternatively, one can generate the operator using an `OpSum` object instead of `cstr_vec` and `fac` using the function [`OperatorFromOpSum`](@ref).

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
hmt = OperatorFromOpSum(bs, bs, ops_hmt ; red_q = 1, sym_q = 1)
```

## Generate the sparse matrix and diagonalise

After specifying the Hamiltonian, we then use the [`OpMat`](@ref) to generate a sparse matrix from the operator. To get the 10 lowest eigenstates and their energies
```julia
@time "Initialise the Hamiltonian matrix" hmtmat = OpMat(hmt)
@show hmtmat.nel
@time "Diagonalise Hamiltonian" enrg, st = GetEigensystem(hmtmat, 10)
@show real(enrg)
```

## Measuring the angular momentum

We can measure the inner product of a final state, an operator or its matrix and an initial state or the action of an operator or its matrix on a state by directly using the [`*`](@ref) operator. To measure the total angular momentum ``L^2`` by definition

```math
\begin{aligned}
L^2&=L^+L^-+(L^z)^2-L^z\\
L^z&=\sum_{\sigma m}mc^\dagger_mc_m\\
L^\pm&=\sum_{\sigma m}\sqrt{(s\mp m)(s\pm m+1)}c^\dagger_{m\pm 1}c_m
\end{aligned}
```

The construction of the operator can be simplified by the (addition)[@ref +(tms1 :: Vector{Term}, tms2 :: Vector{Term})], (multiplication)[@ref *(tms1 :: Vector{Term}, tms2 :: Vector{Term})] and (Hermitian conjugate)[@ref adjoint(tms :: Vector{Term})] of terms. The following code measures the angular momentum of each eigenstate and verify whether ``|T\rangle`` is an eigenstate of ``L^2`` by measuring 

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
tms_l2 = tms_lz * tms_lz - tms_lz + tms_lp * tms_lm
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

## Measuring the density operator

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
@time "Initialise Hamiltonian" hmtmat = OpMat(hmt)
@show hmtmat.nel
@time "Diagonalise Hamiltonian" enrg1, st1 = GetEigensystem(hmtmat, 10)
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

## Use the built-in models

The basis, Hamiltonian and the total angular momentum of the Ising model is built in the package (_cf._ [built-in models](@ref Built-in-models) and [Ising model](@ref Ising-model)). So one can write instead
```julia
nm = 12
cfs = GenerateIsingConfs(nm) # Implement the conserved quantities and generate the confs
bs = GenerateIsingBasis(cfs ; PH = 1, Ry = 1, Z2 = 1) # Implement the discrete symmetries and initialise the basis
tms_hmt = GenerateIsingHamiltonianTerms(nm ; ps_pot = Number[4.75, 1.], fld_h = 3.16)
hmt = Operator(bs, bs, tms_hmt ; red_q = 1, sym_q = 1) # Record the Hamiltonian operator
hmtmat = OpMat(hmt) # Generate the sparse matrix and diagonalise
enrg, st = GetEigensystem(hmtmat, 10)
@show real(enrg)

# Measure the total angular momentum observable
tms_l2 = GenerateL2Terms(nm, 2)
l2 = Operator(bs, bs, tms_l2 ; red_q = 1, sym_q = 1)
l2_mat = OpMat(l2)
l2_val = [ st[:, i]' * l2_mat * st[:, i] for i = 1 : length(enrg)]
@show real(l2_val)
```