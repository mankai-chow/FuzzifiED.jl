# `FuzzifiED` explained in an example

The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculation on the fuzzy sphere. In this example, we will illustrate how to use  `FuzzifiED` to calculate the spectrum of Ising model on fuzzy sphere and how to calculate the OPE coefficient $\lambda_{\sigma\sigma\epsilon}$ by measuring the expectation value of the density operator $n^z$. 

## Compile and header

The compiled `lib_fuzzifi_ed.so` file can be used directly. Alternatively, one can also compile it from the Fortran source files by
```bash
ifort -shared -fPIC -larpack -qopenmp -O3 -o lib_fuzzifi_ed.so ./fuzzifi_ed_fortran/*.f90
```
Include at the start of your Julia script
```julia
LibpathFuzzifiED = "./lib_fuzzifi_ed.so"
include("./fuzzified.jl")
```
where `LibpathFuzzifiED` points to the Path of the `.so` file

## Implement the conserved quantities and generate the configurations

`FuzzifiED` supports conserved quantities in the form of 

$$
Q_i=\sum_{o=1}^{N_o}q_{io}n_o
$$

where $i=1,\dots,N_U$ is the index of conserved quantities, $o$ is the index of orbital, $n_o=c^\dagger_oc_o$, and $q_o$ is a set of coefficients that must be non negative integer valued. 

The function used to implement the conserved quantities and generate all the configurations (_i.e._, direct product states) is 
```julia
function Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Vector{Int64}} ; nor = div(no, 2) :: Int64) :: Confs
```
where the arguments are 
* `no :: Int64` is the number of orbitals $N_o$ ;
* `qnu_s :: Vector{Int64}` is the set of $Q_i$ for the selected configurations ;
* `qnu_o :: Vector{Vector{Int64}}` is the set of $q_{io}$ for each quantum number and for each orbital. It should contain $N_U$ elements and each element should be a vector of length $N_o$. 
* (`nor :: Int64` is the number of less significant bits used to generate the Lin table.)

As an example, in the Ising model, there are two conserved quantities, _viz._ the particle number and the angular momentum

$$
\begin{aligned}
Q_1&=N_e,& q_{1,m\sigma}&=1\\
Q_2&=L_z+N_es,&q_{2,m\sigma}&=m+s
\end{aligned}
$$

where the orbital index $o$ contains both $m$ and $\sigma$. In the code, we store the spin-up orbitals in $o=1,\dots,N_m$ and the spin-down orbitals in $o=N_m+1,\dots,2N_m$. Thus, if we want to look at the $L_z=0$ half-filled sector

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
push!(qnu_o, vcat(fill(collect(0 : nm - 1), nf)...)) # qnu_o[2] = [0,1,...,7,1,2,...,7] to qnu_o
push!(qnu_s, ne * s) 
# Generate the configurations and print the number
@time "Initialise Basis 0" cfs = Confs(no, qnu_s, qnu_o)
@show cfs.ncf
```
The resulting `Confs` object contains several elements :
* `no :: Int64` is the number of orbitals
* `ncf :: Int64` is the number of configurations 
* `conf :: Vector{Int64}` is an array of length `ncf` containing all the configurations. Each configuration is expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th orbital in the `i`-th configuration is occupied ; if the bit is 0, then the orbital is empty. 
* `nor :: Int64`, `lid :: Vector{Int64}` and `rid :: Vector{Int64}` contain the information of Lin table that is used to inversely look up the index `i` from the configuration. 

## Implement the  discrete symmetries and initialise the basis

`FuzzifiED` supports discrete $\mathbb{Z}_n$ symmetries in the form of 

$$
\mathscr{Z}:\ c_o\mapsto \alpha_o^* c^{(p_o)}_{\pi_o},\quad c_o^\dagger\mapsto \alpha_o c^{(1-p_o)}_{\pi_o}
$$

where we use a notation $c^{(1)}=c^\dagger$ and $c^{0}=c$ for convenience, where $\pi_o$ is a permutation of $1,\dots N_o$, $\alpha_o$ is a coefficient, and $p_o$ specified whether or not particle-hole transformation is performed for the orbital. Note that one must guarentee that all these transformations commute with each other and also commute with the conserved quantities. 

After implementing these symmetries, a state in the new basis should look like 

$$
|I\rangle=\lambda_{i_{I1}}|i_{I1}\rangle+\lambda_{i_{I2}}|i_{I2}\rangle+\cdots+\lambda_{i_{Im_I}}|i_{Im_I}\rangle
$$

where the $|i\rangle$'s are configurations in the `Confs` generated in the last section. The $|I\rangle$ is a linear combination, and can be regarded as a grouping of $m_I$ configurationss.

The function used to implement the discrete symmetries is 
```julia
Basis(cfs :: Confs, qnz_s :: Vector{ComplexF64}, cyc :: Vector{Int64}, perm_o :: Vector{Vector{Int64}}, ph_o :: Vector{Vector{Int64}}, fac_o :: Vector{Vector{ComplexF64}}) :: Basis
```
where the arguments are 
* `cfs :: Confs` is the configuration set with only conserved quantities generated in the last step ;
* `qnz_s :: Vector{ComplexF64}` is a vector of length the same as the number of discrete symmetries $N_Z$ that records the eigenvalue of each transformation ;
* `cyc :: Vector{Int64}` records the cycle of each transformation. For $\mathbb{Z}_n$ symmetry, record $n$ ;
* `perm_o :: Vector{Vector{Int64}}` records the permutation $\pi_o$. It has $N_Z$ elements and each of its elements is a vector of length $N_o$. 
* `ph_o :: Vector{Vector{Int64}}` records $p_o$ to determine whether or not to perform a particle-hole transformation. It has $N_Z$ elements and each of its elements is a vector of length $N_o$. 
* `fac_o :: Vector{Vector{ComplexF64}}` records the factor $p_o$ is determine whether or not to perform a particle-hole transformation. It has $N_Z$ elements and each of its elements is a vector of length $N_o$. 

As an example, in the Ising model, there are three $\mathbb{Z}_2$ transformations, _viz._ the particle-hole transformation $\mathscr{P}$, the $\pi$-rotation along the $y$-axis $\mathscr{R}_y$ and the flavour symmetry $\mathscr{Z}$

$$
\begin{aligned}
    \mathscr{P}:c^\dagger_{\sigma m}&\mapsto\sigma c_{-\sigma,m}\\
    \mathscr{Z}:c^\dagger_{\sigma m}&\mapsto c^\dagger_{-\sigma,m}\\
    \mathscr{R}_y:c^\dagger_{\sigma m}&\mapsto c^\dagger_{\sigma,-m}\\
\end{aligned}
$$

Thus, if we want to look at the all-positive sector
```julia
cyc = [ 2, 2, 2 ] # Input three Z_2 symmetries 
qnz_s = ComplexF64[ 1, 1, 1 ] # Quantum numbers are all positive 
# Initialise the vectors
perm_o = []
ph_o = []
fac_o = []
# Record the parity
push!(perm_o, [ collect(nm + 1 : 2 * nm) ; collect(1 : nm) ]) # perm_o[1] = [9,10,...,16,1,2,...,8]
push!(ph_o, fill(1, no)) # ph_o[1] = [1,1,...,1] meaning PH
push!(fac_o, [ fill(ComplexF64(1), nm) ; fill(ComplexF64(-1), nm) ]) # fac_o[1] = [1,1,...,1,-1,-1,...,-1]
# Record the flavour symmetry
push!(perm_o, [ collect(nm + 1 : 2 * nm) ; collect(1 : nm) ]) # perm_o[2] = [9,10,...,16,1,2,...,8]
push!(ph_o, fill(0, no)) # ph_o[2] = [0,0,...,0] meaning no PH
push!(fac_o, fill(ComplexF64(1), no)) # fac_o[2] = [1,1,...,1]
# Record the pi-rotation
push!(perm_o, [ collect(nm : -1 : 1) ; collect(2 * nm : -1 : nm + 1) ]) # perm_o[3] = [8,7,...,1,16,15,...,9]
push!(ph_o, fill(0, no)) # ph_o[3] = [0,0,...,0] meaning no PH
push!(fac_o, fill(ComplexF64(1), no)) # fac_o[3] = [1,1,...,1]
# Generate the basis and print the dimension
@time "Initialise Basis Z" bs = Basis(cfs, qnz_s, cyc, perm_o, ph_o, fac_o)
@show bs.dim 
```
The resulting `Basis` object contains several elements 
* `cfs :: Confs` is the basis with only conserved quantities generated in the last step ;
* `dim :: Int64` is the dimension of the basis ;
* `szz :: Int64` records the maximum size $\max m_g$ of groups;
* `cfgr :: Vector{Int64}` is a vector of length `cfs.ncf` and records which group $|I\rangle$ each configuration $|i\rangle$ belong to ;
* `cffac :: Vector{ComplexF64}` is a vector of length `cfs.ncf` and records the coefficients $\lambda_i$ ;
* `grel :: Matrix{Int64}` is a `szz`$\times$`dim` matrix that records the configurations in each group $|i_{I1}\rangle,\dots,|i_{Im_I}\rangle$
* `grsz :: Vector{Int64}` is a vector of length `dim` that records the size $m_I$ of each group.

Note that if no discrete symmetry is needed, one can simply put instead 
```julia
bs = Basis(conf)
```

## Record the Hamiltonian operator
The operator here refers to the sum of product of $c$ and $c^\dagger$'s in the form 

$$
\Phi=\sum_{t=1}^{N_t}U_tc^{(p_{t1})}_{o_{t1}}c^{(p_{t2})}_{o_{t2}}\dots c^{(p_{tl})}_{o_{tl}}
$$

where $c^{(0)}=c$ and $c^{(1)}=c^\dagger$. Here the operator string sum is recorded together with the basis of the initial state and the basis of the final state
```julia
function Operator(bsd :: Basis, bsf :: Basis, red_q :: Int64, sym_q :: Int64, cstr_vec :: Vector{Vector{Int64}}, fac :: Vector{ComplexF64}) :
```
where the arguments are 
* `bsd :: Basis` is the basis of the initial state ;
* `bsf :: Basis` is the basis of the final state ;
* `red_q :: Int64` is a flag that records whether or not the conversion to a sparse martrix can be simplified : if `bsd` and `bsf` have exactly the same quantum numbers, and all the elements in `bsd.cffac` and `bsf.cffac` has the same absolute value, then `red_q = 1` ; otherwise `red_q = 0` ; 
* `sym_q :: Int64` records the symmetry of the operator : if the matrix is Hermitian, then `sym_q = 1` ; if it is symmetric, then `sym_q = 2` ; otherwise `sym_q = 0`. 
* `cstr_vec :: Vector{Vector{Integer}}` records the $c$ and $c^\dagger$ string of each term. A term with $l$ operators $c^{(p_{t1})}_{o_{t1}}c^{(p_{t2})}_{o_{t2}}\dots c^{(p_{tl})}_{o_{tl}}$ correspond to a length-$2l$ vector $(p_{t1},o_{t1},p_{t2},o_{t2},\dots p_{tl},o_{tl})$. Note that each element can have different length. Also note that this format bears a certain similarity with the `OpSum` in `ITensors` ; 
* `fac :: Vector{ComplexF64}` corresponds to the factor $U_t$ in each term.

For example, the Hamiltonian for the fuzzy sphere Ising model

$$
H=\sum_{m_1m_2m_3m_4}U_{m_1m_2m_3m_4}\delta_{m_1+m_2,m_3+m_4}c^\dagger_{m_1\uparrow}c^\dagger_{m_2\downarrow}c_{m_3\downarrow}c_{m_4\uparrow}-h\sum_m(c^\dagger_{m\uparrow}c_{m\downarrow}+\mathrm{h.c.})
$$

can be recorded with 

``` julia
using WignerSymbols
# Input the parameters of the Hamiltonian
ps_pot = [ 4.75, 1. ] * 2.
h = 3.16
cstr_hmt = []
fac_hmt = Array{ComplexF64, 1}(undef, 0)
# Go through all the m1-up, m2-down, m3-down, m4-up and m4 = m1 + m2 - m3
for m1 = 0 : nm - 1
    f1 = 0
    o1 = m1 + f1 * nm + 1
    m1r = m1 - s
    for m2 = 0 : nm - 1
        f2 = 1
        o2 = m2 + f2 * nm + 1
        m2r = m2 - s
        for m3 = 0 : nm - 1
            f3 = 1
            o3 = m3 + f3 * nm + 1
            m3r = m3 - s
            m4 = m1 + m2 - m3 
            if (m4 < 0 || m4 >= nm) continue end
            f4 = 0
            o4 = m4 + f4 * nm + 1
            m4r = m4 - s
            # Calculate the matrix element val from pseudopotentials
            val = ComplexF64(0)
            for l in 1 : length(ps_pot)
                if (abs(m1r + m2r) > nm - l || abs(m3r + m4r) > nm - l) break end 
                val += ps_pot[l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, m1r, m2r, -m1r - m2r) * wigner3j(s, s, nm - l, m4r, m3r, -m3r - m4r)
            end 
            # Record the interaction term val, "Cdag", o1, "Cdag", o2, "C", o3, "C", o4
            push!(cstr_hmt, [1, o1, 1, o2, 0, o3, 0, o4])
            push!(fac_hmt, val)
        end
    end
    o1x = o1 + nm
    # Record the transverse field term
    push!(cstr_hmt, [1, o1, 0, o1x])
    push!(fac_hmt, -h)
    push!(cstr_hmt, [1, o1x, 0, o1])
    push!(fac_hmt, -h)
end
# Generate the Hamiltonian operator
hmt = Operator(bs, bs, 1, 1, cstr_hmt, fac_hmt)
```

## Generate the sparse matrix and diagonalise

After specifying the Hamiltonian, we then use the function `OpMat` to generate a sparse matrix from the operator.
```julia
OpMat(op :: Operator) :: OpMat
```
The fields of an `OpMat` object are
* `dimd :: Int64` and `dimf :: Int64` are the number of columns and rows of the matrix ;
* `symq :: Int64` records whether or not the matrix is Hermitian or symmetric ;
* `nel :: Int64` records the number of elements ;
* `colptr :: Vector{Int64}`, `rowid :: Vector{Int64}` and `elval :: Vector{ComplexF64}` records the elements of the sparse matrix as in the `SparseMatrixCSC` elements of Julia. 

After that, the function `GetEigensystem` can be used to find the lowest eigenvalues and eigenstates
```julia
OpMat(op :: Operator) :: OpMat
GetEigensystem(mat :: OpMat, nst :: Int64 ; tol = 1E-8 :: Float64 ) :: Tuple{Vector{ComplexF64}, Matrix{ComplexF64}}
```
where the outputs are 
* A length-`nst` array recording the eigenvalues, and 
* A `dimd`$\times$`nst` matrix where every column records an eigenstate. 

For example, in the Ising model, to get the 10 lowest eigenstates and their energies
```julia
@time "Initialise the Hamiltonian matrix" hmtmat = OpMat(hmt)
@show hmtmat.nel
@time "Diagonalise Hamiltonian" enrg, st = GetEigensystem(hmtmat, 10)
@show real(enrg)
```

## Measuring the angular momentum

We can measure the inner product of a final state, an operator or its matrix and an initial state by directly using the `*` operator. 
```julia 
*(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat, st_d :: Vector{ComplexF64}) :: ComplexF64
*(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, op :: Operator, st_d :: Vector{ComplexF64}) :: ComplexF64
```
where `st_d` must be of length `op.bsd.dim` or `mat.dimd` and `st_fp` must be of length `op.bsf.dim` or `mat.dimf`, and `st_fp` must be an adjoint. 
We can similarly measure the action of an operator or its matrix on a state 
```julia
*(mat :: OpMat, st_d :: Vector{ComplexF64}) :: Vector{ComplexF64}
*(op :: Operator, st_d :: Vector{ComplexF64}) :: Vector{ComplexF64}
```
where `st_d` must be of length `op.bsd.dim` or `mat.dimd` and the result is a length-`op.bsf.dim` or `mat.dimf` vector. 

For example, to measure the total angular momentum $L^2$ by definition

$$
\begin{aligned}
L^2&=L^+L^-+(L^z)^2-L^z\\
L^z&=\sum_{\sigma m}mc^\dagger_mc_m\\
L^\pm&=\sum_{\sigma m}\sqrt{(s\mp m)(s\pm m+1)}c^\dagger_{m\pm 1}c_m
\end{aligned}
$$

The following code measures the angular momentum of each eigenstate and verify whether $|T\rangle$ is an eigenstate of $L^2$ by measuring 

$$
|L^2T\rangle=L^2|T\rangle,\quad\frac{|\langle T|L^2T\rangle|^2}{\langle T|T\rangle\langle L^2T|L^2T\rangle}
$$

```julia
cstr_l2 = []
fac_l2 = Array{ComplexF64, 1}(undef, 0)
for o1 = 1 : no 
    m1 = mod(o1 - 1, nm) 
    # record the -Lz term
    push!(cstr_l2, [1, o1, 0, o1])
    push!(fac_l2, -(m1 - s))
    for o2 = 1 : no 
        m2 = mod(o2 - 1, nm)
        # record the Lz^2 term
        push!(cstr_l2, [1, o2, 0, o2, 1, o1, 0, o1])
        push!(fac_l2, (m1 - s) * (m2 - s))
        if m1 == nm - 1 continue end
        if m2 == 0 continue end 
        # record the L+L- term
        push!(cstr_l2, [1, o1 + 1, 0, o1, 1, o2 - 1, 0, o2])
        push!(fac_l2, sqrt(m2 * (nm - m2) * (m1 + 1) * (nm - m1 - 1)))
    end
end
# Initialise the L2 operator
l2 = Operator(bs, bs, 1, 1, cstr_l2, fac_l2)
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

$$ 
\lambda_{\sigma\sigma\epsilon}=\frac{\langle\sigma|n^z_{00}|\epsilon\rangle}{\langle\sigma|n^z_{00}|\mathbb{I}\rangle},\quad n^z_{00}=\frac{1}{N_m}\sum_{\sigma m}\sigma c^\dagger_{\sigma m}c_{\sigma m}
$$

To do that, we need to first repeat the calculation in the $\mathbb{Z}_2$-odd sector
```julia
# Repeat the calculation for the Z_2 odd sector (with subscript 1)
qnz_s1 = ComplexF64[ 1, -1, 1 ] # Change only the discrete quantum numbers and generate the basis
@time "Initialise Basis Z" bs1 = Basis(cfs, qnz_s1, cyc, perm_o, ph_o, fac_o) 
@show bs1.dim 
hmt = Operator(bs1, bs1, 1, 1, cstr_hmt, fac_hmt) # Generate and diagonalise Hamiltonian in the new basis
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
cstr_nz = []
fac_nz = Array{ComplexF64, 1}(undef, 0)
for o1u = 1 : nm
    o1d = o1u + nm
    push!(cstr_nz, [1, o1u, 0, o1u])
    push!(fac_nz, 1 / nm)
    push!(cstr_nz, [1, o1d, 0, o1d])
    push!(fac_nz, -1 / nm)
end
# The nz operator sends a state in bs (+) to bs1 (-)
nz = Operator(bs, bs1, 1, 0, cstr_nz, fac_nz)
# Measuring the finite size OPE
@show abs((st_s' * nz * st_e) / (st_s' * nz * st_I))
```