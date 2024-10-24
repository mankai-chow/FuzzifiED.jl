# An introduction to fuzzy sphere

## Introduction

(Under construction)

## Literature review

(Under construction)

## Model construction 

### Projection onto the lowest Landau level

To build the setup of fuzzy sphere, we consider a sphere with radius $R$ and put a $4\pi s$-monopole at its centre. Consider free electrons moving on the sphere. The monopole will modify the single particle Hamilltonian. 

```math
    H_0=\frac{1}{2MR^2}(\partial^\mu+iA^\mu)^2
```

where $\mu=\theta,\phi$ and we take 

```math
    A_\theta=0, A_\phi=-\frac{s}{R}\operatorname{ctg}\theta
```

The eigenstates of the Hamiltonian are the monopole spherical harmonics

```math
    Y_{lm}^{(s)}(\hat{\mathbf{r}}),\quad l=s,s+1,\dots,\quad m=-s,\dots,s-1,s
```

with the eigenenergies 

```math
    E_l=\frac{1}{2MR^2}(l(l+1)-s^2)
```

Each level, known as a Landau level, has a degeneracy of $(2l+1)$. Specifically, the wavefunctions on the lowest Landau level (LLL) $l=s$ is easy to write out :

```math
    Y_{sm}^{(s)}(\hat{\mathbf{r}})=C_me^{im\phi}\cos^{s+m}\frac{\theta}{2}\sin^{s-m}\frac{\theta}{2},\quad C_m=\frac{1}{{\sqrt{4\pi\Beta(s+m+1,s-m+1)}}}
```

where $C_m$ is the normalising factor. The LLL has a degeneracy $N_m=2s+1$

We now consider $N_f$ flavours of fermions moving on the sphere, characterised by the second-quantised fermion operator $\hat{\psi}_f(\hat{\mathbf{r}})$, with a flavour index $f=1,\dots,N_f$. We partially fill the lowest Landau level and set the single energy gap to be much larger than the scale of interaction $H_0\gg H_\mathrm{int}$, so that the quantum fluctuation can be constrained on the lowest Landau level. In practice, we often fill integer number of flavours $N_e=kN_m$ ($k\in\mathbb{Z}$) so that a quantum Hall ferromagnet (_i.e._, the state where integer number of LLLs are filled) can exist on the phase diagram. 

We then project the system onto the LLL. Technically, this can be done by write the fermion operators in terms of the annihilation operators of the LLL orbitals
```math
    \hat{\psi}_f(\hat{\mathbf{r}})=\sum_{m=-s}^s Y^{(s)}_{sm}(\hat{\mathbf{r}})\hat{c}_{mf}
```
where $\hat{c}^{(\dagger)}_{mf}$ annihilates/creates an electron with $L^z$-quantum number $m$ at the $f$-th flavour of the lowest Landau level. In the following sections, we will omit the hats on the operators.

### Density operator

The simplest building block of an interacting many-body Hamiltonian is density operators, which are local fermion bilinears. 

```math
    n_M(\hat{\mathbf{r}})=\psi_{f'}^\dagger(\hat{\mathbf{r}})M_{f'f}\psi_f(\hat{\mathbf{r}})
```

Here the matrix insertion $M$ put the density operators at a certain representation of the flavour symmetry. Like the fermion operator, the density operator can also be expressed in the orbital space. 

```math
    n_M(\hat{\mathbf{r}})=\sum_{lm}Y_{lm}(\hat{\mathbf{r}})n_{lm}
```

Conversely,

```math
    \begin{aligned}
        n_{lm}&=\int\mathrm{d}^2\hat{\mathbf{r}}\,\bar{Y}_{lm}n_M(\hat{\mathbf{r}})\\
        &=\int\mathrm{d}^2\hat{\mathbf{r}}\,\bar{Y}_{lm}\left(\sum_{m_1}\bar{Y}^{(s)}_{sm_1}c^\dagger_{m_1f_1}\right)M_{f_1f_2}\left(\sum_{m_2}Y^{(s)}_{sm_2}c_{m_1f_2}\right)\\
        &=\sum_{m_1m_2}c^\dagger_{m_1f_1}M_{f_1f_2}c_{m_1f_2}\int\mathrm{d}^2\hat{\mathbf{r}}\,\bar{Y}_{lm}\bar{Y}^{(s)}_{sm_1}Y^{(s)}_{sm_2}\\
        &=\sum_{m_1}c^\dagger_{m_1f_1}M_{f_1f_2}c_{m+m_1,f_2}(-1)^{s+m+m_1}(2s+1)\sqrt{\frac{2l+1}{4\pi}}\begin{pmatrix}s&l&s\\m_1&m&-m_1-m\end{pmatrix}\begin{pmatrix}s&l&s\\m_1&m&-m_1-m\end{pmatrix}
    \end{aligned}
```

where [various properties of the monopole spherical harmonics](https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics) are used. In this way, we have fully expressed the density operator in terms of the operators in the orbital space $c^{(\dagger)}_{mf}$. 

### Density-density interaction

The most straightforward way to construct an interaction term is to add a density-density interaction potential $U(r)$. We note that this is not the simplest construction and we will present the simpler construction in terms of pseudopotentials in the next section. 

```math
    H_\mathrm{int}=\int\mathrm{d}^2\hat{\mathbf{r}}_1\,\mathrm{d}^2\hat{\mathbf{r}}_2\,U(|\hat{\mathbf{r}}_1-\hat{\mathbf{r}}_2|)n_M(\hat{\mathbf{r}}_1)n_M(\hat{\mathbf{r}}_2)
```
The interacting potentials can be expanded in terms of the Legendre polynomials 
```math
    U(\theta_{12})=\sum_lU_lP_l(\cos\theta_{12})=\sum_{lm}\frac{4\pi}{2l+1}\bar{Y}_{lm}(\hat{\mathbf{r}}_1)Y_{lm}(\hat{\mathbf{r}}_2)
```
Conversely
```math
    U_l=\int \sin\theta_{12}\mathrm{d}\theta_{12}\,\frac{2l+1}{2}U(\theta_{12})P_l(\cos\theta_{12})
```
Specifically, for local and super-local interactions
```math
    \begin{aligned}
        U(|\mathbf{r}_{12}|)&=\delta(\mathbf{r}_{12})&U_l&=2l+1\\
        U(|\mathbf{r}_{12}|)&=\nabla^2\delta(\mathbf{r}_{12})&U_l&=-l(l+1)(2l+1)
    \end{aligned}
```
By expanding the density operators into the orbital space and completing the integrals,
```math
    H_\mathrm{int}=\sum_{lm}\frac{4\pi U_l}{2l+1}n^\dagger_{lm}n_{lm}
```

### Interaction in terms of pseudopotentials

(under construction)

### Conformal generators

(under construction)