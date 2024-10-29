# An introduction to fuzzy sphere

## Introduction

(Under construction)

## Review of the existing work

The study of 3d CFTs on fuzzy sphere can mainly be devided into four catagories :

1. Accessing various conformal data, 

2. Realising various 3d CFTs,

3. Studying conformal defects and boundaries, and 

4. Exploring applicable numerical techniques.

### Accessing various conformal data

The first direction is to develop methods to calculate as many data and quantities of 3d CFTs on fuzzy sphere. Typically, these methods are tested on the simplest example of 3d Ising CFT. For many of those CFT data, fuzzy sphere is the first non-perturbative method to access them ; for the others, fuzzy sphere has achieved great consistency with previous methods such as quantum Monte Carlo and conformal bootstrap. So far, the accessible CFT data include operator spectrum, OPE coefficients, correlation functions, entropic $F$-function and conformal generators. 

1. __Operator spectrum__ [[Zhu 2022](#References)] This seminal paper opens a new avenue for studying 3d conformal field theories. It calculates and analyses the energy spectra at the 3d Ising transition, and explicitly demonstrate the state-operator correspondence as a fingerprint of conformal field theory, thus directly elucidates the emergent conformal symmetry of the 3d Ising transition.

2. __OPE coefficients__ [[Hu 2023Mar](#References)] This paper computes 17 OPE coefficients of low-lying CFT primary fields with high accuracy, including 4 that has not being reported before.

3. __Correlation functions__ [[Han 2023Jun](#References)] This paper computes the four-point correlators and verify the crossing symmetry. 

4. __Entropic $F$-function__ [[Hu 2024](#References)] This paper have performed the first non-perturbative computation of the $F$-function for the paradigmatic 3d Ising conformal field theory through entanglement entropy. 

5. __Conformal generators__ [[Fardelli 2024](#References), [Fan 2024](#References)] These papers investigate the conformal generators of translations and special conformal transformations which are emergent in the infrared and construct these generators using the energy momentum tensor.

### Realising various 3d CFTs

The second direction is study various other CFTs beyond 3d Ising. Fuzzy sphere has revealed many new information about these theories ; the previously known data are also consistent with fuzzy sphere results. So far, the accessible CFTs include $\mathrm{SO}(5)$ deconfined criticality, $\mathrm{O}(3)$ Wilson-Fisher and a series of new theories with $\mathrm{Sp}(N)$ symmetry. 

1. __The $\mathrm{SO}(5)$ deconfined criticality__ [[Zhou 2023](#References)] This paper provides clear evidence that the DQCP exhibits approximate conformal symmetry, and demonstrate that the DQCP is more likely pseudo-critical.

2. __The $\mathrm{O}(3)$ Wilson-Fisher__ [[Han 2023Dec](#References)] This paper design a microscopic model of Heisenberg magnet bilayer and study the underlying Wilson-Fisher $\mathrm{O}(3)$ transition through the lens of fuzzy sphere regularization. 

3. __A series of new $\mathrm{Sp}(N)$-symmetric CFTs__ [[Zhou 2024Oct](#References)] This paper discovers a series of new CFTs with global symmetry $\mathrm{Sp}(N)$ in the fuzzy sphere models that are closely related to the SO(5) deconfined phase transition, and are related to non-linear sigma model with a Wess-Zumino-Witten term and Chern-Simons-matter theories. The emergent conformal symmetry is numerically verified by observing the integer-spaced conformal multiplets and the quality of conformal generators. 

### Studying conformal defects and boundaries

Apart from the bulk CFTs, fuzzy sphere can also be used to study their conformal defects and boundaries. So far, the accessible defects/boundaries include the magnetic line defect of 3d Ising CFT, including its defect operator spectrum, correlators, $g$-function, defect changing operators, and its cusp, and the conformal boundaries of 3d Ising CFT.

1. __Conformal magnetic line defect__ [[Hu 2023Aug](#References)] This paper studies the magnetic line defect of 3D Ising CFT and clearly demonstrates that it flows to a conformal defect fixed point. The authors have identified 6 low-lying defect primary operators and extract their scaling dimensions, as well as computing one-point bulk correlators and two-point bulk-defect correlators.

2. __The $g$-function and defect changing operators__ [[Zhou 2024Jan](#References)] This paper have performed the non-perturbative computations of the scaling dimensions of defect-changing, creation operators and the $g$-function for the pinning defect in 3d Ising model. 

3. __Cusp__ [[Cuomo 2024](#References)] This paper study the general properties of the cusp anomalous dimension and in particular calculates the pinning field defects in the 3d Ising model on fuzzy sphere.

4. __Conformal boundaries of 3d Ising CFT__ [[Zhou 2024Jul](#References), [Dedushenko 2024](#References)] These papers demonstrates that conformal field theory (CFT) with a boundary, known as surface CFT in three dimensions, can be studied with the setup of fuzzy sphere, and in particular in the example of surface criticality, proposes two schemes by cutting a boundary in the orbital space or the real space to realise the ordinary and the normal surface CFTs on the fuzzy sphere.

### Exploring applicable numerical techniques

So far, the numerical methods that has been applied to fuzzy sphere to include exact diagonalisation (ED), density matrix renormalisation group (DMRG) and determinant quantum Monte Carlo (DQMC). The former two has been reviewed in previous sections.

1. __Quantum Monte Carlo on fuzzy sphere__ [[Hofmann 2023](#References)] This paper presents a numerical quantum Monte Carlo (QMC) method for simulating the 3D phase transition on the recently proposed fuzzy sphere.

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

## Numerical methods

(under construction)

## References

For a more detailed summary of the background, please visit [« The fuzzified world »](https://www.fuzzified.world/fuzzified-world)

* __[Zhu 2022]__ Uncovering Conformal Symmetry in the 3D Ising Transition: State-Operator Correspondence from a Quantum Fuzzy Sphere Regularization, Wei Zhu, Chao Han, Emilie Huffman, Johannes S. Hofmann, and Yin-Chen He, [Phys. Rev. X 13, 021009 (2023)](https://doi.org/10.1103/PhysRevX.13.021009).
* __[Hu 2023Mar]__ Operator Product Expansion Coefficients of the 3D Ising Criticality via Quantum Fuzzy Sphere, Liangdong Hu, Yin-Chen He, and Wei Zhu, [Phys. Rev. Lett 131, 031601 (2023)](https://doi.org/10.1103/PhysRevLett.131.031601).
* __[Han 2023Jun]__ Conformal four-point correlators of the 3D Ising transition via the quantum fuzzy sphere, Chao Han, Liangdong Hu, Wei Zhu, and Yin-Chen He, [Phys. Rev. B 108, 235123 (2023)](https://doi.org/10.1103/PhysRevB.108.235123).
* __[Zhou 2023]__ The ``\mathrm{SO}(5)`` Deconfined Phase Transition under the Fuzzy Sphere Microscope: Approximate Conformal Symmetry, Pseudo-Criticality, and Operator Spectrum, Zheng Zhou, Liangdong Hu, Wei Zhu, and Yin-Chen He, [Phys. Rev. X 14, 021044 (2024)](https://doi.org/10.1103/PhysRevX.14.021044).
* __[Hu 2023Aug]__ Solving Conformal Defects in 3D Conformal Field Theory using Fuzzy Sphere Regularization, Liangdong Hu, Yin-Chen He, and Wei Zhu, [Nat. Commun. 15, 3659 (2024)](https://doi.org/10.1038/s41467-024-47978-y).
* __[Hofmann 2024]__ Quantum Monte Carlo Simulation of the 3D Ising Transition on the Fuzzy Sphere, Johannes S. Hofmann, Florian Goth, Wei Zhu, Yin-Chen He, and Emilie Huffman, [SciPost Phys. Core 7, 028 (2024)](https://doi.org/10.21468/SciPostPhysCore.7.2.028).
* __[Han 2023Dec]__ Conformal Operator Content of the Wilson-Fisher Transition on Fuzzy Sphere Bilayers, Chao Han, Liangdong Hu, and Wei Zhu, [Phys. Rev. B 110, 115113 (2024)](https://doi.org/10.1103/PhysRevB.110.115113).
* __[Zhou 2024Jan]__ The ``g``-function and Defect Changing Operators from Wavefunction Overlap on a Fuzzy Sphere, Zheng Zhou, Davide Gaiotto, Yin-Chen He, Yijian Zou, [SciPost Phys. 17, 021 (2024)](https://doi.org/10.21468/SciPostPhys.17.1.021).
* __[Hu 2024]__ Entropic ``F``-function of 3D Ising conformal field theory via the fuzzy sphere regularization, Liangdong Hu, Wei Zhu, and Yin-Chen He, [arXiv : 2401.17362](https://arxiv.org/abs/2401.17362).
* __[Cuomo 2024]__ Impurities with a cusp: general theory and 3d Ising, Gabriel Cuomo, Yin-Chen He, Zohar Komargodski, [arXiv : 2406.10186](https://arxiv.org/abs/2406.10186). 
* __[Zhou 2024Jul]__ Studying the 3d Ising surface CFTs on the fuzzy sphere, Zheng Zhou, and Yijian Zou, [arXiv : 2407.15914](https://arxiv.org/abs/2407.15914).
* __[Dedushenko 2024]__ Ising BCFTs from the fuzzy hemisphere, Mykola Dedushenko, [arXiv : 2407.15948](https://arxiv.org/abs/2407.15948).
* __[Fardelli 2024]__ Constructing the Infrared Conformal Generators on the Fuzzy Sphere, Giulia Fardelli, A. Liam Fitzpatrick, and Emanuel Katz, [arXiv : 2409.02998](https://arxiv.org/abs/2409.02998).
* __[Fan 2024]__ Note on explicit construction of conformal generators on the fuzzy sphere, Ruihua Fan, [arXiv : 2409.08257](https://arxiv.org/abs/2409.08257).
* __[Zhou 2024Oct]__ A new series of 3D CFTs with ``\mathrm{Sp}(N)`` global symmetry on fuzzy sphere, Zheng Zhou, and Yin-Chen He, [arXiv : 2410.00087](https://arxiv.org/abs/2410.00087).