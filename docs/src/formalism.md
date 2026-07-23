# Formalism

Our purpose is to build a Hilbert space equipped with full $\mathrm{SO}(3)$ rotation symmetry and part of the non-Abelian flavour symmetry, in addition to the $\mathrm{U}(1)$ quantum numbers of which FuzzifiED is capable, and find the eigen-states and the energy of the Hamiltonian on that basis. Although it is impractical to fully diagonalize the total angular momentum $L^2$ on the full system, such diagonalization is practical on a segment (_e. g._ one of the two flavours in the model for the Ising CFT), and the segment states can assemble into states on the full system through angular-momentum composition
```math
|(l_1l_2)lm\rangle=\sum_{m_1m_2}|l_1m_1\rangle|l_2m_2\rangle\langle l_1m_1,l_2m_2|lm\rangle
```
where $\langle l_1m_1,l_2m_2|lm\rangle$ is the Clebsch-Gordan (CG) coefficient.

Our strategy is therefore to 

* split the system into segments, 
* construct the Hilbert space resolved by $\mathrm{SO}(3)$ spatial rotation and flavour symmetry and the operators on each segment, 
* assemble the segment Hilbert spaces and the segment operators into composed Hilbert space and composed operators,
* find the eigen-system of the Hamiltonian on the composed Hilbert space, and make measurements. 

## Composing the Segment Hilbert Spaces

For simplicity, we first consider bi-partition of the full system. _E. g._, for the Ising model, we split it into two segment Hilbert spaces of spin-up and spin-down. We will build a Hilbert space
```math
    \mathscr{H}(\{Q\},l,\{C_{2,p}\})
```
with an assigned total angular momentum $l$, total $\mathrm{U}(1)$ quantum number $\{Q\}$, and representations under flavour symmetry within each segment specified by the values of the quadratic Casimir $C_{2,p}$ ($p=1,2$). For that purpose, on each of the system we build the basis with a set of $\mathrm{U}(1)$ quantum numbers $\{Q\}_p$ and diagonalize simultaneously the total angular momentum $L^2$ and the quadratic Casimir $C_{2,p}$ of the sub-group of the flavour symmetry on that segment. 
```math
    \mathscr{H}_p(C_{2,p})=\bigoplus_{\{Q\}_p,l_p}\left|\{Q\}_pC_{2,p},l_pm_p,\alpha_p\right\rangle
```
Here, we only need to keep states with a specific value of $C_{2,p}$. On the contrary, although we only need a single value of $\{Q\}$ and $l$ on the composed space, we need to keep all values of $\{Q\}_p$ and $l$ in the segment space, as different sets of $\{Q\}_p$ and $l_p$ can all compose into the same total $\{Q\}$ and $l$. The index $\alpha_p$ denotes the multiplicity. 

Now we build the composed Hilbert space from the segment Hilbert spaces. Taking one state from each segment, the composed state is
```math
    \left|\{Q\}_{\{12\}}C_{2,{\{12\}}},(l_1l_2)lm,\alpha_{\{12\}}\right\rangle=\sum_{m_1m_2}\left|\{Q\}_1C_{2,1},l_1m_1,\alpha_1\right\rangle\left|\{Q\}_2C_{2,p},l_2m_2,\alpha_2\right\rangle\langle l_1m_1,l_2m_2|lm\rangle
```
Here $\{Q\}_{\{12\}}$ is a short-hand notation for $\left\{\{Q\}_1,\{Q\}_2\right\}$, and the same notation applies to the multiplicity $\alpha_{\{12\}}$ and $C_{2,{\{12\}}}$. The segment quantum numbers and angular momenta constraints that $\{Q\}_1+\{Q\}_2=\{Q\}$ and $|l_1-l_2|\leq l\leq l_1+l_2$. 

The composed Hilbert space is then the collection of these states
```math
    \mathscr{H}(\{Q\},l,C_{2,{12}})=\bigoplus_{\{Q\}_{\{12\}},l_{\{12\}}}\left|\{Q\}_{\{12\}}C_{2,{\{12\}}},(l_1l_2)lm,\alpha_{\{12\}}\right\rangle
```
In practice, we pick a state with a representative $m$ within each $\mathrm{SO}(3)$ multiplet, and work in a $L^z$-free treatment. The composition can therefore be denoted by
```math
    \left\|\{Q\}_{\{12\}}C_{2,{\{12\}}},(l_1l_2)l,\alpha_{\{12\}}\right\rangle=\big[\left\|\{Q\}_1C_{2,1},l_1,\alpha_1\right\rangle\otimes\left\|\{Q\}_2C_{2,2},l_2,\alpha_2\right\rangle\big]_l.
```

## Decomposing the Operator

Having contructed the Hilbert space, we now calculate the matrix element of the Hamiltonian. For that purpose, we decompose the Hamiltonian into direct product of -covariant (_i. e._ carrying definite spin) operators acting on the segments that carry definite $\mathrm{SO}(3)$ spin.
```math
    H=\sum_d g_d\Big[[\Phi^{(d)}_1]_{L_1^{(d)}}\otimes[\Phi^{(d)}_2]_{L^{(d)}_2}\Big]_{LM=00}
```
where $d$ is the index of decomposition. As above, this is a short-hand notation for 
```math
    \Big[[\Phi_1]_{L_1}\otimes[\Phi_2]_{L_2}\Big]_{LM}=\sum_{M_{\{12\}}}[\Phi_1]_{L_1M_1}[\Phi_2]_{L_2M_2}\langle L_1M_1,L_2M_2|LM\rangle
```
We take the Ising model as an example and explain how this decomposition is accomplished.

### Recoupling the Pseudo-Potentials

We consider a two-body interaction written in terms of pseudo-potentials
```math
    H=\sum_l\tilde{V}_l\Big[[c^\dagger_1 c^\dagger_2]_l\otimes [c_3c_4]_l\Big]_0
```
_E. g._, in the Ising model, the fermions 1 and 4 are spin-up, the fermions 2 and 3 are spin-down, and $\tilde{V}_{0,1}$ are non-zero. Here, we use the brackets to indicate objects composed by the rule of CG-coefficient, and we adopt a notation where $[\Phi]_{lm}$ _increases_ the $L^z$ by $m$. In this notation, 
```math
    [c^\dagger]_{sm}=c^\dagger_m,\qquad[c]_{sm}=(-1)^{m}c_{-m}.
```
In this convension, the $\tilde{V}_l$ is connected with the ordinary pseudo-potential by $\tilde{V}_l=(-1)^lV_l/\sqrt{2l+1}$. We now want to rewrite the Hamiltonian in a different channel 
```math 
    H=\sum_l\tilde{W}_j\Big[[c^\dagger_1 c_4]_j\otimes [c^\dagger_2c_3]_j\Big]_0.
```
_E. g._, in the Ising model, now $c^\dagger_1 c_4$ is completely in the spin-up segment, and $c^\dagger_2c_3$ is completely in the spin-down segment, so the re-coupling accomplishes the decomposition. The relation of $\tilde{V}$ and $\tilde{W}$ can be expressed in a Dirac bracket notation, which corresponds to a $9j$-symbol that could be further reduced to a $6j$-symbol.
```math
\begin{aligned}
    \tilde{W}_j&=\sum_l\tilde{V}_l\langle((s_1s_2)l,(s_3s_4)l)0|((s_1s_4)j,(s_2s_3)j)0\rangle\\
    &=\sum_l\tilde{V}_l\,(-1)^{s_3+s_4-l}(2l+1)(2j+1)\begin{Bmatrix}s_1&s_2&l\\s_4&s_3&l\\j&j&0\end{Bmatrix}.
\end{aligned}
```

### Translate from a Contact Coupling

We now discuss a second case where the Hamiltonian can be written in terms of a contact coupling 
```math
    H=\int\mathrm{d}^2\mathbf{r}\,\Phi_1(\mathbf{r})\Phi_2(\mathbf{r}).
```
The local operators $\Phi_{1,2}$ can be decomposed into spherical components 
```math
    \Phi^{(s)}(\mathbf{r})=\sum_{lm} (\Phi)_{lm}^{(s)}Y_{lm}^{(s)}.
```
Here we use the parentheses to indicate that these objects compose by the rule of monopole harmonics instead of CG-coefficients.
```math
    (\Phi_1\Phi_2)^{(s)}_{lm}=\sum_{m_1m_2}(\Phi_1)^{(s_1)}_{l_1m_1}(\Phi_2)^{(s_2)}_{l_2m_2}\langle l_1m_1,l_2m_2|lm\rangle\langle l_1(-s_1),l_2(-s_2)|l(-s)\rangle\sqrt{\frac{(2l_1+1)(2l_2+1)}{4\pi(2l+1)}}.
```
_E. g._, for the electron and density operators 
```math
    (c^\dagger)^{(s)}_{sm}=c^\dagger_m,\qquad(c)^{(s)}_{sm}=(-1)^{m-s}c_{-m}\\
    (n)_{lm}=(c^\dagger c)_{lm}=\sum_{m_1m_2}c^\dagger_{m_1} c_{-m_2}\frac{(-1)^{m_2-s}(2s+1)}{\sqrt{4\pi(2l+1)}}\langle sm_1,s(-m_2)|lm\rangle\langle s(-s),ss|l0\rangle.
```
The Hamiltonian can therefore be expressed as 
```math
    H=\sum_l\tilde{V}_l\Big[[\Phi_1]_{l}\otimes[\Phi_2]_{l}\Big]_0
```
where
```math
    [\Phi_p]_{lm}=(\Phi_p)_{lm}\\
    \tilde{V}_l=\langle l(-s_1),l(-s_2)|00\rangle\sqrt{\frac{(2l_1+1)(2l_2+1)}{(2l+1)}}.
```
_E. g._, for the Ising interaction, we take $\Phi_1(\mathbf{r})=n_\uparrow(\mathbf{r})$ and $\Phi_2(\mathbf{r})=n_\downarrow(\mathbf{r})$ or $\nabla^2n_\downarrow(\mathbf{r})$ ; for the transverse field, we take $\Phi_{1,2}(\mathbf{r})=c^{(\dagger)}_{\uparrow,\downarrow}(\mathbf{r})$. They completely lie in one of the two segments ; hence, this expression accomplishes the decomposition.

## Composing the Segment Operators

Having obtained the decomposed Hamiltonian, we first calculate the matrix elements of the segment operator $[H_{d,p}]_{l_d}$. By Wigner-Eckart theorem, the matrix elements with different $L^z$-numbers are connected through a $3j$-symbol. Omitting the irrelevant indices, 
```math
    \langle l'm'|[\Phi]_{LM}|lm\rangle
    =(-1)^{l'-m'}\begin{pmatrix}l'&l&l\\-m'&m&m\end{pmatrix}\langle l'\|[\Phi]_L\|l\rangle.
```
Here, $\langle l'\|[\Phi]_L\|l\rangle$ is independent of the choice of states $m,m'$ in the multiplet and is called the _reduced matrix element_. To calculate the reduced matrix of the segment operator, we need only pick one representative state from each multiplet in the segment Hilber space. We typically pick the representative state as $m'=M=m=0$ except when the $3j$-symbols vanish, _i. e._ when $l'+L+l\in 2\mathbb{Z}+1$, in which cases we pick $m=1$ instead.

We now assemble the segment operator into the composed operator acting on the full Hilbert space. This involves recoupling through the $9j$-symbol.
```math
    \langle(l_1'l_2')l'\|[\Phi_{1,L_1}\otimes\Phi_{2,L_2}]_L\|(l_1 l_2)l\rangle
    =\sqrt{(2l+1)(2l'+1)(2L+1)}\begin{Bmatrix}l_1&l_2&l\\L_1&L_2&L\\l_1'&l_2'&l'\end{Bmatrix}
    \langle l_1'\|\Phi_{1,L_1}\|l_1\rangle\,\langle l_2'\|\Phi_{2,L_2}\|l_2\rangle .
```
Note that when exchanging $|l_1\rangle$ and $\Phi_2$, a factor related to fermion parity may arise. 

## Generalization to Multiple Parts

For multiple parts, the angular-momentum composes successively along a chain, so a state is specified not only by the total angular momentum $l$, but also by the cumulative angular momentum of parts $1\cdots p$ for each $p$
```math
    \|(\cdots((l_1 l_2)l_{12}\,l_3)l_{123}\cdots l_p)l_{1\cdots p}\cdots l_{N_p})l\rangle,
```
Each channel is specified by $2N_p+1$ angular momenta in total — $l_p$ ($p=1,\dots,N_p$), $l_{1\cdots p}$ ($p=2,\dots,N_p-1$) and $l$. The coupling of segment operators are specified in a similar way. 
_E. g_., for the contact coupling
```math
    H=\int\mathrm{d}^2\mathbf{r}\,\Phi_1(\mathbf{r})\Phi_2(\mathbf{r})\cdots\Phi_p(\mathbf{r})\cdots\Phi_{N_p}(\mathbf{r}),
```
it can be rewritten in terms of channels as 
```math
    H=\sum V_{\{L\}}\Big[\dots[[\Phi_1]_{L_1}[\Phi_2]_{L_2}]_{L_{12}}\dots [\Phi_p]_{L_p}\big]_{L_{1\cdots p}} \cdots [\Phi_{N_p}]_{L_{N_p}}\Big]_{L=0}
```
Here we use $\{L\}$ as a short-hand notation for the coupling channel. The potentiaL reads
```math
    V_{\{L\}}=\sqrt{\frac{\prod_p(2L_p+1)}{(4\pi)^{N_p-2}(2L+1)}}\prod_p\langle L_{1\dots (p-1)}(-s_{L\dots(p-1)}),L_p(-s_p)|L_{1\cdots p}(-s_{1\cdots p})\rangle
```
where $s_{1\cdots p}=s_1+\cdots+s_p$. To calculate the matrix element, each re-coupling produces a $9j$-symbol, so the final 
```math
    \langle \{l'\}\|\Phi_{\{L\}}\|\{l\}\rangle
    =\prod_{p=2}^{N_p}\left[\sqrt{(2l_{1\cdots p}+1)(2k_{1\cdots p}+1)(2l'_{1\cdots p}+1)}
    \begin{Bmatrix}l_{1\cdots p-1}&l_p&l_{1\cdots p}\\k_{1\cdots p-1}&L_p&L_{1\cdots p}\\L'_{1\cdots p-1}&l'_p&l'_{1\cdots p}\end{Bmatrix}\right]
    \prod_{p=1}^{N_p}\langle l'_p\|\Phi_{p,L_p}\|l_p\rangle ,
```
