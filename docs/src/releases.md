# Release notes 

## July 2024

### 28 July, 2024 (Version 0.8.2)

- Fix bugs in GetDenIntTerms and multiplication in SphereObs. 

### 26 July, 2024 (Version 0.8.0)

- Improve the performance of SimplifyTerms. 
- Add file operation. 

### 24 July, 2024 (Version 0.7.2)

- Add example and support for surface CFTs. 
- Add example of Ising cusp.

## June 2024

### 11 June, 2024 (Version 0.7.1)

- Add some new interfaces for built-in operators. 
- Add new examples. 
- Fix bugs

### 9 June, 2024 (Version 0.7.0)

- Revise the implementation of diagonal and off-diagonal quantum number. 

### 8 June, 2024 (Version 0.6.3)

- Add support for calculating entanglement spectrum. 
- Add global parameters to control the number of threads, the output and the path of the dynamic library. 
- Fix bugs. 

### 5 June, 2024 (Version 0.6.0)

- Add support for full diagonalisation. 
- Fix bugs and typos.
- Ready for formal release !

### 3 June, 2024 (Version 0.5.8)

- Change the binary dependence to Julia Binary Builder. 

## May 2024

### 30 May, 2024 (Version 0.5.0)

- Enable simplification of terms.
- Add general observables and built-in electrons and density operators. 
- Reorganise the realisations of built-in models.
- Cancel ITensorMPOConstruction dependence. 

### 29 May, 2024 (Version 0.4.3)

- Add QNU truncation for ITensors use.
- Change the Fortran code to be robust against QNU breaking terms
- Add built-in density-density interaction. 
- Add built-in 3-state Potts model.
- Add built-in Ising model with magnetic line defect. 

### 28 May, 2024 (Version 0.4.0)

- Add support for DMRG.
- Add convertion from diagonal QNs to sites. 
- Add the support of ``\mathbb{Z}_n`` diagonal quantum numbers in `Confs`.
- Merge the submodules to the main package. 
- Add Ising model in X basis.  
- Add functions in built-in models to export diagonal QNs. 

### 27 May, 2024 (Version 0.3.0)

- Add the support for Hamiltonians with real elements. 
- Add conversion with SparseMatrixCSC. 
- Add the conversion from terms to OpSum.
- Add the look-up of configurations. 
- For the built-in Ising model, add density operator.
- Add built-in ``\mathrm{Sp}(N)`` model. 

### 26 May, 2024 (Version 0.2.0)

- Add operations of terms.
- Add built-in Ising model. 