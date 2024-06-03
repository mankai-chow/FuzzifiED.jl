The package `FuzzifiED` is designed to do exact diagonalisation (ED) calculation on the fuzzy sphere, and also facilitates the DMRG calculations by ITensors. It can also be used for generic fermion models. 

Documentations can be found at [this link](https://mankai-chow.github.io/FuzzifiED/). 

This package is developped by Zheng Zhou (周正) at Perimeter Institute and collaborators. If this package is helpful in your research, we would appreciate it if you mention in the acknowledgement. If you have any questions, please contact at [zzhou@pitp.ca](mailto:zzhou@pitp.ca).

## Install

To install the package, please first enter Julia by entering in the command line `julia`, and then enter the commands
```julia
julia> using Pkg; 
julia> Pkg.add(url="https://github.com/mankai-chow/FuzzifiED_jll.jl.git")
julia> Pkg.add(url="https://github.com/mankai-chow/FuzzifiED.jl.git")
```
Include at the start of your Julia script
```julia
using FuzzifiED
```

## Compile 

If the built-in dynamic library does not work for you, right now you need to compile with the following steps. Better solutions with Julia Binary Builder will be provided in the next version. 

1. Make sure that [Intel oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) and [Intel HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html) is installed, so that you can use `ifort` and OpenMP. 
2. Make sure that Arpack libraries are installed. On linux, the command is 
```bash
sudo apt-get install libarpack++2-dev
```
3. Download the code from the FuzzifiED Fortran code GitHub
```bash
cd $workdir # Change $workdir here to your work directory
git clone --branch=Intel-compiler https://github.com/mankai-chow/FuzzifiED_Fortran
```
4. Compile the fortran code with the ifort compiler/
```bash
cd FuzzifiED_Fortran
ifort -fPIC -larpack -qopenmp -O3 -c cfs.f90
ifort -fPIC -larpack -qopenmp -O3 -c bs.f90
ifort -fPIC -larpack -qopenmp -O3 -c op.f90
ifort -fPIC -larpack -qopenmp -O3 -c diag.f90
ifort -fPIC -larpack -qopenmp -O3 -c diag_re.f90
```
On linux
```bash
ifort -shared -fPIC -larpack -qopenmp -O3 -o lib_fuzzifi_ed.so ./*.o
```
On MacOS
```bash
ifort -dynamiclib -fPIC -larpack -qopenmp -O3 -o lib_fuzzifi_ed.dylib ./*.o
```
Or equivalent on Windows.

5. When using FuzzifiED, declare at the header 
```julia
FuzzifiED.LibFuzzifiED = "$workdir/FuzzifiED_Fortran/lib_fuzzifi_ed.so"
```
change `.so` to `.dylib` or `.dll` if you are not on Linux. 
