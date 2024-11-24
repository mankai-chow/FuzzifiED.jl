# Other extensions

Apart from ITensors extension, FuzzifiED also provides other extensions, _viz._ HDF5 extension. 

## HDF5 extension 

The HDF5 extension supports writing the types `Confs`, `Basis`, `Vector{Term}`, `Operator`, `OpMat{ComplexF64}` and `OpMat{Float64}` into HDF5 files and reading them from groups and subgroups in HDF5 format. 
```julia
using HDF5 
open(file_name, "cw")
# include the file name as a string 
# Modes : "cw" for write and "r" for read
...
close(f)
```

To write, include in the middle 
```julia
write(f, group_name :: String, cfs :: Confs)
write(f, group_name :: String, bs  :: Basis)
write(f, group_name :: String, tms :: Vector{Term})
write(f, group_name :: String, op  :: Operator)
write(f, group_name :: String, mat :: OpMat{ComplexF64})
write(f, group_name :: String, mat :: OpMat{Float64})
```
To read, include in the middle 
```julia
cfs = read(f, group_name :: String, Confs)
bs  = read(f, group_name :: String, Basis)
tms = read(f, group_name :: String, Vector{Term})
op  = read(f, group_name :: String, Operator)
mat = read(f, group_name :: String, OpMat{ComplexF64})
mat = read(f, group_name :: String, OpMat{Float64})
```
