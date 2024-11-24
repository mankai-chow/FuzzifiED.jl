"""
    TruncateQnu(; qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}, trunc_lth :: Int64, trunc_wt :: Vector{Int64}) 

**We have improved the interface for this function. Please consider using in the future [`TruncateQNDiag`](@ref)**
```julia 
TruncateQNDiag(qnd ; trunc_lth, trunc_wt)
```

truncates the list of ``N_U`` QNU's from to a number ``N'_U`` acceptable by ITensors. The new quantum numbers are 
```math
\\begin{aligned}
    &Q'_1=Q_1,\\ Q'_2=Q_2,\\ …,\\ Q'_{N'_U-1}=Q_{N'_U-1}\\\\
    &Q'_{N'_U}=λ_{N'_U}Q_{N'_U}+λ_{N'_U+1}Q_{N'_U+1}+…+λ_{N_U}Q_{N_U}
\\end{aligned}
```

# Arguments

- `qnu_o :: Vector{Vector{Int64}}` stores the charge of each orbital under each conserved quantity. See [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)))) for detail. 
- `qnu_name :: Vector{String}` stores the name of each quantum number. Facultative, QN1, QN2, ... by default. 
- `modul :: Vector{Int64}` stores the modulus of each quantum number. Store 1 if no modulus. Facultative, all 1 by default. 
- `trunc_lth :: Int64` stores the truncated numbers of QNU. Facultative, 3 by default. 
- `trunc_wt :: Vecotr{Int64}` stores the ``N_U-N'_U+1`` coefficients ``λ``. Facultative, ``1,2,4,8,…`` by default. 

# Output

A named tuple with three elements that can be directly fed into [`SitesFromQnu`](@ref)

- `qnu_o :: Vector{Vector{Int64}}` stores the charge of each orbital under each conserved quantity. See [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)))) for detail.
- `qnu_name :: Vector{String}` stores the name of each quantum number.
- `modul :: Vector{Int64}` stores the modulus of each quantum number, 1 if no modulus. 

"""
function FuzzifiED.TruncateQnu(; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], trunc_lth :: Int64 = 3, trunc_wt :: Vector{Int64} = [ 2 ^ (i - trunc_lth) for i = trunc_lth : length(qnu_o)]) 
    if (!SilentStd)
        @info """
        We have improved the interface for the function `TruncateQnu`. Please consider using in the future
            TruncateQNDiag(qnd ; trunc_lth, trunc_wt)
        For detail please visit http://docs.fuzzified.world/itensors/#Format-conversion. This function may be superceded in the future version. 
        """
    end
    qnu_o1 = [ qnu_o[1 : trunc_lth - 1] ; [sum([ qnu_o[i + trunc_lth - 1] .* trunc_wt[i] for i in eachindex(trunc_wt)])] ] 
    return (qnu_o = qnu_o1, qnu_name = qnu_name[1 : trunc_lth], modul = modul[1 : trunc_lth])
end

"""
    function SitesFromQnu(; qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64})

**We have improved the interface for the function. Please consider using in the future [`SitesFromQNDiag`](@ref)**
```Julia
SitesFromQNDiag(qnd)
```

returns the ITensors Sites object from the information of quantum numbers 

# Arguments 

- `qnu_o :: Vector{Vector{Int64}}` stores the charge of each orbital under each conserved quantity. See [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)))) for detail. 
- `qnu_name :: Vector{String}` stores the name of each quantum number. Facultative, QN1, QN2, ... by default. 
- `modul :: Vector{Int64}` stores the modulus of each quantum number. Store 1 if no modulus. Facultative, all 1 by default. 
"""
function FuzzifiED.SitesFromQnu(; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)])
    if (!SilentStd)
        @info """
        We have improved the interface for the function `SitesFromQnu`. Please consider using in the future 
            SitesFromQNDiag(qnd)
        For detail please visit http://docs.fuzzified.world/itensors/#Easy-sweep. This function may be superceded in the future version. 
        """
    end
    qnd = [ QNDiag(qnu_name[i], qnu_o[i], modul[i]) for i in eachindex(qnu_o)]
    return [ siteind("Fermion" ; o, qnd) for o in eachindex(qnu_o[1]) ]
end



"""
    function GetMPOSites(id :: String, tms :: Union{Vector{Term}, Sum{Scaled{ComplexF64, Prod{Op}}}} ; path :: String, qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}, mpo_method :: Function) :: Tuple{MPO, Vector{<:Index}}

**We have improved the interface for the function. Please consider using in the future**
```Julia
GetMPOSites(id, tms, qnd)
```

# Function 

This function returns the MPO and sites for a given operator and a Hilbert space with given quantum numbers. It first checks the file `op_\$(id).h5` in a specified repository. If the file exists, it will try to read the fields `mpo` and `sites` and return the MPO and Sites. Otherwise it will first generates the sites with the quantum numbers given in `qnu_o`, `qnu_name` and `modul` (these objects are often results of a function named `Get*Qnu`). Then it will generate the MPO with the terms of the operator given in `tms`. The MPO and sites will be written into the file `op_\$(id).h5` in the fields `mpo` and `sites`. 

# Arguments 

- `id :: String` is a string identifying the file to which the results will be accessed and written.
- `tms :: Vector{Term}` or `tms :: Sum{Scaled{ComplexF64, Prod{Op}}}` is either an array of terms or a `OpSum` objects that specifies the expression of the operator. 
- `path :: String` identifies the path where the results will be accessed and stored. Facultative, `./` by default. 
- `qnu_o :: Vector{Vector{Int64}}` is a list where each element specifies a quantum number. Each element is a list that specifies the charges of each orbital under the quantum number. Obligatory. 
- `qnu_name :: Vector{String}` specifies the name of each quantum number. Facultative, QN1, QN2, ... by default. 
- `module :: Vector{Int64}` specifies the modulus of each quantum number. Facultative, all 1 by default.
- `mpo_method :: Function` is a function `mpo_method(os :: OpSum, sites :: Sites) :: MPO` that generates the MPO from OpSum and Sites. Facultative, `MPO` by default. We suggest using `MPO_new` in `ITensorMPOConstruction` package. See [`example_ising_dmrg_easysweep.jl`](https://github.com/mankai-chow/FuzzifiED.jl/blob/main/examples/example_ising_dmrg_easysweep.jl) for example. 

"""
function FuzzifiED.GetMPOSites(id :: String, tms :: Union{Vector{Term}, Sum{Scaled{ComplexF64, Prod{Op}}}} ; path :: String = "./", qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], mpo_method :: Function = MPO)
    if (!SilentStd)
        @info """
        We have improved the interface for the function `GetMPOSites`. Please consider using in the future 
            GetMPOSites(id, tms, qnd)
        For detail please visit http://docs.fuzzified.world/itensors/#Easy-sweep. This function may be superceded in the future version. 
        """
    end
    if (isfile("$(path)op_$(id).h5"))
        f = h5open("$(path)op_$(id).h5","r")
        mpo = read(f, "mpo", MPO)
        sites = read(f, "sites", Vector{<:Index})
        close(f)
        println("FINISHED READING OPERATOR $(id)")
    else
        sites = SitesFromQnu(; qnu_o, qnu_name, modul)
        os = TermsOrOpSum(tms)

        @time "GENERATE OPERATOR $(id) MPO" mpo = mpo_method(os, sites)
        f = h5open("$(path)op_$(id).h5","cw")
        write(f, "mpo", mpo)
        write(f, "sites", sites)
        close(f)
        println("FINISHED GENERATING OPERATOR $(id), BOND DIM $(maxlinkdim(mpo))")
    end 
    return mpo, sites
end