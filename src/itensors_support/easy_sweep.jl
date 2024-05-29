mutable struct EasySweepObserver <: AbstractObserver
    e_tol :: Float64
    e_last :: Float64
    EasySweepObserver(e_tol = 0.0) = new(e_tol, 1000.0)
end

function ITensors.checkdone!(o :: EasySweepObserver; kwargs...)
    sw = kwargs[:sweep]
    energy = kwargs[:energy]
    if abs(energy - o.e_last) < o.e_tol
        println("Stopping DMRG after sweep $sw")
        return true
    end
    # Otherwise, update last_energy and keep going
    o.e_last = energy
    return false
end


"""
    function SweepOne(id :: String, hmt :: MPO, st0 :: MPS, dim1 :: Int64 ; path :: String, cutoff :: Vector{Float64}, maxdim :: Vector{Int64}, nsweeps :: Int64, noise :: Vector{Float64}, proj :: Vector{String}, e_tol :: Float64, weight :: Float64, observer :: AbstractObserver) :: Tuple{Float64, MPS}

# Function 

This function performs one round of `nsweeps` sweeps. It first checks the file `st_\$(id).h5` in a specified repository. If the key `st_d\$(dim1)` exists, it reads the energy and MPS from the file and return the energy and MPS, otherwise it will perform the DMRG process with the maximal bond dimension specified by `maxdim` if it exists, or `dim1`. The projected states will be read from the key `st_d\$(dim1)` if it exists or `st_fin` in the file `st_\$(fi).h5` in the same repository for each string `fi` in the array `proj`. The sweeps will be ended if the energy difference is less than `etol` or whatever criteria is given in `observer`. The resulting energy will be written into the key `E_d\$(dim1)` in the file `st_\$(id).h5`, and the MPS written into `st_d\$(dim1)`. The function returns a tuple of energy and the final MPS. 

# Arguments 

- `id :: String` is a string identifying the file to which the results will be accessed and written.
- `hmt :: MPO` is an MPO specifying the Hamiltonian.
- `st0 :: MPS` is an MPS specifying the initial state. 
- `dim1 :: Int64` is a bond dimension that will be used to identify the result. 
- `path :: String` identifies the path where the results will be accessed and stored. Facultive, `./` by default. 
- `cutoff :: Vector{Float64}` is the cutoff that will be sent into DMRG. Facultive, `[1.E-9]` by default. 
- `maxdim :: Vector{Int64}` specifies the maximal bond dimension of each sweep. Facultive, `[dim1]` by default. 
- `nsweeps :: Int64` specifies the number of sweeps in the round. Facultive, 10 by default. 
- `noise :: Vector{Float64}` specifies the noise of each sweep and will be sent into DMRG. Facultive, `[1E-6,1E-7,0]` by default. 
- `proj :: Vector{String}` specifies the name of the states that will be projected. Facultive, empty by default. 
- `e_tol :: Float64` specifies the energy tolerence as a criteria to end the sweeps. Facultive, `1E-6` by default. 
- `weight :: Float64` specifies the weight of projected states and will be sent into DMRG. Facultive, 100.0 by default.
- `observer :: AbstractObserver` specifies the measurement and cutoff condition for each sweep. Facultive, by default the observer will print the energy and cutoff once the energy difference is less than `e_tol` at each sweep. 

"""
function SweepOne(id :: String, hmt :: MPO, st0 :: MPS, dim1 :: Int64 ; path :: String = "./", cutoff :: Vector{Float64} = [1E-9], maxdim :: Vector{Int64} = [dim1], nsweeps :: Int64 = 10, noise :: Vector{Float64} = [1E-6,1E-7,0], proj :: Vector{String} = String[], e_tol :: Float64 = 1E-6, weight :: Float64 = 100.0, observer :: AbstractObserver = EasySweepObserver(e_tol))
    if (isfile("$(path)st_$(id).h5"))
        f = h5open("$(path)st_$(id).h5","r")
        if (haskey(f, "st_d$(dim1)"))
            st1 = read(f, "st_d$(dim1)", MPS)
            E1 = read(f, "E_d$(dim1)")
            close(f)
            println("FINISHED READING STATE $id, BOND DIM $(maxlinkdim(st1)), ENERGY $E1")
            return E1, st1
        end
        close(f)
    end

    if (isempty(proj))
        E1, st1 = dmrg(hmt, st0 ; nsweeps, maxdim, cutoff, noise, observer, outputlevel = 1)  # ground state
    else
        fs = [ h5open("$(path)st_$(fi).h5", "r") for fi in proj ]
        grs = [ "st_d$(dim1)" for fi in proj ]
        for i = 1 : length(proj)
            if (haskey(fs[i], grs[i])) continue end
            grs[i] = "st_fin"
        end
        sts = [ read(fs[i], grs[i], MPS) for i = 1 : length(proj) ]
        for fi in fs
            close(fi) 
        end 
        # strategy for reading excited states : 
        # first try to read the same bond dimension
        # if not exist, read the final
        E1, st1 = dmrg(hmt, sts, st0 ; nsweeps, maxdim, cutoff, noise, observer, outputlevel = 1, weight)  # ground state
    end

    f = h5open("$(path)st_$(id).h5","cw")
    write(f, "st_d$(dim1)", st1)
    write(f, "E_d$(dim1)", E1)
    close(f)
    println("FINISHED GENERATING STATE $id, BOND DIM $(maxlinkdim(st1)), ENERGY $E1")
    return E1, st1
end


"""
    function EasySweep(id :: String, hmt :: MPO, st00 :: MPS ; path :: String, dim_list :: Vector{Int64}, proj :: Vector{String}, e_tol1 :: Float64, e_tol :: Float64, cutoff :: Vector{Float64}, maxdim0 :: Vector{Float64}, noise0 :: Vector{Float64}, noise :: Vector{Int64}, nsweeps :: Int64, weight :: Float64, observer :: AbstractObserver) :: Tuple{Float64, MPS}

# Function 

This function automatically performs several rounds of DMRG sweeps with increasing bond dimensions. It first checks the file `st_\$(id).h5` in a specified repository. If the key `st_fin` exists, it reads the energy and MPS from the file and return the energy and MPS, otherwise it will perform the DMRG process. For each round, it will try to access the results from the key `st_d\$(dim_i)` in `st_\$(id).h5`, where `dim_i` is either `0` representing the initial round, or an element of array `dim_list`. If the key exist, it will read the result ; otherwise it will perform the sweeps using [`SweepOne`](@ref). For the initial round, it will take the initial state from `st00`, the maximal bond dimensions from `maxdim0`, noise from `noise0` and record the results in the key `E_d0` and `st_d0` in `st_\$(id).h5`. For each of the following round, it will take the result from the previous round as the initial state and perform `nsweeps` sweeps with the bond dimension `dim_list[i]`. Each round will be stopped if the energy difference is less than `e_tol1`. The entire process will be stopped if the energy difference between two rounds is less than `e_tol` or the bond dimension of the result is less than 0.9 times the maximal bond dimension. The projected states will be accessed from the files specified by `proj`. It will try to access first the states with the same bond dimension as the projected states. If such states do not exist, it will then access the final state. The resulting energy will be written into the key `E_fin` in the file `st_\$(id).h5`, and the MPS written into `st_fin`. The function returns a tuple of energy and the final MPS. 

# Arguments 

- `id :: String` is a string identifying the file to which the results will be accessed and written.
- `hmt :: MPO` is an MPO specifying the Hamiltonian.
- `st00 :: MPS` is an MPS specifying the initial state. 
- `path :: String` identifies the path where the results will be accessed and stored. Facultive, `./` by default. 
- `dim_list :: Vector{Int64} :: Int64` is a list that specifies the maximal bond dimensions of each round of sweeps starting from the second round. Facultive, `[1000,2000,3000,4000,5000,6000]` by default
- `proj :: Vector{String}` specifies the name of the states that will be projected. Facultive, empty by default. 
- `e_tol1 :: Float64` specifies the energy tolerence as a criteria to end the round of sweeps for each round of sweeps. Facultive, `1E-6` by default. 
- `e_tol :: Float64` specifies the energy tolerence as a criteria to end the entire process. Facultive, `1E-7` by default. 
- `cutoff :: Vector{Float64}` is the cutoff that will be sent into DMRG. Facultive, `[1.E-9]` by default. 
- `maxdim0 :: Vector{Int64}` specifies the maximal bond dimensions of the first round of sweeps. Facultive, `[1E-4,3E-5,1E-5,3E-6,1E-6,3E-7]` by default. 
- `noise0 :: Vector{Float64}` specifies the noise of each sweep in the initial round and will be sent into DMRG. Facultive, `[1E-4,3E-5,1E-5,3E-6,1E-6,3E-7]` by default. 
- `noise :: Vector{Float64}` specifies the noise of each sweep from the second round and will be sent into DMRG. Facultive, `[1E-6,1E-7,0]` by default. 
- `nsweeps :: Int64` specifies the number of sweeps in each round from the second rounds. Facultive, 10 by default. 
- `weight :: Float64` specifies the weight of projected states and will be sent into DMRG. Facultive, 100.0 by default.
- `observer :: AbstractObserver` specifies the measurement and cutoff condition for each sweep starting from the second round. Facultive, by default the observer will print the energy and cutoff once the energy difference is less than `e_tol` at each sweep. 

"""
function EasySweep(id :: String, hmt :: MPO, st00 :: MPS ; path :: String = "./", dim_list :: Vector{Int64} = [1000,2000,3000,4000,5000,6000], proj :: Vector{String} = String[], e_tol1 :: Float64 = 1E-6, e_tol :: Float64 = 1E-7, cutoff :: Vector{Float64} = [1E-9], maxdim0 :: Vector{Int64} = [10,20,50,100,200,500], noise0 :: Vector{Float64} = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], noise :: Vector{Float64} = [1E-6,2E-7,5E-8,1E-8,0], nsweeps :: Int64 = 10, weight :: Float64 = 100.0, observer :: AbstractObserver = EasySweepObserver(e_tol1))
    if (isfile("$(path)st_$(id).h5"))
        f = h5open("$(path)st_$(id).h5","r")
        if (haskey(f, "st_fin"))
            st1 = read(f, "st_fin", MPS)
            E1 = read(f, "E_fin")
            close(f)
            println("FINISHED READING STATE $id, BOND DIM $(maxlinkdim(st1)), ENERGY $E1")
            return E1, st1
        end
        close(f)
    end

    E0, st0 = SweepOne(id, hmt, st00, 0 ; path, proj, cutoff, maxdim = maxdim0, nsweeps = length(maxdim0), noise = noise0, e_tol = 0., weight)
    for dim_i in dim_list
        E1, st1 = SweepOne(id, hmt, st0, dim_i ; path, proj, e_tol = e_tol1, cutoff, noise, nsweeps, weight, observer)
        if (abs(E1 - E0) < e_tol || maxlinkdim(st1) < .9 * dim_i) break end
        E0 = E1 
        st0 = st1
    end 

    f = h5open("$(path)st_$(id).h5","cw")
    write(f, "st_fin", st1)
    write(f, "E_fin", E1)
    close(f)

    println("FINISHED GENERATING STATE $id, FINAL, BOND DIM $(maxlinkdim(st1)), ENERGY $E1")
    return E1, st1
end

TermsOrOpSum(tms :: Vector{Term}) = OpSumFromTerms(tms)
TermsOrOpSum(os :: Sum{Scaled{ComplexF64, Prod{Op}}}) = os

"""
    function GetMpoSites(id :: String, tms :: Union{Vector{Term}, Sum{Scaled{ComplexF64, Prod{Op}}}} ; path :: String, qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}, old = false) :: Tuple{MPO, Vector{<:Index}}

# Function 

This function returns the MPO and sites for a given operator and a Hilbert space with given quantum numbers. It first checks the file `op_\$(id).h5` in a specified repository. If the file exists, it will try to read the fields `mpo` and `sites` and return the MPO and Sites. Otherwise it will first generates the sites with the quantum numbers given in `qnu_o`, `qnu_name` and `modul` (these objects are often results of a function named `Get*Qnu`). Then it will generate the MPO with the terms of the operator given in `tms`. The MPO and sites will be written into the file `op_\$(id).h5` in the fields `mpo` and `sites`. 

# Arguments 

- `id :: String` is a string identifying the file to which the results will be accessed and written.
- `tms :: Vector{Term}` or `tms :: Sum{Scaled{ComplexF64, Prod{Op}}}` is either an array of terms or a `OpSum` objects that specifies the expression of the operator. 
- `path :: String` identifies the path where the results will be accessed and stored. Facultive, `./` by default. 
- `qnu_o :: Vector{Vector{Int64}}` is a list where each element specifies a quantum number. Each element is a list that specifies the charges of each orbital under the quantum number. Obligatory. 
- `qnu_name :: Vector{String}` specifies the name of each quantum number. Facultive, QN1, QN2, ... by default. 
- `module :: Vector{Int64}` specifies the modulus of each quantum number. Facultive, all 1 by default.
- `old :: Bool` specifies the method of MPO generation. If true, then the `MPO()` method in ITensors will be used ; if false, then `MPO_new()` method in `ITensorMPOConstruction` package will be used. Facultive, `false` by default. 

"""
function GetMpoSites(id :: String, tms :: Union{Vector{Term}, Sum{Scaled{ComplexF64, Prod{Op}}}} ; path :: String = "./", qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], old :: Bool = false)
    if (isfile("$(path)op_$(id).h5"))
        f = h5open("$(path)op_$(id).h5","r")
        mpo = read(f, "mpo", MPO)
        sites = read(f, "sites", Vector{<:Index})
        close(f)
        println("FINISHED READING OPERATOR $(id)")
    else
        sites = SitesFromQnu(; qnu_o, qnu_name, modul)
        os = TermsOrOpSum(tms)

        if (old)
            @time "GENERATE OPERATOR $(id) MPO" mpo = MPO(os, sites)
        else
            operatorNames = [ "I", "C", "Cdag", "N" ]
            opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
            @time "GENERATE OPERATOR $(id) MPO" mpo = MPO_new(os, sites ; basisOpCacheVec = opCacheVec)
        end

        f = h5open("$(path)op_$(id).h5","cw")
        write(f, "mpo", mpo)
        write(f, "sites", sites)
        close(f)
        println("FINISHED GENERATING OPERATOR $(id), BOND DIM $(maxlinkdim(mpo))")
    end 
    return mpo, sites
end


"""
    function GetMpo(id :: String, tms :: Union{Vector{Term}, Sum{Scaled{ComplexF64, Prod{Op}}}} ; path :: String, qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}, old = false) :: Tuple{MPO, Vector{<:Index}}

# Function 

This function returns the MPO for a given operator and a given set of sites. It first checks the file `op_\$(id).h5` in a specified repository. If the file exists, it will try to read the fields `mpo` and return the MPO it has read. Otherwise it will generate the MPO with the terms of the operator given in `tms`. The MPO will be written into the file `op_\$(id).h5` in the fields `mpo`. 

# Arguments 

- `id :: String` is a string identifying the file to which the results will be accessed and written.
- `tms :: Vector{Term}` or `tms :: Sum{Scaled{ComplexF64, Prod{Op}}}` is either an array of terms or a `OpSum` objects that specifies the expression of the operator. 
- `sites :: Vector{<:Index}` specifies the sites that the operator is acting on. 
- `path :: String` identifies the path where the results will be accessed and stored. Facultive, `./` by default. 
- `old :: Bool` specifies the method of MPO generation. If true, then the `MPO()` method in ITensors will be used ; if false, then `MPO_new()` method in `ITensorMPOConstruction` package will be used. Facultive, `false` by default. 

"""
function GetMpo(id :: String, tms :: Union{Vector{Term}, Sum{Scaled{ComplexF64, Prod{Op}}}}, sites :: Vector{<:Index} ; path :: String = "./", old :: Bool = false)
    if (isfile("$(path)op_$(id).h5"))
        f = h5open("$(path)op_$(id).h5","r")
        mpo = read(f, "mpo", MPO)
        close(f)
        println("FINISHED READING OPERATOR $(id)")
    else
        os = TermsOrOpSum(tms)

        if (old)
            @time "GENERATE OPERATOR $(id) MPO" mpo = MPO(os, sites)
        else
            operatorNames = [ "I", "C", "Cdag", "N" ]
            opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
            @show os
            @time "GENERATE OPERATOR $(id) MPO" mpo = MPO_new(os, sites; basisOpCacheVec = opCacheVec)
        end

        f = h5open("$(path)op_$(id).h5","cw")
        write(f, "mpo", mpo)
        close(f)
        println("FINISHED GENERATING OPERATOR $(id), BOND DIM $(maxlinkdim(mpo))")
    end 
    return mpo
end