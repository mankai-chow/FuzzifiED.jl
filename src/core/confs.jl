"""
    mutable struct Confs

This type stores all the configurations that respects the diagonal quantum numbers (QNDiag) and also a table to inversely look up the index from the configuration. 

# Fields

* `no :: Int64` is the number of orbital\\*flavour.
* `ncf :: Int64` is the number of configurations.
* `conf :: Vector{Int64}` is an array of length `ncf` containing all the configurations. Each configuration is expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th orbital in the `i`-th configuration is occupied ; if the bit is 0, then the orbital is empty. 
* `nor :: Int64`, `lid :: Vector{Int64}` and `rid :: Vector{Int64}` contain the information of Lin table that is used to inversely look up the index `i` from the configuration. 
"""
mutable struct Confs
    no :: Int64 
    nor :: Int64
    ncf :: Int64
    conf :: Vector{Int64} 
    lid :: Vector{Int64}
    rid :: Vector{Int64}
end 


"""
    function Confs(no :: Int64, secd :: Vector{Int64}, qnd :: Vector{QNDiag} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64}, num_th :: Int64, disp_std :: Bool) :: Confs

generates the configurations from the list of QNDiags. 

# Arguments

* `no :: Int64` is the number of orbital\\*flavour ``N_o``.
* `secd :: Vector{Int64}` is the set of ``Q_i`` for the selected configurations in the sector.
* `qnd :: Vector{QNDiag}` is the set of [QNDiags](@ref QNDiag). 
* `nor :: Int64` is the number of less significant bits used to generate the Lin table. Facultative, ``N_o/2`` by default.
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* `cfs :: Confs` is a [Confs](@ref Confs) object.

# Note 

If your `qnd` has negative entries, QNDiags must contain the number of electrons.
"""
function Confs(no :: Int64, secd :: Vector{Int64}, qnd :: Vector{QNDiag} ; nor :: Int64 = div(no, 2), num_th :: Int64 = NumThreads, disp_std :: Bool = !SilentStd)
    nqnd = length(secd)
    lid = Vector{Int64}(undef, 2 ^ (no - nor) + 1)
    ref_ncf = Ref{Int64}(0)
    secd1 = Int64[]
    qnd1 = Vector{Int64}[]
    modul = [ qnd_i.modul for qnd_i in qnd]
    id_ne = 0
    for i = 1 : nqnd 
        if (maximum(qnd[i].charge) == 1 && minimum(qnd[i].charge) == 1)
            id_ne = i
            break
        end 
    end
    for i = 1 : nqnd
        if (modul[i] > 1 || minimum(qnd[i].charge) ≥ 0) 
            push!(secd1, secd[i]) 
            push!(qnd1, qnd[i].charge)
            continue
        end
        qm = minimum(qnd[i].charge)
        push!(secd1, secd[i] .- secd[id_ne] * qm)
        push!(qnd1, qnd[i].charge .- qm)
    end
    qnd1_mat = Matrix{Int64}(reduce(hcat, qnd1))

    @ccall Libpath.__cfs_MOD_count_cfs(
        no :: Ref{Int64}, nor :: Ref{Int64}, 
        nqnd :: Ref{Int64}, secd1 :: Ref{Int64}, 
        qnd1_mat :: Ref{Int64}, modul :: Ref{Int64}, 
        ref_ncf :: Ref{Int64}, lid :: Ref{Int64}, 
        num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
    ) :: Nothing
    ncf = ref_ncf[]
    rid = Vector{Int64}(undef, 2 ^ nor + 1)
    conf = Vector{Int64}(undef, ncf)
    @ccall Libpath.__cfs_MOD_generate_cfs(
        no :: Ref{Int64}, nor :: Ref{Int64}, 
        nqnd :: Ref{Int64}, secd1 :: Ref{Int64}, 
        qnd1_mat :: Ref{Int64}, modul :: Ref{Int64}, 
        ncf :: Ref{Int64}, lid :: Ref{Int64}, rid :: Ref{Int64}, conf :: Ref{Int64}, 
        num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}
    ) :: Nothing
    return Confs(no, nor, ncf, conf, lid, rid)
end

"""
    GetConfId(cfs :: Confs, cf :: Int64) :: Int64

inversely look up the index from the configuration

# Arguments

* `cfs :: Confs` stores the configurations.
* `cf :: Int64` stores the configuration to be looked-up expressed in a binary number. If the `o-1`-th bit of `conf[i]` is 1, then the `o`-th orbital in the `i`-th configuration is occupied ; if the bit is 0, then the orbital is empty. 

# Output
* `id :: Int64` is the id of the configuration such that `cfs.conf[id] == cf`

"""
function GetConfId(cfs :: Confs, cf :: Int64)
    cf_l = cf >> cfs.nor
    cf_r = cf & (1 << cfs.nor - 1)
    return cfs.lid[cf_l + 1] + cfs.rid[cf_r + 1]
end 
