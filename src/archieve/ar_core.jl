"""
    function Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Vector{Int64}} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64}, num_th :: Int64, disp_std :: Bool) :: Confs

**We have improved the interface for this function. Please consider using in the future**
```julia
Confs(no :: Int64, secd :: Vector{Int64}, qnd :: Vector{QNDiag}) :: Confs
```

generates the configurations that has the diagonal quantum numbers given by `qnu_s` of certain conserved quantities specified by `qnu_o :: Vector{Vector{Int64}}`
```math
Q_i=∑_{o=1}^{N_o}q_{io}n_o
```
or
```math
Q_i=∑_{o=1}^{N_o}q_{io}n_o\\ \\mathrm{mod}\\ p_i
```
where ``i=1,…,N_U`` is the index of quantum number, ``o`` is the index of orbital, ``n_o=c^†_oc_o``, and ``q_o`` is a set of coefficients that must be non negative integer valued. (A list of ``q_o`` with both positive and negative entries can be adapted by shifting every elements by a same value)

# Arguments

* `no :: Int64` is the number of orbital\\*flavour ``N_o`` ;
* `qnu_s :: Vector{Int64}` is the set of ``Q_i`` for the selected configurations ;
* `qnu_o :: Vector{Vector{Int64}}` is the set of ``q_{io}`` for each quantum number and for each orbital. It should contain ``N_U`` elements and each element should be a vector of length ``N_o``. 
* `nor :: Int64` is the number of less significant bits used to generate the Lin table. Facultative, ``N_o/2`` by default.
* `modul :: Vector{Int64}` is the modulus of each quantum number. Setting it to 1 means there is no modulus. Facultative, all 1 by default. 
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* `cfs :: Confs` is a [`Confs`](@ref) object
"""
function Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)), num_th = NumThreads, disp_std = !SilentStd)
    if (disp_std)
        @info """
        We have improved the interface for the function `Confs`. Please consider using in the future
            Confs(no :: Int64, secd :: Vector{Int64}, qnd :: Vector{QNDiag}) :: Confs
        For detail please visit http://docs.fuzzified.world/core/#Confs. This function may be superceded in the future version. 
        """
    end
    nqnu = length(qnu_s)
    lid = Array{Int64, 1}(undef, 2 ^ (no - nor) + 1)
    ref_ncf = Ref{Int64}(0)
    qnu_o_mat = Matrix{Int64}(reduce(hcat, qnu_o))
    @ccall Libpath.__cfs_MOD_count_cfs(no :: Ref{Int64}, nor :: Ref{Int64}, nqnu :: Ref{Int64}, qnu_s :: Ref{Int64}, qnu_o_mat :: Ref{Int64}, modul :: Ref{Int64}, ref_ncf :: Ref{Int64}, lid :: Ref{Int64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    ncf = ref_ncf[]
    rid = Array{Int64, 1}(undef, 2 ^ nor + 1)
    conf = Array{Int64, 1}(undef, ncf)
    @ccall Libpath.__cfs_MOD_generate_cfs(no :: Ref{Int64}, nor :: Ref{Int64}, nqnu :: Ref{Int64}, qnu_s :: Ref{Int64}, qnu_o_mat :: Ref{Int64}, modul :: Ref{Int64}, ncf :: Ref{Int64}, lid :: Ref{Int64}, rid :: Ref{Int64}, conf :: Ref{Int64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    return Confs(no, nor, ncf, conf, lid, rid)
end 


"""
    function Basis(cfs :: Confs, qnz_s :: Vector{ComplexF64} ; cyc :: Vector{Int64}, perm_o :: Vector{Vector{Int64}}, ph_o :: Vector{Vector{Int64}}, fac_o :: Vector{Vector{ComplexF64}} ; num_th :: Int64, disp_std :: Bool) :: Basis

**We have improved the interface for this function. Please consider using in the future**
```julia
Basis(cfs :: Confs, secf :: Vector{ComplexF64}, qnf :: Vector{QNOffd}) :: Basis
```
        
generates the basis that respects the off-diagonal ``ℤ_n`` quantum numbers (QNZ) from the diagonal QN–preserving configurations. The discrete ``ℤ_n`` symmetries are in the form of 

```math
𝒵:\\ c_o↦ α_o^* c^{(p_o)}_{π_o},  c_o^†↦ α_o c^{(1-p_o)}_{π_o}
```

where we use a notation ``c^{(1)}=c^†`` and ``c^{0}=c`` for convenience, ``π_o`` is a permutation of ``1,…,N_o``, ``α_o`` is a coefficient, and ``p_o`` specified whether or not particle-hole transformation is performed for the orbital. Note that one must guarentee that all these transformations commute with each other and also commute with the diagonal QNs. 

# Arguments 

* `cfs :: Confs` is the diagonal QN–preserving configurations ;
* `qnz_s :: Vector{ComplexF64}` is a vector of length the same as the number of discrete symmetries that records the eigenvalue of each transformation ;
* `cyc :: Vector{Int64}` records the cycle of each transformation. For ``ℤ_n`` symmetry, record ``n`` ;
* `perm_o :: Vector{Vector{Int64}}` records the permutation ``π_o``. It has ``N_Z`` elements and each of its elements is a vector of length ``N_o``. 
* `ph_o :: Vector{Vector{Int64}}` records ``p_o`` to determine whether or not to perform a particle-hole transformation. It has ``N_Z`` elements and each of its elements is a vector of length ``N_o``. 
* `fac_o :: Vector{Vector{ComplexF64}}` records the factor ``α_o`` in the transformation. Each of its elements is a vector of length ``N_o``. 
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* `bs :: Basis` is the resulting `Basis` object
"""
function Basis(cfs :: Confs, qnz_s :: Vector{ComplexF64} ; cyc :: Vector{Int64}, perm_o :: Vector{Any}, ph_o :: Vector{Any}, fac_o :: Vector{Any}, num_th = NumThreads, disp_std = !SilentStd)
    if (disp_std)
        @info """
        We have improved the interface for the function `Basis`. Please consider using in the future
            Basis(cfs :: Confs, secf :: Vector{ComplexF64}, qnf :: Vector{QNOffd}) :: Basis
        For detail please visit http://docs.fuzzified.world/core/#Basis. This function may be superceded in the future version. 
        """
    end
    if (length(qnz_s) == 0) return Basis(cfs) end
    nqnz = length(qnz_s)
    perm_o_mat = reduce(hcat, perm_o)
    ph_o_mat = reduce(hcat, ph_o)
    fac_o_mat = ComplexF64.(reduce(hcat, fac_o))
    dim_ref = Ref{Int64}(0)
    cfgr = Array{Int64, 1}(undef, cfs.ncf)
    cffac = Array{ComplexF64, 1}(undef, cfs.ncf)
    szz = prod([ abs(qnz_s[i]) < 1E-8 ? 1 : cyc[i] for i = 1 : nqnz ])
    @ccall Libpath.__bs_MOD_generate_bs_cfgr(cfs.no :: Ref{Int64}, cfs.nor :: Ref{Int64}, cfs.ncf :: Ref{Int64}, cfs.lid :: Ref{Int64}, cfs.rid :: Ref{Int64}, cfs.conf :: Ref{Int64}, nqnz :: Ref{Int64}, qnz_s :: Ref{ComplexF64}, cyc :: Ref{Int64}, perm_o_mat :: Ref{Int64}, ph_o_mat :: Ref{Int64}, fac_o_mat :: Ref{ComplexF64}, szz :: Ref{Int64}, dim_ref :: Ref{Int64}, cfgr :: Ref{Int64}, cffac :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    dim = dim_ref[]
    grel = Array{Int64, 2}(undef, szz, dim)
    grsz = Array{Int64, 1}(undef, dim)
    @ccall Libpath.__bs_MOD_generate_bs_grel(cfs.ncf :: Ref{Int64}, szz :: Ref{Int64}, dim :: Ref{Int64}, cfgr :: Ref{Int64}, grel :: Ref{Int64}, grsz :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    return Basis(cfs, dim, szz, cfgr, cffac, grel, grsz)
end 

"""
    function GetEntSpec(st :: Vector{<:Number}, bs0 :: Basis, qnu_s_lst :: Vector{Vector{Vector{Int64}}}, qnz_s_lst :: Vector{Vector{Vector{ComplexF64}}} ; qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}, cyc :: Vector{Int64}, perm_o :: Vector{Vector{Int64}}, ph_o :: Vector{Vector{Int64}}, fac_o :: Vector{Vector{ComplexF64}}, amp_oa :: Vector{<:Number}, amp_ob :: Vector{<:Number} = sqrt.(1 .- abs.(amp_oa .^ 2))) :: Dict{@NamedTuple{qnu_sa, qnz_sa, qnu_sb, qnz_sb}, Vector{Float64}}

**We have improved the interface for this function. Please consider using in the future**
```julia
GetEntSpec(st :: Vector{<:Number}, bs0 :: Basis, secd_lst :: Vector{Vector{Vector{Int64}}}, secf_lst :: Vector{Vector{Vector{<:Number}}} ; qnd_a :: Vector{QNDiag} qnf_a :: Vector{QNOffd}, amp_oa :: Vector{<:Number})
```

# Arguments 

- `st :: Vector{<:Number}` is the state to be decomposed into direct-product basis of two subsystems.
- `bs0 :: Basis` is the total basis. 
- `qnu_s_lst :: Vector{Vector{Vector{Int64}}}` gives the list of QNU sectors of subsystems to be calculated. Each of its elements is a two element vector ; the first specifies the QNUs for subsystem A, and the second specifies the QNU for subsystem B. 
- `qnz_s_lst :: Vector{Vector{Vector{ComplexF64}}}` gives the list of QNZ sectors of subsystems to be calculated. Each of its elements is a two element vector ; the first specifies the QNZs for subsystem A, and the second specifies the QNZs for subsystem B. 
- `qnu_o :: Vector{Vector{Int64}}`, `qnu_name :: Vector{String}` and `modul :: Vector{Int64}` specifies the diagonal quantum numbers of the subsystems A and B. 
- `cyc :: Vector{Int64}`, `perm_o :: Vector{Vector{Int64}}`, `ph_o :: Vector{Vector{Int64}}` and `fac_o :: Vector{Vector{ComplexF64}}` specifies the off-diagonal quantum numbers of the subsystems A and B. 
- `amp_oa :: Vector{ComplexF64}` and `amp_ob :: Vector{ComplexF64}` are complex lists of length `no` that specify the amplitute of each orbital in the subsystems A and B. For a non-local basis, we decompose each electron into creation operators in two subsystems ``c^†_o=a_{o,A}c^†_{o,A}+a_{o,B}c^†_{o,B}`` and this list specifies ``a_{o,A}``. This is equivalent to ``√{ℱ_{m,A}}`` in [PRB 85, 125308 (2012)](https://dx.doi.org/10.1103/PhysRevB.85.125308) with an extra phase factor. 

# Output

A dictionary whose keys are named tuples that specify the sector containing entries `qnu_sa`, `qnz_sq`, `qnu_sb`, `qnz_sb` and values are lists of eigenvalues of the density matrix in those sectors. 

"""
function GetEntSpec(st :: Vector{<:Number}, bs0 :: Basis, qnu_s_lst :: Vector{Any}, qnz_s_lst :: Vector{Any} ; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], cyc :: Vector{Int64}, perm_o :: Vector{Any}, ph_o :: Vector{Any}, fac_o :: Vector{Any}, amp_oa :: Vector{<:Number}, amp_ob :: Vector{<:Number} = sqrt.(1 .- abs.(amp_oa .^ 2)))
    if (disp_std)
        @info """
        We have improved the interface for the function `GetEntSpec`. Please consider using in the future
            GetEntSpec(st :: Vector{<:Number}, bs0 :: Basis, secd_lst :: Vector{Vector{Vector{Int64}}}, secf_lst :: Vector{Vector{Vector{<:Number}}} ; qnd_a :: Vector{QNDiag} qnf_a :: Vector{QNOffd}, amp_oa :: Vector{<:Number})
        For detail please visit http://docs.fuzzified.world/core/#Entanglement. This function may be superceded in the future version. 
        """
    end
    no = bs0.cfs.no 
    nor = bs0.cfs.nor
    dictlock = ReentrantLock()
    EntSpec = Dict{@NamedTuple{qnu_sa :: Vector{<:Number}, qnz_sa :: Vector{<:Number}, qnu_sb :: Vector{<:Number}, qnz_sb :: Vector{<:Number}}, Vector{Float64}}()
    Threads.@threads for qnu_s in qnu_s_lst 
        cfsa = Confs(no, qnu_s[1], qnu_o ; nor, modul, num_th = 1, disp_std = false)
        if (cfsa.ncf == 0) continue end 
        cfsb = Confs(no, qnu_s[2], qnu_o ; nor, modul, num_th = 1, disp_std = false)
        if (cfsb.ncf == 0) continue end 
        for qnz_s in qnz_s_lst 
            bsa = Basis(cfsa, ComplexF64.(qnz_s[1]) ; cyc, perm_o, ph_o, fac_o, num_th = 1, disp_std = false)
            if (bsa.dim == 0) continue end 
            bsb = Basis(cfsb, ComplexF64.(qnz_s[2]) ; cyc, perm_o, ph_o, fac_o, num_th = 1, disp_std = false)
            if (bsb.dim == 0) continue end 
            st_dcp = StateDecompMat(st, bs0, bsa, bsb, amp_oa, amp_ob) 
            ent_spec = abs.(svdvals(st_dcp)) .^ 2
            Threads.lock(dictlock) 
            try
                EntSpec[(qnu_sa = qnu_s[1], qnz_sa = qnz_s[1], qnu_sb = qnu_s[2], qnz_sb = qnz_s[2])] = ent_spec
            finally
                Threads.unlock(dictlock)
            end
        end
    end
    return EntSpec
end