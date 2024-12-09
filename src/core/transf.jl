export Transf


"""
    Transf 

The mutable type `Transf` records a transformation in the same form as a [QNOffd](@ref)
```math
𝒵:\\ c_o↦ α_o^* c^{(p_o)}_{π_o},  c_o^†↦ α_o c^{(1-p_o)}_{π_o}
```
together with information about its symmetry and the basis of the state it acts on and the basis of the resulting state.

# Fields
* `bsd :: Basis` is the basis of the initial state.
* `bsf :: Basis` is the basis of the final state.
* `perm :: Vector{Int64}`, `ph :: Vector{Int64}` and `fac :: Vector{ComplexF64}` records the transformation in the same form as a [QNOffd](@ref).
"""
mutable struct Transf 
    bsd :: Basis 
    bsf :: Basis
    perm :: Vector{Int64}
    ph :: Vector{Int64}
    fac :: Vector{ComplexF64}
end


"""
    Transf(bsd :: Basis, bsf :: Basis, qnf :: QNOffd)

generates a transformation object from a QNOffd. 

# Arguments
* `bsd :: Basis` is the basis of the initial state.
* `bsf :: Basis` is the basis of the final state. Facultative, the same as `bsd` by default.
* `qnf :: QNOffd` records the transformation.
"""
function Transf(bsd :: Basis, bsf :: Basis, qnf :: QNOffd)
    return Transf(bsd, bsf, qnf.perm, qnf.ph, qnf.fac)
end
Transf(bsd :: Basis, qnf :: QNOffd) = Transf(bsd, bsd, qnf)


"""
    *(trs :: Transf, st_d :: Vector{ComplexF64} ; num_th = NumThreads) :: Vector{ComplexF64}
    *(trs :: Transf, st_d :: Vector{Float64} ; num_th = NumThreads) :: Vector{Float64}

Act a transformation on a state. `st_d` must be of length `trs.bsd.dim`. Returns a vector of length `trs.bsf.dim` that represents the final state.

# Facultative arguments

* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
"""
function Base.:*(trs :: Transf, st_d :: Vector{ComplexF64} ; num_th = NumThreads)
    st_f = Vector{ComplexF64}(undef, trs.bsf.dim)
    @ccall Libpath.__bs_MOD_action_trs(trs.bsd.cfs.no :: Ref{Int64}, trs.bsd.cfs.nor :: Ref{Int64}, 
    trs.bsd.cfs.ncf :: Ref{Int64}, trs.bsd.dim :: Ref{Int64}, trs.bsd.cfs.conf :: Ref{Int64}, trs.bsd.cfs.lid :: Ref{Int64}, trs.bsd.cfs.rid :: Ref{Int64}, trs.bsd.szz :: Ref{Int64}, trs.bsd.cfgr :: Ref{Int64}, trs.bsd.cffac :: Ref{ComplexF64}, trs.bsd.grel :: Ref{Int64}, trs.bsd.grsz :: Ref{Int64}, 
    trs.bsf.cfs.ncf :: Ref{Int64}, trs.bsf.dim :: Ref{Int64}, trs.bsf.cfs.conf :: Ref{Int64}, trs.bsf.cfs.lid :: Ref{Int64}, trs.bsf.cfs.rid :: Ref{Int64}, trs.bsf.szz :: Ref{Int64}, trs.bsf.cfgr :: Ref{Int64}, trs.bsf.cffac :: Ref{ComplexF64}, trs.bsf.grel :: Ref{Int64}, trs.bsf.grsz :: Ref{Int64}, 
    trs.perm :: Ref{Int64}, trs.ph :: Ref{Int64}, trs.fac :: Ref{ComplexF64}, ComplexF64.(st_d) :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}, num_th :: Ref{Int64}) :: Nothing 
    return st_f
end
function Base.:*(trs :: Transf, st_d :: Vector{Float64} ; num_th = NumThreads)
    st_f = Vector{ComplexF64}(undef, trs.bsf.dim)
    @ccall Libpath.__bs_MOD_action_trs(trs.bsd.cfs.no :: Ref{Int64}, trs.bsd.cfs.nor :: Ref{Int64}, 
    trs.bsd.cfs.ncf :: Ref{Int64}, trs.bsd.dim :: Ref{Int64}, trs.bsd.cfs.conf :: Ref{Int64}, trs.bsd.cfs.lid :: Ref{Int64}, trs.bsd.cfs.rid :: Ref{Int64}, trs.bsd.szz :: Ref{Int64}, trs.bsd.cfgr :: Ref{Int64}, trs.bsd.cffac :: Ref{ComplexF64}, trs.bsd.grel :: Ref{Int64}, trs.bsd.grsz :: Ref{Int64}, 
    trs.bsf.cfs.ncf :: Ref{Int64}, trs.bsf.dim :: Ref{Int64}, trs.bsf.cfs.conf :: Ref{Int64}, trs.bsf.cfs.lid :: Ref{Int64}, trs.bsf.cfs.rid :: Ref{Int64}, trs.bsf.szz :: Ref{Int64}, trs.bsf.cfgr :: Ref{Int64}, trs.bsf.cffac :: Ref{ComplexF64}, trs.bsf.grel :: Ref{Int64}, trs.bsf.grsz :: Ref{Int64}, 
    trs.perm :: Ref{Int64}, trs.ph :: Ref{Int64}, trs.fac :: Ref{ComplexF64}, ComplexF64.(st_d) :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}, num_th :: Ref{Int64}) :: Nothing 
    return real.(st_f)
end
