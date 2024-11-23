"""
    mutable struct Operator

An `Operator` object records the sum of terms together with information about its symmetry and the basis of the state it acts on and the basis of the resulting state.

# Fields
* `bsd :: Basis` is the basis of the initial state ;
* `bsf :: Basis` is the basis of the final state ;
* `red_q :: Int64` is a flag that records whether or not the conversion to a sparse martrix can be simplified : if `bsd` and `bsf` have exactly the same set of quantum numbers, and the operator fully respects the symmetries, and all the elements in `bsd.cffac` and `bsf.cffac` has the same absolute value, then `red_q = 1` ; otherwise `red_q = 0` ; 
* `sym_q :: Int64` records the symmetry of the operator : if the matrix is Hermitian, then `sym_q = 1` ; if it is symmetric, then `sym_q = 2` ; otherwise `sym_q = 0` ;
* `ntm :: Int64` is the number of terms ;
* `nc :: Int64` is the maximum number of operators in an operator string
* `cstrs :: Matrix{Int64}` is a matrix recording the operator string of each term. Each column corresponds to a term and is padded to the maximum length with `-1`'s.
* `coeffs :: Vector{ComplexF64}` corresponds to the coefficients in each term.
"""
mutable struct Operator 
    bsd :: Basis 
    bsf :: Basis
    red_q :: Int64
    sym_q :: Int64
    ntm :: Int64
    nc :: Int64
    cstrs :: Matrix{Int64}
    coeffs :: Vector{ComplexF64}
end 


"""
    function Operator(bsd :: Basis[, bsf :: Basis], terms :: Vector{Term} ; red_q :: Int64, sym_q :: Int64, num_th :: Int64, disp_std :: Bool) :: Operator

generates an operator object from a series of terms. 

# Arguments
* `bsd :: Basis` is the basis of the initial state ;
* `bsf :: Basis` is the basis of the final state. Facultative, the same as `bsd` by default. 
* `terms :: Vector{Term}` records the terms ; 
* `red_q :: Int64` is a flag that records whether or not the conversion to a sparse martrix can be simplified : if `bsd` and `bsf` have exactly the same set of quantum numbers, and the operator fully respects the symmetries, and all the elements in `bsd.cffac` and `bsf.cffac` has the same absolute value, then `red_q = 1` ; otherwise `red_q = 0` ; Facultative, if `bsf` is not given, 1 by default, otherwise 0 by default.
* `sym_q :: Int64` records the symmetry of the operator : if the matrix is Hermitian, then `sym_q = 1` ; if it is symmetric, then `sym_q = 2` ; otherwise `sym_q = 0`. Facultative, if `bsf` is not given, 1 by default, otherwise 0 by default.
* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

"""
function Operator(bsd :: Basis, bsf :: Basis, terms :: Vector{Term} ; red_q :: Int64 = 0, sym_q :: Int64 = 0)
    ntm = length(terms)
    nc = div(maximum([length(tm.cstr) for tm in terms]), 2)
    cstrs_vec = [ [tm.cstr ; fill(-1, 2 * nc - length(tm.cstr))] for tm in terms]
    cstrs = reduce(hcat, cstrs_vec)
    coeffs = [ tm.coeff for tm in terms ]
    return Operator(bsd, bsf, red_q, sym_q, ntm, nc, cstrs, coeffs)
end
Operator(bsd :: Basis, terms :: Vector{Term} ; red_q :: Int64 = 1, sym_q :: Int64 = 1) = Operator(bsd, bsd, terms ; red_q, sym_q)


"""
    function *(op :: Operator, st_d :: Vector{ComplexF64} ; num_th :: Int64, disp_std :: Bool) :: Vector{ComplexF64}
    function *(op :: Operator, st_d :: Vector{Float64} ; num_th :: Int64, disp_std :: Bool) :: Vector{Float64}

Measure the action of an operator on a state. `st_d` must be of length `op.bsd.dim`. Returns a vector of length `op.bsf.dim` that represents the final state.

# Facultative arguments

* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

"""
function *(op :: Operator, st_d :: Vector{ComplexF64} ; num_th = NumThreads, disp_std = !SilentStd)
    st_f = Vector{ComplexF64}(undef, op.bsf.dim)
    @ccall Libpath.__op_MOD_action_op(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
    op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
    op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, st_d :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    return st_f
end
function *(op :: Operator, st_d :: Vector{Float64} ; num_th = NumThreads, disp_std = !SilentStd)
    st_f = Vector{ComplexF64}(undef, op.bsf.dim)
    @ccall Libpath.__op_MOD_action_op(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
    op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
    op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, ComplexF64.(st_d) :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing
    return real.(st_f)
end


"""
    function *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, op :: Operator, st_d :: Vector{ComplexF64} ; num_th :: Int64, disp_std :: Bool) :: ComplexF64
    function *(st_fp :: LinearAlgebra.Adjoint{Float64, Vector{Float64}}, op :: Operator, st_d :: Vector{Float64} ; num_th :: Int64, disp_std :: Bool) :: Float64

Measuring the inner product between two states and an operator. `st_d` must be of length `op.bsd.dim` and `st_fp` must be of length `op.bsf.dim`, and `st_fp` must be an adjoint. 

# Facultative arguments

* `num_th :: Int64`, the number of threads. Facultative, `NumThreads` by default. 
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 
"""
function *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, op :: Operator, st_d :: Vector{ComplexF64} ; num_th = NumThreads, disp_std = !SilentStd)
    ovl_ref = Ref{ComplexF64}(0)
    @ccall Libpath.__op_MOD_overlap_op(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
    op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
    op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, ComplexF64.(st_d) :: Ref{ComplexF64}, ComplexF64.(st_fp') :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing 
    return ovl_ref[]
end
function *(st_fp :: LinearAlgebra.Adjoint{Float64, Vector{Float64}}, op :: Operator, st_d :: Vector{Float64} ; num_th = NumThreads, disp_std = !SilentStd)
    ovl_ref = Ref{ComplexF64}(0)
    @ccall Libpath.__op_MOD_overlap_op(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
    op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
    op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstrs :: Ref{Int64}, op.coeffs :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, ComplexF64.(st_d) :: Ref{ComplexF64}, ComplexF64.(st_fp') :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}, num_th :: Ref{Int64}, (disp_std ? 1 : 0) :: Ref{Int64}) :: Nothing 
    return real(ovl_ref[])
end


function HDF5.write(parent :: Union{HDF5.File, HDF5.Group}, name :: String, op :: Operator)
    grp = create_group(parent, name)
    write(grp, "red_q", op.red_q)
    write(grp, "sym_q", op.sym_q)
    write(grp, "bsd", op.bsd)
    if (op.sym_q == 0) write(grp, "bsf", op.bsf) end
    write(grp, "ntm", op.ntm)
    write(grp, "nc", op.nc)
    write(grp, "cstrs", op.cstrs)
    write(grp, "coeffs", op.coeffs)
    close(grp)
end

function HDF5.read(parent :: Union{HDF5.File, HDF5.Group}, name :: String, :: Type{Operator})
    grp = open_group(parent, name)
    red_q = read(grp, "red_q")
    sym_q = read(grp, "sym_q")
    bsd = read(grp, "bsd", Basis)
    bsf = (sym_q == 0) ? read(grp, "bsf", Basis) : bsd
    ntm = read(grp, "ntm")
    nc = read(grp, "nc")
    cstrs = read(grp, "cstrs")
    coeffs = read(grp, "coeffs")
    close(grp)
    return Operator(bsd, bsf, red_q, sym_q, ntm, nc, cstrs, coeffs)
end