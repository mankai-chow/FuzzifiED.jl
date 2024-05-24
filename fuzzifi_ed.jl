using LinearAlgebra
import Base.:*

mutable struct Confs
    no :: Int64 
    nor :: Int64
    ncf :: Int64
    conf :: Array{Int64,1} 
    lid :: Array{Int64,1}
    rid :: Array{Int64,1}
    function Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2))
        # qnu_o :: Vector{Vector{Int64}}
        nqnu = length(qnu_s)
        lid = Array{Int64, 1}(undef, 2 ^ (no - nor) + 1)
        ref_ncf = Ref{Int64}(0)
        qnu_o_mat = reduce(hcat, qnu_o)
        @ccall LibpathFuzzifiED.cfs_mp_count_cfs_(no :: Ref{Int64}, nor :: Ref{Int64}, nqnu :: Ref{Int64}, qnu_s :: Ref{Int64}, qnu_o_mat :: Ref{Int64}, ref_ncf :: Ref{Int64}, lid :: Ref{Int64}) :: Nothing
        ncf = ref_ncf[]
        rid = Array{Int64, 1}(undef, 2 ^ nor + 1)
        conf = Array{Int64, 1}(undef, ncf)
        @ccall LibpathFuzzifiED.cfs_mp_generate_cfs_(no :: Ref{Int64}, nor :: Ref{Int64}, nqnu :: Ref{Int64}, qnu_s :: Ref{Int64}, qnu_o_mat :: Ref{Int64}, ncf :: Ref{Int64}, lid :: Ref{Int64}, rid :: Ref{Int64}, conf :: Ref{Int64}) :: Nothing
        return new(no, nor, ncf, conf, lid, rid)
    end 
end 

mutable struct Basis
    cfs :: Confs
    dim :: Int64
    szz :: Int64
    cfgr :: Vector{Int64}
    cffac :: Vector{ComplexF64}
    grel :: Array{Int64, 2}
    grsz :: Vector{Int64}

    function Basis(cfs :: Confs, qnz_s :: Vector{ComplexF64}, cyc :: Vector{Int64}, perm_o :: Vector{Any}, ph_o :: Vector{Any}, fac_o :: Vector{Any})
        nqnz = length(qnz_s)
        perm_o_mat = reduce(hcat, perm_o)
        ph_o_mat = reduce(hcat, ph_o)
        fac_o_mat = reduce(hcat, fac_o)
        dim_ref = Ref{Int64}(0)
        cfgr = Array{Int64, 1}(undef, cfs.ncf)
        cffac = Array{ComplexF64, 1}(undef, cfs.ncf)
        szz = prod([ abs(qnz_s[i]) < 1E-8 ? 1 : cyc[i] for i = 1 : nqnz ])
        @ccall LibpathFuzzifiED.bs_mp_generate_bs_cfgr_(cfs.no :: Ref{Int64}, cfs.nor :: Ref{Int64}, cfs.ncf :: Ref{Int64}, cfs.lid :: Ref{Int64}, cfs.rid :: Ref{Int64}, cfs.conf :: Ref{Int64}, nqnz :: Ref{Int64}, qnz_s :: Ref{ComplexF64}, cyc :: Ref{Int64}, perm_o_mat :: Ref{Int64}, ph_o_mat :: Ref{Int64}, fac_o_mat :: Ref{ComplexF64}, szz :: Ref{Int64}, dim_ref :: Ref{Int64}, cfgr :: Ref{Int64}, cffac :: Ref{ComplexF64}) :: Nothing
        dim = dim_ref[]
        grel = Array{Int64, 2}(undef, szz, dim)
        grsz = Array{Int64, 1}(undef, dim)
        @ccall LibpathFuzzifiED.bs_mp_generate_bs_grel_(cfs.ncf :: Ref{Int64}, szz :: Ref{Int64}, dim :: Ref{Int64}, cfgr :: Ref{Int64}, grel :: Ref{Int64}, grsz :: Ref{Int64}) :: Nothing
        return new(cfs, dim, szz, cfgr, cffac, grel, grsz)
    end 

    function Basis(cfs :: Confs)
        nqnz = 1
        dim = cfs.ncf
        szz = 1 
        cfgr = collect(1 : dim)
        cffac = fill(ComplexF64(1), dim)
        grel = reshape(collect(1 : dim), 1, :)
        grsz = fill(1, dim)
        return new(cfs, dim, szz, cfgr, cffac, grel, grsz)
    end 
end 

mutable struct Operator 
    bsd :: Basis 
    bsf :: Basis
    red_q :: Int64
    sym_q :: Int64
    ntm :: Int64
    nc :: Int64
    cstr :: Array{Int64, 2}
    fac :: Vector{ComplexF64}
    function Operator(bsd :: Basis, bsf :: Basis, cstr_vec :: Vector{Any}, fac :: Vector{ComplexF64} ; red_q :: Int64 = 0, sym_q :: Int64 = 0)
        ntm = length(cstr_vec)
        nc = div(maximum(length.(cstr_vec)), 2)
        cstr_vec_eq = [ [tm ; fill(-1, 2 * nc - length(tm))] for tm in cstr_vec]
        cstr = reduce(hcat, cstr_vec_eq)
        return new(bsd, bsf, red_q, sym_q, ntm, nc, cstr, fac)
    end
end 

mutable struct OpMat
    dimd :: Int64
    dimf :: Int64
    sym_q :: Int64
    nel :: Int64 
    colptr :: Vector{Int64}
    rowid :: Vector{Int64}
    elval :: Vector{ComplexF64}
    function OpMat(op :: Operator)
        colptr = Array{Int64, 1}(undef, op.bsd.dim + 1)
        nel_ref = Ref{Int64}(0)
        @ccall LibpathFuzzifiED.op_mp_count_op_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
            op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
            op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
            op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstr :: Ref{Int64}, op.fac :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}) :: Nothing
        nel = nel_ref[]
        rowid = Array{Int64, 1}(undef, nel)
        elval = Array{ComplexF64, 1}(undef, nel)
        @ccall LibpathFuzzifiED.op_mp_generate_op_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
            op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
            op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
            op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstr :: Ref{Int64}, op.fac :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, op.sym_q :: Ref{Int64}, nel_ref :: Ref{Int64}, colptr :: Ref{Int64}, rowid :: Ref{Int64}, elval :: Ref{ComplexF64}) :: Nothing
        return new(op.bsd.dim, op.bsf.dim, op.sym_q, nel, colptr, rowid, elval)
    end
end 

function GetEigensystem(mat :: OpMat, nst :: Int64 ; tol :: Float64 = 1E-8)
    eigval = Array{ComplexF64, 1}(undef, nst + 1)
    eigvec = Array{ComplexF64, 2}(undef, mat.dimd, nst)
    @ccall LibpathFuzzifiED.diag_mp_diagonalisation_(mat.dimd :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, nst :: Ref{Int64}, tol :: Ref{Float64}, eigval :: Ref{ComplexF64}, eigvec :: Ref{ComplexF64}) :: Nothing
    return eigval[1 : end - 1], eigvec
end 

function *(op :: Operator, st_d :: Vector{ComplexF64})
    st_f = Array{ComplexF64, 1}(undef, op.bsf.dim)
    @ccall LibpathFuzzifiED.op_mp_action_op_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
    op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
    op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstr :: Ref{Int64}, op.fac :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, st_d :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}) :: Nothing
    return st_f
end

function *(mat :: OpMat, st_d :: Vector{ComplexF64})
    st_f = Array{ComplexF64, 1}(undef, mat.dimf)
    @ccall LibpathFuzzifiED.diag_mp_vec_prod_(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, st_d :: Ref{ComplexF64}, st_f :: Ref{ComplexF64}) :: Nothing
    return st_f
end

function *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, op :: Operator, st_d :: Vector{ComplexF64})
    ovl_ref = Ref{ComplexF64}(0)
    @ccall LibpathFuzzifiED.op_mp_overlap_op_(op.bsd.cfs.no :: Ref{Int64}, op.bsd.cfs.nor :: Ref{Int64}, 
    op.bsd.cfs.ncf :: Ref{Int64}, op.bsd.dim :: Ref{Int64}, op.bsd.cfs.conf :: Ref{Int64}, op.bsd.cfs.lid :: Ref{Int64}, op.bsd.cfs.rid :: Ref{Int64}, op.bsd.szz :: Ref{Int64}, op.bsd.cfgr :: Ref{Int64}, op.bsd.cffac :: Ref{ComplexF64}, op.bsd.grel :: Ref{Int64}, op.bsd.grsz :: Ref{Int64}, 
    op.bsf.cfs.ncf :: Ref{Int64}, op.bsf.dim :: Ref{Int64}, op.bsf.cfs.conf :: Ref{Int64}, op.bsf.cfs.lid :: Ref{Int64}, op.bsf.cfs.rid :: Ref{Int64}, op.bsf.szz :: Ref{Int64}, op.bsf.cfgr :: Ref{Int64}, op.bsf.cffac :: Ref{ComplexF64}, op.bsf.grel :: Ref{Int64}, op.bsf.grsz :: Ref{Int64}, 
    op.ntm :: Ref{Int64}, op.nc :: Ref{Int64}, op.cstr :: Ref{Int64}, op.fac :: Ref{ComplexF64}, op.red_q :: Ref{Int64}, st_d :: Ref{ComplexF64}, st_fp' :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}) :: Nothing 
    return ovl_ref[]
end

function *(st_fp :: LinearAlgebra.Adjoint{ComplexF64, Vector{ComplexF64}}, mat :: OpMat, st_d :: Vector{ComplexF64})
    ovl_ref = Ref{ComplexF64}(0)
    @ccall LibpathFuzzifiED.diag_mp_scal_prod_(mat.dimd :: Ref{Int64}, mat.dimf :: Ref{Int64}, mat.sym_q :: Ref{Int64}, mat.nel :: Ref{Int64}, mat.colptr :: Ref{Int64}, mat.rowid :: Ref{Int64}, mat.elval :: Ref{ComplexF64}, st_d :: Ref{ComplexF64}, st_fp' :: Ref{ComplexF64}, ovl_ref :: Ref{ComplexF64}) :: Nothing
    return ovl_ref[]
end