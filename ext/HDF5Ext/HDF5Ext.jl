module HDF5Ext 

using HDF5
using FuzzifiED

function HDF5.write(parent :: Union{HDF5.File, HDF5.Group}, name :: String, cfs :: Confs)
    grp = create_group(parent, name)
    write(grp, "no", cfs.no)
    write(grp, "nor", cfs.nor)
    write(grp, "ncf", cfs.ncf)
    write(grp, "conf", cfs.conf)
    write(grp, "lid", cfs.lid)
    write(grp, "rid", cfs.rid)
    close(grp)
end

function HDF5.read(parent :: Union{HDF5.File, HDF5.Group}, name :: String, :: Type{Confs})
    grp = open_group(parent, name)
    no = read(grp, "no")
    nor = read(grp, "nor")
    ncf = read(grp, "ncf")
    conf = read(grp, "conf")
    lid = read(grp, "lid")
    rid = read(grp, "rid")
    close(grp)
    return Confs(no, nor, ncf, conf, lid, rid)
end

function HDF5.write(parent :: Union{HDF5.File, HDF5.Group}, name :: String, bs :: Basis)
    grp = create_group(parent, name)
    write(grp, "cfs", bs.cfs)
    write(grp, "dim", bs.dim)
    write(grp, "szz", bs.szz)
    write(grp, "cfgr", bs.cfgr)
    write(grp, "cffac", bs.cffac)
    write(grp, "grel", bs.grel)
    write(grp, "grsz", bs.grsz)
    close(grp)
end

function HDF5.read(parent :: Union{HDF5.File, HDF5.Group}, name :: String, :: Type{Basis})
    grp = open_group(parent, name)
    cfs = read(grp, "cfs", Confs)
    dim = read(grp, "dim")
    szz = read(grp, "szz")
    cfgr = read(grp, "cfgr")
    cffac = read(grp, "cffac")
    grel = read(grp, "grel")
    grsz = read(grp, "grsz")
    close(grp)
    return Basis(cfs, dim, szz, cfgr, cffac, grel, grsz)
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

function HDF5.write(parent :: Union{HDF5.File, HDF5.Group}, name :: String, tms :: Terms)
    grp = create_group(parent, name)
    coeffs = [ tm.coeff for tm in tms ]
    nc = maximum([length(tm.cstr) for tm in tms]) ÷ 2
    cstrs = hcat([ [tm.cstr ; fill(-1, 2 * nc - length(tm.cstr))] for tm in tms]...)
    write(grp, "coeffs", coeffs)
    write(grp, "cstrs", cstrs)
    close(grp)
end

function HDF5.read(parent :: Union{HDF5.File, HDF5.Group}, name :: String, :: Type{Terms})
    grp = open_group(parent, name)
    coeffs = read(grp, "coeffs")
    cstrs = read(grp, "cstrs")
    close(grp)
    tms = [ Term(coeffs[i], cstrs[:, i]) for i in eachindex(coeffs) ]
    return tms
end

function HDF5.write(parent :: Union{HDF5.File, HDF5.Group}, name :: String, mat :: OpMat{T}) where T <: Union{Float64, ComplexF64}
    grp = create_group(parent, name)
    write(grp, "dimd", mat.dimd)
    write(grp, "dimf", mat.dimf)
    write(grp, "sym_q", mat.sym_q)
    write(grp, "nel", mat.nel)
    write(grp, "colptr", mat.colptr)
    write(grp, "rowid", mat.rowid)
    write(grp, "elval", mat.elval)
    close(grp)
end

function HDF5.read(parent :: Union{HDF5.File, HDF5.Group}, name :: String, :: Type{OpMat{T}}) where T <: Union{Float64, ComplexF64}
    grp = open_group(parent, name)
    dimd = read(grp, "dimd")
    dimf = read(grp, "dimf")
    sym_q = read(grp, "sym_q")
    nel = read(grp, "nel")
    colptr = read(grp, "colptr")
    rowid = read(grp, "rowid")
    elval = read(grp, "elval")
    close(grp)
    return OpMat{T}(dimd, dimf, sym_q, nel, colptr, rowid, elval)
end

end 