export TermsDecomp, CompOperator
import FuzzifiEDFullRotation: GetEigensystem

mutable struct TermsDecomp
    amd :: Vector{AngModes}
    ch :: Vector{Matrix{Int64}}
    coeff :: Vector{ComplexF64}
    sec :: Matrix{Int64} # sec[iqn, p]
end

mutable struct CompOperator{T <: Union{Float64, ComplexF64}}
    cpspd :: CompSpace 
    cpspf :: CompSpace
    nd :: Int64
    ltot :: Int64
    coeff :: Vector{ComplexF64}
    sgop :: Matrix{SegOperator}
    mat9j :: Array{Float64, 3}
end

function CompOperator(cpspd :: CompSpace{T}, cpspf :: CompSpace{T}, tmd :: Vector{TermsDecomp}, ltot :: Int64, modul :: Vector{Int64} = fill(1, length(cpspd.sgsp[1].sec[1]))) where T <: Union{Float64, ComplexF64}
    id_tc = vcat([ fill(i, length(tmd[i].ch)) for i in eachindex(tmd)]...)
    ch = vcat([ tmd[i].ch for i in eachindex(tmd) ]...)
    coeff = vcat([ tmd[i].coeff for i in eachindex(tmd) ]...)
    nd = length(id_tc)
    np = cpspd.np

    sgop = Matrix{SegOperator}(undef, np, nd)
    Threads.@threads :greedy for (p, d) in collect(Iterators.product(1 : np, 1 : nd))
        amd = tmd[id_tc[d]].amd[p]
        secop = tmd[id_tc[d]].sec[:, p]
        ll = ch[d][1, p]
        sgop[p, d] = SegOperator(cpspd.sgsp[p], cpspf.sgsp[p], amd, ll, secop, modul)
    end
    
    mat9j = Array{Float64}(undef, cpspf.nch, cpspd.nch, nd)
    lfac = 1 / √(cpspd.ltot + 1) # UNKNOWN SOURCE OF FACTOR !!! 
    Threads.@threads :greedy for (isec, jsec, d) in collect(Iterators.product(eachindex(cpspf.idsec), eachindex(cpspd.idsec), 1 : nd))
        chh = ch[d]
        pfh = mod.(tmd[id_tc[d]].sec[1, :], 2)
        pfi = [ mod(cpspd.sgsp[p].sec[isec][1], 2) for p = 1 : np]
        pftot = sum([pfi[p] * sum(pfh[p + 1 : end]) for p = 1 : np]) % 2

        jst = cpspd.ptr_ch[jsec] - 1
        ist = cpspf.ptr_ch[isec] - 1
        for j in eachindex(cpspd.chs[jsec])
            chj = cpspd.chs[jsec][j]
            for i in eachindex(cpspd.chs[isec])
                chi = cpspf.chs[isec][i]
                fac = lfac
                for p = 2 : np 
                    fac = fac * float(d9j(chj[2,p-1], chj[1,p], chj[2,p],  chh[2,p-1], chh[1,p], chh[2,p],  chi[2,p-1], chi[1,p], chi[2,p])) * √((chj[2,p]+1) * (chh[2,p]+1) * (chi[2,p]+1))
                end
                (pftot == 1) && (fac = -fac)
                mat9j[i + ist, j + jst, d] = fac
            end
        end
    end

    return CompOperator{T}(cpspd, cpspf, nd, ltot, coeff, sgop, mat9j)
end

function Base.:*(cpop :: CompOperator{T}, std :: Vector{T}) where T <: Union{Float64, ComplexF64}
    th_lock = ReentrantLock()
    stf = zeros(T, cpop.cpspf.dim)
    np = cpop.cpspd.np
    Threads.@threads :greedy for (jsec, d) in collect(Iterators.product(eachindex(cpop.cpspd.idsec), 1 : cpop.nd))
        stf1 = zeros(T, cpop.cpspf.dim)
        idsecj = cpop.cpspd.idsec[jsec]
        coeff = cpop.coeff[d]
        idel_rng = [ cpop.sgop[p, d].colptr[idsecj[p]] : cpop.sgop[p, d].colptr[idsecj[p] + 1] - 1 for p = 1 : np ]

        jrng_sg = Vector{UnitRange{Int64}}(undef, np)
        irng_sg = Vector{UnitRange{Int64}}(undef, np)
        for idel in Iterators.product(idel_rng...)
            idseci = [ cpop.sgop[p, d].rowid[idel[p]] for p = 1 : np]
            isec_rng = searchsorted(cpop.cpspf.idsec, idseci)
            isempty(isec_rng) && continue 
            isec = isec_rng[1]
            for jch in eachindex(cpop.cpspd.chs[jsec])
                jrng = cpop.cpspd.ptr_st[jsec][jch] : cpop.cpspd.ptr_st[jsec][jch + 1] - 1
                for p = 1 : np 
                    lj = cpop.cpspd.chs[jsec][jch][1, p]
                    idlj = cpop.cpspd.sgsp[p].l_lookup[idsecj[p]][lj]
                    jrng_sg[p] = (cpop.cpspd.sgsp[p].ptr_st[idsecj[p]][idlj] + 1 : cpop.cpspd.sgsp[p].ptr_st[idsecj[p]][idlj + 1]) .- cpop.cpspd.sgsp[p].ptr_st[idsecj[p]][1]
                end
                for ich in eachindex(cpop.cpspf.chs[isec])
                    fac9j = cpop.mat9j[cpop.cpspf.ptr_ch[isec] - 1 + ich, cpop.cpspf.ptr_ch[jsec] - 1 + jch, d]
                    abs(fac9j) < √eps(Float64) && continue
                    irng = cpop.cpspf.ptr_st[isec][ich] : cpop.cpspf.ptr_st[isec][ich + 1] - 1
                    for p = 1 : np 
                        li = cpop.cpspf.chs[isec][ich][1, p]
                        idli = cpop.cpspf.sgsp[p].l_lookup[idseci[p]][li]
                        irng_sg[p] = (cpop.cpspf.sgsp[p].ptr_st[idseci[p]][idli] + 1 : cpop.cpspf.sgsp[p].ptr_st[idseci[p]][idli + 1]) .- cpop.cpspf.sgsp[p].ptr_st[idseci[p]][1]
                    end
                    @views stf1[irng] .+= coeff * fac9j * ⊗([cpop.sgop[p, d].elmat[idel[p]][irng_sg[p], jrng_sg[p]] for p = 1 : np]...) * std[jrng]
                end
            end
        end
        lock(th_lock) do 
            stf .+= stf1 
        end
    end
    return stf
end

function FuzzifiED.GetEigensystem(cpop :: CompOperator{T}, nst :: Int64 ; tol :: Float64 = 1E-8, ncv :: Int64 = max(2 * nst, nst + 10), initvec = rand(T, cpop.cpspd.dim), disp_std = !FuzzifiED.SilentStd, kwargs...) where T <: Union{ComplexF64,Float64}
    kwargs1 = haskey(kwargs, :krylovdim) ? kwargs : (kwargs..., krylovdim = ncv)
    eigval, eigvec, info = eigsolve(x -> cpop * x, initvec, nst, :SR ; tol, kwargs1...)
    if (disp_std) print(info) end
    return Vector{T}(eigval), Matrix{T}(hcat(eigvec...))
end