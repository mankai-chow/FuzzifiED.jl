export SSegOperator, BuildSSegOperator, BuildSSegOperators

mutable struct SSegOperator{T <: Union{Float64, ComplexF64}}
    sgspd :: SSegSpace
    sgspf :: SSegSpace
    colptr :: Vector{Int64}
    rowid :: Vector{Int64}
    elmat :: Vector{Matrix{T}}
end

function BuildSSegOperator(sgspd :: SSegSpace, sgspf :: SSegSpace, amd :: SAngModes, ll :: Int64, secop :: Vector{Int64}, modul :: Vector{Int64} = fill(1, length(sgspd.sec[1])) ; eltype = FuzzifiED.ElementType, num_th = 1)
    index = 0
    colptr = zeros(Int64, length(sgspd.sec) + 1)
    colptr[1] = 1
    rowid = Int64[]
    elmat = Matrix{eltype}[]

    for j in eachindex(sgspd.sec)
        secd = sgspd.sec[j]
        bsd = SBasis(sgspd.cfs[j])
        std = sgspd.sts[j]
        md = secd[2]
        if (md == 0)
            bsd1 = SBasis(sgspd.cfs1[j])
            std1 = sgspd.sts1[j]
        end

        for i in eachindex(sgspf.sec)
            secf = sgspf.sec[i]
            EquivSec(secd .+ secop, secf, modul) || continue
            index += 1
            push!(rowid, i)
            mf = secf[2]
            mm = mf - md

            bsf = SBasis(sgspf.cfs[i])
            stf = sgspf.sts[i]
            tms = GetComponent(amd, ll/2, mm/2)
            op = SOperator(bsd, bsf, tms)
            op_mat = Matrix(OpMat(op ; num_th))
            hmt_block = stf' * op_mat * std

            if (md == 0 && mf == 0 && ll > 0) # when 3j could vanish
                tms1 = GetComponent(amd, ll/2, mm/2 - 1)
                op1 = SOperator(bsd1, bsf, tms1)
                op_mat1 = Matrix(OpMat(op1 ; num_th))
                hmt_block1 = stf' * op_mat1 * std1
            end

            for jl in eachindex(sgspd.l_rng[j])
                ld = sgspd.l_rng[j][jl]
                rngj = (sgspd.ptr_st[j][jl] + 1 : sgspd.ptr_st[j][jl + 1]) .- sgspd.ptr_st[j][1]
                for il in eachindex(sgspf.l_rng[i])
                    lf = sgspf.l_rng[i][il]
                    (ll < abs(ld - lf) || ll > ld + lf) && continue
                    rngi = (sgspf.ptr_st[i][il] + 1 : sgspf.ptr_st[i][il + 1]) .- sgspf.ptr_st[i][1]
                    fac3j = wigner3j(lf/2, ll/2, ld/2, -mf/2, mm/2, md/2)
                    ((lf - mf) % 4 == 2) && (fac3j = -fac3j)
                    if (fac3j ≠ 0)
                        hmt_block[rngi, rngj] /= fac3j
                    else
                        fac3j1 = wigner3j(lf/2, ll/2, ld/2, -mf/2, mm/2 - 1, md/2 + 1)
                        ((lf - mf) % 4 == 2) && (fac3j1 = -fac3j1)
                        hmt_block[rngi, rngj] = hmt_block1[rngi, rngj] / √(ld/2 * (ld/2 + 1)) / fac3j1
                    end
                end
            end
            push!(elmat, hmt_block)
        end
        colptr[j + 1] = index + 1
    end
    return SSegOperator{eltype}(sgspd, sgspf, colptr, rowid, elmat)
end

function BuildSSegOperators(sgspd :: Vector{SSegSpace{T}}, sgspf :: Vector{SSegSpace{T}}, cpd :: Vector{SCoupleDecomp}, modul :: Vector{Int64} = fill(1, length(sgspd[1].sec[1]))) where T <: Union{Float64, ComplexF64}
    id_tc = vcat([ fill(i, length(cpd[i].ch)) for i in eachindex(cpd)]...)
    ch = vcat([ cpd[i].ch for i in eachindex(cpd) ]...)
    coeff = vcat([ cpd[i].coeff for i in eachindex(cpd) ]...)
    nd = length(id_tc)
    np = length(sgspd)
    sgop = Matrix{SSegOperator}(undef, np, nd)
    Threads.@threads :greedy for (p, d) in collect(Iterators.product(1 : np, 1 : nd))
        amd = cpd[id_tc[d]].amd[p]
        secop = cpd[id_tc[d]].sec[:, p]
        ll = ch[d][1, p]
        sgop[p, d] = BuildSSegOperator(sgspd[p], sgspf[p], amd, ll, secop, modul)
    end
    return sgop
end
BuildSSegOperators(sgspd :: Vector{SSegSpace{T}}, cpd :: Vector{SCoupleDecomp}, modul :: Vector{Int64} = fill(1, length(sgspd[1].sec[1]))) where T <: Union{Float64, ComplexF64} = BuildSSegOperators(sgspd, sgspd, cpd, modul)
