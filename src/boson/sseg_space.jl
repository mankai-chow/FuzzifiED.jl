export SSegSpace, BuildSSegSpace

mutable struct SSegSpace{T <: Union{Float64, ComplexF64}}
    sec :: Vector{Vector{Int64}}
    l_rng :: Vector{Vector{Int64}}
    l_lookup :: Vector{Dict{Int64, Int64}}
    ptr_st :: Vector{Vector{Int64}}
    cfs :: Vector{SConfs}
    cfs1 :: Vector{SConfs}
    sts :: Vector{Matrix{T}}
    sts1 :: Vector{Matrix{T}}
end

function BuildSSegSpace(nof :: Int64, nob :: Int64, nebm :: Vector{Int64}, sec :: Vector{Vector{Int64}}, qnd :: Vector{SQNDiag}, tms_lzlp :: Tuple{STerms, STerms}, tms_c2 :: STerms = 0 * one(STerms), c2_rng :: Vector{Float64} = [0.0] ; eltype = FuzzifiED.ElementType, num_th = FuzzifiED.NumThreads)
    nsec = length(sec)
    cfs = Vector{SConfs}(undef, nsec)
    cfs1 = Vector{SConfs}(undef, nsec)
    sts = Vector{Matrix{eltype}}(undef, nsec)
    sts1 = Vector{Matrix{eltype}}(undef, nsec)
    l_rng = [ Int64[] for _ ∈ sec ]
    ptr_st = [ Int64[] for _ ∈ sec ]
    l_lookup = [ Dict{Int64, Int64}() for _ ∈ sec ]
    tms_l2 = GetL2STerms(tms_lzlp)
    BLAS.set_num_threads(num_th)
    for isec ∈ eachindex(sec)
        seci = sec[isec]
        cfs[isec] = SConfs(nof, nob, nebm[isec], seci, qnd ; num_th)
        bs = SBasis(cfs[isec])

        l2c2_mat = Matrix(OpMat(SOperator(bs, √2 * tms_l2 + tms_c2) ; num_th))
        l2c2_val, st = eigen(Hermitian(l2c2_mat))

        l2_mat = OpMat(SOperator(bs, tms_l2) ; num_th)
        l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in axes(st, 2)]
        c2_val = l2c2_val .- √2 .* l2_val
        l_val = round.(Int64, sqrt.(real.(4 * l2_val) .+ 1) .- 1)

        ls = sort(unique(l_val))
        index = 1
        i_rng = Int64[]
        for l in ls
            for i = 1 : bs.dim
                (l_val[i] == l) || continue
                flag = false
                for c2 in c2_rng
                    abs(c2_val[i] - c2) > 1E-4 && continue
                    flag = true
                    break
                end
                flag || continue
                if (isempty(l_rng[isec]) || l_rng[isec][end] ≠ l)
                    push!(l_rng[isec], l)
                    push!(ptr_st[isec], index)
                end
                push!(i_rng, i)
                index += 1
            end
        end
        for il in eachindex(l_rng[isec])
            l_lookup[isec][l_rng[isec][il]] = il
        end
        push!(ptr_st[isec], index)
        sts[isec] = st[:, i_rng]

        if (seci[2] == 0)
            seci1 = deepcopy(seci)
            seci1[2] = 2
            cfs1[isec] = SConfs(nof, nob, nebm[isec], seci1, qnd ; num_th)
            bs1 = SBasis(cfs1[isec])
            lp = SOperator(bs, bs1, tms_lzlp[2])
            lp_mat = Matrix(OpMat(lp))
            sts1[isec] = lp_mat * sts[isec]
        end
        #println("SECTOR $(seci), TOTAL DIMENSION $(bs.dim), SELECTED DIMENSION $(index - 1).")
    end
    BLAS.set_num_threads(1)
    ptr_sec = cumsum([ptr_st[isec][end] - 1 for isec ∈ eachindex(sec)])
    for isec = 2 : length(sec)
        ptr_st[isec] .+= ptr_sec[isec - 1]
    end
    return SSegSpace{eltype}(sec, l_rng, l_lookup, ptr_st, cfs, cfs1, sts, sts1)
end
