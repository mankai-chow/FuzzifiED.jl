export SegSpace

mutable struct SegSpace{T <: Union{Float64, ComplexF64}}
    sec :: Vector{Vector{Int64}}
    l_rng :: Vector{Vector{Int64}}
    l_lookup :: Vector{Dict{Int64, Int64}}
    ptr_st :: Vector{Vector{Int64}}
    cfs :: Vector{Confs}
    cfs1 :: Vector{Confs}
    sts :: Vector{Matrix{T}}
    sts1 :: Vector{Matrix{T}}
end

function SegSpace(no :: Int64, sec :: Vector{Vector{Int64}}, qnd :: Vector{QNDiag}, tms_lzlp :: Tuple{Terms, Terms}, tms_c2 :: Terms = 0 * one(Terms), c2_rng :: Vector{Float64} = [0.0] ; eltype = FuzzifiED.ElementType)
    nsec = length(sec)
    cfs = Vector{Confs}(undef, nsec)
    cfs1 = Vector{Confs}(undef, nsec)
    sts = Vector{Matrix{eltype}}(undef, nsec)
    sts1 = Vector{Matrix{eltype}}(undef, nsec)
    l_rng = [ Int64[] for _ ∈ sec ]
    ptr_st = [ Int64[] for _ ∈ sec ]
    l_lookup = [ Dict{Int64, Int64}() for _ ∈ sec ]
    tms_l2 = GetL2Terms(tms_lzlp)
    Threads.@threads for isec ∈ eachindex(sec)
        seci = sec[isec]
        cfs[isec] = Confs(no, seci, qnd ; num_th = 1)

        bs = Basis(cfs[isec])

        l2_mat = Matrix(OpMat(Operator(bs, tms_l2), num_th = 1))
        c2_mat = Matrix(OpMat(Operator(bs, tms_c2), num_th = 1))
        _, st = eigen(Hermitian(l2_mat + √2 * c2_mat))

        l2_val = vec(sum(conj.(st) .* (l2_mat * st) ; dims = 1))
        c2_val = vec(sum(conj.(st) .* (c2_mat * st) ; dims = 1))
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
                !flag && continue
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
            cfs1[isec] = Confs(no, seci1, qnd ; num_th = 1)
            bs1 = Basis(cfs1[isec])
            lp = Operator(bs, bs1, tms_lzlp[2])
            lp_mat = Matrix(OpMat(lp))
            sts1[isec] = lp_mat * sts[isec]
        end
    end
    ptr_sec = cumsum([ptr_st[isec][end] - 1 for isec ∈ eachindex(sec)])
    for isec = 2 : length(sec)
        ptr_st[isec] .+= ptr_sec[isec - 1]
    end
    return SegSpace{eltype}(sec, l_rng, l_lookup, ptr_st, cfs, cfs1, sts, sts1)
end
