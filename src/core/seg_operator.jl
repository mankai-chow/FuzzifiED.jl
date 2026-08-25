export SegOperator, BuildSegOperator, BuildSegOperators, BuildTransfSegOperator


"""
    SegOperator{Float64}
    SegOperator{ComplexF64}

The mutable type `SegOperator` stores the action of a spherical-symmetric operator from an initial [SegSpace](@ref SegSpace) to a final SegSpace. By virtue of the Wigner-Eckart theorem, the full ``m``-dependence is factored out and only the reduced matrix elements ``⟨l_2\\|[Φ]_l\\|l_1⟩`` – independent of ``l^z`` – need to be kept, and 
```math
    ⟨l_2m_2|[Φ]_{lm}|l_1m_1⟩=(-1)^{l_2-m_2}\\begin{pmatrix}l_2&l&l_1\\\\-m_2&m&m_1\\end{pmatrix}⟨l_2\\|[Φ]_l\\|l_1⟩
```
The matrix ``⟨Q_2l_2α_2\\|[Φ]_l\\|Q_1l_1α_1⟩`` is stored in blocks of QNDiag – for given sectors ``Q_{\\{12\\}}``, the matrix elements ``(M_{l_2l_1})_{α_2α_1}`` are stored. 

# Fields

* `sgspd :: SegSpace` and `sgspf :: SegSpace` are the initial and final segment spaces.
* `colptr :: Vector{Int64}` and `rowid :: Vector{Int64}` store the allowed blocks of sectors ``Q_{12}`` in the format of a CSC sparse matrix.
* `elmat :: Vector{Matrix{Matrix{T}}}` stores, for each block of sectors, the reduced matrix elements. It takes five indices `elmat[e][ich, jch][i, j]`, where `e` is the index for the sector block, `ich` and `jch` are the channel index, and `i` and `j` are the state index within each channel. 
"""
mutable struct SegOperator{T <: Union{Float64, ComplexF64}}
    colptr :: Vector{Int64}
    rowid :: Vector{Int64}
    elmat :: Vector{Matrix{Matrix{T}}}
end


"""
    BuildSegOperator(sgspd :: SegSpace{T}[, sgspf :: SegSpace{T}], amd :: AngModes, ll :: Int64, secop :: Vector{Int64} ; full_mat :: Bool, num_th :: Int64) :: SegOperator{T}
    BuildSegOperator(sgspd :: SSegSpace{T}[, sgspf :: SSegSpace{T}], amd :: SAngModes, ll :: Int64, secop :: Vector{Int64} ; full_mat :: Bool, num_th :: Int64) :: SegOperator{T}

constructs a [SegOperator](@ref SegOperator) from the angular modes `amd` acting on a single segment and the angular momentum `ll`. For each allowed pair of sectors and angular momenta, it computes the reduced matrix element from the full matrix element via the Wigner-Eckart theorem by dividing out the phase and the ``3j``-symbol. When the ``3j``-symbol vanishes (for odd ``l'+L+l`` when ``m_1=m_2=0``) the reduced matrix element is instead calculated from the pre-stored``m=1`` components.

# Arguments

* `sgspd :: SegSpace{T}` or `sgspd :: SSegSpace{T}` is the initial segment space. 
* `sgspf :: SegSpace{T}` or `sgspf :: SSegSpace{T}` is the final segment space. Facultative, the same as `sgspd` by default.
* `amd :: AngModes` or `amd :: SAngModes` is the angular modes.
* `ll :: Int64` is twice the angular momentum ``2l`` of the segment operator.
* `secop :: Vector{Int64}` is the QNDiag shift of the operator.
* `full_mat :: Bool` determines whether a full matrix instead of sparse matrix is generated when generating the segment operators. It should be turned on if the flavour symmetry ``C_2`` is not selected when constructing the segment space. Facultative, `false` by default.
* `num_th :: Int64` is the number of threads. Facultative, `NumThreads` by default.

# Output

* `sgop :: SegOperator{T}` is the resulting segment operator.
"""
function BuildSegOperator(sgspd :: AbstractSegSpace{T}, sgspf :: AbstractSegSpace{T}, amd :: Union{AngModes, SAngModes, Symbol}, ll :: Int64, secop :: Vector{Int64} ; full_mat :: Bool = false, num_th = FuzzifiED.NumThreads) where T <: Union{Float64, ComplexF64}
    _SegOperator(:: SegSpace, bsd, bsf, tms) = Operator(bsd, bsf, tms)
    _SegOperator(:: SSegSpace, bsd, bsf, tms) = SOperator(bsd, bsf, tms)
    _SegOneAngModes(:: SegSpace) = one(AngModes)
    _SegOneAngModes(:: SSegSpace) = one(SAngModes)

    (amd === :Identity) && (amd = _SegOneAngModes(sgspd))
    index = 0
    colptr = zeros(Int64, size(sgspd.sec, 2) + 1)
    colptr[1] = 1
    rowid = Int64[]
    elmat = Matrix{Matrix{T}}[]
    modul = sgspd.sec_modul

    for j in axes(sgspd.sec, 2)
        secd = sgspd.sec[:, j]
        std = sgspd.sts[j]
        md = secd[2]
        (md == 0) && (std1 = sgspd.sts1[j])

        for i in axes(sgspf.sec, 2)
            secf = sgspf.sec[:, i]
            EquivSec(secd .+ secop, secf, modul) || continue
            index += 1
            push!(rowid, i)
            mf = secf[2]
            mm = mf - md

            stf = sgspf.sts[i]
            tms = GetComponent(amd, ll/2, mm/2)
            op = _SegOperator(sgspd, sgspd.bs[j], sgspf.bs[i], tms)
            op_mat = OpMat(op ; num_th, disp_std = false)
            if full_mat 
                hmt_block = stf' * Matrix(op_mat) * std
            else
                hmt_block = stf' * *(op_mat, std ; num_th)
            end

            if (md == 0 && mf == 0 && ll > 0) # when 3j could vanish
                tms1 = GetComponent(amd, ll/2, mm/2 - 1)
                op1 = _SegOperator(sgspd, sgspd.bs1[j], sgspf.bs[i], tms1)
                op_mat1 = OpMat(op1 ; num_th, disp_std = false)
                if full_mat 
                    hmt_block1 = stf' * Matrix(op_mat1) * std1
                else
                    hmt_block1 = stf' * *(op_mat1, std1 ; num_th)
                end
            end

            hmt_mat = Matrix{Matrix{T}}(undef, length(sgspf.l_rng[i]), length(sgspd.l_rng[j]))

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
                        hmt_mat[il, jl] = hmt_block[rngi, rngj] / fac3j
                    else
                        fac3j1 = wigner3j(lf/2, ll/2, ld/2, -mf/2, mm/2 - 1, md/2 + 1)
                        ((lf - mf) % 4 == 2) && (fac3j1 = -fac3j1)
                        hmt_mat[il, jl] = hmt_block1[rngi, rngj] / √(ld/2 * (ld/2 + 1)) / fac3j1
                    end
                end
            end
            push!(elmat, hmt_mat)
        end
        colptr[j + 1] = index + 1
    end
    return SegOperator{T}(colptr, rowid, elmat)
end
BuildSegOperator(sgspd :: AbstractSegSpace, amd :: Union{AngModes, SAngModes, Symbol}, ll :: Int64, secop :: Vector{Int64} ; full_mat :: Bool = false, num_th = FuzzifiED.NumThreads) = BuildSegOperator(sgspd, sgspd, amd, ll, secop ; full_mat, num_th)


"""
    BuildSegOperators(sgspd :: Vector{<:AbstractSegSpace{T}}[, sgspf :: Vector{<:AbstractSegSpace{T}}], cpd :: CoupleDecomps ; p_rng :: Vector{Int64}, ident_seg :: Vector{Int64}, full_mat :: Bool, num_th :: Int64) :: Matrix{SegOperator}

constructs, in parallel, all the [SegOperators](@ref SegOperator) required to assemble the composite operators described by the coupling decompositions `cpd`.

# Arguments

* `sgspd :: Vector{<:AbstractSegSpace{T}}` is the list of the initial segment spaces. Its elements are either `SegSpace` or `SSegSpace`. 
* `sgspf :: Vector{<:AbstractSegSpace{T}}` is the list of the final segment spaces. Its elements are either `SegSpace` or `SSegSpace`. Facultative, the same as `sgspd` by default.
* `cpd :: CoupleDecomps` is the list of coupling decompositions, _e. g._, an assembled Hamiltonian.
* `p_rng :: Vector{Int64}`. When specified, only the SegOperators of the specified parts will be generated. It must be of the same length as `sgspd`. An array ``1:N_p`` by default.
* `ident_seg :: Vector{Int64}` labels the identical segments. If given, an array of length ``N_p``, identical segments carry identical index. Facultative, empty by default, marking no identification. 
* `full_mat :: Bool` determines whether a full matrix instead of sparse matrix is generated when generating the segment operators. It should be turned on if the flavour symmetry ``C_2`` is not selected when constructing the segment space. Facultative, `false` by default.
* `num_th :: Int64` is the number of threads. BLAS and Generation of ``OpMat`` are parallelized. Facultative, `NumThreads` by default.
* `disp_std :: Bool`, whether or not the log shall be displayed. Facultative, `!SilentStd` by default. 

# Output

* `sgop :: Matrix{SegOperator}` is a matrix of segment operators of size ``N_p×N_d``, where ``N_p`` is the number of parts and ``N_d`` the total number of channels of `cpd`. It is passed together with the same `cpd` to [BuildCompOperator](@ref BuildCompOperator).
"""
function BuildSegOperators(sgspd :: Vector{<:AbstractSegSpace{T}}, sgspf :: Vector{<:AbstractSegSpace{T}}, cpd :: CoupleDecomps ; p_rng :: Vector{Int64} = collect(eachindex(sgspd)), ident_seg :: Vector{Int64} = Int64[], full_mat :: Bool = false, num_th :: Int64 = FuzzifiED.NumThreads, disp_std = !FuzzifiED.SilentStd) where T <: Union{Float64, ComplexF64}
    nd = length(cpd)
    np = length(p_rng)

    sgop_cnx = Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}()
    if (!isempty(ident_seg))
        _id_coeff(coeff :: ComplexF64) = round(abs(coeff) + real(coeff)/2 + imag(coeff)/2, digits = 8)
        _id_tms(tms :: Union{Terms, STerms}) = isempty(tms) ? 0.0 : sum(_id_coeff.(getproperty.(tms, :coeff)))
        _id_amd(:: Symbol) = [-1.0]
        _id_amd(amd :: Union{AngModes, SAngModes}) = [amd.l2m ; [_id_tms(GetComponent(amd, l, l)) for l = amd.l2m / 2 : -1 : 0]]
        id = [ [ident_seg[p] ; cpdi.ch[1, p] ; cpdi.sec[:, p] ; _id_amd(cpdi.amd[p]) ] for p in p_rng, cpdi in cpd]
        coord_sort = sort(vec(collect(Iterators.product(1 : length(p_rng), 1 : length(cpd)))), by = ij -> id[ij...])

        for i in eachindex(coord_sort)
            coordi = coord_sort[i]
            for j = i - 1 : -1 : 1
                coordj = coord_sort[j]
                id[coordi...] ≠ id[coordj...] && continue 
                sgop_cnx[coordi] = haskey(sgop_cnx, coordj) ? sgop_cnx[coordj] : coordj
            end
        end
    end

    sgop = Matrix{SegOperator}(undef, np, nd)
    for d = 1 : nd, ip = 1 : np
        haskey(sgop_cnx, (ip, d)) && continue
        p = p_rng[ip]
        amd = cpd[d].amd[p]
        secop = cpd[d].sec[:, p]
        ll = cpd[d].ch[1, p]
        sgop[ip, d] = BuildSegOperator(sgspd[ip], sgspf[ip], amd, ll, secop ; full_mat, num_th)
    end
    for ((ip, d), (ip1, d1)) in sgop_cnx
        sgop[ip, d] = sgop[ip1, d1]
    end
    info_str = "FINISH BUILDING $np * $nd SEG OPERATORS"
    isempty(sgop_cnx) || (info_str *= ", $(np * nd - length(sgop_cnx)) INDEPENDENT")
    disp_std && @info info_str
    return sgop
end
BuildSegOperators(sgspd :: Vector{<:AbstractSegSpace{T}}, cpd :: CoupleDecomps ; p_rng :: Vector{Int64} = collect(eachindex(sgspd)), ident_seg :: Vector{Int64} = collect(1 : maximum(p_rng)), full_mat :: Bool = false, num_th :: Int64 = FuzzifiED.NumThreads, disp_std = !FuzzifiED.SilentStd) where T <: Union{Float64, ComplexF64} = BuildSegOperators(sgspd, sgspd, cpd ; p_rng, ident_seg, full_mat, num_th, disp_std)


function BuildTransfSegOperator(sgspd :: AbstractSegSpace{T}, sgspf :: AbstractSegSpace{T}, qnf :: Union{QNOffd,SQNOffd}, sec_cnx :: Vector{Int64} ; full_mat :: Bool = false, num_th = FuzzifiED.NumThreads) where T <: Union{Float64, ComplexF64}
    _SegTransf(:: SegSpace, bsd, bsf, qnf :: QNOffd) = Transf(bsd, bsf, qnf)
    _SegTransf(:: SSegSpace, bsd, bsf, qnf :: SQNOffd ) = STransf(bsd, bsf, qnf)

    colptr = collect(1 : length(sec_cnx) + 1)
    rowid = sec_cnx
    elmat = Matrix{Matrix{T}}[]
    modul = sgspd.sec_modul

    for j in axes(sgspd.sec, 2)
        i = sec_cnx[j]
        std = sgspd.sts[j]
        stf = sgspf.sts[i]
        trs = _SegTransf(sgspd, sgspd.bs[j], sgspf.bs[i], qnf)
        hmt_block = stf' * *(trs, std ; num_th)

        hmt_mat = Matrix{Matrix{T}}(undef, length(sgspf.l_rng[i]), length(sgspd.l_rng[j]))
        for jl in eachindex(sgspd.l_rng[j])
            ld = sgspd.l_rng[j][jl]
            rngj = (sgspd.ptr_st[j][jl] + 1 : sgspd.ptr_st[j][jl + 1]) .- sgspd.ptr_st[j][1]
            for il in eachindex(sgspf.l_rng[i])
                lf = sgspf.l_rng[i][il]
                ld == lf || continue 
                rngi = (sgspf.ptr_st[i][il] + 1 : sgspf.ptr_st[i][il + 1]) .- sgspf.ptr_st[i][1]
                hmt_mat[il, jl] = hmt_block[rngi, rngj] * √(ld+1)
            end
        end
        push!(elmat, hmt_mat)
    end

    return SegOperator{T}(colptr, rowid, elmat)
end
BuildTransfSegOperator(sgspd :: AbstractSegSpace{T}, qnf :: Union{QNOffd}, sec_cnx :: Vector{Int64} ; num_th = FuzzifiED.NumThreads) where T <: Union{Float64, ComplexF64} = BuildTransfSegOperator(sgspd, sgspd, qnf, sec_cnx ; num_th)