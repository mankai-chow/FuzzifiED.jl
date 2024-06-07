"""
    function StateDecompMat(st :: Vector{<:Number}, bs0 :: Basis, bsa :: Basis, bsb :: Basis, amp_oa :: Vector{ComplexF64}, amp_ob :: Vector{ComplexF64}) :: Matrix{ComplexF64}

Decompose a state ``|ψ⟩=v_I|I⟩`` into a direct-product basis of two subsystems ``|ψ⟩=M_{JI}|I_A⟩|J_B⟩``

# Arguments 

- `st :: Vector{<:Number}` is the state to be decomposed into direct-product basis of two subsystems.
- `bs0 :: Basis` is the total basis. 
- `bsa :: Basis` is the basis for the subsystem A.
- `bsb :: Basis` is the basis for the subsystem B. 
- `amp_oa :: Vector{ComplexF64}` is a complex list of length `no` that specifies the amplitute of each orbital in the subsystem A. For a non-local basis, we decompose each electron into creation operators in two subsystems ``c^†_o=a_{o,A}c^†_{o,A}+a_{o,B}c^†_{o,B}`` and this list specifies ``a_{o,A}``. This is equivalent to ``√{ℱ_{m,A}}`` in [PRB 85, 125308 (2012)](https://dx.doi.org/10.1103/PhysRevB.85.125308) with an extra phase factor. 
- `amp_ob :: Vector{ComplexF64}` is a complex list of length `no` that specifies the amplitute of each orbital in the subsystem B. 

# Output

A complex matrix of dimension `bsb.dim * bsa.dim` that corresponds to the state in the decomposed basis ``|ψ⟩=M_{JI}|I_A⟩|J_B⟩``. This is equivalent to ``R_{μν}^A/√p`` in [PRB 85, 125308 (2012)](https://dx.doi.org/10.1103/PhysRevB.85.125308). After calculating all the sectors, the reduced density matrix will be ``ρ_B=⊕\\mathbf{M}\\mathbf{M}^†``

"""
function StateDecompMat(st :: Vector{<:Number}, bs0 :: Basis, bsa :: Basis, bsb :: Basis, amp_oa :: Vector{<:Number}, amp_ob :: Vector{<:Number})
    st_dcp = Matrix{ComplexF64}(undef, bsb.dim, bsa.dim) ;
    @ccall FuzzifiED_jll.LibpathFuzzifiED.__ent_MOD_decomp_basis(bs0.cfs.no :: Ref{Int64}, bs0.cfs.nor :: Ref{Int64}, 
        bs0.cfs.ncf :: Ref{Int64}, bs0.dim :: Ref{Int64}, bs0.cfs.conf :: Ref{Int64}, bs0.cfs.lid :: Ref{Int64}, bs0.cfs.rid :: Ref{Int64}, bs0.szz :: Ref{Int64}, bs0.cfgr :: Ref{Int64}, bs0.cffac :: Ref{ComplexF64}, bs0.grel :: Ref{Int64}, bs0.grsz :: Ref{Int64}, 
        bsa.cfs.ncf :: Ref{Int64}, bsa.dim :: Ref{Int64}, bsa.cfs.conf :: Ref{Int64}, bsa.cfs.lid :: Ref{Int64}, bsa.cfs.rid :: Ref{Int64}, bsa.szz :: Ref{Int64}, bsa.cfgr :: Ref{Int64}, bsa.cffac :: Ref{ComplexF64}, bsa.grel :: Ref{Int64}, bsa.grsz :: Ref{Int64}, 
        bsb.cfs.ncf :: Ref{Int64}, bsb.dim :: Ref{Int64}, bsb.cfs.conf :: Ref{Int64}, bsb.cfs.lid :: Ref{Int64}, bsb.cfs.rid :: Ref{Int64}, bsb.szz :: Ref{Int64}, bsb.cfgr :: Ref{Int64}, bsb.cffac :: Ref{ComplexF64}, bsb.grel :: Ref{Int64}, bsb.grsz :: Ref{Int64}, 
        ComplexF64.(amp_oa) :: Ref{ComplexF64}, ComplexF64.(amp_ob) :: Ref{ComplexF64}, ComplexF64.(st) :: Ref{ComplexF64}, st_dcp :: Ref{ComplexF64}) :: Nothing
    return st_dcp
end

"""
    function GetEntSpec(st :: Vector{<:Number}, bs0 :: Basis, qnu_s_lst :: Vector{Vector{Vector{Int64}}}, qnz_s_lst :: Vector{Vector{Vector{ComplexF64}}} ; qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}, cyc :: Vector{Int64}, perm_o :: Vector{Vector{Int64}}, ph_o :: Vector{Vector{Int64}}, fac_o :: Vector{Vector{ComplexF64}}, amp_oa :: Vector{<:Number}, amp_ob :: Vector{<:Number} = sqrt.(1 .- abs.(amp_oa .^ 2))) :: Dict{@NamedTuple{qnu_sa, qnz_sa, qnu_sb, qnz_sb}, Vector{Float64}}

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
    no = bs0.cfs.no 
    nor = bs0.cfs.nor
    EntSpec = Dict{@NamedTuple{qnu_sa :: Vector{<:Number}, qnz_sa :: Vector{<:Number}, qnu_sb :: Vector{<:Number}, qnz_sb :: Vector{<:Number}}, Vector{Float64}}()
    for qnu_s in qnu_s_lst 
        cfsa = Confs(no, qnu_s[1], qnu_o ; nor, modul)
        if (cfsa.ncf == 0) continue end 
        cfsb = Confs(no, qnu_s[2], qnu_o ; nor, modul)
        if (cfsb.ncf == 0) continue end 
        for qnz_s in qnz_s_lst 
            bsa = Basis(cfsa, ComplexF64.(qnz_s[1]) ; cyc, perm_o, ph_o, fac_o)
            if (bsa.dim == 0) continue end 
            bsb = Basis(cfsb, ComplexF64.(qnz_s[2]) ; cyc, perm_o, ph_o, fac_o)
            if (bsb.dim == 0) continue end 
            st_dcp = StateDecompMat(st, bs0, bsa, bsb, amp_oa, amp_ob) 
            ent_spec = abs.(svdvals(st_dcp)) .^ 2
            EntSpec[(qnu_sa = qnu_s[1], qnz_sa = qnz_s[1], qnu_sb = qnu_s[2], qnz_sb = qnz_s[2])] = ent_spec
        end
    end
    return EntSpec
end