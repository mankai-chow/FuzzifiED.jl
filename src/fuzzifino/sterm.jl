"""
The mutable type `STerm` records a STerm that looks like ``Ua^{(p_1)}_{o_1}a^{(p_2)}_{o_2}… a^{(p_l)}_{o_l}`` in an operator, where positive ``o`` denotes fermions and negative ``o`` denotes bosons

```math
    a^{(0)}_o=f_o, a^{(1)}_o=f_o^†, a^{(0)}_{-o}=b_o, a^{(1)}_{-o}=b_o^†
```

# Fields
- `coeff :: ComplexF64` records the coefficient ``U``
- `cstr :: Vector{Int64}` is a length-``2l`` vector ``(p_1,o_1,p_2,o_2,… p_l,o_l)`` recording the operator string

# Method 

It can be generated by the function

    STerm(coeff :: ComplexF64, cstr :: Vector{Int64}) :: ComplexF64 
"""
mutable struct STerm 
    coeff :: ComplexF64 
    cstr :: Vector{Int64}
end 

""" 

`STerms` is an alias for `Vector{STerm}` for convenience

# Initialisation

    STerms(coeff :: Number, cstr :: Vector{Int64})

Gives a `STerms` with a single `STerm`.

"""
const STerms = Vector{STerm}
Base.Vector{T}(coeff :: Number, cstr :: Vector{Int64}) where T <: STerm = [STerm(coeff, cstr)]
zero( :: Type{STerms}) = STerm[]

"""
    *(fac :: Number, tms :: STerms) :: STerms
    -(tms :: STerms) :: STerms
    *(tms :: STerms, fac :: Number) :: STerms
    /(tms :: STerms, fac :: Number) :: STerms

Return the product of a collection of STerms with a number. 
"""
function *(fac :: Number, tms :: STerms)
    return [ STerm(fac * tm.coeff, tm.cstr) for tm in tms ]
end
function -(tms :: STerms)
    return (-1) * tms
end
function *(tms :: STerms, fac :: Number)
    return fac * tms
end
function /(tms :: STerms, fac :: Number)
    return (1 / fac) * tms
end

"""
    +(tms1 :: STerms, tms2 :: STerms) :: STerms
    -(tms1 :: STerms, tms2 :: STerms) :: STerms

Return the naive sum of two series of STerms by taking their union. 
"""
function +(tms1 :: STerms, tms2 :: STerms)
    return [ tms1 ; tms2 ]
end
function -(tms1 :: STerms, tms2 :: STerms)
    return tms1 + (-tms2)
end
function +(tms1 :: STerms, tms2 :: Vararg{STerms})
    return tms1 + +(tms2...)
end


"""
    *(tms1 :: STerms, tms2 :: STerms) :: STerms
    ^(tms :: STerms, pow :: Int64) :: STerms

Return the naive product of two series of STerms or the power of one STerms. The number of STerms equals the product of the number of STerms in `tms1` and `tms2`. For each STerm in `tms1` ``Ua^{(p_1)}_{o_1}…`` and `tms2` ``U'a^{(p'_1)}_{o'_1}…``, a new STerm is formed by taking ``UU'a^{(p_1)}_{o_1}… a^{(p'_1)}_{o'_1}…``
"""
function *(tms1 :: STerms, tms2 :: STerms)
    return vcat([ STerm(tm1.coeff * tm2.coeff, [tm1.cstr ; tm2.cstr])
        for tm1 in tms1, tm2 in tms2 ]...)
end
function *(tms1 :: STerms, tms2 :: Vararg{STerms})
    return tms1 * *(tms2...)
end
function ^(tms :: STerms, pow :: Int64)
    if pow == 1
        return tms
    else
        return tms * tms ^ (pow - 1)
    end
end

function adjoint(tm :: STerm)
    nc = length(tm.cstr)
    cstr1 = [ isodd(i) ? 1 - tm.cstr[nc - i] : tm.cstr[nc + 2 - i] for i = 1 : nc]
    return STerm(conj(tm.coeff), cstr1)
end
"""
    adjoint(tm :: STerm) :: STerm
    adjoint(tms :: STerms) :: STerms

Return the Hermitian conjugate of a series of STerms. For each STerm ``Ua^{(p_1)}_{o_1}a^{(p_2)}_{o_2}… a^{(p_l)}_{o_l}``, the adjoint is ``\\bar{U}a^{(1-p_l)}_{o_l}… a^{(1-p_2)}_{o_2}a^{(1-p_1)}_{o_1}``
"""
function adjoint(tms :: STerms)
    return adjoint.(tms)
end


"""
    NormalOrder(tm :: STerm) :: STerms

rearrange a STerm such that 
- the creation operators must be commuted in front of the annihilation operator 
- the orbital index of the creation operators are in ascending order and the annihilation operators in descending order. 
return a list of STerms whose result is equal to the original STerm. 
"""
function NormalOrder(tm :: STerm)
    coeff0 = tm.coeff
    cstr0 = tm.cstr
    for i = 1 : 2 : length(cstr0) - 3
        if (cstr0[i] == -1) 
            cstr1 = deepcopy(cstr0)
            deleteat!(cstr1, i : i + 1)
            return(NormalOrder(STerm(coeff0, cstr1)))
        end
        if (cstr0[i] == 0 && cstr0[i + 2] == 1)
            if (cstr0[i + 1] == cstr0[i + 3])
                comm = -sign(cstr0[i + 1])
                cstr_nrm = deepcopy(cstr0)
                cstr_com = deepcopy(cstr0)
                cstr_nrm[i : i + 1], cstr_nrm[i + 2 : i + 3] = cstr_nrm[i + 2 : i + 3], cstr_nrm[i : i + 1]
                deleteat!(cstr_com, i : i + 3)
                return([ NormalOrder(STerm(comm * coeff0, cstr_nrm)) ; 
                    NormalOrder(STerm(coeff0, cstr_com))])
            else
                comm = (cstr0[i + 1] > 0 && cstr0[i + 3] > 0) ? -1 : 1
                cstr_nrm = deepcopy(cstr0)
                cstr_nrm[i : i + 1], cstr_nrm[i + 2 : i + 3] = cstr_nrm[i + 2 : i + 3], cstr_nrm[i : i + 1]
                return(NormalOrder(STerm(comm * coeff0, cstr_nrm)))
            end
        elseif (cstr0[i] == cstr0[i + 2])
            if (cstr0[i + 1] == cstr0[i + 3]) 
                if (cstr0[i + 1] > 0) return STerm[] end 
            elseif ((cstr0[i] == 1) == (cstr0[i + 1] > cstr0[i + 3]))
                comm = (cstr0[i + 1] > 0 && cstr0[i + 3] > 0) ? -1 : 1
                cstr_nrm = deepcopy(cstr0)
                cstr_nrm[i : i + 1], cstr_nrm[i + 2 : i + 3] = cstr_nrm[i + 2 : i + 3], cstr_nrm[i : i + 1]
                return(NormalOrder(STerm(comm * coeff0, cstr_nrm)))
            end 
        end
    end
    if length(cstr0) == 0 return STerms(coeff0, [-1, -1]) end
    if (cstr0[end - 1] == -1 && length(cstr0) > 2) return STerms(coeff0, cstr0[1 : end - 2]) end
    return STerm[tm]
end


"""
    SimplifyTerms(tms :: STerms ; cutoff :: Float64 = eps(Float64)) :: STerms

simplifies the sum of STerms such that 
- each STerm is normal ordered,
- like STerms are combined, and STerms with zero coefficients are removed.

# Argument 

- `cutoff :: Float64` is the cutoff such that STerms with smaller absolute value of coefficients will be neglected. Facultative, `eps(Float64)` by default. 

"""
function SimplifyTerms(tms :: STerms ; cutoff :: Float64 = eps(Float64)) :: STerms
    dictlock = [ ReentrantLock() for i = 1 : 64 ]
    dict_tms = [ Dict{Vector{Int64}, ComplexF64}() for i = 1 : 64 ]
    
    Threads.@threads for tm in tms 
        tm1 = NormalOrder(tm)
        for tmi in tm1 
            id = 1 + tmi.cstr[2] & 63
            Threads.lock(dictlock[id]) 
            try
                if haskey(dict_tms[id], tmi.cstr) 
                    dict_tms[id][tmi.cstr] += tmi.coeff
                else
                    dict_tms[id][tmi.cstr] = tmi.coeff
                end
            finally
                Threads.unlock(dictlock[id])
            end
        end
    end
    tms_f = [ 
        STerm(coeff_i, cstr_i)
        for dict_tms_i in dict_tms 
        for (cstr_i, coeff_i) in dict_tms_i if abs(coeff_i) > cutoff
    ]
    return tms_f
end