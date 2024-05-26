"""
    mutable struct Term

A `Term` object records a term that looks like ``Uc^{(p_1)}_{o_1}c^{(p_2)}_{o_2}\\dots c^{(p_l)}_{o_l}`` in an operator

# Fields
- `coeff :: ComplexF64` records the coefficient ``U``
- `cstr :: Vector{Int64}` is a length-``2l`` vector ``(p_1,o_1,p_2,o_2,\\dots p_l,o_l)`` recording the operator string
"""
mutable struct Term 
    coeff :: ComplexF64 
    cstr :: Vector{Int64}
end 

"""
    function *(fac :: Number, tms :: Vector{Term}) :: Vector{Term}
    function -(tms :: Vector{Term}) :: Vector{Term}
    function *(tms :: Vector{Term}, fac :: Number) :: Vector{Term}
    function /(tms :: Vector{Term}, fac :: Number) :: Vector{Term}

Return the product of a number with terms
"""
function *(fac :: Number, tms :: Vector{Term})
    return [ Term(fac * tm.coeff, tm.cstr) for tm in tms ]
end
function -(tms :: Vector{Term})
    return (-1) * tms
end
function *(tms :: Vector{Term}, fac :: Number)
    return fac * tms
end
function /(tms :: Vector{Term}, fac :: Number)
    return (1 / fac) * tms
end

"""
    function +(tms1 :: Vector{Term}, tms2 :: Vector{Term}) :: Vector{Term}
    function -(tms1 :: Vector{Term}, tms2 :: Vector{Term}) :: Vector{Term}

Return the sum of two series of terms
"""
function +(tms1 :: Vector{Term}, tms2 :: Vector{Term})
    return [ tms1 ; tms2 ]
end
function -(tms1 :: Vector{Term}, tms2 :: Vector{Term})
    return tms1 + (-tms2)
end

"""
    function *(tms1 :: Vector{Term}, tms2 :: Vector{Term}) :: Vector{Term}

Return the product of two series of terms
"""
function *(tms1 :: Vector{Term}, tms2 :: Vector{Term})
    return vcat([ Term(tm1.coeff * tm2.coeff, [tm1.cstr ; tm2.cstr])
        for tm1 in tms1, tm2 in tms2 ]...)
end

function adjoint(tm :: Term)
    nc = length(tm.cstr)
    cstr1 = [ isodd(i) ? 1 - tm.cstr[nc - i] : tm.cstr[nc + 2 - i] for i = 1 : nc]
    return Term(conj(tm.coeff), cstr1)
end
"""
    function adjoint(tms :: Vector{Term}) :: Vector{Term}

Return the Hermitian conjugate of a series of terms
"""
function adjoint(tms :: Vector{Term})
    return adjoint.(tms)
end
