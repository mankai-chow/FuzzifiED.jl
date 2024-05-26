"""
    function GenerateL2Terms(nm :: Int64, nf :: Int64) :: Vector{Term}

Return the terms for the total angular momentum.
"""
function GenerateL2Terms(nm :: Int64, nf :: Int64)
    s = (nm - 1) / 2.0
    no = nm * nf
    tms_lz = 
        [ begin m = div(o - 1, nf)
            Term(m - s, [1, o, 0, o])
        end for o = 1 : no ]
    tms_lp = 
        [ begin m = div(o - 1, nf)
            Term(sqrt(m * (nm - m)), [1, o, 0, o - nf])
        end for o = nf + 1 : no ]
    tms_lm = tms_lp' 
    tms_l2 = tms_lz * tms_lz - tms_lz + tms_lp * tms_lm
    # Initialise the L2 operator
    return tms_l2
end