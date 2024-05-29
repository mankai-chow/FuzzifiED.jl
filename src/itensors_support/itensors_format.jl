"""
    function ConfsFromSites(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, qn_s :: QN) :: Confs

Converts a `Sites` object in the `ITensors` package to the `Confs` object

# Arguments 

* `sites :: Vector{Index{Vector{Pair{QN, Int64}}}}` is a `Sites` object. Only `Fermion` site type is supported, and the quantum numbers of the `0` state must be all zero. Note that this will subject to the limitation in ITensors that the number of conserved quantities must be less than 4. 
* `qn_s :: QN` is a `QN` object that specifies the the quantum number of the selected configuration.

"""
function ConfsFromSites(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, qn_s :: QN)
    no = length(sites)
    qnu_names = []
    modul = Vector{Int64}()
    for site in sites 
        for qn in qn(site, 2)
            if abs(qn.modulus) == 0 continue end
            if qn.name in qnu_names continue end
            push!(qnu_names, qn.name)
            push!(modul, abs(qn.modulus))
        end
    end
    qnu_o = []
    for qnu_name in qnu_names 
        qnu_oi = []
        for site in sites 
            push!(qnu_oi, val(qn(site, 2), qnu_name))
        end
        push!(qnu_o, qnu_oi)
    end
    qnu_s = [ val(qn_s, qnu_name) for qnu_name in qnu_names ]
    return Confs(no, qnu_s, qnu_o ; modul)
end 

"""
    function ConfsFromSites(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, cf_ref :: Vector{Int64}) :: Confs

* `sites :: Vector{Index{Vector{Pair{QN, Int64}}}}` is a `Sites` object, where the modulus-``\\pm 1`` quantum numbers will be taken out, and QNs with other modulus will be discarded automatically. Also note that only `Fermion` site type is supported, and the quantum numbers of the `0` state must be all zero. Note that this will subject to the limitation in ITensors that the number of conserved quantities must be less than 4. 
* `cf_ref :: Vector{Int64})` is a reference configuration composed of `0` and `1` in `ITensors` format.

"""
function ConfsFromSites(sites :: Vector{Index{Vector{Pair{QN, Int64}}}}, cf_ref :: Vector{Int64})
    no = length(sites)
    qnu_names = []
    modul = Vector{Int64}()
    for site in sites 
        for qn in qn(site, 2)
            if abs(qn.modulus) == 0 continue end
            if qn.name in qnu_names continue end
            push!(qnu_names, qn.name)
            push!(modul, abs(qn.modulus))
        end
    end
    qnu_o = []
    for qnu_name in qnu_names 
        qnu_oi = []
        for site in sites 
            push!(qnu_oi, val(qn(site, 2), qnu_name))
        end
        push!(qnu_o, qnu_oi)
    end
    qnu_s = [ cf_ref' * qnu_oi for qnu_oi in qnu_o ]
    return Confs(no, qnu_s, qnu_o ; modul)
end 

function ITensors.space( :: SiteType"Fermion"; no :: Int64 = 1, o :: Int = 1, 
    qnu_o :: Vector{Any} = [fill(1, no)], 
    qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], 
    modul = [1 for qn in eachindex(qnu_o)] )
    return [
        QN(
            [ (qnu_name[i], qnu_o[i][o] * n, modul[i]) for i in eachindex(qnu_o)]...
        ) => 1 for n = 0 : 1
    ]
end

"""
    TruncateQnu(; qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}, trunc_lth :: Int64, trunc_wt :: Vector{Int64}) 

truncates the list of ``N_U`` QNU's from to a number ``N'_U`` acceptable by ITensors. The new quantum numbers are 
```math
    Q'_1=Q_1,\\ Q'_2=Q_2,\\ \\dots,\\ Q'_{N'_U-1}=Q_{N'_U-1}\\ Q'_{N'_U}=\\lambda_{N'_U}Q_{N'_U}+\\lambda_{N'_U+1}Q_{N'_U+1}+\\dots+\\lambda_{N_U}Q_{N_U}
```

# Arguments

- `qnu_o :: Vector{Vector{Int64}}` stores the charge of each orbital under each conserved quantity. See [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)))) for detail. 
- `qnu_name :: Vector{String}` stores the name of each quantum number. Facultive, QN1, QN2, ... by default. 
- `modul :: Vector{Int64}` stores the modulus of each quantum number. Store 1 if no modulus. Facultive, all 1 by default. 
- `trunc_lth :: Int64` stores the truncated numbers of QNU. Facultive, 3 by default. 
- `trunc_wt :: Vecotr{Int64}` stores the ``N_U-N'_U+1`` coefficients ``\\lambda``. Facultive, ``1,2,4,8,\\dots`` by default. 

# Output

A named tuple with three elements that can be directly fed into [`SitesFromQnu`](@ref)

- `qnu_o :: Vector{Vector{Int64}}` stores the charge of each orbital under each conserved quantity. See [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)))) for detail.
- `qnu_name :: Vector{String}` stores the name of each quantum number.
- `modul :: Vector{Int64}` stores the modulus of each quantum number, 1 if no modulus. 

"""
function TruncateQnu(; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], trunc_lth :: Int64 = 3, trunc_wt :: Vector{Int64} = [ 2 ^ (i - trunc_lth) for i = trunc_lth : length(qnu_o)]) 
    qnu_o1 = [ qnu_o[1 : trunc_lth - 1] ; [sum([ qnu_o[i + trunc_lth - 1] .* trunc_wt[i] for i in eachindex(trunc_wt)])] ] 
    return (qnu_o = qnu_o1, qnu_name = qnu_name[1 : trunc_lth], modul = modul[1 : trunc_lth])
end

"""
    function SitesFromQnu(; qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64})

returns the ITensors Sites object from the information of quantum numbers 

# Arguments 

- `qnu_o :: Vector{Vector{Int64}}` stores the charge of each orbital under each conserved quantity. See [`Confs`](@ref Confs(no :: Int64, qnu_s :: Vector{Int64}, qnu_o :: Vector{Any} ; nor :: Int64 = div(no, 2), modul :: Vector{Int64} = fill(1, length(qnu_s)))) for detail. 
- `qnu_name :: Vector{String}` stores the name of each quantum number. Facultive, QN1, QN2, ... by default. 
- `modul :: Vector{Int64}` stores the modulus of each quantum number. Store 1 if no modulus. Facultive, all 1 by default. 
"""
function SitesFromQnu(; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)])
    no = length(qnu_o[1])
    return [ siteind("Fermion" ; no, o, qnu_o, qnu_name, modul) for o in eachindex(qnu_o[1]) ]
end

"""
    function TermsFromOpSum(opsum :: Sum{Scaled{ComplexF64, Prod{Op}}}) :: Vector{Term}

Converts a `OpSum` object in `ITensors` to a series of terms. Note that the only operators supported are `"C"`, `"Cdag"` `"N"` and `"I"`.

"""
function TermsFromOpSum(opsum :: Sum{Scaled{ComplexF64, Prod{Op}}})
    tms = Vector{Term}(undef,0)
    for i = 1 : length(opsum)
        cstr = []
        for j = 1 : length(opsum[i])
            optype = opsum[i][j].which_op
            opsite = opsum[i][j].sites[1]
            if optype == "Cdag"
                append!(cstr, [1, opsite])
            elseif optype == "C"
                append!(cstr, [0, opsite])
            elseif optype == "N"
                append!(cstr, [1, opsite, 0, opsite])
            elseif optype == "I"
                append!(cstr, [-1, -1])
            else 
                print("The only operator supported are C, Cdag and N")
            end
        end 
        if length(opsum[i]) == 0
            append!(cstr, [-1, -1])
        end
        fac = coefficient(opsum[i])
        push!(tms, Term(fac, cstr))
    end
    return tms
end

"""
    function OpSumFromTerms(tms :: Vector{Term}) :: Sum{Scaled{ComplexF64, Prod{Op}}}

Converts a series of terms to `OpSum` object in `ITensors`.

"""
function OpSumFromTerms(tms :: Vector{Term})
    opsum = OpSum()
    for tm in tms
        optm = Any[tm.coeff]
        for i = 1 : 2 : length(tm.cstr)
            push!(optm, tm.cstr[i] == 0 ? "C" : "Cdag")
            push!(optm, tm.cstr[i + 1])
        end
        opsum += Tuple(optm)
    end
    return opsum
end