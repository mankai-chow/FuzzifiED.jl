"""
    GetRmOrbQnu(no :: Int64 ; qnu_o :: Vector{Vector{Int64}}, qnu_name :: Vector{String}, modul :: Vector{Int64}, empt_o :: Vector{Int64}, fill_o :: Vector{Int64})

Pin part of the orbitals by setting the total electron number in the pinned orbitals to be a quantum number. This method is used when calculating defects and orbital space--cut entanglement.  

# Arguments 

- `no :: Int64` is the number of orbitals 
- `qnu_o`, `qnu_name`, `modul` specifies the QNUs before the orbitals are pinned. The latter two are facultative. 
- `empt_o :: Vector{Int64}` specifies the collection of orbitals to be set empty. 
- `fill_o :: Vector{Int64}` specifies the collection of orbitals to be set filled. 

# Output

A named tuple containing updated `qnu_o`, `qnu_name` and `modul`. A quantum number containing the total electron number of the empty orbitals is appended to `qnu_o`. In generating the Confs, this QNU should be set to 0. If fill_o is not empty, another QNU containing the total electron number of the filled orbitals is appended, to be set to `length(fill_o)`. 
"""
function GetRmOrbQnu(no :: Int64 ; qnu_o :: Vector{Any}, qnu_name :: Vector{String} = [ "QN" * string(qn) for qn in eachindex(qnu_o)], modul :: Vector{Int64} = [1 for qn in eachindex(qnu_o)], empt_o :: Vector{Int64}, fill_o :: Vector{Int64} = Int64[])
    qnu_o1 = [ qnu_o ; [[ i in empt_o ? 1 : 0 for i = 1 : no ]] ]
    qnu_name1 = [ qnu_name1 ; "Ne empty" ]
    modul1 = [ modul ; 1 ]
    if (!isempty(fill_o))
        push!(qnu_o1, [ i in fill_o ? 1 : 0 for i = 1 : no ])
        push!(qnu_name1, "Ne filled")
        push!(modul1, 1)
    end
    return (qnu_o = qnu_o1, qnu_name = qnu_name1, modul = modul1)
end