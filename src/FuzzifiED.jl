module FuzzifiED 

using LinearAlgebra
import Base.:*
LibpathFuzzifiED = filter(f -> endswith(f, ".so"), readdir(dirname(@__FILE__), join = true))[1]

include("basics.jl")

export Confs
export Basis
export Operator
export OpMat
export GetEigensystem

module ITensorSupport 
    using ..FuzzifiED
    using ITensors
    include("itensors_support.jl")

    export ConfsFromSites
    export OperatorFromOpSum
end

end