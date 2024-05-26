module FuzzifiED 

using LinearAlgebra
import Base.:*
LibpathFuzzifiED = filter(f -> endswith(f, ".so"), readdir(dirname(@__FILE__), join = true))[1]

include("confs.jl")
include("basis.jl")
include("operator.jl")
include("opmat.jl")

export Confs
export Basis
export Term
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