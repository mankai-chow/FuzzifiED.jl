module FuzzifiED 

using LinearAlgebra
import Base.:+
import Base.:-
import Base.:*
import Base.:/
import Base.adjoint

LibpathFuzzifiED = filter(f -> endswith(f, ".so"), readdir(dirname(@__FILE__), join = true))[1]

include("basics/confs.jl")
include("basics/basis.jl")
include("basics/term.jl")
include("basics/operator.jl")
include("basics/opmat.jl")

export Confs
export Basis
export Term
export Operator
export OpMat
export GetEigensystem

module ITensorsSupport 
    using ..FuzzifiED
    using ITensors
    include("itensors_support.jl")

    export ConfsFromSites
    export OperatorFromOpSum
end

module Models 
    using ..FuzzifiED 
    using WignerSymbols
    
    include("models/threej.jl")
    include("models/l2.jl")
    include("models/ising.jl")

    export GetInteractionMatrix
    export GenerateL2Terms
    export GenerateIsingConfs
    export GenerateIsingBasis
    export GenerateIsingHamiltonianTerms
end

end