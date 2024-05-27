module FuzzifiED 

using LinearAlgebra
using SparseArrays
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
export SparseMatrixCSCFromOpMat
export GetConfId
export GetConfWeight

module ITensorsSupport 
    using ..FuzzifiED
    using ITensors
    include("itensors_support.jl")

    export ConfsFromSites
    export TermsFromOpSum
    export OpSumFromTerms
    export SitesFromQN
end

module Models 
    using ..FuzzifiED 
    using WignerSymbols
    
    include("models/threej.jl")
    include("models/l2.jl")
    include("models/ising.jl")
    include("models/spn.jl")

    export GetIntMatrix
    export GetL2Terms

    export GetIsingConfs
    export GetIsingBasis
    export GetIsingIntTerms
    export GetXPolTerms
    export GetZPolTerms

    export GetSpnConfs 
    export GetSpnBasis 
    export GetIdDenIntTerms 
    export GetSpnPairIntTerms
    export GetSpnC2Terms
end

end