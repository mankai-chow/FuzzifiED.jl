module FuzzifiED 

using LinearAlgebra
using SparseArrays
using Requires
using WignerSymbols
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

include("models/threej.jl")
include("models/l2.jl")
include("models/ising.jl")
include("models/spn.jl")

export GetIntMatrix
export GetL2Terms

export GetIsingQnu
export GetIsingConfs
export GetIsingBasis
export GetIsingIntTerms
export GetXPolTerms
export GetZPolTerms

export GetSpnQnu
export GetSpnConfs 
export GetSpnBasis 
export GetIdDenIntTerms 
export GetSpnPairIntTerms
export GetSpnC2Terms

function __init__()
    @require ITensors = "9136182c-28ba-11e9-034c-db9fb085ebd5" begin
        using ITensors
        include("itensors_support.jl")
        include("easy_sweep.jl")
    
        export ConfsFromSites
        export TermsFromOpSum
        export OpSumFromTerms
        export SitesFromQN
    end
end

end