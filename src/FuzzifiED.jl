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

include("core/confs.jl")
export Confs
export GetConfId

include("core/basis.jl")
export Basis
export GetConfWeight

include("core/term.jl")
export Term

include("core/operator.jl")
export Operator

include("core/opmat.jl")
export OpMat
export GetEigensystem
export SparseMatrixCSCFromOpMat

include("models/threej.jl")
export GetIntMatrix

include("models/l2.jl")
export GetL2Terms

include("models/ising.jl")
export GetIsingQnu
export GetIsingConfs
export GetIsingBasis
export GetIsingIntTerms
export GetXPolTerms
export GetZPolTerms

include("models/ising_x.jl")
export GetIsingXQnu
export GetIsingXConfs
export GetIsingXIntTerms

include("models/ising_def.jl")
export GetIsingDefQnu
export GetIsingDefConfs
export GetIsingDefIntTerms
export GetDefXPolTerms

include("models/spn.jl")
export GetSpnQnu
export GetSpnConfs 
export GetSpnBasis 
export GetIdDenIntTerms 
export GetSpnPairIntTerms
export GetSpnC2Terms

include("models/nn_int.jl")
export GetLzQnu
export GetLzConfs 
export GetS3Basis
export GetDenIntTerms
export GetPolTerms

function __init__()

    @require ITensors = "9136182c-28ba-11e9-034c-db9fb085ebd5" begin
        using ITensors
        using ITensors.HDF5
        using ITensorMPOConstruction

        include("itensors_support/itensors_format.jl")
        export ConfsFromSites
        export TermsFromOpSum
        export OpSumFromTerms
        export SitesFromQnu
        export TruncateQnu

        include("itensors_support/easy_sweep.jl")
        export SweepOne
        export EasySweep
        export GetMpoSites
        export GetMpo
    end
end

end