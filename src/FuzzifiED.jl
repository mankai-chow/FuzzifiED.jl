module FuzzifiED 

using LinearAlgebra
using SparseArrays
using Requires
using WignerSymbols
using SphericalHarmonics
using FuzzifiED_jll
import Base.:+
import Base.:-
import Base.:*
import Base.:/
import Base.zero
import Base.adjoint

"""
    SilentStd :: bool = false 

a flag to determine whether logs of the FuzzifiED functions should be turned off. False by default. If you want to evaluate without log, put `FuzzifiED.SilentStd = true`
"""
SilentStd = false
NumThreads = Threads.nthreads() 
export SilentStd
export NumThreads

include("core/confs.jl")
export Confs
export GetConfId

include("core/basis.jl")
export Basis
export GetConfWeight

include("core/term.jl")
export Term
export NormalOrder
export SimplifyTerms

include("core/operator.jl")
export Operator

include("core/opmat.jl")
export OpMat
export GetEigensystem
export SparseMatrixCSCFromOpMat
export MatrixFromOpMat

include("core/entangle.jl")
export StateDecompMat
export GetEntSpec

include("models/threej.jl")
export GetIntMatrix

include("models/l2.jl")
export GetL2Terms

include("models/nn_int.jl")
export GetSnBasis
export GetDenIntTerms
export GetPairIntTerms
export GetPolTerms

include("models/sphere_obs.jl")
export SphereObs
export StoreComps!
export StoreComps
export Laplacian
export GetComponent
export GetPointValue
export Electron
export Density

include("models/ising.jl")
export GetLzQnu
export GetLzZnQnu
export GetLzConfs 
export GetLzZnConfs
export GetIsingQnz
export GetIsingBasis
export GetIsingIntTerms
export GetXPolTerms
export GetZPolTerms

include("models/spn.jl")
export GetSpnQnu
export GetSpnConfs
export GetSpnBasis
export GetSpnPairIntTerms
export GetSpnC2Terms

include("models/ising_def.jl")
export GetIsingDefQnu
export GetIsingDefConfs
export GetIsingDefIntTerms
export GetDefXPolTerms

function __init__()

    @require ITensors = "9136182c-28ba-11e9-034c-db9fb085ebd5" begin
        using ITensors
        using ITensors.HDF5

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