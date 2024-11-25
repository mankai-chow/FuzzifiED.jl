module FuzzifiED 

using LinearAlgebra
using SparseArrays
using WignerSymbols
using SphericalHarmonics
using FuzzifiED_jll
using KrylovKit
import Base.:+
import Base.:-
import Base.:*
import Base.:/
import Base.:÷
import Base.:^
import Base.zero
import Base.adjoint

include("core/param.jl")
export NumThreads
export SilentStd
export Libpath
export ElementType
export OpenHelp!

include("core/qn.jl")
export QNDiag
export QNOffd

include("core/confs.jl")
export Confs
export GetConfId

include("core/basis.jl")
export Basis
export GetConfWeight

include("core/term.jl")
export Term
export ParticleHole
export NormalOrder
export SimplifyTerms
export SimplifyTermsOld

include("core/operator.jl")
export Operator

include("core/opmat.jl")
export OpMat
export GetEigensystem
export GetEigensystemKrylov
export SparseMatrixCSCFromOpMat
export MatrixFromOpMat

include("core/entangle.jl")
export StateDecompMat
export GetEntSpec

include("models/qndiag.jl")
export GetNeQNDiag
export GetLz2QNDiag
export GetFlavQNDiag
export GetZnfChargeQNDiag
export GetPinOrbQNDiag

include("models/qnoffd.jl")
export GetParityQNOffd
export GetFlavPermQNOffd
export GetRotyQNOffd

include("models/opterms.jl")
export GetIntMatrix
export GetDenIntTerms
export GetPairIntTerms
export GetPolTerms
export GetIsingIntTerms
export GetL2Terms
export GetC2Terms

include("models/sphere_obs.jl")
export SphereObs
export StoreComps!
export StoreComps
export Laplacian
export GetComponent
export GetPointValue
export Electron
export Density
export Pairing
export PairObs

include("models/ang_modes.jl")
export AngModes 
export ElecMod
export PairMod
export DenMod
export FilterComponent
export FilterL2

include("archieve/ar_core.jl")

include("archieve/ar_models.jl")
export GetLzQnu
export GetLzZnQnu
export GetLzConfs 
export GetLzZnConfs
export GetIsingQnz
export GetIsingBasis
export GetSnBasis
export GetXPolTerms
export GetZPolTerms

include("fuzzifino/fuzzifino.jl")
export Fuzzifino

export QNDiagFromSites
export ConfsFromSites
export TermsFromOpSum
export OpSumFromTerms
export SitesFromQNDiag
export TruncateQNDiag
export SweepOne
export EasySweep
export GetMPOSites
export GetMPO
export TruncateQnu
export SitesFromQnu
export GetEigensystemCuda

function QNDiagFromSites end
function ConfsFromSites end
function TermsFromOpSum end
function OpSumFromTerms end
function SitesFromQNDiag end
function TruncateQNDiag end
function SweepOne end
function EasySweep end
function GetMPOSites end
function GetMPO end
function TruncateQnu end
function SitesFromQnu end
function GetEigensystemCuda end

function __init__()
    FuzzifiED.NumThreads = Threads.nthreads()
end

end