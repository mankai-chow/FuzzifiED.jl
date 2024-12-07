module FuzzifiED 

using LinearAlgebra
using WignerSymbols
using SphericalHarmonics
using FuzzifiED_jll
import Base.:+
import Base.:-
import Base.:*
import Base.:/
import Base.:÷
import Base.:^
import Base.zero
import Base.one
import Base.adjoint

include("core/param.jl")

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
export Terms
export ParticleHole
export NormalOrder
export SimplifyTerms

include("core/operator.jl")
export Operator

include("core/opmat.jl")
export OpMat
export GetEigensystem
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
export GetElectronObs
export GetDensityObs
export GetPairingObs

include("models/ang_modes.jl")
export AngModes 
export GetElectronMod
export GetPairingMod
export GetDensityMod
export FilterComponent
export FilterL2

include("fuzzifino/fuzzifino.jl")
export Fuzzifino

export QNDiagFromSites
export ConfsFromSites
export TermsFromOpSum
export OpSumFromTerms
export SitesFromQNDiag
export GetSites
export TruncateQNDiag
export SweepOne
export EasySweep
export GetMPOSites
export GetMPO
export SparseMatrixCSCFromOpMat
export GetEigensystemKrylov
export GetEigensystemCuda

function QNDiagFromSites end
function ConfsFromSites end
function TermsFromOpSum end
function OpSumFromTerms end
function SitesFromQNDiag end
function GetSites end
function TruncateQNDiag end
function SweepOne end
function EasySweep end
function GetMPOSites end
function GetMPO end
function SparseMatrixCSCFromOpMat end
function GetEigensystemKrylov end
function GetEigensystemCuda end

function __init__()
    FuzzifiED.NumThreads = Threads.nthreads()
end

end