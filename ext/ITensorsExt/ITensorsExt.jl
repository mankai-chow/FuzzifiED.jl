module ITensorsExt

using FuzzifiED
using ITensors 
using ITensorMPS
include("itensors_format.jl")
export QNDiagFromSites
export ConfsFromSites
export TermsFromOpSum
export OpSumFromTerms
export SitesFromQNDiag
export TruncateQNDiag
include("easy_sweep.jl")
export SweepOne
export EasySweep
export GetMPOSites
export GetMPO
include("ar_itensor.jl")
export TruncateQnu
export SitesFromQnu

function __init__()
    # Define the space method at runtime
    @eval ITensorMPS.space( :: SiteType"Fermion"; o :: Int, qnd :: Vector{QNDiag}) = [
        QN(
            [ (qndi.name, qndi.charge[o] * n, qndi.modul) for qndi in qnd ]...
        ) => 1 for n = 0 : 1
    ]
end

end 