module ITensorsExt

using FuzzifiED
import ITensors

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
    function ITensors.space( :: SiteType"Fermion"; o :: Int, qnd :: Vector{QNDiag})
        return [
            QN(
                [ (qndi.name, qndi.charge[o] * n, qndi.modul) for qndi in qnd ]...
            ) => 1 for n = 0 : 1
        ]
    end
end

end
