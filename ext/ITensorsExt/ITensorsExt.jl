module ITensorsExt

using ITensors 
using ITensorMPS
using FuzzifiED

import FuzzifiED.QNDiagFromSites
import FuzzifiED.ConfsFromSites
import FuzzifiED.TermsFromOpSum
import FuzzifiED.OpSumFromTerms
import FuzzifiED.SitesFromQNDiag
import FuzzifiED.TruncateQNDiag
import FuzzifiED.SweepOne
import FuzzifiED.EasySweep
import FuzzifiED.GetMPOSites
import FuzzifiED.GetMPO
import FuzzifiED.TruncateQnu
import FuzzifiED.SitesFromQnu

include("itensors_format.jl")
include("easy_sweep.jl")
include("ar_itensor.jl")

function __init__()
    # Define the space method at runtime
    @eval ITensorMPS.space( :: SiteType"Fermion"; o :: Int, qnd :: Vector{QNDiag}) = [
        QN(
            [ (qndi.name, qndi.charge[o] * n, qndi.modul) for qndi in qnd ]...
        ) => 1 for n = 0 : 1
    ]
end

end 