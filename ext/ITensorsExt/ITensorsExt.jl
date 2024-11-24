module ITensorsExt

using FuzzifiED
using ITensors 
using ITensorMPS

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