module FuzzifiEDFullRotation

using FuzzifiED
using LinearAlgebra
using SparseArrays
using WignerSymbols
using CGcoefficient
using Kronecker
using KrylovKit

function __init__()
    BLAS.set_num_threads(1)
    FuzzifiED.ObsMomIncr = true
    FuzzifiED.SilentStd = true
end

include("core/seg_space.jl")
include("core/comp_space.jl")
include("core/seg_operator.jl")
include("core/comp_operator.jl")

include("model/recouple.jl")

end # module FuzzifiEDFullRotation
