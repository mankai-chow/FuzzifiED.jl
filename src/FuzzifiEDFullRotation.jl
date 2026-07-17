module FuzzifiEDFullRotation

using FuzzifiED
using FuzzifiED.Fuzzifino
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
include("core/couple_decomp.jl")
include("core/seg_operator.jl")
include("core/comp_operator.jl")

include("boson/sseg_space.jl")
include("boson/scomp_space.jl")
include("boson/scouple_decomp.jl")
include("boson/sseg_operator.jl")
include("boson/scomp_operator.jl")

end # module FuzzifiEDFullRotation
