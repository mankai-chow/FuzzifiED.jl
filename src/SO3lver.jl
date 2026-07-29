module SO3lver

using FuzzifiED
using FuzzifiED.Fuzzifino
using LinearAlgebra
using WignerSymbols
using CGcoefficient
using KrylovKit

function __init__()
    FuzzifiED.ObsMomIncr = true
end

include("core/seg_space.jl")
include("core/sseg_space.jl")
include("core/comp_space.jl")
include("core/couple_decomp.jl")
include("core/seg_operator.jl")
include("core/comp_operator.jl")

end
