export Fuzzifino

module Fuzzifino

using FuzzifiED_jll
using FuzzifiED
using LinearAlgebra

import FuzzifiED: NumThreads, SilentStd, ElementType

"""
    FuzzifiED.Fuzzifino.Libpathino :: String = FuzzifiED_jll.LibpathFuzzifino

define path of the Fortran library `libfuzzifino.so`. You do not need to modify that by yourself. However, if you compile the Fortran codes by yourself, you need to point this to your compiled library. 
"""
Libpathino :: String = FuzzifiED_jll.LibpathFuzzifino

include("sqn.jl")
include("sconfs.jl")
include("sbasis.jl")
include("sterm.jl")
include("soperator.jl")
include("sopmat.jl")
include("stransf.jl")
include("sentangle.jl")

end
