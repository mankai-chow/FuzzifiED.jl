module Fuzzifino
    using FuzzifiED_jll
    using FuzzifiED
    using LinearAlgebra
    
    import Base.:+
    import Base.:-
    import Base.:*
    import Base.:/
    import Base.:÷
    import Base.:^
    import Base.zero
    import Base.adjoint
    import FuzzifiED.NormalOrder
    import FuzzifiED.SimplifyTerms
    import FuzzifiED.OpMat
    import FuzzifiED.NumThreads
    import FuzzifiED.SilentStd

    """
        Libpathino :: String = FuzzifiED_jll.LibpathFuzzifino

    define where the Fortran library are compiled. You do not need to modify that by yourself. However, if you compile the Fortran codes by yourself, you need to point this to your compiled library. 
    """
    Libpathino :: String = FuzzifiED_jll.LibpathFuzzifino
    export Libpathino

    include("sqn.jl")
    export SQNDiag
    export SQNOffd

    include("sconfs.jl")
    export SConfs

    include("sbasis.jl")
    export SBasis

    include("sterm.jl")
    export STerm

    include("soperator.jl")
    export SOperator

    include("sopmat.jl")

end