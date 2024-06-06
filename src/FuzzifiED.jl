module FuzzifiED 

using LinearAlgebra
using SparseArrays
using Requires
using WignerSymbols
using SphericalHarmonics
using FuzzifiED_jll
import Base.:+
import Base.:-
import Base.:*
import Base.:/
import Base.zero
import Base.adjoint

"""
    SilentStd :: bool = false 

a flag to determine whether logs of the FuzzifiED functions should be turned off. False by default. If you want to evaluate without log, put 

```julia
    SilentStd = True 
```
"""
SilentStd = false
export SilentStd
macro ctrlstd(ex)
    quote
        if SilentStd
            # Save the current stdout and stderr
            original_stdout = stdout
            original_stderr = stderr
            
            # Redirect stdout and stderr to /dev/null
            open("/dev/null", "w") do devnull
                redirect_stdout(devnull)
                redirect_stderr(devnull)
                
                try
                    # Evaluate the expression silently
                    $(esc(ex))
                finally
                    # Restore the original stdout and stderr
                    redirect_stdout(original_stdout)
                    redirect_stderr(original_stderr)
                end
            end
        else
            # Evaluate the expression normally
            $(esc(ex))
        end
    end
end

include("core/confs.jl")
export Confs
export GetConfId

include("core/basis.jl")
export Basis
export GetConfWeight

include("core/term.jl")
export Term
export NormalOrder
export SimplifyTerms

include("core/operator.jl")
export Operator

include("core/opmat.jl")
export OpMat
export GetEigensystem
export SparseMatrixCSCFromOpMat
export MatrixFromOpMat

include("models/threej.jl")
export GetIntMatrix

include("models/l2.jl")
export GetL2Terms

include("models/nn_int.jl")
export GetSnBasis
export GetDenIntTerms
export GetPairIntTerms
export GetPolTerms

include("models/sphere_obs.jl")
export SphereObs
export StoreComps!
export StoreComps
export Laplacian
export GetComponent
export GetPointValue
export Electron
export Density

include("models/ising.jl")
export GetLzQnu
export GetLzZnQnu
export GetLzConfs 
export GetLzZnConfs
export GetIsingBasis
export GetIsingIntTerms
export GetXPolTerms
export GetZPolTerms

include("models/spn.jl")
export GetSpnQnu
export GetSpnConfs
export GetSpnBasis
export GetSpnPairIntTerms
export GetSpnC2Terms

include("models/ising_def.jl")
export GetIsingDefQnu
export GetIsingDefConfs
export GetIsingDefIntTerms
export GetDefXPolTerms

function __init__()

    @require ITensors = "9136182c-28ba-11e9-034c-db9fb085ebd5" begin
        using ITensors
        using ITensors.HDF5

        include("itensors_support/itensors_format.jl")
        export ConfsFromSites
        export TermsFromOpSum
        export OpSumFromTerms
        export SitesFromQnu
        export TruncateQnu

        include("itensors_support/easy_sweep.jl")
        export SweepOne
        export EasySweep
        export GetMpoSites
        export GetMpo
    end
end

end