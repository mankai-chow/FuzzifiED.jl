"""
    mutable struct SphereObs

stores the information of a local observable (or local operator) ``\\Phi`` that can be decomposed into angular components.
```math 
    \\Phi(\\Omega)=\\sum_{lm}\\Phi_{lm}Y^{(s)}_{lm}
```

# Fields

- `s2 :: Int64` is twice the spin ``2s`` of the observable.
- `l2m :: Int64` is twice the maximal angular momentum ``2l_{\\max}`` of the components of the observable. 
- `get_comp :: Function` is a function `get_comp(l2 :: Int64, m2 :: Int64) :: Vector{Term}` that sends the component specified by a tuple of integers ``(2l,2m)`` where ``|s|\\leq l\\leq l_{\\max}, -l\\leq m\\leq l`` to a list of terms that specifies the expression of the component. 
- `stored_q :: Bool` is a boolean that specifies whether or not each component of the observable is stored.
- `comps :: Dict{Tuple{Int64, Int64}, Vector{Term}}` stores each component of the observable in the format of a dictionary whose keys are the tuples of integers ``(2l,2m)`` and values are the lists of terms that specifies the expression of the component. 
"""
mutable struct SphereObs
    s2 :: Int64
    l2m :: Int64
    get_comp :: Function
    stored_q :: Bool 
    comps :: Dict{Tuple{Int64, Int64}, Vector{Term}}
end


"""
    function SphereObs(s2 :: Int64, l2m :: Int64, get_comp :: Function) :: SphereObs

initialises the observable from ``2s``, ``2l_{\\max}`` and the function ``(l,m)\\mapsto\\Phi_{lm}``

# Arguments

- `s2 :: Int64` is twice the spin ``2s`` of the observable.
- `l2m :: Int64` is twice the maximal angular momentum ``2l_{\\max}`` of the components of the observable. 
- `get_comp :: Function` is a function `get_comp(l2 :: Int64, m2 :: Int64) :: Vector{Term}` that sends the component specified by a tuple of integers ``(2l,2m)`` where ``|s|\\leq s\\leq l_{\\max}, -l\\leq m\\leq l`` to a list of terms that specifies the expression of the component. 
"""
function SphereObs(s2 :: Int64, l2m :: Int64, get_comp :: Function)
    return SphereObs(s2, l2m, get_comp, false, Dict{Tuple{Int64, Int64}, Vector{Term}}())
end

"""
    function SphereObs(s2 :: Int64, l2m :: Int64, get_comp :: Function) :: SphereObs

initialises the observable from ``2s``, ``2l_{\\max}`` and a list of ``\\Phi_{lm}`` specified by a dictionary. 

# Arguments

- `s2 :: Int64` is twice the spin ``2s`` of the observable.
- `l2m :: Int64` is twice the maximal angular momentum ``2l_{\\max}`` of the components of the observable. 
- `comps :: Dict{Tuple{Int64, Int64}, Vector{Term}}` stores each component of the observable in the format of a dictionary whose keys are the tuples of integers ``(2l,2m)`` and values are the lists of terms that specifies the expression of the component. 
"""
function SphereObs(s2 :: Int64, l2m :: Int64, cmps :: Dict{Tuple{Int64, Int64}, Vector{Term}})
    return SphereObs(s2, l2m, (l2, m2) -> (l2 ≤ l2m && l2 ≥ abs(s2) && abs(m2) ≤ l2 && haskey(cmps, (l2, m2))) ? cmps[(l2, m2)] : Term[], true, cmps)
end


"""
    function StoreComps!(obs :: SphereObs) :: Nothing 
    
calculates and stores each component of the observable `obs` and replace the function in `obs` by the list of calculated components. 
"""
function StoreComps!(obs :: SphereObs)
    if (obs.stored_q) return end
    cmps = Dict{Tuple{Int64, Int64}, Vector{Term}}()
    s2 = obs.s2
    l2m = obs.l2m
    for l2 = abs(s2) : 2 : l2m
        for m2 = -l2 : 2 : l2 
            cmps[(l2, m2)] = SimplifyTerms(obs.get_comp(l2, m2))
        end
    end
    obs.stored_q = true 
    obs.comps = cmps
    obs.get_comp = (l2, m2) -> (l2 ≤ l2m && l2 ≥ abs(s2) && abs(m2) ≤ l2 && haskey(obs.comps, (l2, m2))) ? obs.comps[(l2, m2)] : Term[]
end


"""
    function StoreComps(obs :: SphereObs) :: SphereObs
    
calculates and stores each component of the observable `obs` and return a new observable with the list of calculated components. 
"""
function StoreComps(obs :: SphereObs)
    if (obs.stored_q) return obs end
    cmps = Dict{Tuple{Int64, Int64}, Vector{Term}}()
    s2 = obs.s2
    l2m = obs.l2m
    for l2 = abs(s2) : 2 : l2m
        for m2 = -l2 : 2 : l2 
            cmps[(l2, m2)] = obs.get_comp(l2, m2)
        end
    end
    obs.stored_q = true 
    obs.comps = cmps
    return SphereObs(s2, l2m, cmps)
end


"""
    function *(fac :: Number, obs :: SphereObs) :: SphereObs
    function *(obs :: SphereObs, fac :: Number) :: SphereObs
    function /(obs :: SphereObs, fac :: Number) :: SphereObs
    function -(obs :: SphereObs) :: SphereObs
    
enables the multiplication of an observable with a number.
"""
function *(fac :: Number, obs :: SphereObs) 
    return SphereObs(obs.s2, obs.l2m, (l2, m2) -> fac * obs.get_comp(l2, m2))
end
function *(obs :: SphereObs, fac :: Number) 
    return fac * obs
end
function /(obs :: SphereObs, fac :: Number) 
    return (1 / fac) * obs
end
function -(obs :: SphereObs) 
    return (-1) * obs
end


"""
    function +(obs1 :: SphereObs, obs2 :: SphereObs) :: SphereObs
    function -(obs1 :: SphereObs, obs2 :: SphereObs) :: SphereObs
    
enables the addition of two observables.
"""
function +(obs1 :: SphereObs, obs2 :: SphereObs) 
    if (obs1.s2 ≠ obs2.s2) 
        print("Additions must have equal S")
        return 
    end
    s2 = obs1.s2
    l2m = max(obs1.l2m, obs2.l2m)
    return SphereObs(s2, l2m, (l2, m2) -> obs1.get_comp(l2, m2) + obs2.get_comp(l2, m2))
end
function -(obs1 :: SphereObs, obs2 :: SphereObs)
    return obs1 + (-1) * obs2
end


"""
    function adjoint(obs :: SphereObs) :: SphereObs
    
enables the Hermitian conjugate of an observable.
```math
\\begin{aligned}
    Φ^†(Ω)&=∑_{lm}(Φ_{lm})^†\\bar{Y}^{(s)}_{lm}(Ω)=∑_{lm}(Φ_{lm})^†(-1)^{s+m)Y^{(-s)}_{l,-m}(Ω)\\\\
    (Φ^†)_{lm}&=(-1)^{s-m}(Φ_{l,-m})^†
\\end{aligned}
```
"""
function adjoint(obs :: SphereObs)
    s2 = obs.s2
    l2m = obs.l2m
    obs1 = SphereObs(-s2, l2m, (l2, m2) -> obs.get_comp(l2, -m2)' * (iseven((s2 - m2) ÷ 2) ? 1 : -1))
    return obs1
end


"""
    function *(obs1 :: SphereObs, obs2 :: SphereObs) :: SphereObs
    
enables the multiplication of two observable by making use of the composition of two monopole harmonics into one. 
"""
function *(obs1 :: SphereObs, obs2 :: SphereObs)
    s21 = obs1.s2 
    s22 = obs2.s2
    s2 = s21 + s22
    l2m1 = obs1.l2m
    l2m2 = obs2.l2m
    l2m = l2m1 + l2m2
    gc = ((l2, m2) -> sum(Vector{Term}[sum(Vector{Term}[sum(Vector{Term}[
            (iseven(s2 + m2) ? 1 : -1) *
            sqrt((l21 + 1) * (l22 + 1) * (l2 + 1) / (4 * π)) *
            wigner3j(l21/2, l22/2, l2/2, -s21/2, -s22/2, s2/2) *
            wigner3j(l21/2, l22/2, l2/2, m21/2, (m2-m21)/2, -m2/2) *
            obs1.get_comp(l21, m21) * obs2.get_comp(l22, m2 - m21)
        for m21 = max(-l21, -l22 + m2) : 2 : min(l21, l22 + m2)])
        for l21 = max(abs(s21), abs(l2 - l22)) : 2 : min(l2m1, l2 + l22)])
        for l22 = abs(s22) : 2 : l2m2]))
    return SphereObs(s2, l2m, gc)
end


"""
    function Laplacian(obs :: SphereObs) :: SphereObs
    
Takes the Laplacian of an observable
```math
    (∇^2Φ)_{lm}=-l(l+1)Φ_{lm}
```
"""
function Laplacian(obs :: SphereObs)
    return SphereObs(obs.s2, obs.l2m, (l2, m2) -> -l2 / 2 * (l2 / 2 + 1) * obs.get_comp(l2, m2))
end


"""
    function GetComponent(obs :: SphereObs, l :: Number, m :: Number) :: Vector{Term}

returns an angular component ``Φ_{lm}`` of an observable in the format of a list of terms.
"""
function GetComponent(obs :: SphereObs, l :: Number, m :: Number)
    return obs.get_comp(Int64(2 * l), Int64(2 * m))
end


"""
    function GetPointValue(obs :: SphereObs, θ :: Float64, ϕ :: Float64) :: Vector{Term}

evaluates an observable at one point ``Φ(θ,ϕ)`` in the format of a list of terms.
"""
function GetPointValue(obs :: SphereObs, θ :: Float64, ϕ :: Float64)
    if (obs.s2 ≠ 0) 
        println("S ≠ 0 not supported")
        return
    end
    lm = obs.l2m ÷ 2
    Ylm = computeYlm(θ, ϕ, lmax = lm)
    tms = Term[]
    for l = 0 : lm 
        for m = -l : l 
            tms += obs.get_comp(2 * l, 2 * m) * Ylm[(l,m)]
        end
    end
    return tms
end

"""
    function Electron(nm :: Int64, nf :: Int64, f :: Int64) :: SphereObs

returns the electron annihilation operator ``ψ_f``

# Arguments

- ``nf :: Int64`` is the number of flavours.
- ``nm :: Int64`` is the number of orbitals.
- ``f :: Int64`` is the index of the orbital to be taken.
"""
function Electron(nm :: Int64, nf :: Int64, f :: Int64)
    gc = (l2, m2) -> (l2 == nm - 1) ? [Term(1.0, [0, f + nf * ((m2 + nm - 1) ÷ 2)])] : Term[]
    return SphereObs(nm - 1, nm - 1, gc)
end



"""
    function Density(nm :: Int64, nf :: Int64 ; mat :: Matrix{<:Number}) :: SphereObs

returns the density operator ``n=∑_{ff'}ψ^†_{f}M_{ff'}ψ_{f'}``

# Arguments

- ``nf :: Int64`` is the number of flavours.
- ``nm :: Int64`` is the number of orbitals.
- ``mat :: Int64`` is the matrix ``M_{ff'}``. Facultive, identity matrix ``\\mathbb{I}`` by default.
"""
function Density(nm :: Int64, nf :: Int64 ; mat :: Matrix{<:Number} = Matrix{Float64}(I, nf, nf))
    el = [ StoreComps(Electron(nm, nf, f)) for f = 1 : nf ]
    obs = SphereObs(0, 0, Dict{Tuple{Int64, Int64}, Vector{Term}}())
    for f1 = 1 : nf, f2 = 1 : nf
        if abs(mat[f1, f2]) < 1E-13 continue end 
        obs += mat[f1, f2] * el'[f1] * el[f2]
    end
    return obs
end