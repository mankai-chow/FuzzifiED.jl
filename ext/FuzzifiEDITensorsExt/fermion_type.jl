
"""
    ITensors.space( :: SiteType"FuzzyFermion" ; o :: Int, qnd :: Vector{QNDiag})

Define a new site type "FuzzyFermion" which inherits all the features of ITensor type "Fermion". It can be initialised by a set of QNDiag's and the site index ``o``.
"""
function ITensors.space( :: SiteType"FuzzyFermion"; o :: Int, qnd :: Vector{QNDiag})
    return [
        QN(
            [ (qndi.name, qndi.charge[o] * n, qndi.modul) for qndi in qnd ]...
        ) => 1 for n = 0 : 1
    ]
end

ITensors.val(::ValName"Emp", ::SiteType"FuzzyFermion") = 1
ITensors.val(::ValName"Occ", ::SiteType"FuzzyFermion") = 2
ITensors.val(::ValName"0", st::SiteType"FuzzyFermion") = val(ValName("Emp"), st)
ITensors.val(::ValName"1", st::SiteType"FuzzyFermion") = val(ValName("Occ"), st)

ITensors.state(::StateName"Emp", ::SiteType"FuzzyFermion") = [1.0 0.0]
ITensors.state(::StateName"Occ", ::SiteType"FuzzyFermion") = [0.0 1.0]
ITensors.state(::StateName"0", st::SiteType"FuzzyFermion") = state(StateName("Emp"), st)
ITensors.state(::StateName"1", st::SiteType"FuzzyFermion") = state(StateName("Occ"), st)

function ITensors.op!(Op::ITensor, ::OpName"N", ::SiteType"FuzzyFermion", s::Index)
  return Op[s' => 2, s => 2] = 1.0
end
function ITensors.op!(Op::ITensor, on::OpName"n", st::SiteType"FuzzyFermion", s::Index)
  return op!(Op, alias(on), st, s)
end

function ITensors.op!(Op::ITensor, ::OpName"C", ::SiteType"FuzzyFermion", s::Index)
  return Op[s' => 1, s => 2] = 1.0
end
function ITensors.op!(Op::ITensor, on::OpName"c", st::SiteType"FuzzyFermion", s::Index)
  return op!(Op, alias(on), st, s)
end

function ITensors.op!(Op::ITensor, ::OpName"Cdag", ::SiteType"FuzzyFermion", s::Index)
  return Op[s' => 2, s => 1] = 1.0
end
function ITensors.op!(Op::ITensor, on::OpName"c†", st::SiteType"FuzzyFermion", s::Index)
  return op!(Op, alias(on), st, s)
end
function ITensors.op!(Op::ITensor, on::OpName"cdag", st::SiteType"FuzzyFermion", s::Index)
  return op!(Op, alias(on), st, s)
end

function ITensors.op!(Op::ITensor, ::OpName"F", ::SiteType"FuzzyFermion", s::Index)
  Op[s' => 1, s => 1] = +1.0
  return Op[s' => 2, s => 2] = -1.0
end

ITensors.has_fermion_string(::OpName"C", ::SiteType"FuzzyFermion") = true
function ITensors.has_fermion_string(on::OpName"c", st::SiteType"FuzzyFermion")
  return has_fermion_string(alias(on), st)
end
ITensors.has_fermion_string(::OpName"Cdag", ::SiteType"FuzzyFermion") = true
function ITensors.has_fermion_string(on::OpName"c†", st::SiteType"FuzzyFermion")
  return has_fermion_string(alias(on), st)
end
function ITensors.has_fermion_string(on::OpName"cdag", st::SiteType"FuzzyFermion")
  return has_fermion_string(alias(on), st)
end