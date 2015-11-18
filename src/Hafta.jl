module Hafta

abstract QuantumBasis

"""
A `ManyBodySystem` is the abstract type that defines an pairwise interacting
quantum system. It should have the following methods:

  - `Base.size(s::T)`
  - `H0(s::T, i, j)`
  - `V(s::T, i, j, k, l)`
"""
abstract ManyBodySystem

function H0 end
export H0

function V end
export V

include("HarmonicOscillator.jl")

end
