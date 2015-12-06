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
function H0 end; export H0
function V end; export V

# Methods of the iterative solvers
function energy end; export energy
function iterate! end; export iterate!
function solve! end; export solve!

include("HarmonicOscillator.jl")

include("HarmonicSystems.jl")
using .HarmonicSystems
export Harmonic1DFermionSystem
export Harmonic2DSystem

include("HartreeFock.jl")
using .HartreeFock
export hartreefock

end
