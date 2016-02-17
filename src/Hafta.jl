module Hafta

abstract QuantumBasis

"""
A `ManyBodySystem` is the abstract type that defines an pairwise interacting
quantum system. It should have the following methods:

  - `Base.length(s::T) -> Int`
  - `H0(s::T, i, j) -> Float64`
  - `V(s::T, i, j, k, l) -> Float64`

Additionally, the following methods are deprecated.

  - `Base.size(s::T) -> Int`
"""
abstract ManyBodySystem
function H0 end; export H0
function V end; export V

# Methods of the iterative solvers
function energy end; export energy
function iterate! end; export iterate!
function solve! end; export solve!

include("Utils.jl")

include("HarmonicOscillator.jl")
include("HarmonicSystems.jl")
using .HarmonicSystems
export Harmonic1DFermionSystem
export Harmonic2DSystem

include("HartreeFock.jl")
using .HartreeFock
export hartreefock

include("HFB.jl")
using .HFB
export hfb

end
