module Hafta

using DocStringExtensions

abstract type QuantumBasis end

"""
$(TYPEDEF)

A `ManyBodySystem` is the abstract type that defines a pairwise interacting quantum system.
It should have the following methods:

  - `Base.length(s::T) -> Int`
  - `H0(s::T, i, j) -> Float64`
  - `V(s::T, i, j, k, l) -> Float64`

Additionally, the following methods are deprecated.

  - `Base.size(s::T) -> Int`
"""
abstract type ManyBodySystem end
function H0 end; export H0
function V end; export V

# Functions to calculate observables
# Two possible signatures currently envisioned:
#   - f(::State) -> expectation value for the state
#   - f(::System, i) -> expectation value for the i-th basis state
function energy end; export energy
function spin end; export spin
function parity end; export parity

# Observables for only many-body states
function particle_number end; export particle_number

# Functions of the iterative solvers
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

include("QRPA.jl")
using .QRPA
export QRPAOperator
export arneigs

end
