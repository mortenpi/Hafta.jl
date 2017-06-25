# The 2D harmonic oscillator
module HarmonicSystems

using Hafta
import Hafta: H0, V, energy, spin, parity

import Hafta: HarmonicOscillator
import Hafta.HarmonicOscillator: Psi_parity

import Base: size, length, getindex

export Harmonic1DFermionSystem
export Harmonic2DSystem

"""
Number of energy shells in a basis.
"""
function nshells end

# ===================================================================
# The 1D harmonic oscillator
# ===================================================================

"""
Creates a `Hafta.ManyBodySystem` that represents a one-dimensional harmonic
oscillator with a delta-interaction.

The constructor is `Harmonic1DFermionSystem(N,vcoef)`.
"""
mutable struct Harmonic1DFermionSystem <: Hafta.ManyBodySystem
    shells::Int
    wmatrix::Hafta.HarmonicOscillator.WMatrix
    vcoef::Float64

    function Harmonic1DFermionSystem(N,vcoef)
        wmatrix = Hafta.HarmonicOscillator.generate_wmatrix(N)
        new(N, wmatrix, vcoef)
    end
end

# size() is deprecated
function Base.size(s::Harmonic1DFermionSystem)
    2*s.shells
end

function Base.length(s::Harmonic1DFermionSystem)
    2*s.shells
end

nshells(s::Harmonic1DFermionSystem) = s.shells

function H0(s::Harmonic1DFermionSystem, i, j)
    if i==j
        float(div(i-1,2))+0.5
    else
        0.0
    end
end

function V(s::Harmonic1DFermionSystem, i, j, k, l)
    i,is = divrem(i-1, 2)
    j,js = divrem(j-1, 2)
    k,ks = divrem(k-1, 2)
    l,ls = divrem(l-1, 2)
    if is == ks && js == ls
        s.vcoef*s.wmatrix[i,j,k,l]
    else
        0.0
    end
end

spin(::Harmonic1DFermionSystem, i) = (i%2 == 0) ? -.5 : .5
parity(::Harmonic1DFermionSystem, i) = Psi_parity(div(i-1, 2))
energy(::Harmonic1DFermionSystem, i) = float(div(i-1,2))+0.5


# ===================================================================
# The 2D harmonic oscillator
# ===================================================================

# Number of basis functions for N shells
basis_size(N::Integer) = N*(N+1)

mutable struct BasisIdentifier
    nx::UInt16
    ny::UInt16
    s::Bool
end

mutable struct HarmonicBasis
    shells::Int64
    basislookup::Array{BasisIdentifier,1}
    W::HarmonicOscillator.WMatrix
end

function Base.length(b::HarmonicBasis)
    basis_size(b.shells)
end
function Base.getindex(b::HarmonicBasis, idx)
    b.basislookup[idx]
end
function arridx(nx::Integer, ny::Integer, s::Bool)
    shell = nx+ny
    basis_size(shell) + 2*ny + (s?1:0) + 1
end


function harmonicbasis(N)
    # Generate the 2D hosc. basis with N shells
    W = HarmonicOscillator.generate_wmatrix(N)

    basislookup = Array{BasisIdentifier}(basis_size(N))
    idx = 1
    for shell=0:(N-1)
        for ny=0:shell
            nx = shell - ny
            basislookup[idx] = BasisIdentifier(nx,ny,false)
            basislookup[idx+1] = BasisIdentifier(nx,ny,true)
            idx += 2
        end
    end

    HarmonicBasis(N, basislookup, W)
end

function _V(b::HarmonicBasis, i,j,k,l)
    i,j,k,l = b[i],b[j],b[k],b[l]
    #i::BasisIdentifier=b[i]
    #j=b[j]; k=b[k]; l=b[l]
    if (i.s!=k.s) || (j.s!=l.s)
        return 0.0
    end
    Wx = b.W[i.nx,j.nx,k.nx,l.nx]
    Wy = b.W[i.ny,j.ny,k.ny,l.ny]
    float(Wx*Wy)
end

# Create the Hafta.ManyBodySystem type
mutable struct Harmonic2DSystem <: Hafta.ManyBodySystem
    vcoef::Float64
    basis::HarmonicBasis

    function Harmonic2DSystem(N,vcoef)
        basis = harmonicbasis(N)
        new(vcoef, basis)
    end
end

# size() is deprecated
function Base.size(s::Harmonic2DSystem)
    length(s.basis)
end

function Base.length(s::Harmonic2DSystem)
    length(s.basis)
end

nshells(s::Harmonic2DSystem) = s.basis.shells

function H0(s::Harmonic2DSystem, i, j)
    if i==j
        be = s.basis[i]
        be.nx+be.ny+1
    else
        0.0
    end
end

function V(s::Harmonic2DSystem, i, j, k, l)
    s.vcoef*_V(s.basis,i,j,k,l)
end

spin(system::Harmonic2DSystem, i) = system.basis[i].s ? -.5 : .5
parity(system::Harmonic2DSystem, i) = Psi_parity(div(system.basis[i].nx-1, 2)) * Psi_parity(div(system.basis[i].ny-1, 2)) # TODO: does this make sense?
energy(system::Harmonic2DSystem, i) = float(system.basis[i].nx + system.basis[i].ny)+0.5

end
