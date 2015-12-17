module HFB

using Gadfly
using Formatting
using Hafta
#export hfb

"""
`HFBState` stores a Bogoliubov transformation.

The transformation is defined by the `U` and `V` matrices. The object also stores
the `rho` and `kappa` matrices that are calculated from the `U` and `V`.
Calculations should use the `rho` and `kappa` directly and `U`/`V` are stored
only for debugging purposes.

Sometimes it might be that the the `rho`/`kappa` are defined directly, and then
`U`/`V` will be undefined.
"""
type HFBState{T <: Hafta.ManyBodySystem}
    # The many body system related to this state
    system::T

    # Intermediate objects constructed from the transformations
    rho::Matrix{Float64}
    kappa::Matrix{Float64}

    # The defining matrices of the Bogoliubov transformation
    U::Matrix{Float64}
    V::Matrix{Float64}

    """Constructs a `HFBState` from the `U` and `V` matrices."""
    function HFBState(system::T, U::Matrix{Float64}, V::Matrix{Float64})
        N = size(system)
        @assert (N,N) == size(U) && (N,N) == size(V)

        VT = transpose(V)
        rho = V*VT
        kappa = -U*VT
        new(system,rho,kappa,U,V)
    end
end
# TODO: this weird hack is necessary to append to the documentation of a type.
@doc Docs.catdoc((@doc HFBState), Docs.typesummary(HFBState)) HFBState

HFBState{T<:Hafta.ManyBodySystem}(system::T, U, V) = HFBState{T}(system, U, V)

import Base: size
function size(state::HFBState)
    size(state.system)
end

import Base: writemime
function writemime(io, ::MIME"text/html", state::HFBState)
    width,height = 10cm,8cm

    write(io, "<table>")
    write(io, "<tr><th colspan=\"2\" style=\"text-align: center;\">HFBState $(size(state.U))</th></tr>")

    # Table of energies and other values
    E, Ek, Ei, Ep = energy(state)
    Aest = trace(state.rho)
    html = """
    <tr><td colspan="2">
    <table style="width: 100%; border: 1px solid gray;">
    <tr>
        <th>Energy</th>
        <th>Free particles E</th>
        <th>Interaction energy (Γ)</th>
        <th>Pairing energy (Δ)</th>
        <th>Particle number</th>
    </tr>
    <tr>
        <td>{:.5f}</td>
        <td>{:.5f}</td>
        <td>{:.5f}</td>
        <td>{:.5f}</td>
        <td>{:.5f}</td>
    </tr>
    </table>
    """
    write(io,format(html, E, Ek, Ei, Ep, Aest))

    # rho and kappa matrices
    rho,kappa = state.rho, abs(state.kappa)
    minz = min(minimum(rho), minimum(kappa))
    maxz = max(maximum(rho), maximum(kappa))
    scale = Scale.color_continuous(minvalue=minz, maxvalue=maxz)

    write(io, "<tr>")
    write(io, "<td>")
    p = spy(rho, Guide.title("rho matrix"), scale)
    draw(SVG(io,width,height,false), p)
    write(io, "</td>")
    write(io, "<td>")
    p = spy(kappa, Guide.title("kappa matrix"), scale)
    draw(SVG(io,width,height,false), p)
    write(io, "</td>")
    write(io, "</tr>")

    # U and V matrices
    minz = min(minimum(state.U), minimum(state.V))
    maxz = max(maximum(state.U), maximum(state.V))
    scale = Scale.color_continuous(minvalue=minz, maxvalue=maxz)

    write(io, "<tr>")
    write(io, "<td>")
    p = spy(abs(state.U), Guide.title("U matrix"), scale)
    draw(SVG(io,width,height,false), p)
    write(io, "</td>")
    write(io, "<td>")
    p = spy(abs(state.V), Guide.title("V matrix"), scale)
    draw(SVG(io,width,height,false), p)
    write(io, "</td>")
    write(io, "</tr>")

    # the end
    write(io, "</table>")
end


"""
The `HFBIterator` is what stores the state of the iteration.

A `HFBIterator` object is the basis, which then can be iterated to solve
the equations. The object should be constructed with `HFB.hfb`.
"""
type HFBIterator{T <: Hafta.ManyBodySystem}
    # setup
    system::T
    A::Int64
    # iteration variables
    states::Vector{HFBState{T}}
    lambdas::Vector{Float64}
    es::Vector{Float64}
    eigenvalues::Vector{Vector{Float64}}
end
# TODO: this weird hack is necessary to append to the documentation of a type.
@doc Docs.catdoc((@doc HFBIterator), Docs.typesummary(HFBIterator)) HFBIterator

import Base: length
length(hfbi::HFBIterator) = length(hfbi.es)

import Hafta: energy
import Hafta: iterate!
import Hafta: solve!
end
