module HFB

using Gadfly
using Formatting

using Hafta
import Hafta: ManyBodySystem
import Hafta.Utils: find_value

export hfb
export particle_number

"""
`HFBSystem{T}` is a wrapper around a `ManyBodySystem`.

It provides some temporary variables useful for the HFB solver
that are constructed from the systems `H0(i,j)` and `V(i,j,k,l)`
functions.

The fields are:

    - `system::T` -- the many body system related to this state
    - `Tij::Matrix{Float64}` -- a matrix form of the `H0(i,j)`

The automatically provided methods are `H0`, `V` and `length`,
and they just dispatch to the corresponding methods for the type `T`.

An additional method provided is the `Vbar(::HFBSystem, i,j,k,l)`,
which is the antisymmetrized version of `V`, which makes expressions
for HFB simpler.
"""
type HFBSystem{T <: ManyBodySystem}
    system::T
    Tij::Matrix{Float64}

    function HFBSystem(s::T)
        M = length(s)
        Tij = zeros(Float64, (M,M))
        for i=1:M, j=1:M
            Tij[i,j] = H0(s, i,j)
        end
        new(s,Tij)
    end
end
HFBSystem{T}(s::T) = HFBSystem{T}(s)

import Hafta: H0, V
import Base: length
H0(s::HFBSystem, i, j) = H0(s.system, i,j)
V(s::HFBSystem, i, j, k, l) = V(s.system, i,j,k,l)
length(s::HFBSystem) = length(s.system) # TODO: should it be "static" as well?

Vbar(s::HFBSystem, i,j,k,l) = V(s.system, i,j,k,l) - V(s.system, i,j,l,k)

"""
`HFBState` stores a Bogoliubov transformation.

The transformation is defined by the `U` and `V` matrices. The object also stores
the `rho` and `kappa` matrices that are calculated from the `U` and `V`.
Calculations should use the `rho` and `kappa` directly and `U`/`V` are stored
only for debugging purposes.

Sometimes it might be that the the `rho`/`kappa` are defined directly, and then
`U`/`V` will be undefined.
"""
type HFBState{T <: ManyBodySystem}
    # The many body system related to this state
    system::HFBSystem{T}
    lambda::Float64
    energies::Vector{Float64}

    # potentials
    gamma::Matrix{Float64}
    delta::Matrix{Float64}

    # Intermediate objects constructed from the transformations
    rho::Matrix{Float64}
    kappa::Matrix{Float64}

    # The defining matrices of the Bogoliubov transformation
    U::Matrix{Float64}
    V::Matrix{Float64}

    """Constructs a `HFBState` from the `U` and `V` matrices."""
    function HFBState(system::HFBSystem{T}, lambda, gamma, delta, energies, U, V)
        M = length(system)
        @assert (M,M) == size(U) && (M,M) == size(V)
        @assert (M,M) == size(gamma) && (M,M) == size(delta)
        rho = conj(V)*transpose(V)
        kappa = -U*ctranspose(V)
        new(system, lambda, energies, gamma, delta, rho, kappa, U, V)
    end

    """Constructs a `HFBState` from the `U` and `V` matrices."""
    function HFBState(system::HFBSystem{T}, lambda, energies, U, V)
        M = length(system)
        @assert (M,M) == size(U) && (M,M) == size(V)
        rho = conj(V)*transpose(V)
        kappa = -U*ctranspose(V)
        gamma, delta = gamma_delta(system, rho, kappa)
        new(system, lambda, energies, gamma, delta, rho, kappa, U, V)
    end
end

HFBState{T<:ManyBodySystem}(system::HFBSystem{T}, lambda, energies, U, V) = HFBState{T}(system, lambda, energies, U, V)
HFBState{T<:ManyBodySystem}(system::T, lambda, energies, U, V) = HFBState{T}(HFBSystem(system), lambda, energies, U, V)

import Base: size
size(state::HFBState) = size(state.system)
particle_number(s::HFBState) = trace(s.rho)

"""
The `HFBIterator` is what stores the state of the iteration.

A `HFBIterator` object is the basis, which then can be iterated to solve
the equations. The object should be constructed with `HFB.hfb`.
"""
type HFBIterator{T <: ManyBodySystem}
    # setup
    system::HFBSystem{T}
    A::Int64
    # iteration variables
    states::Vector{HFBState{T}}
    es::Vector{Float64}
    eigenvalues::Vector{Vector{Float64}}
end

import Base: length
length(hfbi::HFBIterator) = length(hfbi.es)


"""
`hfb(system, A; maxkappa)` constructs a `HFBIterator` object.

Arguments:

    - `system` is the quantum many-body system (`<: ManyBodySystem`)
    - `A` is the number of particles
    - `maxkappa` -- ???
"""
function hfb(system, A; maxkappa=1)
    M = length(system)
    hfbi = HFBIterator{typeof(system)}(HFBSystem(system), A, [],[],[])

    state = HFBState(system, NaN, Vector{Float64}(), zeros(Float64, (M,M)), zeros(Float64, (M,M)))
    for i=1:A
        state.rho[i,i] = 1.0
    end

    for d=2:1+maxkappa, i=1:div(M,d)
        m = d*(i-1)+1
        n = m+d-1
        state.kappa[m,n] = 0.2
        state.kappa[n,m] = -0.2
    end

    push!(hfbi.states, state)
    push!(hfbi.es, energy(state)[1])

    hfbi
end


"""
`gamma_delta(system, rho, kappa)` calculates the `gamma` and `delta`
matrices from the `rho` and `kappa`. It also needs a system, since
the `gamma` and `delta` also include the interaction `V(i,j,k,l)`.
"""
function gamma_delta(system::HFBSystem, rho::Matrix, kappa::Matrix)
    M = length(system)
    gamma = zeros(Float64, (M,M))
    delta = zeros(Float64, (M,M))
    for i=1:M,j=1:M,k=1:M,l=1:M
        gamma[i,j] += rho[k,l] * Vbar(system, i,k,j,l)
        delta[i,j] += 0.5*kappa[k,l] * Vbar(system, i,j,k,l)
    end
    gamma,delta
end

"""`gamma_delta(::HFBState)` is a convenience wrapper it directly from a `HFBState`"""
gamma_delta(state::HFBState) = gamma_delta(state.system, state.rho, state.kappa)

import Hafta: energy
function energy(state::HFBState)
    N = length(state.system)
    gamma, delta = state.gamma, state.delta

    T = zeros(Float64, (N,N))
    for i=1:N, j=1:N
        T[i,j] = H0(state.system, i,j)
    end

    Ef = trace(T*state.rho)
    Ei = 0.5*trace(gamma*state.rho)
    Ep = -0.5*trace(delta*state.kappa)

    Ef+Ei+Ep, Ef, Ei, Ep
end

function solve_state(system,T,gamma,delta,lambda)
    M = length(system)
    h = T + gamma - lambda*eye(M)
    equation = zeros(Float64, (2*M, 2*M))
    equation[1:M, 1:M] = h
    equation[M+1:2M, M+1:2M] = -h
    equation[1:M, M+1:2M] = delta
    equation[M+1:2M, 1:M] = -delta

    if !ishermitian(equation)
        maxdiff = maximum(abs(equation-transpose(equation)))
        if maxdiff > 1e-14
            warn("Equation not hermitian (|diff| = $maxdiff)")
        end
        equation = 0.5*(equation+transpose(equation))
    end

    efact = eigfact(equation)

    perms = sortperm(efact[:values])[M+1:2M]
    energies = efact[:values][perms]
    U = efact[:vectors][1:M, perms]
    V = efact[:vectors][M+1:2M, perms]

    state = HFBState(system,lambda,energies,U,V)
    state, trace(state.rho), efact
end

# deprecated solve_state
solve_state(system,N,lambda,T,gamma,delta) = solve_state(system,T,gamma,delta,lambda)

"""
`iterate_lambda` uses binary search to find the next solution.

The underlying assumption is that the particle number of the HFB solution
is a monotonically increasing function of lambda.
"""
function iterate_lambda(system::HFBSystem,rho,kappa,A; lambdaepsilon=1e-12, nepsilon=lambdaepsilon, maxiters=100, verbose=false)
    if verbose @show lambdaepsilon nepsilon end
    M = length(system)
    gamma, delta = gamma_delta(system, rho, kappa)

    function f(λ)
        _,n = solve_state(system,system.Tij,gamma,delta,λ)
        n
    end

    λmin,λmax = find_value(f,A; xeps=lambdaepsilon, yeps=nepsilon, maxiters=maxiters, verbose=verbose)

    # for the choice of lambda see:
    #   https://lund.mortenpi.eu/hafta/2016/02/17/hfb-lambda-fix.html
    nextstate,_ = solve_state(system,system.Tij,gamma,delta,0.5*(λmin+λmax))
    nextstate
end

import Hafta: iterate!
function iterate!(hfbi::HFBIterator; mixing=0.0, maxiters=100, nepsilon=1e-10, lambdaepsilon=nepsilon/1e3, verbose=false)
    lambdaepsilon = min(nepsilon, lambdaepsilon)

    if !(0.0 <= mixing < 1.0)
        error("Invalid value for mixing in iterate!() ($mixing). Must be 0.0 <= mixing < 1.0.")
    end

    A,N = hfbi.A, length(hfbi.system)

    T = zeros(Float64, (N,N))
    for i=1:N, j=1:N
        T[i,j] = H0(hfbi.system, i,j)
    end

    rho,kappa = if mixing != 0.0 && length(hfbi.states) > 1
        state = hfbi.states[end]
        oldstate = hfbi.states[end-1]
        rho = (1.0 - mixing)*state.rho + mixing*oldstate.rho
        kappa = (1.0 - mixing)*state.kappa + mixing*oldstate.kappa
        rho, kappa
    else
        state = hfbi.states[end]
        state.rho, state.kappa
    end

    nextstate = iterate_lambda(hfbi.system,rho,kappa,A; lambdaepsilon=lambdaepsilon,nepsilon=nepsilon,maxiters=maxiters,verbose=verbose)

    E,_ = energy(nextstate)
    push!(hfbi.states, nextstate)
    push!(hfbi.es, E)

    #@show E
    hfbi
end

function issolved(hfbi::HFBIterator, epsilon)
    const mindeltas = 5
    if length(hfbi.es) < mindeltas+1
        false
    else
        maximum(abs(diff(hfbi.es[end-mindeltas:end]))) < epsilon
    end
end

import Hafta: solve!
function solve!(hfbi::HFBIterator; epsilon=1e-10, maxiters=20, lambdaiters=50, args...)
    #args = Dict{Symbol, Any}(args)
    efact = nothing
    while !issolved(hfbi, epsilon)
        iterate!(hfbi; maxiters=lambdaiters, nepsilon=epsilon/10, args...)

        maxiters -= 1
        if maxiters == 0
            error("Max iterations reached")
            return nothing
        end
    end
    hfbi.es[end]
end


# writemime functions and other decorative things
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
        <th>Chem.pot.</th>
    </tr>
    <tr>
        <td>{:.5f}</td>
        <td>{:.5f}</td>
        <td>{:.5f}</td>
        <td>{:.5f}</td>
        <td>{:.5f}</td>
        <td>{:.5f}</td>
    </tr>
    </table>
    """
    write(io,format(html, E, Ek, Ei, Ep, Aest, state.lambda))

    # rho and kappa matrices
    rho,kappa = abs(state.rho), abs(state.kappa)
    maxz = max(maximum(rho), maximum(kappa))
    scale = Scale.color_continuous(minvalue=0, maxvalue=maxz)

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
    U,V = abs(state.U), abs(state.V)
    maxz = max(maximum(U), maximum(V))
    scale = Scale.color_continuous(minvalue=0, maxvalue=maxz)

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

function writemime(io, ::MIME"text/html", hfbi::HFBIterator)
    width,height = 20cm,6cm

    write(io, "<table>")
    write(io, """
    <tr>
        <th style="text-align: center;">
            HFBIterator ($(length(hfbi)) iterations)
        </th>
    </tr>""")

    write(io, "<tr><td>")
    p = plot(
        x=1:length(hfbi.es), y=hfbi.es,
        Geom.line,# Geom.point,
        Guide.title("Absolute value of energy"),
        Guide.xlabel("Iteration"), Guide.ylabel("E")
    )
    draw(SVG(io,width, height,false), p)
    #write(io, "</td></tr>")
    write(io, "<br /")

    #write(io, "<tr><td>")
    logdiffs = log10(abs(diff(hfbi.es)))
    p = plot(
        x=1:length(logdiffs), y=logdiffs,
        yintercept=[-10, -15], Geom.hline(color="red"),
        Geom.line,# Geom.point,
        Guide.title("Differences of consecutive energies"),
        Guide.xlabel("Iteration"), Guide.ylabel("log10(ΔE)")
    )
    draw(SVG(io,width, height,false), p)
    write(io, "</td></tr>")

    # Output the first and the last state
    if length(hfbi.states) > 1
        write(io, """
        <tr>
        <th style="text-align: center; border-top: 3px solid black;">
        Final state
            </th>
        </tr>""")
        write(io, "<tr><td>")
        writemime(io, "text/html", hfbi.states[end])
        write(io, "</td></tr>")
    end

    if false
        write(io, """
        <tr>
            <th style="text-align: center; border-top: 3px solid black;">
                Initial state
            </th>
        </tr>""")
        write(io, "<tr><td>")
        writemime(io, "text/html", hfbi.states[1])
        write(io, "</td></tr>")
    end

    write(io, "</table>")
end

end
