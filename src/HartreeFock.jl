module HartreeFock

using Hafta
import Hafta: energy, iterate!, solve!
import Base: size

export hartreefock

# Hartree-Fock iteration
"""
The `HartreeFockIterator` is what stores the state of the iteration.

A `HartreeFockIterator` object is the basis, which then can be iterated to solve
the equations. The object should be constructed with `HartreeFock.hartreefock`.
"""
type HartreeFockIterator{T <: Hafta.ManyBodySystem}
    # setup
    system::T
    A::Int64
    # iteration variables
    phis::Vector{Matrix{Float64}}
    es::Vector{Float64}
    eigenvalues::Vector{Vector{Float64}}
end

"""
`hartreefock(system, A)` constructs a `HartreeFockIterator` object.

Arguments:

- `system` is the 2D system (`<: ManyBodySystem`)
- `A` is the number of particles
"""
function hartreefock(system, A)
    basis_size = size(system)
    hf = HartreeFockIterator{typeof(system)}(system, A, [],[],[])

    phiguess = eye(A,basis_size)
    push!(hf.phis, phiguess)
    push!(hf.es, energy(hf, phiguess))

    hf
end

function Base.size(hf::HartreeFockIterator)
    hf.A, size(hf.system)
end

function energy(hf::HartreeFockIterator, phis::Matrix)
    A,N = size(hf)
    #@show A,N

    rho = zeros(Float64, (N,N))
    for i=1:N, j=1:N, k=1:A
        rho[i,j] += phis[k,i]*phis[k,j]
    end

    # TODO: This part should be optimized for the case where H0
    # is diagonal (i.e. the basis we use to expand is the eigenbasis
    # of the single particle part of the hamiltonian).
    # In that case we can get rid of one of the indices to loop over.
    Et = 0.0
    for i=1:N, j=1:N
        Et += H0(hf.system, i, j)*rho[i,j]
    end


    E1, E2 = 0.0, 0.0
    for i=1:N, j=1:N, k=1:N, l=1:N
        Vijkl = V(hf.system,i,j,k,l)
        E1 += Vijkl * rho[i,k] * rho[j,l]
        E2 += Vijkl * rho[i,l] * rho[j,k]
    end
    Et + 0.5*(E1 - E2)
end

function iterate!(hf::HartreeFockIterator, mixing=0.0)
    if !(0.0 <= mixing < 1.0)
        throw(ArgumentError("Invalid value for mixing in iterate!() ($mixing). Must be 0.0 <= mixing < 1.0."))
    end

    A,N = size(hf)
    phis = hf.phis[end]

    # We precalculate the density matrix
    rho = zeros(Float64, (N,N))
    for i=1:N, j=1:N, k=1:A
        rho[i,j] += phis[k,i]*phis[k,j]
    end

    # Mix in the states from the last iteration
    if mixing != 0.0 && length(hf.phis) > 1
        oldphis = hf.phis[end-1]
        oldrho = zeros(Float64, (N,N))
        for i=1:N, j=1:N, k=1:A
            oldrho[i,j] += oldphis[k,i]*oldphis[k,j]
        end
        rho = (1.0 - mixing)*rho + mixing*oldrho
    end

    # The final operators: T, VH and VF
    T = zeros(Float64, (N,N))
    VH = zeros(Float64, (N,N))
    VF = zeros(Float64, (N,N))

    for i=1:N, j=1:N
        T[i,j] = H0(hf.system, i, j)
    end

    for j=1:N, i=1:N
        for m=1:N, n=1:N
            VH[i,j] += rho[m,n] * V(hf.system, i,m,j,n)
            VF[i,j] += (-1) * rho[m,n] * V(hf.system, i,m,n,j)
        end
    end

    # Solve the eigenvalue problem
    if !ishermitian(T)
        @show maximum(VF-transpose(VF))
    end
    if !ishermitian(VH)
        @show maximum(VF-transpose(VF))
    end
    if !ishermitian(VF)
        deviation = maximum(VF-transpose(VF))
        #warn("VF not hermitian: maxdev $deviation. Fixing..")
        VF = 0.5*(VF+transpose(VF))
    end
    #@show T+VH+VF
    eigs = eigfact(T+VH+VF)

    # pick the first A eigenvalues (lowest eigenvalues/energies)
    eigphis = zeros(phis)
    perms = sortperm(eigs[:values])
    for i=1:A
        eigphis[i,:] = eigs[:vectors][:,perms[i]]
    end

    # Calculate the energy
    e = energy(hf, eigphis)

    # Push the generated values
    push!(hf.phis, eigphis)
    push!(hf.es, e)
    push!(hf.eigenvalues, eigs[:values])

    # Return the matrices
    (e,T,VH,VF,eigs)
end

function issolved(hf::HartreeFockIterator, epsilon)
    const mindeltas = 5
    if length(hf.es) < mindeltas+1
        false
    else
        maximum(abs(diff(hf.es[end-mindeltas:end]))) < epsilon
    end
end

function solve!(hf::HartreeFockIterator; mixing=0.0, epsilon=1e-13, maxiters=20)
    for i=1:maxiters
        if issolved(hf, epsilon)
            break
        end
        iterate!(hf, mixing)
    end
    hf.es[end]
end

end
