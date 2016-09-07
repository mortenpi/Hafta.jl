module QRPA
using Hafta
import Hafta.HFB: HFBState

export QRPAOperator
export arneigs

"""
The `QRPAOperator`.

`M` is the size of the single particle basis. The QRPA operator therefore operatates
on a `M(M-1)`-dimensional (complex) vector space.
The `X` and `Y` matrices are antisymmetric and therefore have only `M(M-1)/2`
independent components.

`E`, `U` and `V` are `MxM` matrices. `E` is a diagonal matrix with the HFB spectrum
on the diagonal. `U` and `V` are the respective matrices from the Bogoliubov
transformation.
"""
type QRPAOperator{T <: Hafta.ManyBodySystem}
    system::T
    M::Int
    E::Diagonal{Float64}
    U::Matrix{Float64}
    V::Matrix{Float64}
end

"""
Constructs a QRPA operator from a `HFBState` and `LinAlg.Eigen` eigenvalue
and eigenvector factorization.

The `efact` is used to get the HFB eigenvalues for the `E` matrix.
"""
function QRPAOperator{T}(state::HFBState{T})
    M = length(state.system)
    E = Diagonal(state.energies)
    QRPAOperator{T}(state.system.system, M, E, state.U, state.V)
end

import Base: size, issym, eltype, *
function size(op::QRPAOperator)
    N = (op.M-1)*op.M
    N,N
end
issym(op::QRPAOperator) = false
eltype{T}(::Type{QRPAOperator{T}}) = Float64

function _gamma(system, rho)
    N = size(system)
    gamma = zeros(Float64, (N,N))
    for i=1:N,j=1:N,k=1:N,l=1:N
        gamma[i,j] += rho[k,l]*( V(system, i,k,j,l)-V(system, i,k,l,j) )
    end
    gamma
end

function _delta(system, kappa)
    N = size(system)
    delta = zeros(Float64, (N,N))
    for i=1:N,j=1:N,k=1:N,l=1:N
        delta[i,j] += 0.5*kappa[k,l]*( V(system, i,j,k,l)-V(system, i,j,l,k) )
    end
    delta
end

"""Unpacks a the vector `vZ` into an antisymmetric matrix."""
function antisymmetric_unpack{T}(vZ::Vector{T})
    Mz = length(vZ)
    M = round(Int, (1 + sqrt(1+8*Mz))/2)
    @assert Mz == div(M*(M-1),2)

    Z = zeros(T, (M,M))
    ptr = 1
    for i=1:M-1
        delta = M-i-1
        Z[i,i+1:M] = vZ[ptr:ptr+delta]
        Z[i+1:M,i] = -vZ[ptr:ptr+delta]
        ptr += delta+1
    end
    Z
end

"""
Packs the antisymmetric matrix `Z` into a vector.

It tries to cancel out any numeric instabilities or biases
by summing the upper and lower components together.
"""
function antisymmetric_pack{T}(Z::Matrix{T})
    M,_ = size(Z)
    @assert _ == M
    @assert maximum(abs(Z + transpose(Z))) < 1e-10
    Mz = div(M*(M-1),2)

    vZ = zeros(T, Mz)
    ptr = 1
    for i=1:M-1
        delta = M-i-1
        vZ[ptr:ptr+delta] += transpose(Z[i,i+1:M])
        vZ[ptr:ptr+delta] -= Z[i+1:M,i]
        vZ[ptr:ptr+delta] *= 0.5
        ptr += delta+1
    end
    vZ
end

"""This implements the action of the linear operator."""
function *(op::QRPAOperator, zs::Vector)
    M = op.M
    E,U,V = op.E, op.U, op.V

    Mz = div(M*(M-1),2)

    X = antisymmetric_unpack(zs[1:Mz])
    Y = antisymmetric_unpack(zs[Mz+1:2Mz])

    rho = U*X*transpose(V) + conj(V)*Y*ctranspose(U)
    #_kappa(X,Y) = U*X*transpose(U) + conj(V)*Y*ctranspose(V)
    #kappa1 = _kappa(X,Y)
    #kappa2 = ctranspose(_kappa(ctranspose(Y),ctranspose(X)))
    kappa1 = U*X*transpose(U) + conj(V)*Y*ctranspose(V)
    kappa2 = V*X*transpose(V) + conj(U)*Y*ctranspose(U)

    gamma = _gamma(op.system, rho)
    delta1 = _delta(op.system, kappa1)
    delta2 = ctranspose(_delta(op.system, ctranspose(kappa2)))

    _W1(Γ,Δ1,Δ2) = ctranspose(U)*Γ*conj(V) + ctranspose(U)*Δ1*conj(U) + ctranspose(V)*Δ2*conj(V) - ctranspose(V)*transpose(Γ)*conj(U)
    _W2(Γ,Δ1,Δ2) =  transpose(V)*Γ*U       +  transpose(V)*Δ1*V       +  transpose(U)*Δ2*U       -  transpose(U)*transpose(Γ)*V
    W1 = _W1(gamma,delta1,delta2)
    W2 = _W2(gamma,delta1,delta2)
    #W2 = _W(gamma',delta2',delta1')'

    Xn = E*X + X*E + W1
    Yn = -(E*Y + Y*E + W2)

    Z1 = antisymmetric_pack(Xn)
    Z2 = antisymmetric_pack(Yn)
    vcat(Z1,Z2)
end


# Storing QRPA solutions
"""
`QRPASolution` stores a particular solution found by the QRPA routine.
"""
type QRPASolution{T}
    hfb::HFBState{T}
    energy::Complex128
    X::Matrix{Complex128}
    Y::Matrix{Complex128}
    function QRPASolution(hfb::HFBState{T}, e, zs::Vector)
        Mz = div(length(zs),2)
        M = round(Int, (1 + sqrt(1+8*Mz))/2)
        @assert Mz == div(M*(M-1),2)
        X = antisymmetric_unpack(zs[1:Mz])
        Y = antisymmetric_unpack(zs[Mz+1:2Mz])
        scaling = sqrt(abs(sum(abs(X).^2) - sum(abs(Y).^2)))
        new(hfb,e,X/scaling,Y/scaling)
    end
end
QRPASolution{T}(hfb::HFBState{T},e,zs) = QRPASolution{T}(hfb,e,zs)

"""
`rho_kappa(::QRPASolution)` calculates the `rho`, `kappa1` and `kappa1` matrices
corresponding to the specified `QRPASolutions`.

It returns a tuple of matrices `(rho, kappa1, kappa2)`.
"""
function rho_kappa(qrpa::QRPASolution)
    U,V = qrpa.hfb.U, qrpa.hfb.V
    X,Y = qrpa.X,qrpa.Y
    rho = U*X*transpose(V) + conj(V)*Y*ctranspose(U)
    kappa1 = U*X*transpose(U) + conj(V)*Y*ctranspose(V)
    kappa2 = V*X*transpose(V) + conj(U)*Y*ctranspose(U)
    rho, kappa1, kappa2
end

@enum QRPASolutionClass MIXED PH PP HH UKW

"""Uses the rho and kappa matrices to classify the QRPA solutions."""
function classify(s::QRPASolution)
    rho,kappa1,kappa2 = QRPA.rho_kappa(s)
    S1,S2,S3 = sum(abs(rho).^2), sum(abs(kappa1).^2), sum(abs(kappa2).^2)

    r,m1,m2 = S1 > 1e-8, S2 > 1e-8, S3 > 1e-8
    if r && (m1 || m2)
        MIXED
    elseif  r && !m1 && !m2
        PH
    elseif !r &&  m1 && !m2
        PP
    elseif !r && !m1 &&  m2
        HH
    else
        UKW
    end
end

# Some UGLY code to wrap the QRPA solving
"""Solves the full QRPA problem."""
function get_qrpa_values{T}(::Type{T}, N,A,c; nev=6, verbose=false, maxiters=200)
    s = T(N,c)
    hfbi = hfb(s,A)
    e0 = solve!(hfbi, maxiters=maxiters)
    state = hfbi.states[end]
    if verbose println("> HFB solved (E=$e0)") end

    # QRPA
    op = QRPAOperator(state)
    if verbose println("> QRPA constructed (size=$(size(op)))") end

    evs,vecs = arneigs(op; nev=nev, which=:SM)
    if verbose println("> QRPA solved.") end
    #evs = filter(ev -> abs(ev).>1e-10 && !isnan(ev), evs)
    #evs = map(abs, evs)
    #@show size(vecs)
    solutions = QRPASolution[]
    for (i,ev)=enumerate(evs)
        if abs(imag(ev))/abs(ev) > 1e-12
            if verbose warn("Complex eigenvalue: $ev") end
        end
        push!(solutions, QRPASolution(state, ev, vecs[:,i]))
    end
    hfbi,op,solutions
end

function qrpa_attempts{T}(::Type{T}, N,A,g; nev=40)
    hfbi,solutions = get_qrpa_values(T, N,A,g; nev=nev, maxiters=1500)
    sort!(solutions, lt=(s1,s2) -> abs(s1.energy) < abs(s2.energy))

    ipos,ineg = 0,0
    nevs, Es, Eims = Vector{Int}(), Vector{Float64}(), Vector{Float64}()
    classes = Vector{ASCIIString}()
    for s=solutions
        class = classify(s)

        E = s.energy
        iscomplex = (abs(imag(E)) > 1e-10 && abs(imag(E)/abs(E)) > 1e-10)

        E = real(E)
        nev = if (class == MIXED || class == PH) && !iscomplex
            (E >= 0.0) ? ipos+=1 : -(ineg += 1)
        else
            0
        end

        push!(Es, E)
        push!(Eims, imag(s.energy))
        push!(nevs, nev)
        push!(classes, string(class))
    end

    DF = DataFrame(nev=nevs, E=Es, Eim=Eims, class=classes)
    DF[:N] = N
    DF[:A] = A
    DF[:g] = g
    DF[:Egs] = hfbi.es[end]
    DF[:λ] = hfbi.states[end].lambda
    DF
end


# Methods related to the Arnoli method for solving the eigenvalue problem
import Base.LinAlg: ARPACK, BlasInt, chksquare
"""
This is a modified version of the `eigs` function from `Base.LinAlg`.
"""
function arneigs(A;
              nev::Integer=6, ncv::Integer=max(20,2*nev+1), which=:LM,
              tol=0.0, maxiter::Integer=300, sigma=nothing, v0::Vector=zeros(eltype(A),(0,)),
              ritzvec::Bool=true)
    B = I # we do not need the general case
    n = chksquare(A)

    T = eltype(A)
    iscmplx = T <: Complex
    isgeneral = B !== I
    sym = issym(A) && !iscmplx
    nevmax=sym ? n-1 : n-2
    if nevmax <= 0
        throw(ArgumentError("Input matrix A is too small. Use eigfact instead."))
    end
    if nev > nevmax
        warn("Adjusting nev from $nev to $nevmax")
        nev = nevmax
    end
    if nev <= 0
        throw(ArgumentError("requested number of eigenvalues (nev) must be ≥ 1, got $nev"))
    end
    ncvmin = nev + (sym ? 1 : 2)
    if ncv < ncvmin
        warn("Adjusting ncv from $ncv to $ncvmin")
        ncv = ncvmin
    end
    ncv = BlasInt(min(ncv, n))
    if isgeneral && !isposdef(B)
        throw(PosDefException(0))
    end
    bmat = isgeneral ? "G" : "I"
    isshift = sigma !== nothing

    if isa(which,AbstractString)
        warn("Use symbols instead of strings for specifying which eigenvalues to compute")
        which=symbol(which)
    end
    if (which != :LM && which != :SM && which != :LR && which != :SR &&
        which != :LI && which != :SI && which != :BE)
        throw(ArgumentError("which must be :LM, :SM, :LR, :SR, :LI, :SI, or :BE, got $(repr(which))"))
    end
    if which == :BE && !sym
        throw(ArgumentError("which=:BE only possible for real symmetric problem"))
    end
    isshift && which == :SM && warn("use of :SM in shift-and-invert mode is not recommended, use :LM to find eigenvalues closest to sigma")

    # This screws up the calculation. When using :SM it tries to be clever or something
    # and it will requre factorization of the operator. Disablin this passes :SM
    # directly to ARPACK and makes it work nicely.
    #if which==:SM && !isshift # transform into shift-and-invert method with sigma = 0
    #    isshift=true
    #    sigma=zero(T)
    #    which=:LM
    #end

    if sigma !== nothing && !iscmplx && isa(sigma,Complex)
        throw(ArgumentError("complex shifts for real problems are not yet supported"))
    end
    sigma = isshift ? convert(T,sigma) : zero(T)

    if !isempty(v0)
        if length(v0) != n
            throw(DimensionMismatch())
        end
        if eltype(v0) != T
            throw(ArgumentError("starting vector must have element type $T, got $(eltype(v0))"))
        end
    end

    whichstr = "LM"
    if which == :SM
        whichstr = "SM"
    end
    if which == :BE
        whichstr = "BE"
    end
    if which == :LR
        whichstr = (!sym ? "LR" : "LA")
    end
    if which == :SR
        whichstr = (!sym ? "SR" : "SA")
    end
    if which == :LI
        if !sym
            whichstr = "LI"
        else
            throw(ArgumentError("largest imaginary is meaningless for symmetric eigenvalue problems"))
        end
    end
    if which == :SI
        if !sym
            whichstr = "SI"
        else
            throw(ArgumentError("smallest imaginary is meaningless for symmetric eigenvalue problems"))
        end
    end

    # Refer to ex-*.doc files in ARPACK/DOCUMENTS for calling sequence
    matvecA(x) = A * x
    if !isgeneral           # Standard problem
        matvecB(x) = x
        if !isshift         #    Regular mode
            mode       = 1
            solveSI(x) = x
        else                #    Shift-invert mode
            mode       = 3
            F = factorize(sigma==zero(T) ? A : A - UniformScaling(sigma))
            solveSI(x) = F \ x
        end
    else                    # Generalized eigenproblem
        matvecB(x) = B * x
        if !isshift         #    Regular inverse mode
            mode       = 2
            F = factorize(B)
            solveSI(x) = F \ x
        else                #    Shift-invert mode
            mode       = 3
            F = factorize(sigma==zero(T) ? A : A-sigma*B)
            solveSI(x) = F \ x
        end
    end

    # Compute the Ritz values and Ritz vectors
    (resid, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, TOL) =
       ARPACK.aupd_wrapper(T, matvecA, matvecB, solveSI, n, sym, iscmplx, bmat, nev, ncv, whichstr, tol, maxiter, mode, v0)

    # Postprocessing to get eigenvalues and eigenvectors
    output = ARPACK.eupd_wrapper(T, n, sym, iscmplx, bmat, nev, whichstr, ritzvec, TOL,
                                 resid, ncv, v, ldv, sigma, iparam, ipntr, workd, workl, lworkl, rwork)

    # Issue 10495, 10701: Check that all eigenvalues are converged
    nev = length(output[1])
    nconv = output[ritzvec ? 3 : 2]
    nev ≤ nconv || warn("not all wanted Ritz pairs converged. Requested: $nev, converged: $nconv")

    return output
end

end # module
