# Functions and data structured related to the 1D harmonic oscillator
module HarmonicOscillator

import Iterators
import Polynomials
import Orthopolys
import Cubature

import Base.start
import Base.done
import Base.next
import Base.length
import Base.getindex

export Psi
export WMatrix
export generate_wmatrix

# ===================================================================
# Generic helper functions and macros
# ===================================================================

# Returns the expression where elements of xs have been chained
# together with op.
function chainoperator(op,xs)
  ret = :()
  for i=1:length(xs)-1
    e = :( $op($(xs[i]), $(xs[i+1])) )
    ret = (i==1) ? e : :( $ret && $e )
  end
  ret
end

# Calls the function f with sorted arguments.
macro ordered_fcall(f, key...)
  ret::Expr = :()
  isfirstpermutation = true
  for p=permutations(key)
    chk = chainoperator(<=, p)
    call = :( $f($(p...)) )
    if !isfirstpermutation
      ret = quote
        if $chk
          $call
        else
          $ret
        end
      end
    else
      isfirstpermutation=false
      ret = call
    end
  end
  ret
end

# ===================================================================
# General harmonic oscillator functions
# ===================================================================

# Returns a 1D harmonic oscillator basis function.
function Psi(n::Integer, x0::AbstractFloat=1.0)
  p=Orthopolys.phhermite(n+1, Float64)
  C = 1/((pi^.25)*sqrt(x0*(2^n)*factorial(n)))
  f(x) = begin
    q = x/x0
    C*Polynomials.polyval(p,q).*exp(-(q.^2)/2)
  end
end

# ===================================================================
# The 1-dimensional W matrix
# ===================================================================

# Generates an element of the W matrix
function W(i,j,k,l)
  # Odd n means odd Psi. Therefore if the sum of the indices
  # is odd then the product is odd as well.
  if isodd(i+j+k+l)
    return (0.0,0.0)
  end
  A = Psi(i)
  B = Psi(j)
  C = Psi(k)
  D = Psi(l)
  quadgk(x -> A(x)*B(x)*C(x)*D(x), -Inf, Inf)
end

# WMatrix stores the W matrix elements as single array.
# We define the getindex method and do some index magic to actually
# get the W: (i,j,k,l) -> W[idx(i,j,k,l)] mapping.
type WMatrix
  values::Array{Float64, 1}
  errors::Array{Float64, 1}
end

# The W(n,m) and Wm(n) calculate how many elements are in
# the (N,m) block where m is the number of indices and N the maximum
# value of the indices.
# The function can be defined as a recurrence relation which can be
# solved. W(n,m) uses the recursive definition directly whereas
# Wm(n) use the polynomial representation for the first four m
# values. Using the polynomials gives a noticeable speed increase.
function W(n,m)
  if n < 0; 0
  elseif n == 0; 1
  elseif m == 1; n+1
  else; W(n-1, m) + W(n, m-1); end
end

W1(n) = n+1
W2(n) = div( (n+3)*n, 2) + 1
W3(n) = div( ((n+6)*n+11)*n, 6) + 1
W4(n) = div( (((n+10)*n+35)*n+50)*n, 24) + 1

# Calculates the array index of a key in the W matrix.
wmatrix_arridx_fromsorted_recursive(i,j,k,l) = W(i-1,1)+W(j-1,2)+W(k-1,3)+W(l-1,4)+1
wmatrix_arridx_fromsorted(i,j,k,l) = W1(i-1)+W2(j-1)+W3(k-1)+W4(l-1)+1
function wmatrix_arridx(i,j,k,l)
  @ordered_fcall wmatrix_arridx_fromsorted i j k l
end

# Return the corresponding matrix element
function Base.getindex(W::WMatrix, i,j,k,l)
  idx = wmatrix_arridx(i,j,k,l)
  W.values[idx]
end

# SortedWMatrixKeys is an helper iterator which generates the sorted keys
# that the use when filling the WMatrix.
type SortedWMatrixKeys
  N
end

function Base.start(b::SortedWMatrixKeys)
  (0,0,0,0)
end
function Base.done(b::SortedWMatrixKeys,state)
  state == (b.N+1,0,0,0)
end
function Base.next(b::SortedWMatrixKeys,state)
  ret = [0,0,0,0]
  i = length(state)
  plusone = true
  while i > 1
    if plusone
      if state[i] < state[i-1]
        ret[i] = state[i]+1
        plusone = false
      end
    else
      ret[i] = state[i]
    end
    i-=1
  end
  ret[1] = state[1]+ (plusone ? 1 : 0)
  state,tuple(ret...)
end
function Base.length(b::SortedWMatrixKeys)
  W4(b.N)
end

# Generates a lookup table (type WMatrix) of the matrix elements of
# the one-dimensional delta interaction for the first `n` basis states.
function generate_wmatrix(n::Integer)
  keys = SortedWMatrixKeys(n)
  values = Array(Float64, length(keys))
  errors = Array(Float64, length(keys))
  for k in keys
    idx = wmatrix_arridx(k...)
    v,e = W(k...)
    values[idx] = v
    errors[idx] = e
  end
  WMatrix(values, errors)
end

end
