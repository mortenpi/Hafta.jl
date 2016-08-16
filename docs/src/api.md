Hafta API overview
==================

A high level overview of structures and functions in the Hafta.jl package.

## HarmonicOscillator
This module defines things necessary to calculate a 1D harmonic oscillator.

It defines two main things:
 - `Psi` function which returns an oscillator basis function.
 - `generate_wmatrix` functions, which generates the W-matrix, which is
   basically a lookup table for the delta interaction in 1D.
   It returns a `WMatrix` object and the different elements can be accessed
   using the usual index notation: `w[i,j,k,l]`.
