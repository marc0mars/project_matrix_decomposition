using Test
using LinearAlgebra 
include("bunchkaufman_impl.jl")
display(bunch(Matrix{Float64}([2 3 0; 3 1 2; 0 2 1])))