using BenchmarkTools
using LinearAlgebra
using DecompositionProject

A = (Symmetric(rand(100,100)))                                # make symmetric

println("This is my bunch")
b = @benchmark bunchfast(A)
A = Matrix(A)
println("This is the LAPACK (?) bunch")
b2 = @benchmark bunchkaufman(A)

display(b)
display(b2)