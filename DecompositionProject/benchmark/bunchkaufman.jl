using BenchmarkTools
using LinearAlgebra
using DecompositionProject

A = (Symmetric(rand(1000,1000)))                                # make symmetric
A = SymTridiagonal(A)                                           # make tridiagonal
println("This is my bunch")
b = @benchmark bunchfast(A)
A = Matrix(A)
println("This is the LAPACK (?) bunch")
b2 = @benchmark bunchkaufman(A)

display(b)
display(b2)