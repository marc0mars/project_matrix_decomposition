using Test
using LinearAlgebra 
using BenchmarkTools
using DecompositionProject

@testset "outputs correct factorization" begin
    for i in 1:10
        A = (Symmetric(rand(10,10)))                                  # make symmetric
        A = [(abs(i-j) <= 1) ? A[i,j] : 0.0 for i in 1:10, j in 1:10] # make tridiagonal
        L, D = bunch_explicit(A)
        @test isapprox(L*D*L', A)
    end
end
# A = (Symmetric(rand(10,10)))                                  # make symmetric
# A = [(abs(i-j) <= 1) ? A[i,j] : 0.0 for i in 1:10, j in 1:10] # make tridiagonal

# b = @benchmark bunch_explicit(A)
# b2 = @benchmark bunchkaufman(A)
# display(b)
# display(b2)

