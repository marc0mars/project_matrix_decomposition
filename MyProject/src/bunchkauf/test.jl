using Test
using LinearAlgebra 
include("bunchkaufman_impl.jl")

@testset "outputs correct factorization" begin
    for i in 1:100
        A = (Symmetric(rand(10,10)))                                  # make symmetric
        A = [(abs(i-j) <= 1) ? A[i,j] : 0.0 for i in 1:10, j in 1:10] # make tridiagonal
        L, D = bunch_explicit(A)
        @test isapprox(L*D*L', A)
    end
end
