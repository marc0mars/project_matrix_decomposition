using LinearAlgebra
include("./gram_schmidt_orth.jl")

@testset "Q is unitary" begin
    one_mat = [(i == j) ? 1 : 0 for i in 1:7, j in 1:7]
    for i in 1:500
        A = rand(10,7)
        Q = gram_schmidt(A)[1]
        @test isapprox(transpose(Q)*Q, one_mat)
    end
end

@testset "function detects when A does not have full rank" begin

end

@testset "function computes correct QR decomposition" begin

end