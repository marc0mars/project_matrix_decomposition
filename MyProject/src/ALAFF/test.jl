using LinearAlgebra
include("./gram_schmidt_orth.jl")

@testset "Q is unitary" begin
    one_mat = [(i == j) ? 1 : 0 for i in 1:7, j in 1:7]
    for i in 1:5
        A = rand(10,7)
        Q = gram_schmidt2(A)[1]
        @test isapprox(transpose(Q)*Q, one_mat)
    end
end

@testset "function detects when A does not have full rank" begin
    mat = [1 2 3; 1 2 5; 2 4 19]
    try 
        gram_schmidt2(mat)
    catch err
        @test isa(err, Exception)
        @test err.msg == "The matrix is singular!"
    end
end

@testset "function computes correct QR decomposition" begin
    for i in 1:5
        A = rand(10,7)
        Q, R = gram_schmidt2(A)   
        @test isapprox(Q*R,A)
    end
end