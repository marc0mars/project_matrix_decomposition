using LinearAlgebra
using Test

include("./householder.jl")
include("../Gram-Schmidt/gram_schmidt_orth.jl")

@testset "HQR does not accept singular matrices" begin
    A = rand(5,3)*[1 1 0;0 0 0; 0 0 1]
    try
        HQR(A)
    catch err
        @test err.msg == "This matrix is singular!"
    end
end

A = rand(10, 7)
display(transpose(Q_part(HQR(A)[1],HQR(A)[2]))*(Q_part(HQR(A)[1],HQR(A)[2])))