using LinearAlgebra
using Test
using DecompositionProject

@testset "Testing HQR" begin
    A = rand(5,3)*[1 1 0;0 0 0; 0 0 1]
    try
        HQR(A)
    catch err
        @test err.msg == "This matrix is singular!"
    end
end

# A = [(i == j) ? 1.0 : 0.0 for i in 1:10, j in 1:10]
# display(applyQ_adj(HQR(A)..., ones(10)))
# display(HQR(A)[1])