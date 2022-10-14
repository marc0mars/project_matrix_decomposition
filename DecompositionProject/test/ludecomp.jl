using LinearAlgebra
using Random
using Test
using Permutations
using DecompositionProject

# Band = Matrix{Float64}([(j < i-5 || j > i+3) ? 0.0 : rand(1:10) for i in 1:10, j in 1:10])

# Skykline-Matrix nach dem Muster unten
#=
1 0 0 0 0 0 0
0 1 0 0 1 0 0
0 0 1 0 1 0 0 
1 1 1 1 1 1 0 
0 0 0 0 1 1 0 
0 1 1 1 1 1 0
=#
#=
S = [
1 0 0 3 0;
1 2 0 3 0;
0 0 1 2 0;
0 1 1 3 0;
0 0 0 0 1; 
]
l,r = lu(S)
display(l)
display(r)
=#

#= M = Matrix{Float64}([(i >= j) ? rand(1:10) : 0 for i in 1:5, j in 1:5])
display(M)
M = M*transpose(M)
d = cholesky(M)
display(d)
=# 

#=
A = rand(1:10, 5, 5)
display(A[2:-1,:])
display(view(A,1:2,4:5))
=#

#=
A = rand(1:10, 5, 5)
display(A)
println(findmax(A[2:2,1:2]))
=#
# A = [2 1 1; -2 -2 1; 4 4 1]
# display(left_right_decomposition_version_2(A)[1])
# display(left_right_decomposition_version_2(A)[2])

@testset "lu_decomp_pivot" begin
    for i in 1:1
        n = 5  # test size
        M = rand(1:10, n, n)
        b = randperm(n)
        l,r, perm = lu_decomp_pivot(M)
        P = Matrix(perm)
        display(P)
        display(perm)
        x = M\b
        solution = r\(l\(P*b))
        @test isapprox(x, solution)
    end
end
