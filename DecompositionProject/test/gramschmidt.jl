using LinearAlgebra
using Plots
using Test
using DecompositionProject

@testset "Q is unitary" begin
    one_mat = [(i == j) ? 1 : 0 for i in 1:7, j in 1:7]
    for i in 1:5
        A = rand(10,7)
        Q = gram_schmidt3(A)[1]
        @test isapprox(transpose(Q)*Q, one_mat)
    end
end

@testset "function detects when A does not have full rank" begin
    mat = [1 2 3; 1 2 5; 2 4 19]
    try 
        gram_schmidt3(mat)
    catch err
        @test isa(err, Exception)
        @test err.msg == "The matrix is singular!"
    end
end

@testset "function computes correct QR decomposition" begin
    for i in 1:5
        A = rand(10,7)
        Q, R = gram_schmidt3(A)   
        @test isapprox(Q*R,A)
    end
end

function benchmark_gs1_vs_gs2()
    # comparison between gram_schmidt3 and gram_schmidt1
    num_runs = 10
    one_mat = [(i == j) ? 1.0 : 0.0 for i in 1:70, j in 1:70]
    runs = []
    results = []
    for j in 1:num_runs
        push!(runs, j)
        A = rand(100,70)
        Q1 = gram_schmidt(A)[1]
        Q3 = gram_schmidt3(A)[1]
        ratio_1_to_3 = norm(adj(Q1)*Q1-one_mat)/norm(adj(Q3)*Q3-one_mat)
        push!(results, ratio_1_to_3)
    end
    plot(runs, results; title = "ratio of different algs for QR-decomp", label = "results", xlabel = "run#", ylabel = "norm(Q1*adj(Q1)) divided by norm(Q3*adj(Q3))", lw = 3, seriestype = :scatter)
    plot!(runs, ones(num_runs), label = "ratio one", lw = 3)
end

# benchmark_gs1_vs_gs2()
