using LinearAlgebra

n = 50
A = UpperTriangular(rand(n,n))
x = rand(n)
b = A * x
x_hat = A \ b
minimum_diag = findmin(abs, diag(A))
println("Min: $(minimum_diag) and condition number: $(cond(A))")
println(norm(x_hat-x))
println(norm(b - A * x_hat))