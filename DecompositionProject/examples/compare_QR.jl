using LinearAlgebra
using Plots
using DecompositionProject

ϵ = eps()                           # machine epsilon
A = [1 1 1; 10ϵ 0 0; 0 10ϵ 0; 0 0 10ϵ]
println("This is the condition number of the used matrix: $(cond(A)) \n")
CGS = gram_schmidt(A)[1]            # classical gram schmidt
MGS = gram_schmidt3(A)[1]           # modified gram schmidt
Householder = Q_part(HQR(A)...)             # Housholder variant
methods_Q = [CGS, MGS, Householder]
diff = [norm(transpose(Q)*Q-[1 0 0; 0 1 0; 0 0 1]) for Q in methods_Q]
display(diff)
plot(["CGS", "MGS", "Householder"], diff; seriestype=:scatter, lw=4, title="Comparison of different QR-Factorization")