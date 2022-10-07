using Test
using LinearAlgebra

function gram_schmidt(A::Matrix)
    m, n = size(A)
    rank(A) == n || throw(SingularException)
    Q = fill(0.0, m, n)
    R = fill(0.0, n, n)
    for k in 1:n
        R[1:(k-1), k] = conj(transpose(Q[:,1:k-1]))*A[:,k] 
        q = A[:,k] - (Q*R[:, k])
        R[k,k] = norm(q)
        Q[:,k] = q/R[k,k]
    end
    return(Q, R)
end


