using LinearAlgebra

# implementation of qr-factorization via the householder transformations

"""
    housev(x)

    x: real valued vector (I did not test with complex numbers)

    Housev gives you all the necessary parameters to determine the Householder-matrix
"""
function housev(x)
    m = size(x,1)
    ᵨ = -sign(x[1])*norm(x)
    v₁ = x[1] - ᵨ
    u₂ = x[2:m]/v₁
    tau = (1 + transpose(u₂)*(u₂))*0.5
    return(ᵨ,u₂,tau)                      # all parameters in order to determine the Householder-matrix H₀
end

# main function

"""
    HQR(A)

    A: any mxn matrix (but should be real)

    returns the first step of the Householder QR-factorization. 
"""
function HQR(A)
    M = copy(A) # this is misleading because we would rather want to implement HQR!(A) but for testing this is more convenient
    m, n = size(M)
    rank(M) == n || throw(error("This matrix is singular!")) 
    t = fill(0.0, n)
    for k in 1:n 
        ᵨ, u₂, tau = housev(M[k:m,k]) 
        u₂ = reshape(u₂, m-k, 1)
        M[k,k] = ᵨ
        M[k+1:m,k] .= u₂
        w = (reshape(M[k,k+1:n], 1, n-k)+transpose(u₂)*M[k+1:m,k+1:n])/tau
        M[k,k+1:n] -= reshape(w,n-k)
        M[k+1:m,k+1:n] -= u₂*w      # rank-1 update
        t[k] = tau
    end
    return(M, t)
end

# make Q explicit

"""
    Q_part(A, t)

    A: Matrix that was created by HQR(M), therefore A must be HQR(M)[1]
    t: column of lenght parameters that was created by HQR(M), therefore t must be HQR(M)[2]

    Q_part returns the Q matrix of the Housholder QR-decomposition. Q is implicitly stored in the lower part of the first tupelentry of HQR()
    and is made explicit here.

    Funnily it basically does the same as HQR but in reverse with a mxn "unit matrix".
"""
function Q_part(A, t)           # this is really costly and requires O(mn^2) computations 
    m, n = size(A)
    M = copy(A)        
    for k in 1:n
        l = (n+1-k)             # as we want to count backwards
        u = M[l+1:m,l]          # current defining householder-vector
        M[l,l] = 1 - 1/t[l] 
        M[l,l+1:n] = (-1)*transpose(u)*M[l+1:m,l+1:n]/t[l]
        M[l+1:m,l+1:n] += reshape(u,size(u)[1],1) * reshape(M[l,l+1:n],1,k-1)
        M[l+1:m,l] = -u/t[l]
    end
    return M
end

# the following are smaller helper functions

function R_part(A)
    n = size(A,2)
    return(UpperTriangular(A[1:n,1:n]))
end

function householder(A)
    Q = Q_part(HQR(A)...)
    R = R_part(HQR(A)[1])
    return(Q,R)
end

function applyQ_adj(A, t, y) # y is a Vector (mx1) and Q is given within HQR and therefore not explicit, t is also a return value of HQR
    m = size(A, 1)
    for k in 1:m
        # ω is only a helper to make it more readable
        ω = (y[k]+transpose(A[k+1:m,k])*y[k+1:m])/t[k]
        # update y
        y[k] -= ω
        y[k+1:m] -= A[k+1:m,k]*ω
    end
    return y
end