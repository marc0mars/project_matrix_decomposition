using LinearAlgebra

# the numerical stability of gram_schmidt is not as good as the numerical stability of gram_schmidt3 which is due 
# to the fact that 3 basically corrects roundoff errors that were made by previous qₖ in the direction of the current qₖ

# shortcut-function
function adj(x)
    return transpose(conj(x))
end


"""
    gram_schmidt(A)

    A: Matrix

    This utilizes the gram schmidt algorithm in order to compute the qr factorization! 
    This method is not as good as those following bellow.
"""
function gram_schmidt(A::Matrix)
    m, n = size(A)
    # rank(A) == n || throw(error("The matrix is singular!"))
    Q = fill(0.0, m, n)
    R = fill(0.0, n, n)
    for k in 1:n
        R[1:(k-1), k] = adj(Q[:,1:k-1])*A[:,k] 
        q = A[:,k] - (Q*R[:, k])
        R[k,k] = norm(q)
        Q[:,k] = q/R[k,k]
    end
    return(Q, R)
end


"""
    gram_schmidt2(A)

    A: Matrix

    This algorithm utilizes the fact that we dont need the j-th column of A after we have already computed the j-th column of Q.
    Therefore Q can be saved inside of A.

    This is also known of modified gram-schmidt (MGS)
"""
###### CURRENTLY NOT WORKING PROPERLY
#=
function gram_schmidt2(A::Matrix)
    M = copy(A) # necessary because this time we save the unitary elements into A but dont want to change matrix A 
    n = size(A, 2)
    rank(M) == n || throw(error("The matrix is singular!"))
    R = fill(0.0, n, n)
    for j in 1:n
        for k in 1:j-1
            R[k,j] = adj(M[:,k])*M[:,j]
            M[:,j] = M[:,j] - R[k,j]*M[:,k]
        end
        R[j,j] = (norm(A[:,j]))
        M[:,j] = M[:,j]/R[j,j]
    end
    return(M, R)
end
=#

"""
    gram_schmidt3(A)

    A: Matrix

    The nicest of the above algorithms as it replaces A with Q and keeps the remaining columns of A orthogonal to those of Q all the time.
    It computes R from left to right in each row and then moves to the next row. (Whereas gram_schmidt1 moves from column to column)
"""
function gram_schmidt3(A::Matrix)
    M = copy(A)
    n = size(M, 2)
    R = fill(0.0, n, n)
    # rank(M) == n || throw(error("The matrix is singular!"))
    for k in 1:n
        R[k,k] = norm(M[:,k])
        M[:,k] = M[:,k]/R[k,k]
        R[k,k+1:n] = adj(M[:,k])*M[:,k+1:n]
        M[:,k+1:n] = M[:,k+1:n]-M[:,k]*transpose(R[k,k+1:n])
    end
    return(M,R)
end