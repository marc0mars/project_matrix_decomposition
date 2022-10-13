using LinearAlgebra

# let T be a tridiagonal, symmetric matrix. We aim for a decomposition where T = LDL^t such that L is unitary.
# T = [B    T_12
#      T_21 T_22]

# first primitive recursive implementation of Bunch

"""
    bunch(A)

    A: tridiagonal and symmetric

    Currently changing A
"""
function bunch(A)
    n = size(A,1)
    k = 1
    while k < n
        s = choosing_pivot_size(A[k:k+1,k:k+1])   # set pivot size
        if s == 1
            if A[k,k] == 0
                k += 1
                continue        # skip to the next round
            end
            # the following lines can also be pressed into one if necessary
            d₁ = A[k,k]
            b₁ = A[k,k+1]
            # update A_k+1,k+1
            A[k+1,k+1] -= b₁^2*(1/d₁) # maybe you want to use a more abstract form of the inverse 
            # saving L in the lower Part of A overwriting the original matrix
            A[k,k+1] = b₁*(1/d₁)
            k += 1 # as we only had a one-block
        else # s == 2
            k != n - 2 || continue # when the last block is chosen to be 2 by 2 we can stop
            d₁, d₂, b₁, b₂ = A[k,k], A[k+1,k+1], A[k+1,k], A[k+2,k+1]
            Δ = d₁*d₂-b₁^2  # det of the 2x2 block
            # update A[k+2,k+2]
            A[k+2,k+2] -= d₁*b₁^2/Δ
            # overwriting A with L
            A[k, k+1] = 0
            A[k, k+2] = -b₂*b₁/Δ
            A[k+1, k+2] = b₂*d₁/Δ
            k += 2
        end
    end
    return A
end

# we can try different pivoting strategies

"""
    choosing_pivot_size(B)

    chooses the pivot size. 
"""
function choosing_pivot_size(B)
    sigma = max(abs.(B)...)                     # get maximum
    alpha = (sqrt(5)-1)/2                       # golden ratio 
    if(sigma*abs(B[1,1]) >= alpha*B[2,2,]^2)    # this also is the case when the one and two block are singular
        return 1                                # we choose a block of size 1x1
    else
        return 2                                # we choose a block of size 2x2
    end
end

 