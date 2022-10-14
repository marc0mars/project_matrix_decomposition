using LinearAlgebra

# let T be a tridiagonal, symmetric matrix. We aim for a decomposition where T = LDL^t such that L is a unit lower triangular matrix.
# T = [B    T_12
#      T_21 T_22]

"""
    bunch_explicit(A)

    A: tridiagonal and symmetric

    Just like the above but outputting the Matrix L directly.
"""
function bunch_explicit(A)
    A = copy(A)
    n = size(A,1)
    k = 1
    L = [(i == j) ? 1.0 : 0.0 for i in 1:n, j in 1:n]
    while k < n
        s = choosing_pivot_size(A[k:n,k:n])   # set pivot size
        if s == 1
            if A[k,k] == 0
                k += 1
                continue        # skip to the next round
            end
            d₁, b₁ = A[k,k], A[k,k+1]
            # update A_k+1,k+1
            A[k+1,k+1] -= b₁^2*(1/d₁) 
            L[k+1,k] = b₁*(1/d₁)
            A[k+1,k] = 0
            A[k, k+1] = 0
            k += 1 # as we only had a one-block
        else # s == 2
            # when the last block is chosen to be 2 by 2 we can stop
            if(k == n-1)
                k += 1 
                continue
            end    
            d₁, d₂, b₁, b₂ = A[k,k], A[k+1,k+1], A[k+1,k], A[k+2,k+1]
            Δ = d₁*d₂-b₁^2  # det of the 2x2 block
            # update A[k+2,k+2]
            A[k+2,k+2] -= d₁*b₂^2/Δ
            # updating L
            L[k+2, k] = -b₂*b₁/Δ
            L[k+2, k+1] = b₂*d₁/Δ
            # deleting the zero corners
            A[k+2,k+1] = 0
            A[k+1,k+2] = 0
            k += 2
        end
    end
    return (L,A)
end

# we can try different pivoting strategies

"""
    choosing_pivot_size(B)

    chooses the pivot size according to the original bunchkaufman
"""
function choosing_pivot_size(B)  
    n = size(B,2)                        
    sigma = max([(i == 1 && j == 1) ? 0.0 : abs(B[i,j]) for i in 1:n, j in 1:n]...)                          # get maximum
    alpha = (sqrt(5)-1)/2                       # golden ratio 
    if(sigma*abs(B[1,1]) >= alpha*B[2,2,]^2)    # this also is the case when the one and two block are singular
        return 1                                # we choose a block of size 1x1
    else
        return 2                                # we choose a block of size 2x2
    end
end

# getting the different matrices

function get_L(A)

end

# first primitive recursive implementation of Bunch

# """
#     bunch(A)

#     A: tridiagonal and symmetric

#     Currently changing A
# """
# function bunch(A)
#     n = size(A,1)
#     k = 1
#     while k < n
#         s = choosing_pivot_size(A[k:k+1,k:k+1])   # set pivot size
#         if s == 1
#             if A[k,k] == 0
#                 k += 1
#                 continue        # skip to the next round
#             end
#             # the following lines can also be pressed into one if necessary
#             d₁ = A[k,k]
#             b₁ = A[k,k+1]
#             # update A_k+1,k+1
#             A[k+1,k+1] -= b₁^2*(1/d₁) # maybe you want to use a more abstract form of the inverse 
#             # saving L in the lower Part of A overwriting the original matrix
#             A[k+1,k] = b₁*(1/d₁)
#             k += 1 # as we only had a one-block
#         else # s == 2
#             k != n - 2 || continue # when the last block is chosen to be 2 by 2 we can stop
#             d₁, d₂, b₁, b₂ = A[k,k], A[k+1,k+1], A[k+1,k], A[k+2,k+1]
#             Δ = d₁*d₂-b₁^2  # det of the 2x2 block
#             # update A[k+2,k+2]
#             A[k+2,k+2] -= d₁*b₁^2/Δ
#             # overwriting A with L
#             A[k+2, k] = -b₂*b₁/Δ
#             A[k+2, k+1] = b₂*d₁/Δ
#             k += 2
#         end
#     end
#     return A
# end


 