# VERSION 1-4 DONT UTILIZE PIVOTING
# version 1 seems to me to be the most natural as a first try as it works similar to the proof
using LinearAlgebra
using Permutations

function clear_row_s(M, s)
    # M is a square matrix
    n = size(M,1)
    # s < n necessary
    if s >= n
        return fill(0.0, n, n)
    end
    # matrix with factors
    L = fill(0.0, n, n)
    for k in (s+1):n
        factor_l = M[k,s]/M[s,s]
        L[k,s] = factor_l
        for j in s:n
            M[k,j] = M[k,j] - M[s,j] * factor_l
        end
    end
    return L
end


function left_right_decomposition_version_1(A::Matrix)
    M = map(Float64, copy(A))
    n = size(M)[1]
    L = [(i == j) ? 1 : 0 for i in 1:n, j in 1:n]

    for s in 1:size(M)[1]
        L = L + clear_row_s(M, s)
    end
    return (L,M)
end

# version 2 is more sophisticated

function left_right_decomposition_version_2(M::Matrix)
    # assuming that M is square
    N = map(Float64, copy(M))
    n = size(M)[1]
    L = [(i == j) ? 1.0 : 0.0 for i in 1:n, j in 1:n]
    for j in 1:(n-1)
        for k in (j+1):n
            L[k,j] = N[k,j]/N[j,j]
            for l in j:n
                N[k,l] = N[k,l] - L[k,j]*N[j,l]
            end
        end
    end
    return(L, N)
end

# version 3 tries to omit the additional space needed for L as it is clear, that all entries on the diagonal are 1

function left_right_decomposition_version_3(M::Matrix)
    N = map(Float64, copy(M))
    n = size(M)[1]
    for j in 1:(n-1)
        for k in (j+1):n
            N[k,j] = N[k,j]/N[j,j]
            for l in (j+1):n
                N[k,l] = N[k,l] - N[j,l]*N[k,j]
            end
        end
    end
    return N
end

# version 4 Doolittle

function left_right_decomposition_version_4(M::Matrix)
    N = Matrix{Float64}(copy(M))
    n = size(M, 1)
    L = [(i == j) ? 1.0 : 0.0 for i in 1:n, j in 1:n]
    R = fill(0.0, n, n)
    for i in 1:n
        # calculating the i-th row of R
        for j in i:n
            #### calculate sum (can also be empty), can this be written more compact?
            sum = 0
            for l in 1:i-1
                sum = sum + L[i,l]*R[l,j]
            end
            #### finish calculating sum
            R[i,j] = N[i,j] - sum
        end
        # calculating the i-th column of R
        for j in i:n
            ##### calculate sum (can also be empty), can this be written more compact?
            sum = 0
            for l in 1:i-1
                sum = sum + L[j,l]*R[l,i]
            end
            #### finish calculating sum
            L[j,i] = (N[j,i] - sum)/R[i,i]
        end
    end
    return (L,R)
end

# version for band matrices, obviously also square. There is no need of computing the zeros

function left_right_decomposition_version_band(M::Matrix, p::Integer, q::Integer) # p lower diagonals, q upper diagonals
    N = Matrix{Float64}(copy(M))
    n = size(N, 1)
    for k in 1:(n-1)
        for i in k+1:(min(n,k+p))
            N[i,k] = N[i,k]/N[k,k]
            for j in k+1:(min(n,k+q))
                N[i,j] = N[i,j] - N[k,j] * N[i,k]
            end
        end
    end
    return N
end

# cholesky-decomposition

function cholesky(M::Matrix)
    # returns Cholesky factor and needs a SPD matrix
    N = copy(M)
    n = size(M,1)
    L = fill(0.0,n,n)
    for i in 1:n
        ### sum 
        sum = 0
        for k in 1:(i-1)
            sum += L[i,k]^2
        end
        ### end sum
        L[i,i] = sqrt(N[i,i]-sum)
        for j in (i+1):n
            ### sum 
            sum = 0
            for k in 1:(i-1)
                sum += L[i,k]*L[j,k]
            end
            ### end sum
            L[j,i] = (1/L[i,i])*(N[j,i]-sum)
        end
    end
    return L
end


# the following is currently not working properly
"""
    swap_row!(M::Matrix, i, j)

Swaps the rows i and j of a given Matrix M
"""
function swap_row!(M::Matrix, i, j)
    M[i,:], M[j,:] = M[j,:], M[i,:]
end

"""
    array_product(A, n)

    A: Array
    n: In this special case the length of the permutations although this doesnt seem to be best practice

Simply multiplies all elemnts inside the array starting with the first and ending with the last
Currently restricted to Permutationarrays as I dont know how to implement the neutral element properly
"""
function array_product(A, n)
    res = Permutation(n)
    for a in A
        res = res*a
    end
    return res    
end
"""
    lu_decomp_pivot(M::Matrix)

Returns the decomposed matrix using column pivoting which works for all square matrices with full rank.
This seems to be numerically stable as long as the matrix and the vector for which the Matrix has to be solved are of a "similar dimension".
Be aware that due to the permutations L is not the left part of the lu decomposition.

https://people.math.ethz.ch/~grsam/Num_Meth_SS06/skript210406.pdf describes nicely how L can be found through permutations_combined (see below)
"""
function lu_decomp_pivot(M::Matrix)
    N = Matrix{Float64}(copy(M))
    n = size(M, 1)
    L = [(i == j) ? 1.0 : 0.0 for i in 1:n, j in 1:n]
    permutations = Vector{Permutation}(undef, n-1)
    for i in 1:(n-1)
        m = findmax(N[i:n, i])[2] + (i-1) # + (i-1) is needed due to the findmax function with start counting from 1 onwards along the shortened column
        if m == i
            permutations[i] = Permutation(n)
        else 
            swap_row!(N,m,i)
            permutations[i] = Transposition(n,m,i)
        end
        for k in (i+1):n
            L[k,i] = N[k,i]/N[i,i]
            for j in i:n
                N[k,j] = N[k,j] - L[k,i]*N[i,j]
            end
        end
    end
    permutations_combined = array_product(reverse(permutations),n)
    return L, N, permutations_combined
end
# lu-decomposition with pivoting  
# Permutationen in Julia, Permutationsmatrizen
