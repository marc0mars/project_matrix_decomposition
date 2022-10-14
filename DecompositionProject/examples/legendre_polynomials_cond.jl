#= 
In the following    n-1 refers to the degree of the polynomial starting with degree 0
                    m refers to the number of points
m = 5000 which means we look at 5000 thousand equally spaced points in [0,1] starting with 0 and ending at 5000
=#
using LinearAlgebra
using BenchmarkTools


"""
    creating_mat_legendre_pol

    m: number of x-values
    n-1: degree of Polynomial

This function is way more efficient than those in vandermonde_matrices_cond. 
"""
function creating_mat_legendre_pol(m, n)
    p = [k/(m-1) for k in 0:m-1]        # points between 0 and 
    L = fill(0.0, m, n)                 # creates a zero matrix
    L[:,1] = ones(m)                    # first row are simply ones
    (n == 1) && return L                # if n == 1 we are already finished
    L[:,2] = 2p .- 1                    # the second row is by definition the former
    for k in 3:n                        # range is empty if n == 2 and L will be returned
        L[:,k] = (((2(k-2)+1) .* ((2*p).-1) .* L[:,k-1]) .- (k-2) .* L[:,k-2])./(k-1)
        # create the base for recursion (P₀(x) = 1, P₁(x) = 2x - 1, ..., P_n+1(x) = ((2n+1)(2x+1)Pₙ(x)-nP_n-1(x))/(n+1)
        # be aware of some subtle indexshifts due to Julia counting from 1 onwards
    end
    return L
end

function test_run(num_x_vals, max_n)
    test_set = [i for i in 1:max_n]
    test_result = [cond(creating_mat_legendre_pol(num_x_vals, k)) for k in 1:max_n]
    plot(test_set, test_result; yaxis=:log, xlabel = "degree of polynomial + 1", ylabel = "condition_number")
    scatter!(test_set, test_result)

    # the columns are not far from orthogonal
    # X = creating_mat_legendre_pol(5000,5)
    # display(transpose(X)*X)
    #d isplay(opnorm(transpose(X)*X-Diagonal(X)))
end

# test_run(5000,20)
# just for fun
# display(@benchmark creating_mat_legendre_pol(5000,10))