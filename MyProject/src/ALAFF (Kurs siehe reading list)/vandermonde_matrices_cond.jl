
using Plots
using LinearAlgebra
gr()

"""
    vandermonde_matrix(m, n)

    m: number of rows (number of points)
    n: n-1 is the degree of the polynomial

    This is a recursive function but slower than vandermonde_matrix2 and should not be used.
"""
function vandermonde_matrix(m,n) # where m is the number of points between 0 and 1; and n-1 is the degree of the Polynomial
    if n == 1
        return [1 for i in 1:m, j in 1:1]
    end
    increase_matrix = [(i == j) ? (1/m)*i : 0 for i in 0:(m-1), j in 0:(m-1)]
    return hcat(fill(1.0, m), increase_matrix * vandermonde_matrix(m,n-1))
end

"""
    vandermonde_matrix2(m,n)

    m: number of rows
    n: n - 1 is the degree of the polynomial

    This is currently my fastest function for creating a m x n Van.-Matrix 
"""
function vandermonde_matrix2(m,n) # where m is the number of points between 0 and 1; and n-1 is the degree of the Polynomial
    if n == 1
        return [1.0 for i in 1:m, j in 1:1]
    end
    last_result = vandermonde_matrix2(m,n-1)
    for k in 1:m
        view(last_result, k, :) .= view(last_result, k, :) .* (k-1)/m
    end
    return hcat(fill(1.0, m), last_result)
end

"""
    testrun

    This function creates the data for the final plot
"""
function testrun(max)
    m = 5000  # that means there are m different points
    test_result = []
    for n in 1:max
        push!(test_result, cond(vandermonde_matrix2(m,n)))
    end
    return test_result
end

"""
    result_plot

    This function plots the final results, where the number n is increased for a steady m (=max)
"""
function result_plot(max)
    test_result = testrun(max)
    increasing_n = [i for i in 1:max]
    plot(increasing_n, test_result; ylabel = "condition number of polynomial", xlabel = "degree of polynomial", yaxis=:log)
    scatter!(increasing_n, test_result)
end

result_plot(20)