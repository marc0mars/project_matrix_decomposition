using MyProject
using Test

SlowSum(a) = sum([i for i in 1:a])

@testset "ArithmeticSum" begin 
    @test ArithmeticSum(10) == SlowSum(10)
    @test ArithmeticSum(5) == SlowSum(5)
    @test ArithmeticSum(100) == SlowSum(100)
end