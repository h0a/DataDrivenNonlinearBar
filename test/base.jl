using Test, SparseArrays, LinearAlgebra

using DataDrivenNonlinearBar



@testset "LagrangePolynomials" begin
    # linear Lagrange polynomials
    L0, L1 = linearLagrangePolynomials()

    @test L0(-1) == 1
    @test L0(0) == 0.5
    @test L0(1) == 0

    @test L1(-1) == 0
    @test L1(0) == 0.5
    @test L1(1) == 1

    
    # 1st derivative of linear Lagrange polynomials
    dL0, dL1 = compute1stDeriv4linearLagrangePolynomials()
    xx = LinRange(-1,1,10)
    @test dL0.(xx) == -0.5 .* ones(length(xx))
    @test dL1.(xx) == 0.5 .* ones(length(xx))


    # constant functions
    L = constantFunctions()
    xx = LinRange(-1,1,10)
    @test L.(xx) == ones(length(xx))
end