using Test

tests = [
            "base"
        ]

@testset "DataDrivenNonlinearBar Unit-testing" begin

    # run all unit-tests
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end

end # @testset
