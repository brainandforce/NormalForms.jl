using NormalForms
using Test

@testset "NormalForms.jl" begin
    # Write your tests here.
    @testset "Known square matrix" begin
        M = [-2 1 1; 2 -1 1; 2 1 -1]
        Fr = hnfr(M)
        Fc = hnfc(M)
        @test M * Fc.U == Fc.H
        @test Fr.U * M == Fr.H
    end
end
