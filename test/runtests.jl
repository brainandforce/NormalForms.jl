using Test
using NormalForms
using StaticArrays

@testset "NormalForms.jl" begin
    # Write your tests here.
    @testset "Known square matrix" begin
        M = [-2 1 1; 2 -1 1; 2 1 -1]
        Fr = hnfr(M)
        Fc = hnfc(M)
        S = snf(M)
        @test M * Fc.U == Fc.H
        @test Fr.U * M == Fr.H
        @test S.U * M * S.V == S.S
    end
    @testset "SMatrix" begin
        S = SMatrix{3,3,Int}([-2 1 1; 2 -1 1; 2 1 -1])
        Fr = hnfr(S)
        Fc = hnfc(S)
        @test isbits(Fr)
        @test S * Fc.U == Fc.H
        @test Fr.U * S == Fr.H
    end
end

#= Known problem matrices:
[ 9  7 -7  7 -3  7
 -1  5  0  7 -4 -1
  2  2  9  7 -5  3
 -2  6  3 -9 -3  5] (hnfc fails - also in HermiteNormalForm.jl)

 Anything with [1 1; 0 1] seems to break snf()
=#
