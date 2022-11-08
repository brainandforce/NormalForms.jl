using Test
using Aqua
using NormalForms
using LinearAlgebra
using StaticArrays

Aqua.test_all(NormalForms; project_toml_formatting=false)

@testset "NormalForms.jl" verbose=true begin
    @testset "Factorization assertions" begin
        # Diagonal checks for Smith normal form
        @test_throws AssertionError Smith([0 1; 1 0], diagm([1,1]), diagm([1,1]))
        # Unimodular checks
        @test_throws AssertionError Smith(diagm([1,2,3]), diagm([1,2,3]), diagm([1,1,1]))
        @test_throws AssertionError Smith(diagm([1,2,3]), diagm([1,1,1]), diagm([1,2,3]))
        @test_throws AssertionError RowHermite(diagm([1,2,3]), diagm([1,2,3]))
        @test_throws AssertionError ColumnHermite(diagm([1,2,3]), diagm([1,2,3]))
    end
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
        M = SMatrix{3,3,Int}([-2 1 1; 2 -1 1; 2 1 -1])
        Fr = hnfr(M)
        Fc = hnfc(M)
        S = snf(M)
        @test isbits(Fr)
        @test M * Fc.U == Fc.H
        @test Fr.U * M == Fr.H
        @test S.U * M * S.V == S.S
    end
    @testset "Transposes" begin
        M = [-2 1 1; 2 -1 1; 2 1 -1]
        Fr = hnfr(transpose(M))
        Fc = hnfc(transpose(M))
        S = snf(transpose(M))
        @test transpose(M) * Fc.U == Fc.H
        @test Fr.U * transpose(M) == Fr.H
        @test S.U * transpose(M) * S.V == S.S
    end
    @testset "Adjoints" begin
        M = [-2 1 1; 2 -1 1; 2 1 -1]
        Fr = hnfr(M')
        Fc = hnfc(M')
        S = snf(M')
        @test M' * Fc.U == Fc.H
        @test Fr.U * M' == Fr.H
        @test S.U * M' * S.V == S.S
    end
end

#= Known problem matrices:
[ 9  7 -7  7 -3  7
 -1  5  0  7 -4 -1
  2  2  9  7 -5  3
 -2  6  3 -9 -3  5] (hnfc fails - also in HermiteNormalForm.jl)

 Anything with [1 1; 0 1] seems to break snf()
=#
