using Test
using Aqua
using NormalForms
using LinearAlgebra
using StaticArrays

Aqua.test_all(NormalForms; project_toml_formatting=false)

@testset "NormalForms.jl" verbose=true begin
    @testset "Algorithms" begin
        # If this fails, Smith normal form calculation may not terminate:
        # The critical part is that for the return value (r,p,q), abs(p) > abs(q)
        @test NormalForms.gcd_kb(2,2) === (2,1,0)
        @test isunimodular(Float64[1 1; 0 1]) === true
        @test isunimodular([1/2 0; 0 2]) == false
    end
    @testset "Factorization assertions" begin
        # Diagonal checks for Smith normal form
        @test_throws AssertionError Smith([0 1; 1 0], diagm([1,1]), diagm([1,1]))
        # Unimodular checks
        @test_throws AssertionError Smith(diagm([1,2,3]), diagm([1,2,3]), diagm([1,1,1]))
        @test_throws AssertionError Smith(diagm([1,2,3]), diagm([1,1,1]), diagm([1,2,3]))
        @test_throws AssertionError RowHermite(diagm([1,2,3]), diagm([1,2,3]))
        @test_throws AssertionError ColumnHermite(diagm([1,2,3]), diagm([1,2,3]))
    end
    @testset "Component destructuring" begin
        M = [-2 1 1; 2 -1 1; 2 1 -1]
        Hc = hnfc(M)
        S = snf(M)
        @test iterate(Hc, Val{:H}()) === (Hc.H, Val{:U}())
        @test iterate(Hc, Val{:U}()) === (Hc.U, nothing)
        @test isnothing(iterate(Hc, nothing))
        @test iterate(S, Val{:S}()) === (S.S, Val{:U}())
        @test iterate(S, Val{:U}()) === (S.U, Val{:V}())
        @test iterate(S, Val{:V}()) === (S.V, nothing)
        @test isnothing(iterate(S, nothing))
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
    @testset "Non-square matrix" begin
        M = [1 2 -2 -1 -6 -4; -3 -3 -8 -5 -1 2; 4 9 7 7 -1 0; -9 2 2 9 -3 -6]
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
 -2  6  3 -9 -3  5] (unimodularity failure in hnfc - also in HermiteNormalForm.jl)
=#
