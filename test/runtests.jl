using Test
using Aqua
using NormalForms
using LinearAlgebra
using StaticArrays

import NormalForms: eye, gcd_kb

Aqua.test_all(NormalForms; project_toml_formatting=false)

@testset "NormalForms.jl" verbose=true begin
    @testset "Algorithms" begin
        # If this fails, Smith normal form calculation may not terminate:
        # The critical part is that for the return value (r,p,q), abs(p) > abs(q)
        @test gcd_kb(2,2) === (2,1,0)
        @test isunimodular(Float64[1 1; 0 1]) === true
        @test isunimodular([1/2 0; 0 2]) == false
        @test eye(transpose(zeros(Bool, 3,4)), 2) == transpose(collect(LinearAlgebra.I(3)))
        @test eye(adjoint(zeros(Bool, 3,4)), 2) == adjoint(collect(LinearAlgebra.I(3)))
        @test eye(Diagonal([1,2,3]), 2) == Diagonal([1,1,1])
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
    @testset "Square matrix with floats" begin
        M = Float64[-2 1 1; 2 -1 1; 2 1 -1]
        @test eltype(hnfr(M)) <: Integer
        @test eltype(hnfc(M)) <: Integer
        @test eltype(snf(M)) <: Integer
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
        Hrt = hnfr(transpose(M))
        Hct = hnfc(transpose(M))
        St = snf(transpose(M))
        # Hermite normal form tests
        # Transpose of the HNF of a transpose should give the opposite form (row vs. column)
        @test transpose(M) * Hct.U == Hct.H
        @test Hrt.U * transpose(M) == Hrt.H
        @test transpose(Hrt) == hnfc(M)
        @test transpose(Hct) == hnfr(M)
        # Smith normal form tests
        #= NOTE:
            The Smith normal form factorization doesn't return the exact same results when
            performed with a transpose vs. when it's transposed. This is probably due to the
            rounding modes of Euclidean divisions being irrelevant. Therefore, we don't always
            expect transpose(snf(transpose(M))) == snf(M).
        =#
        @test diag(snf(transpose(M)).S) == diag(snf(M).S)
        @test St.S == St.U * transpose(M) * St.V
        @test St.S == transpose(snf(M).V) * transpose(M) * transpose(snf(M).U)
        @test snf(M).S == transpose(St.V) * M * transpose(St.U)
    end
    @testset "Adjoints" begin
        M = [-2 1 1; 2 -1 1; 2 1 -1]
        Hra = hnfr(M')
        Hca = hnfc(M')
        Sa = snf(M')
        # Hermite normal form tests
        @test M' * Hca.U == Hca.H
        @test Hra.U * M' == Hra.H
        @test Hra' == hnfc(M)
        @test Hca' == hnfr(M)
        # Smith normal form tests
        # As above, snf(M')' != snf(M) in general
        @test diag(snf(M').S) == diag(snf(M).S)
        @test Sa.S == Sa.U * M' * Sa.V
        @test Sa.S == snf(M).V' * M' * snf(M).U'
        @test snf(M).S == Sa.V' * M * Sa.U'
    end
    @testset "Diagonal matrices" begin
        M = Diagonal([7,12,6])
        Hr = hnfr(M)
        Hc = hnfc(M)
        S = snf(M)
        @test S isa Smith{Int, Matrix{Int}}
        @test Hr isa RowHermite{Int, Diagonal{Int, Vector{Int}}}
        @test Hc isa ColumnHermite{Int, Diagonal{Int, Vector{Int}}}
        @test diag(S) == [1,6,84]
        @test Diagonal(S) == Diagonal([1,6,84])
        @test Diagonal(Hr) == M
        @test Diagonal(Hc) == M
    end
    @testset "Known problematic matrices" begin
        M = [2 0 0; 1 4 0; 2 0 8]
        S = snf(M)
        @test S.S == diagm([1, 8, 8])
        @test M == Int.(S.U \ S.S / S.V)
        for k in 1:3
            N = copy(M)
            NormalForms.zero_row_and_col!(N, NormalForms.eye(N,1), NormalForms.eye(N,2), k)
            @test NormalForms.is_row_zero_after(N, k)
            @test NormalForms.is_col_zero_after(N, k)
            @test NormalForms.detb(M) === NormalForms.detb(N)
        end
    end
end

#= Known problem matrices:
[ 9  7 -7  7 -3  7
 -1  5  0  7 -4 -1
  2  2  9  7 -5  3
 -2  6  3 -9 -3  5] (unimodularity failure in hnfc - also in HermiteNormalForm.jl)
=#
