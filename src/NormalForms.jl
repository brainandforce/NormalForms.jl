module NormalForms

using LinearAlgebra
using Requires

include("algorithms.jl")
export isunimodular
include("hermite.jl")
export AbstractHermite, RowHermite, ColumnHermite
export NegativeOffDiagonal, PositiveOffDiagonal
export hnfc!, hnfc, hnfr!, hnfr
include("smith.jl")
export Smith
export snf!, snf

function __init__()
    # Contains methods specific to SMatrix, which can't be mutated in-place
    @require StaticArrays="90137ffa-7385-5640-81b9-e52037218182" begin
        include("staticarrays.jl")
    end
end

end
