using .StaticArrays

# Matrix appears slightly faster, but MMatrix allocates less memory
# TODO: perhaps extend this for the sake of StaticArrays?
detb(M::StaticMatrix) = detb!(convert(MMatrix, M))

function hnfc(M::StaticMatrix{D1,D2,<:Integer}, R::RoundingMode = NegativeOffDiagonal) where {D1,D2}
    (H, U, info) = hnf_ma!(MMatrix(M), R)
    return ColumnHermite(SMatrix{D1,D2}(H), SMatrix{D2,D2}(U), info)
end

function hnfr(M::StaticMatrix{D1,D2,<:Integer}, R::RoundingMode = NegativeOffDiagonal) where {D1,D2}
    (H, U, info) = hnf_ma!(MMatrix(M'), R)
    return RowHermite(SMatrix{D1,D2}(H'), SMatrix{D1,D1}(U'), info)
end

function snf(M::StaticMatrix{D1,D2,<:Integer}) where {D1,D2}
    result = snf_ma!(MMatrix(M))
    return Smith(
        SMatrix{D1,D2}(result[1]),
        SMatrix{D1,D1}(result[2]),
        SMatrix{D2,D2}(result[3]),
        result[4]
    )
end
