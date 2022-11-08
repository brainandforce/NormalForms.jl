using .StaticArrays

# Matrix appears slightly faster, but MMatrix allocates less memory
# TODO: perhaps extend this for the sake of StaticArrays?
detb(M::SMatrix) = detb!(convert(MMatrix, M))

function hnfc(M::SMatrix{D1,D2,<:Integer}, R::RoundingMode = NegativeOffDiagonal) where {D1,D2}
    (H, U, info) = hnf_ma!(MMatrix(M), R)
    return ColumnHermite(SMatrix{D1,D2}(H), SMatrix{D1,D2}(U), info)
end

function hnfr(M::SMatrix{D1,D2,<:Integer}, R::RoundingMode = NegativeOffDiagonal) where {D1,D2}
    (H, U, info) = hnf_ma!(MMatrix(M'), R)
    return RowHermite(SMatrix{D1,D2}(H'), SMatrix{D1,D2}(U'), info)
end

snf(M::SMatrix{D1,D2,<:Integer}) where {D1,D2} = Smith(snf_ma!(MMatrix(M))...)
