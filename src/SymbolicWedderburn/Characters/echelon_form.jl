Base.@propagate_inbounds function _swap_rows!(A::AbstractMatrix, i, j)
    @boundscheck @assert i in axes(A, 1)
    @boundscheck @assert j in axes(A, 1)

    i == j && return A

    @inbounds for col_idx in axes(A, 2)
        A[i, col_idx], A[j, col_idx] = A[j, col_idx], A[i, col_idx]
    end
    return A
end

Base.@propagate_inbounds function _mul_row!(
    A::AbstractMatrix,
    val,
    row_idx;
    starting_at = 1,
)
    @boundscheck @assert starting_at in axes(A, 2)

    @inbounds for col_idx in starting_at:last(axes(A, 2))
        A[row_idx, col_idx] *= val
    end
    return A
end

# Version over exact field:

Base.@propagate_inbounds function _find_pivot(
    A::AbstractMatrix,
    col_idx;
    starting_at = 1,
)
    k = findnext(!iszero, @view(A[:, col_idx]), starting_at)
    isnothing(k) && return false, col_idx
    return true, k
end

Base.@propagate_inbounds function _reduce_column_by_pivot!(
    A::AbstractMatrix,
    row_idx,
    col_idx;
    starting_at = 1,
)
    @boundscheck checkbounds(A, row_idx, col_idx)
    @boundscheck starting_at in axes(A, 2)

    for ridx in axes(A, 1)
        ridx == row_idx && continue
        v = A[ridx, col_idx]
        iszero(v) && continue
        @inbounds for cidx in starting_at:last(axes(A, 2))
            A[ridx, cidx] -= v * A[row_idx, cidx]
        end
    end
    return A
end

_finalize_pivot_reduce!(A::AbstractMatrix, pivot) = A

# version over AbstractFloat

const FloatOrComplex = Union{AbstractFloat,Complex{<:AbstractFloat}}

Base.@propagate_inbounds function _find_pivot(
    A::AbstractMatrix{T},
    col_idx;
    starting_at = 1,
    atol = eps(real(eltype(A))) * size(A, 1),
) where {T<:FloatOrComplex}
    isempty(starting_at:size(A, 1)) && return false, starting_at
    @boundscheck @assert starting_at in axes(A, 1)
    # find the largest entry in the column below the last pivot
    @boundscheck @assert col_idx in axes(A, 2)

    mval, midx = findmax(abs, @view A[starting_at:end, col_idx])
    if mval < atol # the largest entry is below threshold so we zero everything in the column
        @view(A[starting_at:end, col_idx]) .= zero(T)
        return false, starting_at
    end
    return true, oftype(starting_at, starting_at + midx - 1)
end

function _finalize_pivot_reduce!(
    A::AbstractSparseMatrix{T},
    pivot,
) where {T<:FloatOrComplex}
    m = T <: Complex ? 2abs(eps(T)) : eps(T)
    droptol!(A, max(size(A)...) * m)
    return A
end

Base.@propagate_inbounds function _reduce_column_by_pivot!(
    A::AbstractMatrix{T},
    row_idx,
    col_idx;
    starting_at = 1,
    atol = eps(real(eltype(A))) * size(A, 1),
) where {T<:FloatOrComplex}
    @boundscheck checkbounds(A, row_idx, col_idx)
    @boundscheck starting_at in axes(A, 2)

    for ridx in axes(A, 1)
        ridx == row_idx && continue
        @inbounds v = A[ridx, col_idx]
        if abs(v) < atol
            @inbounds A[ridx, col_idx] = zero(T)
            continue
        end
        @inbounds for cidx in starting_at:last(axes(A, 2))
            A[ridx, cidx] -= v * A[row_idx, cidx]
        end
    end
    return A
end

Base.@propagate_inbounds function row_echelon_form!(A::AbstractMatrix)
    pivots = Int[]
    row_idx = firstindex(axes(A, 1)) - 1 # also: current_rank
    @inbounds for col_idx in axes(A, 2)
        found, j = _find_pivot(A, col_idx; starting_at = row_idx + 1)
        found || continue
        row_idx += 1
        push!(pivots, col_idx)

        # swap the rows so that  A[row_idx, :] is the row with leading nz in
        # i-th column
        A = _swap_rows!(A, row_idx, j)

        w = inv(A[row_idx, col_idx])

        # multiply A[row_idx, :] by w
        A = _mul_row!(A, w, row_idx; starting_at = col_idx)
        # A[current_rank, col_idx] is 1 now

        # zero the whole col_idx-th column above and below pivot:
        # to the left of col_idx everything is zero
        A = _reduce_column_by_pivot!(A, row_idx, col_idx; starting_at = col_idx)
        A = _finalize_pivot_reduce!(A, col_idx)
    end
    return A, pivots
end

row_echelon_form(A::AbstractMatrix) = row_echelon_form!(deepcopy(A))
