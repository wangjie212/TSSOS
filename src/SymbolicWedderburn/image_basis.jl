## Finding basis of the row-space (right image) of an AbstractMatrix

"""
    image_basis(A::AbstractMatrix)
    image_basis([hom::InducedActionHomomorphism, ]χ::Character[, rank])
    image_basis([hom::InducedActionHomomorphism, ]α::AlgebraElement[, rank])
Return basis of the row-space of a matrix.

For characters or algebra elements return a basis of the row-space of
`matrix_projection([hom, ]χ)` or `matrix_representation([hom, ]α)`.

Two methods are employed to achieve the goal.
* By default (symbolic, exact) row echelon form is produced, and therefore there
is no guarantee on the orthogonality of the returned basis vectors (rows).
* If `eltype(A) <: AbstractFloat` a thin `svd` (or `qr` in the sparse case)
decomposition is computed and the appropriate rows of its orthogonal factor are
returned, thus the basis is (numerically) orthonormal.

If the dimension (`=rank` of the matrix) is known beforehand (because of e.g.
algebra) it can be passed to `image_basis` to use more efficient methods.

# Examples:
```julia
julia> a = rand(-3:3, 3,3).//rand(1:4, 3,3);

julia> a[3, :] .= 2a[1, :] .- 1a[2, :]; a
3×3 Matrix{Rational{Int64}}:
  3//2   0//1  1//1
 -1//2  -2//1  1//2
  7//2   2//1  3//2

julia> ib = SymbolicWedderburn.image_basis(a)
2×3 Matrix{Rational{Int64}}:
 1//1  0//1   2//3
 0//1  1//1  -5//12

julia> ibf = SymbolicWedderburn.image_basis(float.(a))
2×3 Matrix{Float64}:
 -0.666651  0.416657  -0.618041
  0.529999  0.847998  -8.85356e-17

```
"""
image_basis(A::AbstractMatrix, rank = nothing) = image_basis!(deepcopy(A), rank)

function image_basis(α::Union{Character,AlgebraElement}, rank = nothing)
    mpr = matrix_representation(α)
    return image_basis!(mpr, rank)
end

function image_basis(
    hom::InducedActionHomomorphism,
    α::Union{Character,AlgebraElement},
    rank = nothing,
)
    mpr = matrix_representation(hom, α)
    return image_basis!(mpr, rank)
end

##
# the internal, (possibly) modifying versions

_eps(T::Type{<:AbstractFloat}) = eps(T)
_eps(::Type{Complex{T}}) where {T} = 2 * _eps(real(T))
_eps(::Type{<:Cyclotomic{T}}) where {T} = 2_eps(T)
_eps(::Any) = 0.0

function image_basis!(A::AbstractMatrix, rank = nothing)
    img, pivots = row_echelon_form!(A)
    if img isa AbstractSparseArray
        ε = _eps(eltype(img))
        if iszero(ε)
            dropzeros!(img)
        else
            droptol!(img, _eps(eltype(img)) * max(size(img)...))
        end
    end
    # TODO: orthogonalize the result
    return img[pivots, :]
end

function image_basis!(
    A::AbstractMatrix{T},
    ::Nothing,
) where {T<:Union{AbstractFloat,Complex}}
    F = svd!(A)
    tol = _eps(T) * first(F.S)
    rk = count(x -> x > tol, F.S)
    return Matrix((@view F.U[1:rk, :])')
end

function image_basis!(
    A::AbstractSparseMatrix{T},
    ::Nothing,
) where {T<:Union{AbstractFloat,Complex}}
    F = qr(A)
    rk = rank(F)
    img = let tmp = F.Q * Matrix(I, size(A, 1), rk)
        pinv = getfield(F, :rpivinv)
        sparse((@view tmp[pinv, :])')
    end
    return droptol!(img, _eps(T) * max(size(img)...))
end

# to disambiguate
function image_basis!(
    M::AbstractMatrix{T},
    rank::Integer,
) where {T<:Union{AbstractFloat,Complex}}
    img = image_basis!(M, nothing)
    if size(img, 1) ≠ rank
        @warn "Possibly wrong numerical rank: (numerical) $(size(img, 1)) ≠ $rank (expected)"
    end
    return img
end

function image_basis!(
    M::AbstractSparseMatrix{T},
    rank::Integer,
) where {T<:Union{AbstractFloat,Complex}}
    # return _image_basis!(M, rank)
    img = image_basis!(M, nothing)
    if size(img, 1) ≠ rank
        @warn "Possibly wrong rank estimate? (numerical) $(size(img, 1)) ≠ $rank (expected)"
    end
    return img
end
