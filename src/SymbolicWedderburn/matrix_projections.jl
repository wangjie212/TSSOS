## preallocation

function _preallocate(::Type{T}, sizes::Tuple, sizehint) where {T}
    res = spzeros(T, sizes...)
    @static if VERSION >= v"1.7"
        sizehint!(res, first(sizes) * sizehint)
    end

    return res
end

__hint(χ::Character) = length(conjugacy_classes(χ))
__hint(α::AlgebraElement) = count(!iszero, SA.coeffs(α))

_projection_size(m::AbstractMatrix) = size(m)
_projection_size(::Nothing, χ::Character) = (d = degree(parent(χ)); (d, d))
function _projection_size(::Nothing, α::AlgebraElement)
    return (d = degree(parent(parent(α))); (d, d))
end
_projection_size(hom::InducedActionHomomorphism, _) = (d = degree(hom); (d, d))

function preallocate(::Type{T}, χ::Union{Character,AlgebraElement}) where {T}
    return preallocate(T, nothing, χ)
end

function preallocate(
    hom::InducedActionHomomorphism,
    χ::Union{Character,AlgebraElement},
)
    T = Base._return_type(*, Tuple{coeff_type(hom),eltype(χ)})
    return preallocate(T, hom, χ)
end

function preallocate(
    ::Type{T},
    a::Union{Nothing,InducedActionHomomorphism},
    χ::Union{Character,AlgebraElement},
) where {T}
    sizes = _projection_size(a, χ)
    return _preallocate(T, sizes, __hint(χ))
end

## matrix projection [irreducible]
"""
    matrix_projection([hom::InducedActionHomomorphism, ]χ::Character{T})
Compute matrix projection associated to character `χ`.

Returned `M<:AbstractMatrix{T}` of size `(d, d)` where the degree `d` of the
projecion is determined by elements in `conjugacy_classes(χ)`. E.g. `d` could
be equal to the `degree` when conjugacy classes consist of `AbstractPermutation`s.
If the homomorphism is passed, the dimension will be derived in similar
manner from the elements of the image of the homomorphism.

The precise type of `M` can be altered by overloading

```
preallocate(::Type{T}, [hom::InducedActionHomomorphism, ]χ::Character)
```
"""
function matrix_projection(χ::Character{T}) where {T}
    return _mproj_outsT!(preallocate(T, χ), χ)
end

function matrix_projection(
    hom::InducedActionHomomorphism,
    χ::Character{T},
) where {T}
    return _mproj_outsT!(preallocate(hom, χ), hom, χ)
end

function matrix_projection(χ::Character{T}) where {T<:Rational}
    if all(isone ∘ Cyclotomics.conductor, table(χ))
        return _mproj_fitsT!(preallocate(T, χ), χ)
    else
        return _mproj_outsT!(preallocate(T, χ), χ)
    end
end

function matrix_projection(χ::Character{T}) where {T<:Union{Cyclotomic,Complex}}
    return _mproj_fitsT!(preallocate(T, χ), χ)
end
function matrix_projection(
    hom::InducedActionHomomorphism,
    χ::Character{T},
) where {T<:Union{Cyclotomic,Complex}}
    return _mproj_fitsT!(preallocate(hom, χ), hom, χ)
end

function matrix_projection(χ::Character{T}) where {T<:AbstractFloat}
    if all(isreal, table(χ))
        return _mproj_fitsT!(preallocate(T, χ), χ)
    else
        return _mproj_outsT!(preallocate(T, χ), χ)
    end
end

function matrix_projection(
    hom::InducedActionHomomorphism,
    χ::Character{T},
) where {T<:AbstractFloat}
    if all(isreal, table(χ))
        return _mproj_fitsT!(preallocate(hom, χ), hom, χ)
    else
        return _mproj_outsT!(preallocate(hom, χ), hom, χ)
    end
end

_mproj_fitsT!(args...) = _mproj_outsT!(args...)

function _mproj_outsT!(mproj::AbstractMatrix{T}, χ::Character) where {T}
    mproj .= sum(
        c .* matrix_projection_irr(Character(table(χ), i)) for
        (i, c) in pairs(multiplicities(χ)) if !iszero(c)
    )
    return mproj
end

function _mproj_outsT!(
    mproj::AbstractMatrix{T},
    hom::InducedActionHomomorphism,
    χ::Character,
) where {T}
    mproj .= sum(
        c .* matrix_projection_irr(hom, Character(table(χ), i)) for
        (i, c) in pairs(multiplicities(χ)) if !iszero(c)
    )
    return mproj
end

"""
    matrix_projection_irr([::Type, ]χ::Character)
    matrix_projection_irr(hom::InducedActionHomomorphism, χ::Character)
Compute matrix projection associated to an *irreducible* character `χ`.

The returned matrix defines so called *isotypical* projection.

See also [matrix_projection](@ref).
"""
matrix_projection_irr(χ::Character{T}) where {T} = matrix_projection_irr(T, χ)
function matrix_projection_irr(::Type{T}, χ::Character) where {T}
    return matrix_projection_irr_acc!(preallocate(T, χ), χ, 1)
end

function matrix_projection_irr(hom::InducedActionHomomorphism, χ::Character)
    return matrix_projection_irr_acc!(preallocate(hom, χ), hom, χ, 1)
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    χ::Character,
    weight,
)
    @assert isirreducible(χ)
    LinearAlgebra.checksquare(result)
    ccls = conjugacy_classes(χ)
    vals = collect(values(χ))
    weight *= degree(χ) // sum(length, ccls)
    result = matrix_projection_irr_acc!(result, vals, ccls, weight)
    return result
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    vals,
    ccls::AbstractVector{<:PG.AbstractOrbit{<:AP.AbstractPermutation}},
    weight,
)
    iszero(weight) && return result
    I = UInt32[]
    J = UInt32[]
    V = eltype(result)[]
    for (val, cc) in zip(vals, ccls)
        iszero(val) && continue
        w = weight * val
        for g in cc
            for i in 1:size(result, 1)
                push!(I, i)
                push!(J, i^g)
                push!(V, w)
            end
        end
    end
    result += sparse(I, J, V)
    return result
end

function matrix_projection_irr_acc!(
    res::AbstractMatrix,
    vals,
    ccls::AbstractVector{<:PG.AbstractOrbit{<:AbstractMatrix}},
    weight,
)
    # TODO: call to inv(Matrix(g)) is a dirty hack, since if `g`
    # is given by a sparse matrix `inv(g)` will fail.
    for (val, cc) in zip(vals, ccls)
        iszero(val) && continue
        for g in cc
            res .+= (weight * val) .* inv(convert(Matrix, g))
        end
    end
    return res
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism,
    χ::Character,
    weight,
)
    @assert isirreducible(χ)
    LinearAlgebra.checksquare(result)
    vals = collect(values(χ))
    ccls = conjugacy_classes(χ)
    weight *= degree(χ) // sum(length, ccls)
    result = matrix_projection_irr_acc!(result, hom, vals, ccls, weight)
    return result
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism{<:ByPermutations},
    class_values,
    conjugacy_cls,
    weight,
)
    iszero(weight) && return result
    I = UInt32[]
    J = UInt32[]
    V = eltype(result)[]
    for (val, ccl) in zip(class_values, conjugacy_cls)
        iszero(val) && continue
        w = weight * val
        for g in ccl
            h = induce(hom, g)
            @assert h isa PG.Perm
            for i in 1:size(result, 1)
                push!(I, i)
                push!(J, i^h)
                push!(V, w)
            end
        end
    end
    result += sparse(I, J, V)
    return result
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    class_values,
    conjugacy_cls,
    weight,
)
    iszero(weight) && return result
    for (val, cc) in zip(class_values, conjugacy_cls)
        iszero(val) && continue
        w = weight * val
        for g in cc
            result .+= w .* induce(hom, inv(g))
        end
    end
    return result
end

"""
    matrix_representation([::Type{T}, ]α::AlgebraElement)
    matrix_representation([::Type{T}, ]hom::InducedActionHomomorphism, α::AlgebraElement)
Compute matrix representative of the given algebra element `α`.

If `α` is an element of the group (or monoid) algebra `RG` of `G`, then every
linear representation `π : G → GL(V)` gives rise to an algebra homomorphism
`π : RG → B(V)` (bounded operators on `V`). This function evaluates this
representation given either by
 * the natural action of `G` on some `V` (e.g. every permutation group of degree
 `d` acts on `d`-dimensional `V` by basis vector permutation), or
 * by action homomorphism `hom`.

See also [matrix_projection](@ref).
"""
matrix_representation(α::AlgebraElement) = matrix_representation(eltype(α), α)
function matrix_representation(::Type{T}, α::AlgebraElement) where {T}
    return matrix_representation_acc!(preallocate(T, α), α)
end

function matrix_representation(
    hom::InducedActionHomomorphism,
    α::AlgebraElement,
)
    return matrix_representation_acc!(preallocate(hom, α), hom, α)
end

function matrix_representation_acc!(
    result::AbstractMatrix,
    α::AlgebraElement{
        <:StarAlgebra{<:PG.AbstractPermutationGroup},
    },
)
    b = basis(parent(α))
    I = UInt32[]
    J = UInt32[]
    V = eltype(result)[]
    for (idx, val) in SA.nonzero_pairs(SA.coeffs(α))
        g = b[idx]
        iszero(val) && continue
        for i in 1:size(result, 1)
            push!(I, i)
            push!(J, i^g)
            push!(V, val)
        end
    end
    result += sparse(I, J, V)
    return result
end

function matrix_representation_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism{<:ByPermutations},
    α::AlgebraElement,
)
    b = basis(parent(α))
    I = UInt32[]
    J = UInt32[]
    V = eltype(result)[]
    for (idx, val) in SA.nonzero_pairs(SA.coeffs(α))
        iszero(val) && continue
        g = induce(hom, b[idx])
        @assert g isa PG.Perm
        for i in 1:size(result, 1)
            push!(I, i)
            push!(J, i^g)
            push!(V, val)
        end
    end
    result += sparse(I, J, V)
    return result
end

function matrix_representation_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    α::AlgebraElement,
)
    b = basis(parent(α))
    for (idx, val) in SA.nonzero_pairs(SA.coeffs(α))
        iszero(val) && continue
        result .+= val .* induce(hom, inv(b[idx]))
    end

    return result
end

# convenience for uniform calls below
function matrix_representation(χ::Character)
    return isirreducible(χ) ? matrix_projection_irr(χ) : matrix_projection(χ)
end

function matrix_representation(hom::InducedActionHomomorphism, χ::Character)
    return isirreducible(χ) ? matrix_projection_irr(hom, χ) :
           matrix_projection(hom, χ)
end
