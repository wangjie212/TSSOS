"""
    AbstractClassFunction{T}
Abstract type representing functions constant on conjugacy classes of a group.

The following functionality is required for an `AbstractClassFunction`:
 * `parent(χ)`: the underlying group
 * `conjugacy_classes(χ)`: iterator over conjugacy classes of `χ`.
 * `values(χ)`: a _typed_ iterator over values
 * `getindex(χ, i::Integer)`: the value on `i`-th conjugacy class.
 Note: Indexing with negative integers should return values on the class which
 contains inverses of the `i`-th class.
"""
abstract type AbstractClassFunction{T} end # <:AbstractVector{T}?

Base.eltype(::AbstractClassFunction{T}) where {T} = T

function LinearAlgebra.dot(χ::AbstractClassFunction, ψ::AbstractClassFunction)
    @assert conjugacy_classes(χ) === conjugacy_classes(ψ)
    CC = conjugacy_classes(χ)
    val = sum(length(cc) * χ[i] * ψ[-i] for (i, cc) in enumerate(CC))
    orderG = sum(length, conjugacy_classes(χ))
    val = _div(val, orderG)
    return val
end

function (χ::AbstractClassFunction)(g::GroupElement)
    for (i, cc) in enumerate(conjugacy_classes(χ))
        g ∈ cc && return χ[i]
    end
    throw(DomainError(g, "element does not belong to conjugacy classes of $χ"))
end

_div(val, orderG) = div(val, orderG)
_div(val::AbstractFloat, orderG) = val / orderG
_div(val::Complex{<:AbstractFloat}, orderG) = val / orderG

## Arbitrary ClassFunctions without decomposition into irreps
struct ClassFunction{T,CCl} <: AbstractClassFunction{T}
    vals::Vector{T}
    conjugacy_classes::CCl
    inv_of::Vector{Int}
end

ClassFunction(vals, cclasses) = ClassFunction(vals, cclasses, _inv_of(cclasses))

## AbstractClassFunction api
Base.parent(χ::ClassFunction) = parent(first(first(conjugacy_classes(χ))))
conjugacy_classes(χ::ClassFunction) = χ.conjugacy_classes
Base.values(χ::ClassFunction) = χ.vals

Base.@propagate_inbounds function Base.getindex(χ::ClassFunction, i::Integer)
    i = i < 0 ? χ.inv_of[-i] : i
    @boundscheck 1 ≤ i ≤ length(conjugacy_classes(χ))
    return values(χ)[i]
end

