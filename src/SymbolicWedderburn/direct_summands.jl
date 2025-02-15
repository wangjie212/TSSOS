struct DirectSummand{T,M<:AbstractMatrix{T},Ch} <: AbstractMatrix{T}
    basis::M
    multiplicity::Int
    character::Ch

    function DirectSummand(
        basis::M,
        multiplicity::Integer,
        character::Characters.Character,
    ) where {M<:AbstractMatrix}
        @assert size(basis, 1) < size(basis, 2)
        let (pr_rank, r) = divrem(size(basis, 1), multiplicity)
            @assert 1 ≤ pr_rank ≤ AP.degree(character)
            @assert r == 0
        end
        return new{eltype(M),M,typeof(character)}(
            basis,
            multiplicity,
            character,
        )
    end
end

# AbstractArray API
Base.size(ds::DirectSummand) = size(ds.basis)
Base.@propagate_inbounds function Base.getindex(ds::DirectSummand, i...)
    return image_basis(ds)[i...]
end

# Accessors
image_basis(ds::DirectSummand) = ds.basis
character(ds::DirectSummand) = ds.character
multiplicity(ds::DirectSummand) = ds.multiplicity

AP.degree(ds::DirectSummand) = AP.degree(character(ds))
function projection_rank(ds::DirectSummand)
    return div(size(image_basis(ds), 1), multiplicity(ds))
end
issimple(ds::DirectSummand) = isone(projection_rank(ds))

function SparseArrays.droptol!(
    ds::DirectSummand{T,<:AbstractSparseArray},
    tol,
) where {T}
    droptol!(image_basis(ds), tol)
    return ds
end
