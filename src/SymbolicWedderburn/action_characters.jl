## characters defined by actions/homomorphisms
function nfixedpoints(p::AP.AbstractPermutation, deg = degree(p))
    return count(i -> i^p == i, Base.OneTo(deg))
end

function _action_class_fun(
    conjugacy_cls::AbstractVector{CCl},
) where {CCl<:PG.AbstractOrbit{<:AP.AbstractPermutation}}
    deg = mapreduce(cc -> maximum(degree, cc), max, conjugacy_cls)
    vals = map(conjugacy_cls) do cc
        repr = first(cc)
        return nfixedpoints(repr, deg)
    end
    # in general:
    # vals = [tr(matrix_representative(first(cc))) for cc in conjugacy_cls]
    return Characters.ClassFunction(vals, conjugacy_cls)
end

function _action_class_fun(
    conjugacy_cls::AbstractVector{CCl},
) where {CCl<:PG.AbstractOrbit{<:AbstractMatrix}}
    vals = [tr(first(cc)) for cc in conjugacy_cls]
    return Characters.ClassFunction(vals, conjugacy_cls)
end

function _action_class_fun(
    hom::InducedActionHomomorphism{<:ByPermutations},
    conjugacy_cls,
)
    deg = AP.degree(hom)
    vals = map(conjugacy_cls) do cc
        repr = induce(hom, first(cc))
        return nfixedpoints(repr, deg)
    end
    # in general:
    # vals = [tr(matrix_representative(first(cc))) for cc in conjugacy_cls]
    return Characters.ClassFunction(vals, conjugacy_cls)
end

function _action_class_fun(
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    conjugacy_cls,
)
    vals = [tr(induce(hom, first(cc))) for cc in conjugacy_cls]
    return Characters.ClassFunction(vals, conjugacy_cls)
end

"""
    action_character([::Type{T}, ]conjugacy_cls, tbl::CharacterTable)
Return the character of the representation given by the elements in the
conjugacy classes `conjugacy_cls`.

This corresponds to the classical definition of characters as a traces of the
representation matrices.
"""
function action_character(conjugacy_clss, tbl::CharacterTable)
    return action_character(eltype(tbl), conjugacy_clss, tbl)
end

function action_character(
    ::Type{T},
    conjugacy_clss,
    tbl::CharacterTable,
) where {T}
    ac_char = _action_class_fun(conjugacy_clss)
    vals = Characters.decompose(Int, ac_char, tbl)
    return Character{T}(tbl, vals)
end

function action_character(
    ::Type{T},
    conjugacy_cls,
    tbl::CharacterTable,
) where {T<:Union{AbstractFloat,ComplexF64}}
    ac_char = _action_class_fun(conjugacy_cls)
    vals = Characters.decompose(T, ac_char, tbl)
    return Character{T}(tbl, round.(Int, abs.(vals)))
end

"""
    action_character([::Type{T}, ]hom::InducedActionHomomorphism, tbl::CharacterTable)
Return the character of the representation given by the images under `hom` of
elements in `conjugacy_classes(tbl)`.

This corresponds to the classical definition of characters as a traces of the
representation matrices.
"""
function action_character(hom::InducedActionHomomorphism, tbl::CharacterTable)
    return action_character(eltype(tbl), hom, tbl)
end

function action_character(
    ::Type{T},
    hom::InducedActionHomomorphism,
    tbl::CharacterTable,
) where {T}
    act_character = _action_class_fun(hom, conjugacy_classes(tbl))
    vals = Characters.decompose(Int, act_character, tbl)
    return Character{T}(tbl, vals)
end
