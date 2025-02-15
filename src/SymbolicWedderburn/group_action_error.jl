### Basic checks if the action defined by the user is really a group action.
"""
    GroupActionError
`GroupActionError` is an indicator that the group action defined by
`(G, action, basis)` might be incorrect.

The checks in `SymbolicWedderburn` only make sure that the group action is well
defined. To verify the correctness of the definition of a group you may
additionally run
```julia
using GroupsCore
include(joinpath(pathof(GroupsCore), "..", "..", "test", "conformance_test.jl"))
let G = ...
    test_Group_interface(G)
    test_GroupElement_interface(rand(G, 2)...)
end
```
"""
struct GroupActionError{A} <: Exception
    action::A
    msg::String
end

function Base.showerror(io::IO, err::GroupActionError)
    print(io, "Group action $(err.action) seems to be ill defined: ")
    return print(io, err.msg)
end

__isexact(::Type{<:Integer}) = true
__isexact(::Type{<:Rational{T}}) where {T} = __isexact(T)
__isexact(::Type{<:AbstractFloat}) = false
__isexact(::Type{<:Complex{T}}) where {T} = __isexact(T)
__isexact(::Type{<:Cyclotomic{T}}) where {T} = __isexact(T)

function __group_action_id(
    hom::InducedActionHomomorphism{<:ByPermutations},
    id::GroupElement,
    x,
)
    @assert isone(id)
    xᵉ = action(action(hom), id, x)
    return x == xᵉ, xᵉ
end

function __group_action_id(
    hom::InducedActionHomomorphism{<:BySignedPermutations},
    id::GroupElement,
    x,
)
    @assert isone(id)
    xᵉ, s = action(action(hom), id, x)
    return x == xᵉ && isone(s), xᵉ
end

function __group_action_id(
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    id::GroupElement,
    x,
)
    @assert isone(id)
    v1 = sparsevec(decompose(x, hom)..., length(basis(hom)))

    xᵉ = action(action(hom), id, x)
    v2 = sparsevec(decompose(xᵉ, hom)..., length(basis(hom)))

    T = coeff_type(hom)
    passed = __isexact(T) ? v1 == v2 : v1 ≈ v2

    return passed, xᵉ
end

function __group_action_right_assoc(
    hom::InducedActionHomomorphism{<:ByPermutations},
    g,
    h,
    x,
)
    xᵍʰ = action(action(hom), g * h, x)
    xᵍ = action(action(hom), g, x)
    xᵍꜝʰ = action(action(hom), h, xᵍ)
    return xᵍʰ == xᵍꜝʰ, xᵍʰ, xᵍꜝʰ
end

function __group_action_right_assoc(
    hom::InducedActionHomomorphism{<:BySignedPermutations},
    g,
    h,
    x,
)
    xᵍʰ, s_gh = action(action(hom), g * h, x)
    xᵍ, s_g = action(action(hom), g, x)
    xᵍꜝʰ, s_g_h = action(action(hom), h, xᵍ)

    passed = s_gh == s_g * s_g_h && xᵍʰ == xᵍꜝʰ
    return passed, (s_gh,  xᵍʰ), (s_g * s_g_h, xᵍꜝʰ)
end

function __group_action_right_assoc(
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    g,
    h,
    x,
)
    xᵍʰ = action(action(hom), g * h, x)
    v1 = sparsevec(decompose(xᵍʰ, hom)..., length(basis(hom)))

    xᵍ⁾ʰ = let xᵍ = action(action(hom), g, x)
        e, c = decompose(xᵍ, hom)
        sum(c[i] * action(action(hom), h, hom[e[i]]) for i in eachindex(c, e))
    end
    v2 = sparsevec(decompose(xᵍ⁾ʰ, hom)..., length(basis(hom)))

    T = coeff_type(hom)
    passed = __isexact(T) ? v1 == v2 : v1 ≈ v2
    return passed, xᵍʰ, xᵍ⁾ʰ
end

function __check_group_action_axioms(
    itr,
    hom::InducedActionHomomorphism,
    elts_idcs,
)
    id = one(first(itr))
    for idx in elts_idcs # checking that id really acts as identity
        x = basis(hom)[idx]
        passed, xᵉ = __group_action_id(hom, id, x)
        if !(passed)
            throw(
                GroupActionError(
                    action(hom),
                    "group identity `id = $(id)` fails to act as identity at `x = $(x)`:\n" *
                    "`x = $(x) ≠ $(xᵉ) = action(act, id, x)`",
                ),
            )
        end
    end

    for idx in elts_idcs
        x = basis(hom)[idx]
        for g in itr, h in itr
            h = inv(h)
            passed, xᵍʰ, xᵍ⁾ʰ = __group_action_right_assoc(hom, g, h, x)
            if !(passed)
                throw(
                    GroupActionError(
                        action(hom),
                        "`g = $(g)` and `h = $(h)` fail right-associativity at `x = $(x)` :\n " *
                        "`action(act, g * h, x) = $(xᵍʰ) ≠ " *
                        "$(xᵍ⁾ʰ) = action(act, h, action(act, g, x))`",
                    ),
                )
            end
        end
    end
    return true
end

function check_group_action(
    G::Group,
    hom::InducedActionHomomorphism;
    full_check = false,
)
    return if full_check
        __check_group_action_axioms(G, hom, 1:length(basis(hom)))
    else
        __check_group_action_axioms(
            gens(G),
            hom,
            1:min(length(basis(hom)), 16), # arbitrary, smallish number
        )
    end
end

function check_group_action(
    G::Group,
    act::Action,
    basis::SA.ExplicitBasis;
    full_check = false,
)
    ehom = SchreierExtensionHomomorphism(G, act, basis; memoize = false)
    return check_group_action(G, ehom; full_check = full_check)
end
