# abstract type Action is defined in ext_homomorphisms

"""
    ByPermutations <: Action
A type of action where a group acts through permutations on a set.

It means that for every group element `g ∈ G` and every `s ∈ S`
```
action(act::ByPermutations, g, s) = s′
```
where `s′ ∈ S` is a (potentially different) element of `S`. By fixing an order
of elements in `S` one can then represent `g` as a permutation of degree `|S|`.
"""
abstract type ByPermutations <: Action end

"""
    ByLinearTransformation <: Action
A type of action where a group acts through linear transformations on an
(implicit) Euclidean space `ℝᴮ` with basis `B = (e₁, … eₙ)`.

That means that for every group element `g ∈ G` and every `eₖ ∈ B`
```
action(act::ByLinearTransformation, g, eₖ) = v
```
where `v = a₁e₁ + ⋯ + aₙeₙ`. Additionally [`decompose`](@ref) must be
implemented to return a (possibly sparse) decomposition of `v` in `B`.
"""
abstract type ByLinearTransformation <: Action end

"""
    decompose(v, hom::InducedActionHomomorphism)
Decompose element `v` in basis `basis(hom)` provided by `hom`.

Let `B = basis(hom)::StarAlgebras.AbstractBasis`. Then `v` should be decomposed
as a unique linear combination `v = a₁b₁ + ⋯ aₙbₙ` and the indices of `bᵢ`s in `B`
and a vector of coefficients `A` returned, so that

```
b, c = decompose(v, hom)
@assert sparsevec(b, c) == v
```

!!! note
    For performance reasons it is best to drop zeros in `c`, i.e. return a
    possibly sparse representation.

See also [`ByLinearTransformation`](@ref).
"""
function decompose(x, hom::InducedActionHomomorphism)
    throw("""No fallback is provided for $(typeof(x)). You need to implement
          `decompose(::$(typeof(x)), ::$(typeof(hom)))`.""")
end

"""
    BySignedPermutations <: ByLinearTransformation
A type of action where a group acts through permutations _with sign_ on an
(implicit) Euclidean space `ℝᴮ` with basis `B = (e₁, …, eₙ)`.

It means that for every group element `g ∈ G` and every `eₖ ∈ B` the result of
action of `g` on `eₖ` is `u·eₗ` for some `1≤l≤n`, where `u` is a root of unity
(i.e. `±1` in the real case). To accomplish this it is required that
```
action(act::BySignedPermutations, g, eₖ) == (eₗ, u)
```

!!! warning
    Only support for the real case (i.e. `u = ±1`) is implemented at the moment.

    Please consult the docstring of [ByLinearTransformation](@ref) for the necessary methods.
"""
abstract type BySignedPermutations <: ByLinearTransformation end

## coeff_types
coeff_type(::ByPermutations) = Int
function coeff_type(ac::ByLinearTransformation)
    throw("No fallback is provided for $(typeof(ac)). You need to implement
          `coeff_type(::$(typeof(ac)))`.")
end
coeff_type(::BySignedPermutations) = Int # lets not worry about roots of unity

## actions on AbstractVectors

function action(::ByPermutations, g::AP.AbstractPermutation, v::AbstractVector)
    # permuting coordinates
    return [v[i^g] for i in eachindex(v)]
end

"""
    action(hom::InducedActionHomomorphism, g::GroupElement, x)
Return the result of `g` acting on `x` through action homomorphism `hom`.

This can be understood as first evaluating the homomorphism: `h = hom(g)` and
then computing `x^h`, the action of the result on `x`.
"""
function action(
    hom::InducedActionHomomorphism{<:ByPermutations},
    g::GroupElement,
    v::AbstractVector,
)
    return action(action(hom), induce(hom, g), v)
end

function action(
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    g::GroupElement,
    v::AbstractVector,
)
    return induce(hom, g) * v
end

# Inducing action via action homomorphism

function induce(hom::InducedActionHomomorphism, g::GroupElement)
    return induce(action(hom), hom, g)
end

function induce(ac::Action, hom::InducedActionHomomorphism, g::GroupElement)
    return throw(
        """No fallback is provided for $(typeof(ac)).
        You need to implement
        `induce(::$(typeof(ac)), ::$(typeof(hom)), ::$(typeof(g)))`.""",
    )
end

function induce(
    ::ByPermutations,
    hom::H,
    g::GroupElement,
) where {H<:InducedActionHomomorphism}
    I = _int_type(hom)
    v = vec(I[hom[action(action(hom), g, f)] for f in basis(hom)])
    return PG.Perm{I}(v; check = false)
end

function induce(
    ac::ByLinearTransformation,
    hom::H,
    g::GroupElement,
) where {H<:InducedActionHomomorphism}
    I = Int[]
    J = Int[]
    V = coeff_type(ac)[]

    for (i, f) in enumerate(basis(hom))
        k = action(action(hom), g, f)
        idcs, vals = decompose(k, hom)
        append!(I, fill(i, length(idcs)))
        append!(J, idcs)
        append!(V, vals)
    end
    n = length(basis(hom))
    return sparse(I, J, V, n, n)
end

function induce(
    ac::BySignedPermutations,
    hom::H,
    g::GroupElement,
) where {H<:InducedActionHomomorphism}
    I = Int[]
    J = Int[]
    V = coeff_type(ac)[]

    for (i, f) in enumerate(basis(hom))
        k, s = action(action(hom), g, f)
        push!(I, i)
        push!(J, basis(hom)[k])
        push!(V, s)
    end
    n = length(basis(hom))
    return sparse(I, J, V, n, n)
end

# disabmiguation methods for custom implementations of InducedActionHomomorphism
function induce(
    ac::ByPermutations,
    hom::CachedExtensionHomomorphism,
    g::GroupElement,
)
    return _induce(ac, hom, g)
end

function induce(
    ac::ByLinearTransformation,
    hom::CachedExtensionHomomorphism,
    g::GroupElement,
)
    return _induce(ac, hom, g)
end

# disabmiguation
function induce(
    ac::BySignedPermutations,
    hom::CachedExtensionHomomorphism,
    g::GroupElement,
)
    return _induce(ac, hom, g)
end

function induce(
    ac::ByPermutations,
    shom::SchreierExtensionHomomorphism,
    g::GroupElement,
)
    return _induce(ac, shom, g)
end

function induce(
    ac::ByLinearTransformation,
    shom::SchreierExtensionHomomorphism,
    g::GroupElement,
)
    return _induce(ac, shom, g)
end

function induce(
    ac::BySignedPermutations,
    shom::SchreierExtensionHomomorphism,
    g::GroupElement,
)
    return _induce(ac, shom, g)
end
