abstract type Action end

abstract type InducedActionHomomorphism{A,T} end

#=
implements:
* basis(action_hom)
* action(hom::InducedActionHomomorphism) → Action
=#

Base.getindex(hom::InducedActionHomomorphism, i::Integer) = basis(hom)[i]
function Base.getindex(hom::InducedActionHomomorphism{A,T}, f::T) where {A,T}
    return basis(hom)[f]
end

AP.degree(hom::InducedActionHomomorphism) = length(basis(hom))

coeff_type(hom::InducedActionHomomorphism) = coeff_type(action(hom))
_int_type(basis::SA.AbstractBasis) = SA.key_type(basis)
_int_type(hom::InducedActionHomomorphism) = _int_type(basis(hom))

# Exceeding typemax(UInt32) here would mean e.g. that you're trying to block-diagonalize
# an SDP constraint of size 4_294_967_295 × 4_294_967_295, which is highly unlikely ;)
_int_type(::Type{<:Action}) = UInt32
_int_type(ac::Action) = _int_type(typeof(ac))

struct ExtensionHomomorphism{A<:Action,T,B<:SA.ExplicitBasis{T}} <:
       InducedActionHomomorphism{A,T}
    action::A
    basis::B
end

# interface:
SA.basis(hom::ExtensionHomomorphism) = hom.basis
action(hom::ExtensionHomomorphism) = hom.action

struct CachedExtensionHomomorphism{A,T,G,H,E<:InducedActionHomomorphism{A,T}} <:
       InducedActionHomomorphism{A,T}
    ehom::E
    cache::Dict{G,H}
    lock::Base.Threads.SpinLock
end

function CachedExtensionHomomorphism{G,H}(
    hom::InducedActionHomomorphism,
) where {G,H}
    return CachedExtensionHomomorphism(hom, Dict{G,H}(), Threads.SpinLock())
end

SA.basis(h::CachedExtensionHomomorphism) = basis(h.ehom)
action(h::CachedExtensionHomomorphism) = action(h.ehom)

function CachedExtensionHomomorphism(
    G::Group,
    action::Action,
    basis;
    precompute = false,
)
    hom = ExtensionHomomorphism(action, basis)
    S = typeof(induce(hom, one(G)))
    chom = CachedExtensionHomomorphism{eltype(G),S}(hom)
    @sync if precompute
        for g in G
            Threads.@spawn begin
                induce(action, chom, g)
            end
        end
    end
    return chom
end

function induce(ac::Action, chom::CachedExtensionHomomorphism, g::GroupElement)
    return _induce(ac, chom, g)
end

function _induce(
    action::Action,
    chom::CachedExtensionHomomorphism,
    g::GroupElement,
)
    if !haskey(chom.cache, g)
        val = induce(action, chom.ehom, g)
        lock(chom.lock) do
            return chom.cache[g] = val
        end
    end
    return chom.cache[g]
end
