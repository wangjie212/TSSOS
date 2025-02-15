struct SchreierExtensionHomomorphism{
    A,
    T,
    E<:InducedActionHomomorphism{A,T},
    Ch,
    St,
} <: InducedActionHomomorphism{A,T}
    ehom::E
    cache::Ch
    schreier_tree::St
    memoize::Bool
    lock::Base.Threads.SpinLock
end

SA.basis(h::SchreierExtensionHomomorphism) = basis(h.ehom)
action(h::SchreierExtensionHomomorphism) = action(h.ehom)

function SchreierExtensionHomomorphism(
    G::Group,
    action::Action,
    basis;
    memoize = false,
)
    hom = ExtensionHomomorphism(action, basis)
    cache = Dict(s => induce(hom, s) for s in gens(G))
    cache[one(G)] = induce(hom, one(G))
    sehom = SchreierExtensionHomomorphism(
        hom,
        cache,
        PG.SchreierTransversal(one(G), gens(G), *),
        memoize,
        Threads.SpinLock(),
    )

    return sehom
end

function memoize!(sehom::SchreierExtensionHomomorphism, val, g, lck = true)
    if lck
        lock(sehom.lock) do
            return sehom.cache[g] = val
        end
    else
        sehom.cache[g] = val
    end
    return sehom
end

function _induce(
    ac::Action,
    sehom::SchreierExtensionHomomorphism,
    g::GroupElement,
)
    if g in keys(sehom.cache)
        return sehom.cache[g]
    else
        s = sehom.schreier_tree.representatives[g]
        sI = induce(ac, sehom, s)
        hI = induce(ac, sehom, g * inv(s))
        gI = hI * sI

        if sehom.memoize
            memoize!(sehom, gI, g)
        end

        return gI
    end
end
