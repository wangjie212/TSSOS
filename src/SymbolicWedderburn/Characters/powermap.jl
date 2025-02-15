struct PowerMap{T<:PG.AbstractOrbit} <: AbstractMatrix{Int}
    cc::Vector{T}
    cache::Matrix{Int}

    function PowerMap(
        ccG::AbstractVector{T},
        exp = exponent(ccG),
    ) where {T<:PG.AbstractOrbit}
        cache = zeros(Int, length(ccG), exp)
        id = one(first(first(ccG)))
        idx = findfirst(cc -> id ∈ cc, ccG)
        cache[:, 1] .= idx
        if size(cache, 2) > 1
            cache[:, 2] .= 1:length(ccG)
        end
        return new{T}(ccG, cache)
    end
end

Base.axes(pm::PowerMap) = (Base.OneTo(length(pm.cc)), 0:size(pm.cache, 2)-1)
Base.size(pm::PowerMap) = size(pm.cache)

function Base.getindex(p::PowerMap, i::Integer, j::Integer)
    j = mod(j, size(p, 2))
    @boundscheck 1 ≤ i ≤ size(p, 1) || throw(Bounds)

    if iszero(p.cache[i, j+1])
        gj = first(p.cc[i])^j
        p.cache[i, j+1] = findfirst(cc -> gj ∈ cc, p.cc)
    end
    return p.cache[i, j+1]
end
