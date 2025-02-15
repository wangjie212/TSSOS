struct CMMatrix{T,C} <: AbstractMatrix{T}
    cc::Vector{C} # vector of conjugacy classes to fix the order
    r::Int # the index of conjugacy class
    m::Matrix{T} # cache of class coefficients

    function CMMatrix(
        cc::A,
        r::Int,
        T::Type = Int,
    ) where {C,A<:AbstractVector{C}}
        M = fill(-one(T), length(cc), length(cc))
        return new{T,C}(cc, r, M)
    end
end

Base.size(M::CMMatrix) = size(M.m)
Base.IndexStyle(::Type{<:CMMatrix}) = IndexCartesian()

function Base.getindex(M::CMMatrix, s::Integer, t::Integer)
    if isone(-M.m[s, t])
        out = one(first(first(M.cc)))

        M.m[s, :] .= 0
        r = M.r

        for g in M.cc[r]
            for h in M.cc[s]
                out = GroupsCore.mul!(out, g, h)
                for t in 1:size(M, 2)
                    if out == first(M.cc[t])
                        M.m[s, t] += 1
                        break
                    end
                end
            end
        end
    end
    return M.m[s, t]
end

conjugacy_classes(G::GroupsCore.Group) = conjugacy_classes_orbit(G)

function conjugacy_classes_orbit(G::GroupsCore.Group)
    id = one(G)
    S = gens(G)
    ordG = order(Int, G)

    cclasses = [PG.Orbit([id])]
    elts_counted = 1

    for g in G
        any(ccl -> g ∈ ccl, cclasses) && continue
        ccl_g = PG.Orbit(g, S, conj)
        elts_counted += length(ccl_g)
        push!(cclasses, ccl_g)
        elts_counted == ordG && break
    end

    @assert elts_counted == ordG
    return cclasses
end
