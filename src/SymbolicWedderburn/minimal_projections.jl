function _group_algebra(G::Group)
    @assert isfinite(G)
    l = if isbitstype(eltype(G))
        UInt16(0)
    else
        convert(UInt16, min(order(Int, G), typemax(UInt16) >> 2))
    end

    fb = SA.FixedBasis(vec(collect(G)), SA.DiracMStructure(*), (l, l))
    SA.complete!(fb.table)
    return StarAlgebra(G, fb)
end

Base.parent(A::StarAlgebra{<:Group}) = A.object
SA.star(g::GroupElement) = inv(g)

function SA.AlgebraElement(
    χ::AbstractClassFunction,
    RG::StarAlgebra{<:Group},
)
    G = parent(RG)
    @assert G === parent(χ)
    b = basis(RG)

    dim = degree(χ)
    ord = order(Int, G)
    c = dim // ord

    T = Base._return_type(*, Tuple{typeof(c),eltype(χ)})

    x = AlgebraElement(zeros(T, length(b)), RG)

    for (v, cc) in zip(values(χ), conjugacy_classes(χ))
        for g in cc
            x[inv(g)] = c * v
        end
    end
    return x
end

function algebra_elt_from_support(
    support,
    RG::StarAlgebra;
    val = 1 // length(support),
)
    b = basis(RG)
    I = [b[s] for s in support]
    V = fill(val, length(support))
    return AlgebraElement(sparsevec(I, V, length(b)), RG)
end
function _small_idem(RG::StarAlgebra{<:Group}, H)
    return algebra_elt_from_support(H, RG; val = 1 // length(H))
end

struct CyclicSubgroups{Gr,GEl}
    group::Gr
    seen::Dict{Int,Set{GEl}}
    min_order::Int
    max_order::Int
end

function CyclicSubgroups(G::Group; min_order = 1, max_order = order(Int, G))
    seen = Dict{Int,Set{eltype(G)}}()
    return CyclicSubgroups{typeof(G),eltype(G)}(G, seen, min_order, max_order)
end

function Base.deepcopy_internal(
    citr::CyclicSubgroups{Gr,E},
    iddict::IdDict,
) where {Gr,E}
    G = citr.group
    seen = Base.deepcopy_internal(citr.seen, iddict)
    return CyclicSubgroups{Gr,E}(G, seen, citr.min_order, citr.max_order)
end

Base.eltype(citr::CyclicSubgroups) = valtype(citr.seen)
Base.IteratorSize(::CyclicSubgroups) = Base.SizeUnknown()

function Base.iterate(citr::CyclicSubgroups)
    g, state = iterate(citr.group) # g is identity here
    @assert isone(g)
    if citr.min_order ≤ 1 ≤ citr.max_order
        citr.seen[1] = Set([g])
        return Set([g]), state
    end
    return iterate(citr, state)
end

function Base.iterate(citr::CyclicSubgroups, state)
    k = iterate(citr.group, state)
    if k === nothing
        # reset the iterator to the initial state
        empty!(citr.seen)
        return nothing
    end
    g, state = k
    ord = order(Int, g)
    if citr.min_order ≤ ord ≤ citr.max_order
        if !haskey(citr.seen, ord)
            citr.seen[ord] = Set([g])
            return Set(g^i for i in 1:ord), state
        else
            if any(g^i in citr.seen[ord] for i in 1:ord-1)
                return iterate(citr, state)
            else
                push!(citr.seen[ord], g)
                return Set(g^i for i in 1:ord), state
            end
        end
    end
    return iterate(citr, state)
end

function (χ::AbstractClassFunction)(α::AlgebraElement{<:StarAlgebra{<:Group}})
    @assert parent(χ) === parent(parent(α))
    return sum(α(g) * χ(g) for g in SA.supp(α))
end

function minimal_rank_projection(
    χ::Character,
    RG::StarAlgebra{<:Group};
    subgroups = CyclicSubgroups(parent(RG); min_order = 2),
    iters = 3,
)
    degree(χ) == 1 && return one(Rational{Int}, RG), 1

    π = one(Rational{Int}, RG)
    rank = degree(χ)

    for H in subgroups
        µ = _small_idem(RG, H)
        d = χ(µ)
        isreal(d) || continue
        rd = float(d)
        isinteger(rd) || continue
        if 0 < rd < rank
            π = µ
            rank = Int(rd)
            isone(rank) && return (π, rank)
        end
    end

    _id = [deepcopy(subgroups), deepcopy(subgroups)]

    for _ in 2:iters
        for sbgrps in Iterators.product(_id...)
            µ = *((_small_idem(RG, H) for H in sbgrps)...)
            π^2 == π || continue
            d = χ(µ)
            isreal(d) || continue
            rd = float(d)
            isinteger(rd) || continue
            if 0 < rd < rank
                π = µ
                rank = Int(rd)
                isone(rank) && return (π, rank)
            end
        end
        push!(_id, deepcopy(subgroups))
    end
    @debug "Could not find minimal projection for $χ"
    return (π, rank)
end

function minimal_projection_system(
    chars::AbstractVector{<:AbstractClassFunction},
    RG::StarAlgebra{<:Group},
)
    res = fetch.([Threads.@spawn minimal_rank_projection(χ, RG) for χ in chars])

    r1p, ranks = first.(res), last.(res) # rp1 are sparse storage

    mps = [
        isone(µ) ? AlgebraElement(χ, RG) : µ * AlgebraElement(χ, RG) for
        (µ, χ) in zip(r1p, chars)
    ] # dense storage

    return mps, ranks
end
