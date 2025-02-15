Base.exponent(G::GroupsCore.Group) = exponent(conjugacy_classes(G))
Base.exponent(cclasses::AbstractArray) = lcm(order.(first.(cclasses)))
dixon_prime(G::GroupsCore.Group) = dixon_prime(order(G), exponent(G))

function dixon_prime(cclasses::AbstractArray)
    ordG = sum(length, cclasses)
    m = exponent(cclasses)
    return dixon_prime(ordG, m)
end

function dixon_prime(ordG::Integer, exponent::Integer)
    # `FiniteFields.GF{p}` needs `p` to be `Int` so no need to do the
    # computations of this function over `BigInt` if `ordG` is so we convert
    # to `Int`.
    p = 2 * convert(Int, isqrt(ordG))
    if isone(ordG)
        @assert isone(exponent)
        return 3one(p)
    end

    while true
        p = nextprime(p + 1)
        isone(p % exponent) && break # we need -1 to be in the field
    end
    return p
end

function common_esd(Ns, F::Type{FiniteFields.GF{q}}) where {q}
    @assert isprime(q)
    itr = iterate(Ns)
    @assert itr !== nothing
    esd = EigenSpaceDecomposition(F.(first(itr)))
    for N in Iterators.rest(Ns, last(itr))
        esd = refine(esd, F.(N, false))
        @debug N esd.eigspace_ptrs
        isdiag(esd) && return esd
    end
    return esd
end

function _multiplicities(
    chars::AbstractVector{<:Character{F}},
    cclasses = conjugacy_classes(first(chars)),
) where {F<:FiniteFields.GF}
    e = Int(exponent(cclasses))
    ie = inv(F(e))

    ω⁻¹ = inv(FiniteFields.rootofunity(F, e))
    ωs = Matrix{typeof(ω⁻¹)}(undef, e, e)
    Threads.@threads for k in 0:e-1
        ω⁻ᵏ = ω⁻¹^k
        for l in 0:e-1
            ωs[l+1, k+1] = ω⁻ᵏ^l
        end
    end

    multiplicities = zeros(Int, length(chars), length(cclasses), e)
    pmap = powermap(table(first(chars)))
    # pmap is a lazy object, computed on demand; to ensure its thread-safety
    # we precompute all its values here so in the loops below it is accessed read-only
    collect(pmap)

    for (i, χ) in enumerate(chars)
        Threads.@threads for j in 1:length(cclasses)
            for k in 0:e-1
                val = Int(ie * sum(χ[pmap[j, l]] * ωs[l+1, k+1] for l in 0:e-1))
                multiplicities[i, j, k+1] = val
            end
        end
    end

    return multiplicities
end
