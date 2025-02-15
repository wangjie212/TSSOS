struct CharacterTable{Gr,T,O} <: AbstractMatrix{T}
    group::Gr
    conjugacy_classes::Vector{O}
    inv_of::Vector{Int}
    pmap::PowerMap
    values::Matrix{T}
    # rows indexed by irreducible characters,
    # columns indexed by conjugacy classes
end

## Array Interface:
Base.size(chtbl::CharacterTable) = size(chtbl.values)
function Base.getindex(chtbl::CharacterTable, i::Integer, j::Integer)
    return chtbl.values[i, j]
end

## Accessors
powermap(chtbl::CharacterTable) = chtbl.pmap
Base.parent(chtbl::CharacterTable) = chtbl.group

conjugacy_classes(chtbl::CharacterTable) = chtbl.conjugacy_classes
nconjugacy_classes(chtbl::CharacterTable) = size(chtbl.values, 2)
nirreps(chtbl::CharacterTable) = size(chtbl.values, 1)

## irreps
function irreducible_characters(chtbl::CharacterTable)
    return irreducible_characters(eltype(chtbl), chtbl)
end

function irreducible_characters(T::Type, chtbl::CharacterTable)
    return [Character{T}(chtbl, i) for i in axes(chtbl, 1)]
end

function irreducible_characters(G::Group, cclasses = conjugacy_classes(G))
    return irreducible_characters(Rational{Int}, G, cclasses)
end

function irreducible_characters(
    R::Type{<:Rational},
    G::Group,
    cclasses = conjugacy_classes(G),
)
    return irreducible_characters(CharacterTable(R, G, cclasses))
end

trivial_character(chtbl::CharacterTable) = Character(chtbl, 1)

## construcing tables

function CharacterTable(G::Group, cclasses = conjugacy_classes(G))
    return CharacterTable(Rational{Int}, G, cclasses)
end

function CharacterTable(
    Fp::Type{<:FiniteFields.GF},
    G::Group,
    cclasses = conjugacy_classes(G),
)
    # make sure that the first class contains the indentity
    k = findfirst(cl -> one(G) in cl, cclasses)
    cclasses[k], cclasses[1] = cclasses[1], cclasses[k]

    Ns = [CMMatrix(cclasses, i) for i in 1:length(cclasses)]
    esd = common_esd(Ns, Fp)
    @assert isdiag(esd)

    tbl = CharacterTable(
        G,
        cclasses,
        _inv_of(cclasses),
        PowerMap(cclasses),
        esd.basis,
    )

    tbl = normalize!(tbl)

    let vals = tbl.values
        # make order of characters deterministic
        vals .= sortslices(vals; dims = 1)
        # and that the trivial character is first
        k = findfirst(i -> all(isone, @views vals[i, :]), axes(vals, 1))
        # doesn't work on julia-1.6
        # k = findfirst(r -> all(isone, r), eachrow(vals))
        if k ≠ 1
            _swap_rows!(vals, 1, k)
        end
    end

    return tbl
end

function CharacterTable(
    R::Type{<:Rational},
    G::Group,
    cclasses = conjugacy_classes(G),
)
    Fp = FiniteFields.GF{dixon_prime(cclasses)}
    tblFp = CharacterTable(Fp, G, cclasses)
    return complex_character_table(R, tblFp)
end

function complex_character_table(
    R::Type{<:Rational},
    tblFp::CharacterTable{<:Group,<:FiniteFields.GF},
)
    charsFp = irreducible_characters(tblFp)
    mult_c = _multiplicities(charsFp)

    e = size(mult_c, 3) # the exponent

    C = Cyclotomics.Cyclotomic{R,Cyclotomics.SparseVector{R,Int}}
    values = Matrix{C}(undef, size(tblFp))

    Es = [E(e, k) for k in 0:e-1]
    Threads.@threads for j in 1:size(tblFp, 2) # conjugacy_classes
        for i in 1:size(tblFp, 1) # characters
            # reduced_embedding may prevent overflow sometimes
            values[i, j] = Cyclotomics.reduced_embedding(
                sum(mult_c[i, j, k+1] * Es[k+1] for k in 0:e-1),
            )
        end
    end

    return CharacterTable(
        parent(tblFp),
        conjugacy_classes(tblFp),
        tblFp.inv_of,
        powermap(tblFp),
        values,
    )
end

function _inv_of(cc::AbstractVector{<:PG.AbstractOrbit})
    inv_of = zeros(Int, size(cc))
    for (i, c) in enumerate(cc)
        g = inv(first(c))
        inv_of[i] = something(findfirst(k -> g in k, cc), 0)
    end
    any(iszero, inv_of) && throw(
        ArgumentError(
            "Could not find the conjugacy class for inverse of $(first(cc[findfirst(iszero, inv_of)])).",
        ),
    )
    return inv_of
end

function normalize!(chtbl::CharacterTable{<:Group,<:FiniteFields.GF})
    id = one(parent(chtbl))
    Threads.@threads for i in axes(chtbl, 1)
        χ = Character(chtbl, i)
        k = χ(id)
        if !isone(k)
            chtbl.values[i, :] .*= inv(k)
        end
        # ⟨χ, χ⟩ = 1/d²

        deg = sqrt(inv(dot(χ, χ)))
        @debug "normalizing with" dot(χ, χ) χ(id) χ

        # normalizing χ
        chtbl.values[i, :] .*= deg
    end
    return chtbl
end
