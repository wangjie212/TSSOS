"""
    Character <: AbstractClassFunction
Struct representing (possibly virtual) character of a group.

Characters are backed by `table(χ)::CharacterTable` which actually stores the
character values. The multiplicities (decomposition into the irreducible
summands) of a given character can be obtained by calling `multiplicities(χ)`
which returns a vector of coefficients of `χ` in the basis of
`irreducible_characters(table(χ))`.

It is assumed that equal class functions on the same group will have
**identical** (ie. `===`) character tables.
"""
struct Character{T,S,ChT<:CharacterTable} <: AbstractClassFunction{T}
    table::ChT
    multips::Vector{S}
end

function Character{R}(
    chtbl::CharacterTable{Gr,T},
    multips::AbstractVector{S},
) where {R,Gr,T,S}
    return Character{R,S,typeof(chtbl)}(chtbl, multips)
end

function Character(chtbl::CharacterTable, multips::AbstractVector)
    R = Base._return_type(*, Tuple{eltype(chtbl),eltype(multips)})
    @assert R ≠ Any
    return Character{R}(chtbl, multips)
end

function Character(chtbl::CharacterTable, i::Integer)
    return Character{eltype(chtbl)}(chtbl, i)
end

function Character{T}(chtbl::CharacterTable, i::Integer) where {T}
    v = zeros(Int, nconjugacy_classes(chtbl))
    v[i] = 1
    return Character{T,Int,typeof(chtbl)}(chtbl, v)
end

function Character{T}(χ::Character) where {T}
    S = eltype(multiplicities(χ))
    ChT = typeof(table(χ))
    return Character{T,S,ChT}(table(χ), multiplicities(χ))
end

## Accessors
table(χ::Character) = χ.table
multiplicities(χ::Character) = χ.multips

## AbstractClassFunction api
Base.parent(χ::Character) = parent(table(χ))
conjugacy_classes(χ::Character) = conjugacy_classes(table(χ))
function Base.values(χ::Character{T}) where {T}
    return T[χ[i] for i in 1:nconjugacy_classes(table(χ))]
end

Base.@propagate_inbounds function Base.getindex(
    χ::Character{T},
    i::Integer,
) where {T}
    i = i < 0 ? table(χ).inv_of[abs(i)] : i
    @boundscheck 1 ≤ i ≤ nconjugacy_classes(table(χ))

    return convert(
        T,
        sum(
            c * table(χ)[idx, i] for
            (idx, c) in enumerate(multiplicities(χ)) if !iszero(c);
            init = zero(T),
        ),
    )
end

## Basic functionality

function Base.:(==)(χ::Character, ψ::Character)
    return table(χ) === table(ψ) && multiplicities(χ) == multiplicities(ψ)
end
Base.hash(χ::Character, h::UInt) = hash(table(χ), hash(multiplicities(χ), h))

function Base.deepcopy_internal(χ::Character{T}, d::IdDict) where {T}
    haskey(d, χ) && return d[χ]
    return Character{T}(table(χ), copy(multiplicities(χ)))
end

## Character arithmetic

for f in (:+, :-)
    @eval begin
        function Base.$f(χ::Character, ψ::Character)
            @assert table(χ) === table(ψ)
            return Character(table(χ), $f(multiplicities(χ), multiplicities(ψ)))
        end
    end
end

Base.:*(χ::Character, c::Number) = Character(table(χ), c .* multiplicities(χ))
Base.:*(c::Number, χ::Character) = χ * c
Base.:/(χ::Character, c::Number) = Character(table(χ), multiplicities(χ) ./ c)

Base.zero(χ::Character) = 0 * χ

function __decompose(T::Type, values::AbstractVector, tbl::CharacterTable)
    ψ = ClassFunction(values, conjugacy_classes(tbl), tbl.inv_of)
    return decompose(T, ψ, tbl)
end

function decompose(cfun::AbstractClassFunction, tbl::CharacterTable)
    return decompose(eltype(cfun), cfun, tbl)
end

function decompose(T::Type, cfun::AbstractClassFunction, tbl::CharacterTable)
    vals = Vector{T}(undef, length(conjugacy_classes(tbl)))
    return decompose!(vals, cfun, tbl)
end

function decompose!(
    vals::AbstractVector,
    cfun::AbstractClassFunction,
    tbl::CharacterTable,
)
    @assert length(vals) == length(conjugacy_classes(tbl))
    @assert conjugacy_classes(tbl) === conjugacy_classes(cfun)

    for (i, idx) in enumerate(eachindex(vals))
        χ = Character(tbl, i)
        vals[idx] = dot(χ, cfun)
    end
    return vals
end

function Base.:*(χ::Character, ψ::Character)
    @assert table(χ) == table(ψ)
    values = Characters.values(χ) .* Characters.values(ψ)
    return Character(table(χ), __decompose(Int, values, table(χ)))
end

Base.:^(χ::Character, n::Integer) = Base.power_by_squaring(χ, n)

## Group-theoretic functions:

AP.degree(χ::Character) = Int(χ(one(parent(χ))))
function AP.degree(
    χ::Character{T,CCl},
) where {T,CCl<:PG.AbstractOrbit{<:AbstractMatrix}}
    return Int(χ[1])
end

function Base.conj(χ::Character{T,S}) where {T,S}
    vals = collect(values(χ))
    all(isreal, vals) && return Character{T}(χ)
    tbl = table(χ)
    ψ = ClassFunction(vals[tbl.inv_of], conjugacy_classes(tbl), tbl.inv_of)
    multips = S[dot(ψ, χ) for χ in irreducible_characters(tbl)]
    return Character{T,eltype(multips),typeof(tbl)}(tbl, multips)
end

function isvirtual(χ::Character)
    return any(<(0), multiplicities(χ)) || any(!isinteger, multiplicities(χ))
end

function isirreducible(χ::Character)
    C = multiplicities(χ)
    k = findfirst(!iszero, C)
    k !== nothing || return false # χ is zero
    isone(C[k]) || return false # muliplicity is ≠ 1
    kn = findnext(!iszero, C, k + 1)
    kn === nothing && return true # there is only one ≠ 0 entry
    return false
end

"""
    affordable_real!(χ::Character)
Return either `χ` or `2re(χ)` depending whether `χ` is afforded by a real
representation, modifying `χ` in place.
"""
function affordable_real!(χ::Character)
    ι = frobenius_schur(χ)
    if ι <= 0 # i.e. χ is complex or quaternionic
        χ.multips .+= multiplicities(conj(χ))
    end
    return χ
end

"""
    frobenius_schur(χ::AbstractClassFunction[, pmap::PowerMap])
Return Frobenius-Schur indicator of `χ`, i.e. `Σχ(g²)` where sum is taken over
the whole group.

If χ is an irreducible `Character`, Frobenius-Schur indicator takes values in
`{1, 0, -1}` which correspond to the following situations:
 1. `χ` is real-valued and is afforded by an irreducible real representation,
 2. `χ` is a complex character which is not afforded by a real representation, and
 3. `χ` is quaternionic character, i.e. it is real valued, but is not afforded by a
 real representation.

In cases 2. and 3. `2re(χ) = χ + conj(χ)` corresponds to an irreducible character
afforded by a real representation.
"""
function frobenius_schur(χ::Character)
    @assert isirreducible(χ)

    pmap = powermap(table(χ))
    ι = sum(
        length(c) * χ[pmap[i, 2]] for (i, c) in enumerate(conjugacy_classes(χ))
    )

    ι_int = Int(ι)
    ordG = sum(length, conjugacy_classes(χ))
    d, r = divrem(ι_int, ordG)
    @assert r == 0 "Non integral Frobenius Schur Indicator: $(ι_int) = $d * $ordG + $r"
    return d
end

Base.isreal(χ::Character) = frobenius_schur(χ) > 0
