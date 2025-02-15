import PrettyTables
import PrettyTables.Tables

Tables.istable(::Type{<:CharacterTable}) = true
Tables.rowaccess(::Type{<:CharacterTable}) = true
Tables.rows(chtbl::CharacterTable) = irreducible_characters(chtbl)
Tables.getcolumn(row::Character, i::Int) = row[i]

# ugly hack to get this going
function Tables.columnnames(char::Character)
    return [Symbol(i) for i in axes(conjugacy_classes(char), 1)]
end
function Tables.getcolumn(row::Character, nm::Symbol)
    i = parse(Int, string(nm))
    return Tables.getcolumn(row, i)
end

__coefs_desc(T::Type) = T
__coefs_desc(::Type{<:Rational{T}}) where {T} = "rationals ($T)"
__coefs_desc(::Type{<:Cyclotomics.Cyclotomic{T}}) where {T} = "cyclotomics ($T)"

function Base.show(io::IO, ::MIME"text/plain", chtbl::CharacterTable)
    hl_odd = PrettyTables.Highlighter(;
        f = (rule, i, j) -> i % 2 == 0,
        crayon = PrettyTables.Crayon(;
            foreground = :dark_gray,
            negative = true,
        ),
    )

    fmt = (v, args...) -> sprint(show, MIME"text/plain"(), v)

    return PrettyTables.pretty_table(
        io,
        chtbl;
        title = "Character table of $(parent(chtbl)) over $(__coefs_desc(eltype(chtbl)))",
        header = ["$(first(cc))^G" for cc in conjugacy_classes(chtbl)],
        row_labels = [
            Symbol('χ', FiniteFields.subscriptify(i)) for i in axes(chtbl, 1)
        ],
        row_label_column_title = "",
        # hlines = [:header, :end],
        # vlines = [1],
        formatters = fmt,
        autowrap = true,
        linebreaks = true,
        columns_width = displaysize(io)[2] ÷ size(chtbl, 2) - 8,
        reserved_display_lines = 3[],
        # vcrop_mode = :middle,
        # equal_columns_width = true,
        # crop = :vertical,
        ellipsis_line_skip = 1,
        # alignment = [:r, :l],
        highlighters = hl_odd,
    )
end

function Base.show(io::IO, ::MIME"text/plain", χ::Character)
    println(io, "Character over ", __coefs_desc(eltype(χ)))
    return _print_char(io, χ)
end

Base.show(io::IO, χ::Character) = _print_char(io, χ)

function _print_char(io::IO, χ::Character)
    first = true
    for (i, c) in enumerate(multiplicities(χ))
        iszero(c) && continue
        first || print(io, " ")
        print(io, ((c < 0 || first) ? "" : '+'))
        !isone(c) && print(io, c, '·')
        print(io, 'χ', FiniteFields.subscriptify(i))
        first = false
    end
end
