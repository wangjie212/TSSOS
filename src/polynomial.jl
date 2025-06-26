import Base: conj
const Mono = DP.Monomial{DP.Commutative{DP.CreationOrder}, Graded{LexOrder}}
const Poly{T} = DP.Polynomial{DP.Commutative{DP.CreationOrder}, Graded{LexOrder}, T}
const PolyLike = Union{Mono,Term,Poly}
export Mono,Poly

mutable struct poly{T}
    supp::Vector{Vector{UInt16}}
    coe::Vector{T}
end

mutable struct poly_matrix
    m::Int # size of the polynomial matrix
    polys::Vector{poly} # store the upper triangular part by colomn
end

mutable struct cpoly{T}
    supp::Vector{Tuple{Vector{UInt16},Vector{UInt16}}}
    coe::Vector{T}
end

function maxdeg(p::poly)
    return maximum(length.(p.supp))
end

function maxdeg(p::cpoly)
    return maximum([sum(length.(item)) for item in p.supp])
end

function maxcdeg(p::cpoly)
    return maximum([max(length(item[1]), length(item[2])) for item in p.supp])
end

function poly(p::Poly{T}, x) where {T<:Union{Number,AffExpr,GenericAffExpr}}
    n = length(x)
    coe = MP.coefficients(p)
    mons = MP.monomials(p)
    supp = [UInt16[] for i=1:length(mons)]
    for (i,mon) in enumerate(mons)
        ind = mon.z .> 0
        vars = mon.vars[ind]
        exp = mon.z[ind]
        for j in eachindex(vars)
            append!(supp[i], ncbfind(x, n, vars[j])*ones(UInt16, exp[j]))
        end
    end
    return poly{T}(supp,coe)
end

function cpoly(p::Poly{T}, x) where {T<:Union{Number,AffExpr,GenericAffExpr}}
    n = length(x)
    cx = MP.conj.(x)
    coe = MP.coefficients(p)
    mons = MP.monomials(p)
    supp = [tuple(UInt16[],UInt16[]) for i=1:length(mons)]
    for (i,mon) in enumerate(mons)
        ind = mon.z .> 0
        vars = mon.vars[ind]
        exp = mon.z[ind]
        for j in eachindex(vars)
            if isconj(vars[j])
                append!(supp[i][2], ncbfind(cx, n, vars[j])*ones(UInt16, exp[j]))
            else
                append!(supp[i][1], ncbfind(x, n, vars[j])*ones(UInt16, exp[j]))
            end
        end
    end
    return cpoly{T}(supp,coe)
end

function conj(a::Tuple{Vector{UInt16},Vector{UInt16}})
    return tuple(a[2], a[1])
end

function conj(p::cpoly)
    return cpoly(conj.(p.supp), conj.(p.coe))
end

function sadd(a::Vector{UInt16}, b::Vector{UInt16}; nb=0)
    c = sort!([a; b])
    if nb > 0
        i = 1
        while i < length(c)
            if c[i] <= nb
                if c[i] == c[i+1]
                    deleteat!(c, i:i+1)
                else
                    i += 1
                end
            else
                break
            end
        end
    end
    return c
end

function sadd(a::Tuple{Vector{UInt16},Vector{UInt16}}, b::Tuple{Vector{UInt16},Vector{UInt16}})
    return tuple(sort!([a[1]; b[1]]), sort!([a[2]; b[2]]))
end

function Base.isless(p::cpoly, q::cpoly)
    return p.supp < q.supp 
end

function Base.:(==)(p::cpoly, q::cpoly)
    return p.supp == q.supp 
end

function supp_multi(a::poly, b::poly, group, action; g=poly([UInt16[]], [1]))
    return unique!(vec([normalform(sadd(sadd(item1, item2), item3), group, action) for item1 in a.supp, item2 in b.supp, item3 in g.supp]))
end

function supp_multi(a::cpoly, b::cpoly, group, action; g=cpoly([tuple(UInt16[],UInt16[])], [1]))
    return unique!(vec([normalform(sadd(sadd(item1, item2), item3), group, action) for item1 in a.supp, item2 in b.supp, item3 in g.supp]))
end

function poly_multi(p::poly{T}, q::poly{T}, group, action; g=poly([UInt16[]], [1])) where {T<:Union{Number,AffExpr,GenericAffExpr}}
    supp = supp_multi(p, q, group, action, g=g)
    sort!(supp)
    coe = zeros(T, length(supp))
    for (i,item1) in enumerate(p.supp), (j,item2) in enumerate(q.supp), (k,item3) in enumerate(g.supp)
        ind = bfind(supp, length(supp), normalform(sadd(sadd(item1, item2), item3), group, action))
        @inbounds coe[ind] += p.coe[i]*q.coe[j]*g.coe[k]
    end
    return poly{T}(supp,coe)
end
