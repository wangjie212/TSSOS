function right_nullspace(M::AbstractMatrix{T}) where {T}
    A, pivots = row_echelon_form(M)
    ncolsA = size(A, 2)
    length(pivots) == ncolsA && return zeros(T, ncolsA, 0)
    W = zeros(T, ncolsA, ncolsA - length(pivots))
    for (i, el) in enumerate(setdiff(1:ncolsA, pivots))
        W[el, i] += 1
        for (j, k) in enumerate(pivots)
            if j < el
                W[k, i] -= A[j, el]
            end
        end
    end
    return W
end

function left_nullspace(M::AbstractMatrix)
    return transpose(right_nullspace(transpose(M)))
end

function left_eigen(M::AbstractMatrix{T}) where {T<:FiniteFields.GF}
    @assert ==(size(M)...)
    eigen = Dict{T,typeof(M)}()
    cumdim = 0
    for i in T # brute force; TODO: use factorization of characteristic/minimal polynomials
        cumdim >= size(M, 1) && break
        #do left eigenspaces!
        basis = first(row_echelon_form!(left_nullspace(M - i * I)))
        if size(basis, 1) > 0 # nullity
            cumdim += size(basis, 1)
            eigen[i] = basis
        end
    end
    return eigen
end

function _find_l(M::AbstractMatrix)
    # this function should be redundant when defining a better structure for echelonized subspaces
    l = Int[]
    for i in 1:size(M, 2)
        j = findfirst(isone, @view M[length(l)+1:end, i])
        if j !== nothing
            push!(l, i)
        end
    end
    return l
end

# EigenSpaceDecomposition

function eigen_decomposition!(M::Matrix{T}) where {T<:FiniteFields.GF}
    eigspace_ptrs = Vector{Int}()
    eigen = left_eigen(M)
    sizehint!(eigspace_ptrs, length(eigen) + 1)
    push!(eigspace_ptrs, 1)
    for val in sort!(collect(keys(eigen))) #to get deterministic behaviour
        basis = eigen[val]
        dim = size(basis, 1)
        cd = eigspace_ptrs[end]
        ran = cd:cd+dim-1
        M[ran, :] = basis
        push!(eigspace_ptrs, cd + dim)
    end
    @assert eigspace_ptrs[end] == size(M, 1) + 1 "Matrix does not split over $T"
    return M, eigspace_ptrs
end

mutable struct EigenSpaceDecomposition{T<:FiniteFields.GF}
    basis::Matrix{T}
    eigspace_ptrs::Vector{Int}

    function EigenSpaceDecomposition(
        basis::Matrix{T},
        eigspace_ptrs::AbstractVector{<:Integer},
    ) where {T<:FiniteFields.GF}
        @assert eigspace_ptrs[1] == 1
        @assert eigspace_ptrs[end] == size(basis, 1) + 1
        return new{T}(basis, eigspace_ptrs)
    end
end

function EigenSpaceDecomposition(M::Matrix{T}) where {T<:FiniteFields.GF}
    return EigenSpaceDecomposition(eigen_decomposition!(deepcopy(M))...)
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    esd::EigenSpaceDecomposition{T},
) where {T}
    println(io, tuple(diff(esd.eigspace_ptrs)...), "-splitting over ", T)
    return print(io, esd.basis)
end

function Base.show(io::IO, esd::EigenSpaceDecomposition{T}) where {T}
    return print(io, tuple(diff(esd.eigspace_ptrs)...), "-splitting over ", T)
end

Base.length(esd::EigenSpaceDecomposition) = length(esd.eigspace_ptrs) - 1

function Base.getindex(esd::EigenSpaceDecomposition, i::Int)
    @boundscheck 1 <= i <= length(esd)
    return esd.basis[esd.eigspace_ptrs[i]:esd.eigspace_ptrs[i+1]-1, :]
end

function Base.iterate(esd::EigenSpaceDecomposition, s = 1)
    s > length(esd) && return nothing
    first_last = esd.eigspace_ptrs[s]:esd.eigspace_ptrs[s+1]-1
    return (esd.basis[first_last, :], s + 1)
end

Base.eltype(::EigenSpaceDecomposition{T}) where {T} = Matrix{T}

function LinearAlgebra.isdiag(esd::EigenSpaceDecomposition)
    return esd.eigspace_ptrs == 1:length(esd)+1
end

function refine(esd::EigenSpaceDecomposition{T}, M::Matrix{T}) where {T}
    nbasis = Array{T}(undef, 0, size(first(esd), 2))
    nptrs = [1]
    for eigspace in esd
        if size(eigspace, 1) > 1
            esd2, ptrs =
                eigen_decomposition!(eigspace * @view M[:, _find_l(eigspace)])
            nbasis = vcat(nbasis, esd2 * eigspace)
            append!(nptrs, ptrs .+ (pop!(nptrs) - 1))
        else
            nbasis = vcat(nbasis, eigspace)
            push!(nptrs, nptrs[end] + 1)
        end
    end
    return EigenSpaceDecomposition(nbasis, nptrs)
end
