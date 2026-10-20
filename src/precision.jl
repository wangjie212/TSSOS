# Delay default solver selection so that a supplied model needs no default solver.
function _sdp_model(model, optimizer, dualize)
    if model !== nothing
        model isa GenericModel || throw(ArgumentError("model must be a JuMP.GenericModel"))
        optimizer === nothing || throw(ArgumentError("Pass either model or optimizer, not both."))
        dualize && throw(ArgumentError("With model, configure dualization on the supplied model instead of using dualize=true."))
        # Optimizer attributes can mark a fresh model dirty, so inspect its contents.
        if num_variables(model) != 0 || num_constraints(model; count_variable_in_set_constraints=true) != 0 ||
            objective_sense(model) != MOI.FEASIBILITY_SENSE
            throw(ArgumentError("Pass an empty model for each solver call."))
        end
        return model
    end
    optimizer = optimizer === nothing ? default_optimizer() : optimizer
    return Model(dualize ? dual_optimizer(optimizer) : optimizer)
end

# LAPACK/Ipopt extraction is a separate heuristic; never overwrite SDP data.
extraction_matrix(m) = eltype(m) <: Real ? Symmetric(Float64.(m)) : Hermitian(ComplexF64.(m))
extraction_moments(moment) = [extraction_matrix(m) for m in moment]
extraction_cmoments(moment) = [extraction_moments(clique) for clique in moment]

# SymbolicWedderburn 0.4 uses SuiteSparse QR for sparse numerical projections,
# which only supports Float64/ComplexF64. Build the same exact projections and
# orthonormalize their row spaces with Julia's generic dense pivoted QR instead.
# This avoids both a Float64 round-trip and methods added to dependency types.
function precision_wedderburn(::Type{T}, group, action, full, half; semisimple=false) where {T}
    if T <: LinearAlgebra.BlasFloat
        return WedderburnDecomposition(T, group, action, full, half; semisimple)
    end
    SW = SymbolicWedderburn
    full_basis = SW.StarAlgebras.Basis{UInt32}(full)
    half_basis = SW.StarAlgebras.Basis{UInt32}(half)
    table = SW.CharacterTable(Rational{BigInt}, group)
    invariants = SW.invariant_vectors(table, action, full_basis)
    hom = SW.CachedExtensionHomomorphism(group, action, half_basis; precompute=true)
    SW.check_group_action(group, hom; full_check=false)
    chars, mults = SW._constituents_decomposition(SW.action_character(hom, table), table)
    if T <: Real
        chars, mults = SW.affordable_real(chars, mults)
    end
    if semisimple || all(isone ∘ SW.degree, chars)
        projections = chars
        ranks = mults .* SW.degree.(chars)
    else
        projections, factors = SW.minimal_projection_system(chars, SW._group_algebra(SW.parent(table)))
        ranks = mults .* factors
    end
    summands = map(projections, ranks) do projection, rank
        exact = SW.matrix_representation(hom, projection)
        numeric = map(v -> projection_coefficient(T, v), Matrix(exact))
        F = qr(Matrix(numeric'), ColumnNorm())
        # rank comes from the character decomposition, not a Float64 threshold.
        Q = Matrix(F.Q)[:, 1:rank]
        rows = SW.sparse(Matrix(Q'))
        SW.droptol!(rows, eps(real(T)) * max(size(rows)...))
        (; basis=rows)
    end
    return (; Uπs=summands, invariants)
end

function projection_coefficient(::Type{T}, value) where {T<:Real}
    isreal(value) || throw(ArgumentError("A real symmetry basis requires real projection coefficients."))
    return real(Complex{T}(value))
end
projection_coefficient(::Type{T}, value) where {T<:Complex} = T(value)
