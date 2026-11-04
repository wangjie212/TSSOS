"""
    weighted_indices(degrees, order)

Exponent tuples of weighted degree at most `order`, with the constant first.
Generators must have positive integral degrees. These are exponents of invariant
generators, not the variable-index supports used elsewhere in TSSOS.
"""
function weighted_indices(degrees, order::Integer)
    !isempty(degrees) && all(d -> d isa Integer && d > 0, degrees) ||
        throw(ArgumentError("invariant degrees must be positive integers"))
    indices = Tuple[]
    function visit(prefix, j, budget)
        if j > length(degrees)
            push!(indices, Tuple(prefix))
        else
            for k in 0:div(budget, degrees[j])
                visit([prefix; k], j+1, budget-k*degrees[j])
            end
        end
    end
    order >= 0 && visit(Int[], 1, order)
    sort!(indices; by=a -> (sum(degrees .* a), reverse(a)))
    return indices
end

_index_sum(a, b) = Tuple(a .+ b)
_index_shift(a, j) = ntuple(i -> a[i] + (i == j), length(a))

"""Assemble an invariant moment matrix from a dictionary keyed by exponent tuples."""
function invariant_moment_matrix(moments::AbstractDict, basis)
    return [moments[_index_sum(a,b)] for a in basis, b in basis]
end

"""
    invariant_moments(atoms, weights, degrees, t)

Generate all invariant moments through weighted degree `2t`. Rational fixtures
remain exact, making this useful as an extraction oracle independent of an SDP.
"""
function invariant_moments(atoms::AbstractVector, weights::AbstractVector, degrees, t::Integer)
    length(atoms) == length(weights) > 0 || throw(ArgumentError("atoms and weights must have equal nonzero lengths"))
    all(v -> length(v) == length(degrees), atoms) || throw(ArgumentError("atom dimension does not match degrees"))
    return Dict(a => sum(w * prod(v[j]^a[j] for j in eachindex(degrees)) for (v,w) in zip(atoms,weights))
        for a in weighted_indices(degrees, 2t))
end

"""
    invariant_moments(data, generators, x, t; degrees=MP.maxdegree.(generators))

Read weighted invariant moments from the coefficient-constraint dual of a real
symmetry relaxation. Missing SDP rows remain missing; they are never filled with
zeros. `generators` must be algebraically independent and generate the invariant
ring (a mathematical assumption, not inferred by this routine).
"""
function invariant_moments(data::polybasis_data, generators, x, t::Integer; degrees=MP.maxdegree.(generators))
    data.obj isa poly && data.moments !== nothing || throw(ArgumentError("real symmetry SDP dual moments are required"))
    length(generators) == length(degrees) && all(MP.maxdegree.(generators) .== degrees) ||
        throw(ArgumentError("generator degrees are inconsistent"))
    for generator in generators
        p = poly(MP.polynomial(generator),x)
        original = Dict(Tuple(mon)=>c for (mon,c) in zip(p.supp,p.coe))
        support_action = data.action isa SignSymmetryAction ? data.action.permutation : data.action
        permutation_group = data.action isa SignSymmetryAction ? data.action.group : data.group
        for g in permutation_group
            image = Dict(Tuple(action_mon(support_action,g,mon))=>c for (mon,c) in zip(p.supp,p.coe))
            image == original || throw(ArgumentError("a generator is not invariant under the permutation action"))
        end
        if data.action isa SignSymmetryAction
            for mon in p.supp, s in eachcol(data.action.signs)
                iseven(sum(s[i] for i in mon;init=0)) || throw(ArgumentError("a generator is not invariant under sign symmetry"))
            end
        end
    end
    moments = Dict{Tuple,data.coe_type}()
    for a in weighted_indices(degrees, 2t)
        p = poly(MP.polynomial(prod(generators[j]^a[j] for j in eachindex(a))), x)
        value = zero(data.coe_type)
        available = true
        for (mon, c) in zip(p.supp, p.coe)
            nf = normalform(mon, data.group, data.action)
            nf === nothing && continue
            if !haskey(data.moments, Tuple(nf))
                available = false
                break
            end
            value += c * data.moments[Tuple(nf)]
        end
        available && (moments[a] = value)
    end
    return moments
end

# Generic Jacobi diagonalization keeps BigFloat computations at their input
# precision without depending on a LAPACK Float64 conversion.
function _invariant_eigen(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    T <: LinearAlgebra.BlasFloat && return eigen(Symmetric(A))
    B = Matrix(Symmetric(A))
    n = size(B,1)
    Q = Matrix{T}(I,n,n)
    threshold = eps(T) * max(one(T), norm(B)) * n
    for sweep in 1:max(30, 10n)
        off = maximum((abs(B[i,j]) for i in 1:n for j in i+1:n); init=zero(T))
        if off <= threshold
            p = sortperm(diag(B))
            return (; values=diag(B)[p], vectors=Q[:,p])
        end
        for p in 1:n-1, q in p+1:n
            abs(B[p,q]) <= threshold && continue
            theta = (B[q,q]-B[p,p]) / (2B[p,q])
            u = (theta >= 0 ? one(T) : -one(T)) / (abs(theta)+sqrt(1+theta^2))
            c = inv(sqrt(1+u^2)); s = u*c
            app, aqq, apq = B[p,p], B[q,q], B[p,q]
            for k in 1:n
                if k != p && k != q
                    bp, bq = B[k,p], B[k,q]
                    B[k,p] = B[p,k] = c*bp-s*bq
                    B[k,q] = B[q,k] = s*bp+c*bq
                end
                qp, qq = Q[k,p], Q[k,q]
                Q[k,p], Q[k,q] = c*qp-s*qq, s*qp+c*qq
            end
            B[p,p], B[q,q] = app-u*apq, aqq+u*apq
            B[p,q] = B[q,p] = zero(T)
        end
    end
    throw(ErrorException("generic symmetric eigensolver did not converge"))
end

function _exact_pivots(A)
    U = Rational{BigInt}.(A)
    row = 1
    pivots = Int[]
    for col in axes(U,2)
        k = findfirst(i -> !iszero(U[i,col]), row:size(U,1))
        k === nothing && continue
        pivot = row+k-1
        U[row,:], U[pivot,:] = copy(U[pivot,:]), copy(U[row,:])
        U[row,:] ./= U[row,col]
        for i in row+1:size(U,1)
            U[i,:] -= U[i,col]*U[row,:]
        end
        push!(pivots,col)
        row += 1
        row > size(U,1) && break
    end
    return pivots
end

function _invariant_basis(M, r)
    # Pivoted Cholesky with the constant pinned, then largest Schur complement.
    selected = [1]
    while length(selected) < r
        remaining = setdiff(collect(axes(M,1)), selected)
        F = cholesky(Symmetric(M[selected,selected]))
        residuals = [M[j,j] - dot(M[selected,j], F \ M[selected,j]) for j in remaining]
        push!(selected, remaining[argmax(residuals)])
    end
    return selected
end

"""
    extract_invariant_atoms(moments, degrees, t; arithmetic=:numeric, ...)
    extract_invariant_atoms(data, generators, x, t; ...)

Extract a positive atomic measure using only the invariant functional. Require
PSD of `M_D` and equal ranks of `M_C` and `M_D`, where `C=D_(t-max(degrees))`.
The constant is pinned in the quotient basis. Cholesky solves and a seeded
symmetric joint eigensolver recover atoms and spectral weights.

`moments` maps nonnegative exponent tuples to real values. All `D+D` indices are
required. Extra supplied moments are checked separately; `full_representation`
also requires coverage of *every* index through weighted degree `2t`.
`arithmetic=:exact` computes ranks, H, K and multiplication matrices over
Rational{BigInt}; atom coordinates are still recovered numerically. Pass
BigFloat moments for higher precision spectral recovery. Algebraic independence
of the invariant generators is a standing assumption.

Returns a named tuple with `status`, `atoms`, `weights`, `rank`, `basis`, `H`,
`multiplication`, `symmetric_multiplication`, and `diagnostics`. Mathematical or
numerical gate failures return an empty measure and a diagnostic status; invalid
inputs throw ArgumentError. No feasibility or optimality is inferred here.
"""
function extract_invariant_atoms(moments::AbstractDict, degrees, t::Integer;
    arithmetic=:numeric, basis=nothing, seed=0, max_attempts=8,
    rtol=nothing, atol=0, validation_tol=nothing, gap_tol=nothing, probability=true)
    started = time_ns()
    D = weighted_indices(degrees,t)
    t >= maximum(degrees) || throw(ArgumentError("order must be at least the largest invariant degree"))
    arithmetic in (:numeric,:exact) || throw(ArgumentError("arithmetic must be :numeric or :exact"))
    max_attempts >= 1 || throw(ArgumentError("max_attempts must be positive"))
    m = length(degrees)
    F = Set(keys(moments))
    all(a -> a isa Tuple && length(a) == m && all(i -> i isa Integer && i >= 0, a) && sum(degrees .* a) <= 2t, F) ||
        throw(ArgumentError("moment keys must be exponent tuples of weighted degree at most 2t"))
    all(v -> v isa Real && isfinite(v), values(moments)) || throw(ArgumentError("moments must be finite real numbers"))
    C = weighted_indices(degrees,t-maximum(degrees))
    required = Set(_index_sum(a,b) for a in D for b in D)
    missing = setdiff(required,F)
    isempty(missing) || throw(ArgumentError("missing required D+D moments: $(sort!(collect(missing)))"))
    FT = float(foldl(promote_type,(typeof(v) for v in values(moments))))
    T = FT <: AbstractFloat ? FT : Float64
    tol = rtol === nothing ? sqrt(eps(T)) : T(rtol)
    vtol = validation_tol === nothing ? 100tol : T(validation_tol)
    gtol = gap_tol === nothing ? 10sqrt(eps(T)) : T(gap_tol)
    tol >= 0 && atol >= 0 && vtol > 0 && gtol >= 0 || throw(ArgumentError("invalid tolerance"))
    mass = T(moments[ntuple(_ -> 0,m)])
    mass > 0 || throw(ArgumentError("L(1) must be positive"))
    probability && abs(mass-1) > vtol && throw(ArgumentError("probability measures require L(1)=1"))
    exact = arithmetic == :exact
    exact && !all(v -> v isa Union{Integer,Rational}, values(moments)) &&
        throw(ArgumentError("exact arithmetic requires integer or rational moments"))
    MC, MD = invariant_moment_matrix(moments,C), invariant_moment_matrix(moments,D)
    assembled = time_ns()
    ec, ed = _invariant_eigen(T.(MC)).values, _invariant_eigen(T.(MD)).values
    threshold_c = T(atol)+tol*maximum(abs,ec)
    threshold_d = T(atol)+tol*maximum(abs,ed)
    rc = exact ? length(_exact_pivots(MC)) : count(>(threshold_c),ec)
    rd = exact ? length(_exact_pivots(MD)) : count(>(threshold_d),ed)
    rank_gap(e, r) = r == 0 ? zero(T) : sort(abs.(e);rev=true)[r] / max(r < length(e) ? sort(abs.(e);rev=true)[r+1] : zero(T), eps(T)*maximum(abs,e))
    diagnostics = Dict{Symbol,Any}(:degrees=>collect(degrees), :order=>t, :C=>C, :D=>D,
        :required_indices=>sort!(collect(required)), :additional_indices=>sort!(collect(setdiff(F,required))),
        :available_indices=>sort!(collect(F)), :rank_C=>rc, :rank_D=>rd,
        :rank_rtol=>tol, :rank_atol=>atol, :rank_thresholds=>(threshold_c,threshold_d),
        :spectra=>(ec,ed), :rank_gaps=>(rank_gap(ec,rc),rank_gap(ed,rd)),
        :minimum_eigenvalue=>minimum(ed), :required_moments_consistent=>false,
        :full_moments_consistent=>false, :full_degree_coverage=>Set(weighted_indices(degrees,2t)) ⊆ F,
        :full_representation=>false, :arithmetic=>arithmetic)
    discarded(e, threshold) = filter(v -> 10length(e)*eps(T)*maximum(abs,e) < v <= threshold,e)
    diagnostics[:subthreshold_positive_eigenvalues] = (discarded(ec,threshold_c),discarded(ed,threshold_d))
    diagnostics[:rank_resolution_warning] = any(!isempty,diagnostics[:subthreshold_positive_eigenvalues])
    empty_result(status; basis=Tuple[], H=nothing, multiplication=Matrix{T}[], symmetric_multiplication=Matrix{T}[]) =
        (; status, atoms=Vector{T}[], weights=T[], rank=rc, basis, H, multiplication, symmetric_multiplication, diagnostics)
    minimum(ed) >= -threshold_d || return empty_result(:not_psd)
    rc == rd && rc > 0 || return empty_result(:not_flat)
    if exact
        # Exact PSD oracle: nonsingular principal rank block plus its Schur complement.
        piv = _exact_pivots(MD)
        P = Rational{BigInt}.(MD[piv,piv])
        for k in 1:rd
            det(P[1:k,1:k]) > 0 || return empty_result(:not_psd)
        end
        Rational{BigInt}.(MD) == Rational{BigInt}.(MD[:,piv]) * (P \ Rational{BigInt}.(MD[piv,:])) || return empty_result(:not_psd)
    end
    ranked = time_ns()
    selected = if basis === nothing
        try
            exact ? _exact_pivots(MC) : _invariant_basis(T.(MC),rc)
        catch e
            e isa PosDefException || rethrow()
            return empty_result(:ill_conditioned_basis)
        end
    else
        length(basis) == rc && first(basis) == first(C) && length(unique(basis)) == rc && all(b -> b in C,basis) ||
            throw(ArgumentError("basis must contain r distinct indices from C with the constant first"))
        [findfirst(==(b),C) for b in basis]
    end
    selected_at = time_ns()
    B = C[selected]
    H = exact ? Rational{BigInt}.(MC[selected,selected]) : T.(MC[selected,selected])
    K = [exact ? Rational{BigInt}.([moments[_index_shift(_index_sum(a,b),j)] for a in B,b in B]) :
        T.([moments[_index_shift(_index_sum(a,b),j)] for a in B,b in B]) for j in 1:m]
    R = try
        cholesky(Symmetric(T.(H))).U
    catch e
        e isa PosDefException || rethrow()
        return empty_result(:ill_conditioned_basis;basis=B,H)
    end
    multiplication = [H \ Kj for Kj in K]
    if exact
        diagnostics[:exact_quotient_consistent] = all(A*B == B*A for A in multiplication,B in multiplication) &&
            all(H*A == A'*H for A in multiplication)
    end
    S = [Matrix(Symmetric((R' \ T.(Kj)) / R)) for Kj in K]
    comm = maximum((norm(A*B-B*A) for A in S,B in S);init=zero(T))
    comm_scaled = maximum((norm(A*B-B*A)/(1+norm(A)*norm(B)) for A in S,B in S);init=zero(T))
    diagnostics[:commutator_raw], diagnostics[:commutator_scaled] = comm, comm_scaled
    rng = MersenneTwister(seed)
    Q = nothing
    for attempt in 1:max_attempts
        coeff = T.(randn(rng,m)); coeff ./= norm(coeff)
        E = _invariant_eigen(sum(coeff[j]*S[j] for j in 1:m))
        gap = rc == 1 ? T(Inf) : minimum(diff(E.values))
        diagnostics[:spectral_attempts] = attempt
        diagnostics[:spectral_gap] = gap
        diagnostics[:spectral_coefficients] = coeff
        gap > gtol*max(one(T),maximum(abs,E.values)) || continue
        Q = E.vectors
        break
    end
    Q === nothing && return empty_result(:spectral_collision;basis=B,H,multiplication,symmetric_multiplication=S)
    joint_residual = maximum(norm(Sj*Q-Q*Diagonal(diag(Q'*Sj*Q)))/(1+norm(Sj)) for Sj in S)
    atoms = [[dot(q,Sj*q) for Sj in S] for q in eachcol(Q)]
    weights = [dot(q,R[:,1])^2 for q in eachcol(Q)]
    order = sortperm(atoms;by=Tuple)
    atoms, weights = atoms[order], weights[order]
    extracted = time_ns()
    residual(a) = abs(sum(w*prod(v[j]^a[j] for j in 1:m) for (v,w) in zip(atoms,weights))-T(moments[a]))
    residuals = Dict(a => residual(a) for a in F)
    scaled = Dict(a => residuals[a]/(1+abs(T(moments[a]))) for a in F)
    required_error = maximum(scaled[a] for a in required)
    additional = setdiff(F,required)
    additional_error = maximum((scaled[a] for a in additional);init=zero(T))
    valid = all(>(zero(T)),weights) && abs(sum(weights)-mass) <= vtol*(1+mass) &&
        required_error <= vtol && comm_scaled <= vtol && joint_residual <= vtol
    if exact
        valid &= diagnostics[:exact_quotient_consistent]
        rational_atoms = [[rationalize(BigInt,v;tol=vtol) for v in atom] for atom in atoms]
        rational_weights = [rationalize(BigInt,w;tol=vtol) for w in weights]
        exact_identity(a) = sum(w*prod(v[j]^a[j] for j in 1:m) for (v,w) in zip(rational_atoms,rational_weights)) == moments[a]
        rational_measure = all(>(0),rational_weights) && all(exact_identity,required)
        diagnostics[:exact_atoms] = rational_measure ? rational_atoms : nothing
        diagnostics[:exact_weights] = rational_measure ? rational_weights : nothing
        diagnostics[:exact_full_moments_consistent] = rational_measure && all(exact_identity,F)
    end
    merge!(diagnostics, Dict(:joint_eigenvector_residual=>joint_residual,
        :moment_residuals_raw=>residuals, :moment_residuals_scaled=>scaled,
        :required_residual=>required_error, :additional_residual=>additional_error,
        :weight_sum_residual=>abs(sum(weights)-mass), :minimum_weight=>minimum(weights),
        :required_moments_consistent=>valid, :full_moments_consistent=>valid && additional_error <= vtol,
        :full_representation=>valid && additional_error <= vtol && diagnostics[:full_degree_coverage]))
    diagnostics[:timings_seconds] = (; moment_assembly=(assembled-started)/1e9,
        rank_checks=(ranked-assembled)/1e9, basis_selection=(selected_at-ranked)/1e9,
        extraction=(extracted-selected_at)/1e9, validation=(time_ns()-extracted)/1e9)
    status = !valid ? :invalid_representation : additional_error > vtol ? :inconsistent_full_moments : :success
    return (; status, atoms, weights, rank=rc, basis=B, H, multiplication,
        symmetric_multiplication=S, diagnostics)
end

function extract_invariant_atoms(data::polybasis_data, generators, x, t::Integer;
    degrees=MP.maxdegree.(generators), kwargs...)
    return extract_invariant_atoms(invariant_moments(data,generators,x,t;degrees),degrees,t;kwargs...)
end

function _elementary(v)
    e = [one(eltype(v)); zeros(eltype(v),length(v))]
    for (i,x) in enumerate(v), j in i:-1:1
        e[j+1] += x*e[j]
    end
    return e[2:end]
end

function _real_fiber(v, symmetry, tol)
    symmetry == :signs && return minimum(v) >= -tol ? sqrt.(max.(v,0)) : nothing
    n = length(v)
    companion = zeros(Float64,n,n)
    for i in 2:n
        companion[i,i-1] = 1
    end
    companion[:,end] = [(-1)^(n-i)*v[n-i+1] for i in 1:n]
    roots = eigvals(companion)
    maximum(abs,imag.(roots)) <= sqrt(tol)*(1+maximum(abs,roots)) || return nothing
    x = sort!(real.(roots))
    # Multiple roots are ill-conditioned in generator coordinates. Resolve a
    # cluster at the declared fiber tolerance, then recheck all invariants below.
    root_tol = sqrt(tol)*(1+maximum(abs,roots))
    start = 1
    for stop in 2:n+1
        if stop == n+1 || x[stop]-x[stop-1] > root_tol
            x[start:stop-1] .= sum(x[start:stop-1])/(stop-start)
            start = stop
        end
    end
    if symmetry == :signed_permutations
        minimum(x) >= -tol || return nothing
        x = sqrt.(max.(x,0))
    end
    return x
end

function _orbit_size(x, symmetry, tol)
    symmetry == :signs && return big(2)^count(v -> abs(v)>tol,x)
    values = sort(symmetry == :signed_permutations ? abs.(x) : x)
    counts = Int[]
    for i in eachindex(values)
        if i == 1 || !isapprox(values[i],values[i-1];atol=tol,rtol=tol)
            push!(counts,1)
        else
            counts[end] += 1
        end
    end
    permutations = factorial(big(length(x))) ÷ prod(factorial(big(c)) for c in counts)
    return symmetry == :signed_permutations ? permutations*big(2)^count(v -> abs(v)>tol,x) : permutations
end

function _builtin_orbit(x, symmetry, tol)
    orbit = [x]
    cursor = 1
    while cursor <= length(orbit)
        y = orbit[cursor]
        images = Vector{eltype(x)}[]
        if symmetry != :signs
            for i in 1:length(x)-1
                z = copy(y); z[i],z[i+1] = z[i+1],z[i]; push!(images,z)
            end
        end
        if symmetry != :permutations
            for i in eachindex(x)
                z = copy(y); z[i] = -z[i]; push!(images,z)
            end
        end
        for z in images
            any(w -> isapprox(w,z;atol=tol,rtol=tol),orbit) || push!(orbit,z)
        end
        cursor += 1
    end
    return orbit
end

"""
    lift_invariant_atoms(result; symmetry=:signs, enumerate=false, ...)

Lift independent-sign, full-permutation (elementary generators), or full signed-
permutation (elementary generators in squared variables) atoms. Use `symmetry=
:custom`, `lift(v)`, `invariants(x)`, `group`, and `action(g,x)` for other finite
actions. Built-in orbit sizes are computed without enumerating points by default.
With `enumerate=true`, each orbit point has mass `weight/orbit_size`.

Optional `inequalities`, `equalities`, and `objective` are functions on original
points and are checked on every orbit point when supplied. Certification requires
`lower_bound_certified=true`, a matching `lower_bound`, a full represented moment
sequence and successful explicit feasibility checks. The caller is responsible
for the validity of that lower bound. Without these, labels remain candidates.
"""
function lift_invariant_atoms(result; symmetry=:signs, lift=nothing, invariants=nothing,
    group=nothing, action=nothing, enumerate=false, inequalities=nothing, equalities=nothing,
    objective=nothing, lower_bound=nothing, lower_bound_certified=false, tol=1e-6)
    tol > 0 || throw(ArgumentError("lifting tolerance must be positive"))
    symmetry in (:signs,:permutations,:signed_permutations,:custom) || throw(ArgumentError("unknown symmetry"))
    result.diagnostics[:required_moments_consistent] || throw(ArgumentError("a valid extracted measure is required"))
    if symmetry == :custom
        all(!isnothing,(lift,invariants,group,action)) || throw(ArgumentError("custom lifting requires lift, invariants, group and action"))
    else
        invariants = symmetry == :signs ? x -> x.^2 : symmetry == :permutations ? _elementary : x -> _elementary(x.^2)
        lift = v -> _real_fiber(v,symmetry,tol)
    end
    orbits = NamedTuple[]
    for (v,weight) in zip(result.atoms,result.weights)
        x = lift(v)
        if x === nothing || !all(isfinite,x) || !all(isreal,x) || !isapprox(invariants(x),v;atol=tol,rtol=tol)
            push!(orbits,(; status=:not_realizable, representative=x, weight, size=nothing,
                stabilizer=nothing, points=nothing, point_weight=nothing, feasible=nothing, label=:candidate))
            continue
        end
        points = nothing
        stabilizer = nothing
        if symmetry == :custom
            images = [action(g,x) for g in group]
            points = Vector{eltype(x)}[]
            for y in images
                all(isfinite,y) && all(isreal,y) || throw(ArgumentError("group action must preserve real finite points"))
                isapprox(invariants(y),v;atol=tol,rtol=tol) || throw(ArgumentError("group action leaves the supplied invariant fiber"))
                any(z -> isapprox(y,z;atol=tol,rtol=tol),points) || push!(points,y)
            end
            size = length(points)
            stabilizer = count(y -> isapprox(y,x;atol=tol,rtol=tol),images)
            size*stabilizer == length(images) || throw(ArgumentError("orbit-stabilizer check failed"))
        else
            size = _orbit_size(x,symmetry,tol)
            group_order = symmetry == :signs ? big(2)^length(x) : factorial(big(length(x))) * (symmetry == :signed_permutations ? big(2)^length(x) : 1)
            stabilizer = group_order ÷ size
            enumerate && (points = _builtin_orbit(x,symmetry,tol))
        end
        points !== nothing && length(points) != size && throw(ArgumentError("orbit enumeration disagrees with orbit size"))
        # Arbitrary user constraints need not be invariant. Materialize the orbit
        # only when these explicit downstream checks were requested.
        checking = inequalities !== nothing || equalities !== nothing || objective !== nothing
        checked_points = checking && points === nothing ? _builtin_orbit(x,symmetry,tol) : points
        feasible = inequalities === nothing || equalities === nothing ? nothing :
            all(y -> all(g -> g(y) >= -tol,inequalities) && all(h -> abs(h(y)) <= tol,equalities),checked_points)
        optimum_matches = objective !== nothing && lower_bound !== nothing &&
            all(y -> abs(objective(y)-lower_bound) <= tol*(1+abs(lower_bound)),checked_points)
        certified = feasible === true && optimum_matches && lower_bound_certified && result.diagnostics[:full_representation]
        push!(orbits,(; status=:lifted, representative=x, weight, size, stabilizer,
            points=enumerate ? points : nothing, point_weight=weight/size, feasible,
            label=certified ? :certified_global_minimizers : :candidate))
    end
    return orbits
end
