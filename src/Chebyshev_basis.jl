mutable struct poly_basis
    cost # cost function
    ineq_cons # inequality constraints
    eq_cons # equality constraints
    coe_type # type of coefficients
    group
    action
    basis # polynomial bases
    ebasis # polynomial bases for equality constraints
    ksupp # extended support at the k-th step
    blocksize # size of blocks
    blocks # block structrue
    eblocks # block structrue for equality constraints
    GramMat # Gram matrices
    multiplier # multipliers for equality constraints
    SDP_status
end

"""
    info = add_psatz_cheby!(model, nonneg, vars, ineq_cons, eq_cons, order; TS="block", SO=1, QUIET=false)

Add a Putinar's style SOS representation of the polynomial `nonneg` in the Chebyshev basis to the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `nonneg`: a nonnegative polynomial constrained to be a Putinar's style SOS on a semialgebraic set
- `vars`: the set of POP variables
- `ineq_cons`: inequality constraints
- `eq_cons`: equality constraints
- `order`: relaxation order
- `TS`: type of term sparsity (`"block"`, `"MD"`, `"MF"`, `false`)
- `SO`: sparse order
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `info`: auxiliary data
"""
function add_psatz_cheby!(model, nonneg::Poly{T}, vars, ineq_cons, eq_cons, order; TS="block", SO=1, merge=false, md=3, QUIET=false) where {T<:Union{Number,AffExpr}}
    m = length(ineq_cons)
    l = length(eq_cons)
    basis = Vector{ChebyshevBasisFirstKind{Poly{Float64}}}(undef, m+l+1)
    basis[1] = basis_covering_monomials(ChebyshevBasis, MP.monomials(vars, 0:order))
    basis[2:m+1] = [basis_covering_monomials(ChebyshevBasis, MP.monomials(vars, 0:order-ceil(Int, MP.maxdegree(g)/2))) for g in ineq_cons]
    basis[m+2:m+1+l] = [basis_covering_monomials(ChebyshevBasis, MP.monomials(vars, 0:2*order-MP.maxdegree(h))) for h in eq_cons]
    tsupp = basis_covering_monomials(ChebyshevBasis, unique([MP.monomials(nonneg); MP.monomials.(ineq_cons)...; MP.monomials.(eq_cons)...]))
    tsupp = [item for item in tsupp]
    sort!(tsupp)
    blocks,cl,blocksize,eblocks = get_blocks(m, l, tsupp, ineq_cons, eq_cons, basis, TS=TS, SO=SO, merge=merge, md=md, QUIET=QUIET)
    pol = nonneg
    pos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, 1+m)
    pos[1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1])
    for i = 1:cl[1]
        bs = blocksize[1][i]
        if bs == 1
            pos[1][i] = @variable(model, lower_bound=0)
            pol -= basis[1][blocks[1][i][1]]^2 * pos[1][i]
        else
            pos[1][i] = @variable(model, [1:bs, 1:bs], PSD)
            pol -= sum(basis[1][blocks[1][i][j]] * pos[1][i][j,k] * basis[1][blocks[1][i][k]] for j = 1:bs, k = 1:bs)
        end
    end
    for k = 1:m
        pos[k+1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k+1])
        for i = 1:cl[k+1]
            bs = blocksize[k+1][i]
            if bs == 1
                pos[k+1][i] = @variable(model, lower_bound=0)
                pol -= basis[k+1][blocks[k+1][i][1]]^2 * pos[k+1][i] * ineq_cons[k]
            else
                pos[k+1][i] = @variable(model, [1:bs, 1:bs], PSD)
                pol -= sum(basis[k+1][blocks[k+1][i][j]] * pos[k+1][i][j,r] * basis[k+1][blocks[k+1][i][r]] for j = 1:bs, r = 1:bs) * ineq_cons[k]
            end
        end
    end
    mul = nothing
    if l > 0
        mul = Vector{Vector{VariableRef}}(undef, l)
        for k = 1:l
            mul[k] = @variable(model, [1:length(eblocks[k])])
            pol -= sum(basis[k+m+1][eblocks[k][j]] * mul[k][j] for j = 1:length(eblocks[k])) * eq_cons[k]
        end
    end
    coefs = MP.coefficients(pol, basis_covering_monomials(ChebyshevBasis, MP.monomials(pol)))
    remove_nearly_zero_terms!(model, coefs)
    drop_zeros!.(coefs)
    @constraint(model, coefs .== 0)
    info = poly_basis(nonneg, ineq_cons, eq_cons, nothing, nothing, nothing, basis, nothing, nothing, blocksize, blocks, eblocks, pos, mul, nothing)
    return info
end

function get_blocks(m::Int, l::Int, tsupp, ineq_cons, eq_cons, basis::Vector{ChebyshevBasisFirstKind{Poly{Float64}}}; TS="block", SO=1, merge=false, md=3, QUIET=false)
    blocks = Vector{Vector{Vector{Int}}}(undef, m+1)
    eblocks = Vector{Vector{Int}}(undef, l)
    blocksize = Vector{Vector{Int}}(undef, m+1)
    cl = Vector{Int}(undef, m+1)
    if TS == false
        for k = 1:m+1
            blocks[k],blocksize[k],cl[k] = [Vector(1:length(basis[k]))],[length(basis[k])],1       
        end
        for k = 1:l
            eblocks[k] = Vector(1:length(basis[k+m+1]))
        end
    else
        for i = 1:SO
            if i > 1
                oblocksize = deepcopy(blocksize)
                oeblocks = deepcopy(eblocks)
            end
            for k = 1:m+1
                if k == 1
                    G = get_graph(tsupp, basis[1])
                else
                    G = get_graph(tsupp, basis[k], g=ineq_cons[k-1])
                end
                if TS == "block"
                    blocks[k] = connected_components(G)
                    blocksize[k] = length.(blocks[k])
                    cl[k] = length(blocksize[k])
                else
                    blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS)
                    if merge == true
                        blocks[k],cl[k],blocksize[k] = clique_merge!(blocks[k], d=md, QUIET=true)
                    end
                end
            end
            for k = 1:l
                eblocks[k] = get_eblock(tsupp, eq_cons[k], basis[k+m+1])
            end
            if i > 1 && blocksize == oblocksize && eblocks == oeblocks
                println("No higher TS step of the TSSOS hierarchy!")
                break
            end
            if i < SO
                tsupp = Poly{Float64}[]
                for t = 1:length(blocks[1]), j = 1:blocksize[1][t], r = j:blocksize[1][t]
                    append!(tsupp, basis_covering_monomials(ChebyshevBasis, MP.monomials(basis[1][blocks[1][t][j]] * basis[1][blocks[1][t][r]])))
                end
                unique!(tsupp)
                sort!(tsupp)
            end
        end
    end
    return blocks,cl,blocksize,eblocks
end

function get_graph(tsupp, basis::ChebyshevBasisFirstKind{Poly{Float64}}; g=1)
    lb = length(basis)
    G = SimpleGraph(lb)
    ltsupp = length(tsupp)
    for i = 1:lb, j = i+1:lb
        flag = 0
        for item in basis_covering_monomials(ChebyshevBasis, MP.monomials(basis[i] * basis[j] * g))
            if lbfind(tsupp, ltsupp, item) !== nothing
                flag = 1
                break
            end
        end
        if flag == 1
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_eblock(tsupp, h, basis::ChebyshevBasisFirstKind{Poly{Float64}})
    ltsupp = length(tsupp)
    eblock = Int[]
    for (i, ba) in enumerate(basis)
        flag = 0
        for item in basis_covering_monomials(ChebyshevBasis, MP.monomials(ba * h))
            if lbfind(tsupp, ltsupp, item) !== nothing
                flag = 1
                break
            end
        end
        if flag == 1
            push!(eblock, i)
        end
    end
    return eblock
end

function Base.isless(P::T, Q::T) where {T<:AbstractPolynomial}
    return MP.monomials(P) < MP.monomials(Q) 
end

function remove_nearly_zero_terms!(model, ps; tol=1e-8)
    vars = all_variables(model)
    for (i, p) in enumerate(ps)
        np = AffExpr(0)
        for x in vars
            c = JuMP.coefficient(p, x)
            add_to_expression!(p, -c, x)
            if abs(c) > tol
                add_to_expression!(np, c, x)
            end
        end
        ps[i] = np + p
    end
    return ps
end
