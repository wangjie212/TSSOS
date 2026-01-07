mutable struct sos_data
    cliquesize # size of cliques
    cliques # clique structrue
    basis # monomial basis
    ebasis # monomial bases for equality constraints
    blocksize # size of blocks
    blocks # block structrue
    eblocks # block structrue for equality constraints
    tsupp # total support
    I # index sets of inequality constraints
    J # index sets of equality constraints
    GramMat # Gram matrices
    multiplier # multipliers for equality constraints
    constrs # constraint name
end

"""
    info = add_psatz!(model, nonneg, x, ineq_cons, eq_cons, order; CS=false, cliques=[], TS="block", eqTS=TS, 
    SO=1, GroebnerBasis=false, QUIET=false, constrs=nothing)

Add a Putinar's style SOS representation of the polynomial `nonneg` to the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `nonneg`: a nonnegative polynomial constrained to be a Putinar's style SOS on a semialgebraic set
- `x`: the set of POP variables
- `ineq_cons`: inequality constraints
- `eq_cons`: equality constraints
- `order`: relaxation order
- `CS`: method of chordal extension for correlative sparsity (`"MF"`, `"MD"`, `"NC"`, `false`)
- `cliques`: the set of cliques used in correlative sparsity
- `TS`: type of term sparsity (`"block"`, `"signsymmetry"`, `"MD"`, `"MF"`, `false`)
- `eqTS`: type of term sparsity for equality constraints (by default the same as `TS`, `false`)
- `SO`: sparse order
- `GroebnerBasis`: exploit the quotient ring structure or not (`true`, `false`)
- `QUIET`: run in the quiet mode (`true`, `false`)
- `constrs`: the constraint name used in the JuMP model

# Output arguments
- `info`: auxiliary data
"""
function add_psatz!(model, nonneg::Poly{T}, x, ineq_cons, eq_cons, order; CS=false, cliques=[], blocks=[], TS="block", eqTS=TS, SO=1, 
    GroebnerBasis=false, QUIET=false, constrs=nothing) where {T<:Union{Number,AffExpr}}
    n = length(x)
    g = [poly([UInt16[]], Float64[1]); poly[poly(p, x) for p in ineq_cons]]
    if GroebnerBasis == true && !isempty(eq_cons)
        h = poly[]
        gb = groebner(convert.(Poly{Float64}, eq_cons), ordering=DegRevLex())
        # f = poly(Groebner.normalform(gb, nonneg, ordering=DegRevLex()), x)
        f = poly(rem(nonneg, gb), x)
        lead = leading_ideal(gb, ordering=DegRevLex())
        leadsupp = [UInt16[] for i=1:length(lead)]
        for (i, mon) in enumerate(lead)
            mon = convert(DP.Monomial, mon)
            ind = mon.z .> 0
            vars = mon.vars[ind]
            exp = mon.z[ind]
            for j in eachindex(vars)
                append!(leadsupp[i], bfind_rev(x, vars[j])*ones(UInt16, exp[j]))
            end
        end
    else
        f = poly(nonneg, x)
        h = poly[poly(p, x) for p in eq_cons]
        gb = []
    end
    if !isempty(cliques)
        cql = length(cliques)
        cliquesize = length.(cliques)
    else
        CS = CS == true ? "MF" : CS
        cliques,cql,cliquesize = clique_decomp([f; g[2:end]; h], n, length(h), order=order+1, alg=CS, QUIET=QUIET)
    end
    ss = nothing
    if TS == "signsymmetry" || eqTS == "signsymmetry"
        ss = get_signsymmetry([f; g[2:end]; h], n)
    end
    I,J,_,_ = assign_constraint(g, h, cliques, cql)
    basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
    ebasis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
    for i = 1:cql
        basis[i] = Vector{Vector{Vector{UInt16}}}(undef, length(I[i]))
        ebasis[i] = Vector{Vector{Vector{UInt16}}}(undef, length(J[i]))
        for (s, k) in enumerate(I[i])
            basis[i][s] = get_basis(cliques[i], order-ceil(Int, maxdeg(g[k])/2))
        end
        for (s, k) in enumerate(J[i])
            ebasis[i][s] = get_basis(cliques[i], 2*order-maxdeg(h[k]))
        end
    end
    if isempty(blocks)
        blocks,cl,blocksize,eblocks = get_pblocks(f, g, h, I, J, cliques, cql, basis, ebasis, TS=TS, eqTS=eqTS, SO=SO, signsymmetry=ss)
    else
        eblocks = nothing
        blocksize = [[length.(block) for block in blocks[1]]]
        cl = [length.(blocksize[1])]
    end
    tsupp = Vector{UInt16}[]
    for i = 1:cql, j = 1:cl[i][1], k = 1:blocksize[i][1][j], r = k:blocksize[i][1][j]
        @inbounds bi = sadd(basis[i][1][blocks[i][1][j][k]], basis[i][1][blocks[i][1][j][r]])
        push!(tsupp, bi)
    end
    if TS != false && TS != "signsymmetry"
        for i = 1:cql
            for (j, w) in enumerate(I[i][2:end]), l = 1:cl[i][j+1], t = 1:blocksize[i][j+1][l], r = t:blocksize[i][j+1][l], item in g[w].supp
                @inbounds bi = sadd(basis[i][j+1][blocks[i][j+1][l][t]], item, basis[i][j+1][blocks[i][j+1][l][r]])
                push!(tsupp, bi)
            end
            for (j, w) in enumerate(J[i]), k in eblocks[i][j], item in h[w].supp
                @inbounds bi = sadd(ebasis[i][j][k], item)
                push!(tsupp, bi)
            end
        end
    end
    if !isempty(gb)
        unique!(tsupp)
        nsupp = Vector{UInt16}[]
        for item in tsupp
            if divide(item, leadsupp)
                append!(nsupp, reminder(item, x, gb).supp)
            else
                push!(nsupp, item)
            end
        end
        tsupp = nsupp
    end
    sort!(tsupp)
    unique!(tsupp)
    cons = [AffExpr(0) for i=1:length(tsupp)]
    pos = Vector{Vector{Vector{Symmetric{VariableRef}}}}(undef, cql)
    for i = 1:cql
        pos[i] = Vector{Vector{Symmetric{VariableRef}}}(undef, length(I[i]))
        for (j, w) in enumerate(I[i])
            pos[i][j] = Vector{Symmetric{VariableRef}}(undef, cl[i][j])
            for l = 1:cl[i][j]
                pos[i][j][l] = @variable(model, [1:blocksize[i][j][l], 1:blocksize[i][j][l]], PSD)
                for t = 1:blocksize[i][j][l], r = t:blocksize[i][j][l], (s, it) in enumerate(g[w].supp)
                    @inbounds bi = sadd(basis[i][j][blocks[i][j][l][t]], it, basis[i][j][blocks[i][j][l][r]])
                    if !isempty(gb) && divide(bi, leadsupp)
                        rem = reminder(bi, x, gb)
                        for (k, item) in enumerate(rem.supp)
                            Locb = bfind(tsupp, item)
                            if t == r
                                @inbounds add_to_expression!(cons[Locb], g[w].coe[s]*rem.coe[k], pos[i][j][l][t,r])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*g[w].coe[s]*rem.coe[k], pos[i][j][l][t,r])
                            end
                        end
                    else
                        Locb = bfind(tsupp, bi)
                        if t == r
                            @inbounds add_to_expression!(cons[Locb], g[w].coe[s], pos[i][j][l][t,r])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2*g[w].coe[s], pos[i][j][l][t,r])
                        end
                    end
                end
            end
        end
    end
    free = Vector{Vector{Vector{VariableRef}}}(undef, cql)
    for i = 1:cql
        free[i] = Vector{Vector{VariableRef}}(undef, length(J[i]))
        for (j, w) in enumerate(J[i])
            free[i][j] = @variable(model, [1:length(eblocks[i][j])])
            for (u, k) in enumerate(eblocks[i][j]), (s, item) in enumerate(h[w].supp)
                @inbounds bi = sadd(ebasis[i][j][k], item)
                Locb = bfind(tsupp, bi)
                @inbounds add_to_expression!(cons[Locb], h[w].coe[s], free[i][j][u])
            end
        end
    end
    for (i, item) in enumerate(f.supp)
        Locb = bfind(tsupp, item)
        if Locb === nothing
            @error "The monomial basis is not enough!"
        else
            cons[Locb] -= f.coe[i]
        end
    end
    if constrs !== nothing
        @constraint(model, cons==zeros(length(cons)), base_name=constrs)
    else
        @constraint(model, cons==zeros(length(cons)))
    end
    info = sos_data(cliquesize, cliques, basis, ebasis, blocksize, blocks, eblocks, tsupp, I, J, pos, free, constrs)
    return info
end

"""
    info = add_complex_psatz!(model, nonneg, x, ineq_cons, eq_cons, order; ipart=true, CS=false, cliques=[], TS="block", eqTS=TS, 
    SO=1, ConjugateBasis=false, normality=!ConjugateBasis, QUIET=false)

Add a complex Putinar's style Hermitian SOS representation of the complex polynomial `nonneg` to the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `nonneg`: a nonnegative complex polynomial constrained to be a complex Putinar's style Hermitian SOS on a semialgebraic set
- `x`: the set of POP variables
- `ineq_cons`: inequality constraints
- `eq_cons`: equality constraints
- `order`: relaxation order
- `ipart`: involving complex coefficients 
- `CS`: method of chordal extension for correlative sparsity (`"MF"`, `"MD"`, `"NC"`, `false`)
- `cliques`: the set of cliques used in correlative sparsity
- `TS`: type of term sparsity (`"block"`, `"MD"`, `"MF"`, `false`)
- `SO`: sparse order
- `ConjugateBasis`: including conjugate variables in monomial bases.
- `normality`: normal order
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `info`: auxiliary data
"""
function add_complex_psatz!(model, nonneg::Poly{T}, x, ineq_cons, eq_cons, order; ipart=true, CS=false, cliques=[], TS="block", eqTS=TS, SO=1, 
    ConjugateBasis=false, normality=!ConjugateBasis, QUIET=false) where {T<:Union{Number,AffExpr,GenericAffExpr}}
    n = length(x)
    if !isempty(ineq_cons)
        g = [cpoly([tuple(UInt16[], UInt16[])], [1]); cpoly[cpoly(p, x) for p in ineq_cons]]
    else
        g = [cpoly([tuple(UInt16[], UInt16[])], [1])]
    end
    if !isempty(eq_cons)
        h = [cpoly(p, x) for p in eq_cons]
    else
        h = cpoly[]
    end
    f = cpoly(nonneg, x)
    if !isempty(cliques)
        cql = length(cliques)
        cliquesize = length.(cliques)
    else
        time = @elapsed begin
        CS = CS == true ? "MF" : CS
        cliques,cql,cliquesize = clique_decomp([f; g[2:end]; h], n, order=order, alg=CS, QUIET=QUIET, ReducedCS=false)
        end
        if CS != false && QUIET == false
            mc = maximum(cliquesize)
            println("Obtained the variable cliques in $time seconds.\nThe maximal size of cliques is $mc.")
        end
    end
    I,J,_,_ = assign_constraint(g, h, cliques, cql)
    ebasis = Vector{Vector{Vector{Tuple{Vector{UInt16},Vector{UInt16}}}}}(undef, cql)
    if ConjugateBasis == false
        if normality == 0
            basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        else
            basis = Vector{Vector{Vector{Union{Vector{UInt16}, Tuple{Vector{UInt16},Vector{UInt16}}}}}}(undef, cql)
        end
    else
        basis = Vector{Vector{Vector{Tuple{Vector{UInt16},Vector{UInt16}}}}}(undef, cql)
    end
    for i = 1:cql
        ebasis[i] = Vector{Vector{Tuple{Vector{UInt16},Vector{UInt16}}}}(undef, length(J[i]))
        if ConjugateBasis == false
            if normality == 0
                basis[i] = [get_basis(cliques[i], order-maxcdeg(g[j])) for j in I[i]]
            else
                basis[i] = Vector{Vector{Union{Vector{UInt16}, Tuple{Vector{UInt16},Vector{UInt16}}}}}(undef, length(I[i])+cliquesize[i])
                basis[i][1] = get_basis(cliques[i], order)
                temp = get_basis(cliques[i], Int(normality))
                for s = 1:cliquesize[i]
                    basis[i][s+1] = [[tuple(item, UInt16[]) for item in temp]; [tuple(item, UInt16[cliques[i][s]]) for item in temp]]
                end
                for s = 1:length(I[i])-1
                    basis[i][s+1+cliquesize[i]] = get_basis(cliques[i], order-maxcdeg(g[I[i][s+1]]))
                end
                I[i] = [ones(Int, cliquesize[i]); I[i]]
            end
            for s = 1:length(J[i])
                if order < maxcdeg(h[J[i][s]])
                    @error "The relaxation order is too small!"
                end
                temp = get_basis(cliques[i], order-maxcdeg(h[J[i][s]]))
                ebasis[i][s] = vec([tuple(item1, item2) for item1 in temp, item2 in temp])
                sort!(ebasis[i][s])
            end
        else
            if normality < d
                basis[i] = Vector{Vector{Tuple{Vector{UInt16},Vector{UInt16}}}}(undef, length(I[i]))
                basis[i][1] = get_conjugate_basis(cliques[i], order)
                for s = 1:length(I[i])
                    basis[i][s] = get_conjugate_basis(cliques[i], order-Int(ceil(maxdeg(g[I[i][s]])/2)))
                end
            else
                basis[i] = Vector{Tuple{Vector{UInt16},Vector{UInt16}}}(undef, length(I[i])+cliquesize[i])
                basis[i][1] = get_conjugate_basis(cliques[i], order)
                temp = get_basis(cliques[i], normality)
                for s = 1:cliquesize[i]
                    basis[i][s+1] = [[tuple(item, UInt16[]) for item in temp]; [tuple(item, UInt16[cliques[i][s]]) for item in temp]]
                end
                for s = 1:length(I[i])-1
                    basis[i][s+1+cliquesize[i]] = get_conjugate_basis(cliques[i], order-Int(ceil(maxdeg(g[I[i][s+1]])/2)))
                end
                I[i] = [ones(Int, cliquesize[i]); I[i]]
            end
            for s = 1:length(J[i])
                ebasis[i][s] = get_conjugate_basis(cliques[i], 2*order-maxdeg(h[J[i][s]]))
            end
        end
    end
    blocks,cl,blocksize,eblocks = get_pblocks(f, g, h, I, J, order, cliques, cql, cliquesize, basis, ebasis, TS=TS, eqTS=eqTS, SO=SO, ConjugateBasis=ConjugateBasis, normality=normality)
    tsupp = Tuple{Vector{UInt16},Vector{UInt16}}[]
    for i = 1:cql
        if ConjugateBasis == false
            a = normality > 0 ? 1 + cliquesize[i] : 1
        else
            a = normality >= order ? 1 + cliquesize[i] : 1
        end
        for s = 1:a, j = 1:cl[i][s], k = 1:blocksize[i][s][j], r = k:blocksize[i][s][j]
            if ConjugateBasis == false && s == 1
                @inbounds bi = tuple(basis[i][s][blocks[i][s][j][k]], basis[i][s][blocks[i][s][j][r]])
            else
                @inbounds bi = sadd(basis[i][s][blocks[i][s][j][k]], conj(basis[i][s][blocks[i][s][j][r]]))
            end
            bi[1] <= bi[2] ? push!(tsupp, bi) : push!(tsupp, conj(bi))
        end
    end
    if TS != false
        csupp = get_csupp(order*ones(Int, cql), basis, ebasis, g, h, I, J, [], [], blocks, eblocks, cl, blocksize, cql, cliquesize, ConjugateBasis=ConjugateBasis, normality=normality)
        append!(tsupp, csupp)
    end
    sort!(tsupp)
    unique!(tsupp)
    rcons = [AffExpr(0) for i=1:length(tsupp)]
    if ipart == true
        icons = [AffExpr(0) for i=1:length(tsupp)]
    end
    pos = Vector{Vector{Vector{Symmetric{VariableRef}}}}(undef, cql)
    free = Vector{Vector{Vector{VariableRef}}}(undef, cql)
    for i = 1:cql
        pos[i] = Vector{Vector{Symmetric{VariableRef}}}(undef, length(I[i]))
        for (j, p) in enumerate(g[I[i]])
            pos[i][j] = Vector{Symmetric{VariableRef}}(undef, cl[i][j])
            for l = 1:cl[i][j]
                bs = blocksize[i][j][l]
                if ipart == true && bs != 1
                    pos[i][j][l] = @variable(model, [1:2bs, 1:2bs], PSD)
                else
                    pos[i][j][l] = @variable(model, [1:bs, 1:bs], PSD)
                end
                for t = 1:bs, r = 1:bs, (s, item) in enumerate(p.supp)
                    if typeof(basis[i][j][1]) == Vector{UInt16}
                        @inbounds bi = tuple(sadd(basis[i][j][blocks[i][j][l][t]], item[1]), sadd(basis[i][j][blocks[i][j][l][r]], item[2]))
                    else
                        @inbounds bi = tuple(sadd(basis[i][j][blocks[i][j][l][t]][1], item[1], basis[i][j][blocks[i][j][l][r]][2]), sadd(basis[i][j][blocks[i][j][l][t]][2], item[2], basis[i][j][blocks[i][j][l][r]][1]))
                    end
                    if bi[1] <= bi[2]
                        Locb = bfind(tsupp, bi)
                        if ipart == true
                            if bs != 1
                                @inbounds add_to_expression!(rcons[Locb], real(p.coe[s]), pos[i][j][l][t,r]+pos[i][j][l][t+bs,r+bs])
                                @inbounds add_to_expression!(rcons[Locb], -imag(p.coe[s]), pos[i][j][l][t,r+bs]-pos[i][j][l][r,t+bs])
                                if bi[1] != bi[2]
                                    @inbounds add_to_expression!(icons[Locb], imag(p.coe[s]), pos[i][j][l][t,r]+pos[i][j][l][t+bs,r+bs])
                                    @inbounds add_to_expression!(icons[Locb], real(p.coe[s]), pos[i][j][l][t,r+bs]-pos[i][j][l][r,t+bs])
                                end
                            else
                                @inbounds add_to_expression!(rcons[Locb], real(p.coe[s]), pos[i][j][l][t,r])
                                if bi[1] != bi[2]
                                    @inbounds add_to_expression!(icons[Locb], imag(p.coe[s]), pos[i][j][l][t,r])
                                end
                            end
                        else
                            @inbounds add_to_expression!(rcons[Locb], p.coe[s], pos[i][j][l][t,r])
                        end
                    end
                end
            end
        end
        free[i] = Vector{Vector{VariableRef}}(undef, length(J[i]))
        for (j, p) in enumerate(h[J[i]])
            mons = ebasis[i][j][eblocks[i][j]]
            temp = mons[[item[1] <= item[2] for item in mons]]
            lb = length(temp)
            if ipart == true
                free[i][j] = @variable(model, [1:2*lb])
            else
                free[i][j] = @variable(model, [1:lb])
            end
            for k in eblocks[i][j], (s, item) in enumerate(p.supp)
                @inbounds bi = tuple(sadd(ebasis[i][j][k][1], item[1]), sadd(ebasis[i][j][k][2], item[2]))
                if bi[1] <= bi[2]
                    Locb = bfind(tsupp, bi)
                    if ebasis[i][j][k][1] <= ebasis[i][j][k][2]
                        loc = bfind(temp, ebasis[i][j][k])
                        tag = 1
                        if ebasis[i][j][k][1] == ebasis[i][j][k][2]
                            tag = 0
                        end
                    else
                        loc = bfind(temp, conj(ebasis[i][j][k]))
                        tag = -1
                    end
                    if ipart == true
                        @inbounds add_to_expression!(rcons[Locb], real(p.coe[s])*free[i][j][loc]-tag*imag(p.coe[s])*free[i][j][loc+lb])
                        if bi[1] != bi[2]
                            @inbounds add_to_expression!(icons[Locb], tag*real(p.coe[s])*free[i][j][loc+lb]+imag(p.coe[s])*free[i][j][loc])
                        end
                    else
                        @inbounds add_to_expression!(rcons[Locb], p.coe[s], free[i][j][loc])
                    end
                end
            end
        end
    end
    if ipart == true
        ind = [item[1] != item[2] for item in tsupp]
        itsupp = tsupp[ind]
        icons = icons[ind]
    end
    ind = [item[1] <= item[2] for item in f.supp]
    nsupp,ncoe = f.supp[ind],f.coe[ind]
    for (i, item) in enumerate(nsupp)
        Locb = bfind(tsupp, item)
        if Locb === nothing
            @error "The monomial basis is not enough!"
        else
            rcons[Locb] -= real(ncoe[i])
            if ipart == true && item[1] != item[2]
                Locb = bfind(itsupp, item)
                icons[Locb] -= imag(ncoe[i])
            end
        end
    end
    @constraint(model, rcons .== 0)
    if ipart == true
        @constraint(model, icons .== 0)
    end
    info = sos_data(cliquesize, cliques, basis, ebasis, blocksize, blocks, eblocks, tsupp, I, J, pos, free, nothing)
    return info
end

function get_pblocks(f, g, h, I, J, cliques, cql, basis, ebasis; TS="block", eqTS=TS, SO=1, signsymmetry=nothing)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
    status = ones(Int, cql)
    if TS != false && TS != "signsymmetry"
        tsupp = vcat(f.supp, [p.supp for p in g]..., [p.supp for p in h]...)
        sort!(tsupp)
        unique!(tsupp)
    end
    for i = 1:cql
        supp = nothing
        if TS != false && TS != "signsymmetry"
            supp = tsupp[[issubset(item, cliques[i]) for item in tsupp]]
            for item in basis[i][1]
                push!(supp, sadd(item, item))
            end
            sort!(supp)
            unique!(supp)
        end
        blocks[i],cl[i],blocksize[i],eblocks[i],status[i] = get_pblocks(g[I[i]], h[J[i]], supp, basis[i], ebasis[i], TS=TS, eqTS=eqTS, SO=SO, signsymmetry=signsymmetry)
    end
    if minimum(status) == 1
        println("No higher TS step of the CS-TSSOS hierarchy!")
    end
    return blocks,cl,blocksize,eblocks
end

function get_pblocks(ineq_cons, eq_cons, tsupp, basis, ebasis; TS="block", eqTS=TS, SO=1, merge=false, md=3, signsymmetry=nothing)
    blocks = Vector{Vector{Vector{Int}}}(undef, length(ineq_cons))
    blocksize = Vector{Vector{Int}}(undef, length(ineq_cons))
    cl = Vector{Int}(undef, length(ineq_cons))
    status = 0
    if TS == false
        for k = 1:length(ineq_cons)
            blocks[k],blocksize[k],cl[k] = [Vector(1:length(basis[k]))],[length(basis[k])],1
        end
    else
        for i = 1:SO
            if i > 1
                oblocksize = deepcopy(blocksize)
                oeblocks = deepcopy(eblocks)
            end
            for k = 1:length(ineq_cons)
                G = get_graph(tsupp, ineq_cons[k].supp, basis[k], TS=TS, signsymmetry=signsymmetry)
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
            if eqTS != false
                eblocks = Vector{Vector{Int}}(undef, length(eq_cons))
                for k = 1:length(eq_cons)
                    eblocks[k] = get_eblock(tsupp, eq_cons[k].supp, ebasis[k], signsymmetry=signsymmetry)
                end
            end
            if i > 1 && blocksize == oblocksize && eblocks == oeblocks
                status = 1
                break
            end
            if i < SO
                tsupp = Vector{UInt16}[]
                for block in blocks[1], j = 1:length(block), r = j:length(block)
                    @inbounds bi = sadd(basis[1][block[j]], basis[1][block[r]])
                    push!(tsupp, bi)
                end
                sort!(tsupp)
                unique!(tsupp)
            end
        end
    end
    if eqTS == false
        eblocks = [Vector(1:length(ebasis[k])) for k = 1:length(eq_cons)]
    end
    return blocks,cl,blocksize,eblocks,status
end

function get_pblocks(f, g, h, I, J, order, cliques, cql, cliquesize, basis, ebasis; TS="block", eqTS=TS, SO=1, ConjugateBasis=false, normality=1)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
    status = ones(Int, cql)
    if TS != false
        tsupp = filter(item -> item[1] <= item[2], vcat([p.supp for p in [f; g; h]]...))
        sort!(tsupp)
        unique!(tsupp)
    end
    for i = 1:cql
        ksupp = TS == false ? nothing : tsupp[[issubset(union(item[1], item[2]), cliques[i]) for item in tsupp]]
        blocks[i],cl[i],blocksize[i],eblocks[i],status[i] = get_pblocks(ksupp, order, cliquesize[i], g[I[i]], h[J[i]], basis[i], ebasis[i], TS=TS, eqTS=eqTS, SO=SO, ConjugateBasis=ConjugateBasis, normality=normality)
    end
    if minimum(status) == 1
        println("No higher TS step of the CS-TSSOS hierarchy!")
    end
    return blocks,cl,blocksize,eblocks
end

function get_pblocks(tsupp, order, n, ineq_cons::Vector{T1}, eq_cons::Vector{T2}, basis, ebasis; TS="block", eqTS=TS, SO=1, ConjugateBasis=false, normality=1) where {T1,T2<:cpoly}
    blocks = Vector{Vector{Vector{Int}}}(undef, length(ineq_cons))
    blocksize = Vector{Vector{Int}}(undef,length(ineq_cons))
    cl = Vector{Int}(undef, length(ineq_cons))
    status = 0
    if TS == false
        for k = 1:length(ineq_cons) 
            blocks[k],blocksize[k],cl[k] = [Vector(1:length(basis[k]))],[length(basis[k])],1
        end
    else
        for i = 1:SO
            if i > 1
                oblocksize = deepcopy(blocksize)
                oeblocks = deepcopy(eblocks)
            end
            for (k, p) in enumerate(ineq_cons)
                G = get_graph(tsupp, p.supp, basis[k])
                if TS == "block"
                    blocks[k] = connected_components(G)
                    blocksize[k] = length.(blocks[k])
                    cl[k] = length(blocksize[k])
                else
                    blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS)
                end
            end
            if eqTS != false
                eblocks = [get_eblock(tsupp, p.supp, ebasis[k]) for (k, p) in enumerate(eq_cons)]
            end
            if i > 1 && blocksize == oblocksize && eblocks == oeblocks
                status = 1
                break
            end
            if i < SO
                tsupp = Tuple{Vector{UInt16},Vector{UInt16}}[]
                if ConjugateBasis == false
                    a = normality > 0 ? 1 + n : 1
                else
                    a = normality >= order ? 1 + n : 1
                end
                for s = 1:a, j = 1:cl[s], k = 1:blocksize[s][j], r = k:blocksize[s][j]
                    if ConjugateBasis == false && s == 1
                        @inbounds bi = tuple(basis[s][blocks[s][j][k]], basis[s][blocks[s][j][r]])
                    else
                        @inbounds bi = sadd(basis[s][blocks[s][j][k]], conj(basis[s][blocks[s][j][r]]))
                    end
                    bi[1] <= bi[2] ? push!(tsupp, bi) : push!(tsupp, conj(bi))
                end
            end
            csupp = get_csupp([order], [basis], [ebasis], ineq_cons, eq_cons, [Vector(1:length(ineq_cons))], [Vector(1:length(eq_cons))], [], [], [blocks], [eblocks], [cl], [blocksize], 1, [n], ConjugateBasis=ConjugateBasis, normality=normality)
            append!(tsupp, csupp)
            sort!(tsupp)
            unique!(tsupp)
        end
    end
    if eqTS == false
        eblocks = [Vector(1:length(ebasis[k])) for k = 1:length(eq_cons)]
    end
    return blocks,cl,blocksize,eblocks,status
end

function get_moment(supp::Vector{Vector{UInt16}}, lb, ub)
    v = zeros(length(supp))
    for (i, item) in enumerate(supp)
        temp = zeros(Int, length(lb))
        for k in item
            temp[k] += 1
        end
        v[i] = prod([(ub[j]^(temp[j]+1)-lb[j]^(temp[j]+1))/(temp[j]+1) for j=1:length(lb)])
    end
    return v
end

function get_moment(mons::Vector{DP.Monomial{V, M}}, x, lb, ub) where {V, M}
    return get_moment([exps(mon, x) for mon in mons], lb, ub)
end

function get_moment_matrix(moment, info)
    MomMat = Vector{Symmetric{Float64}}(undef, info.cql)
    for i = 1:info.cql
        lb = length(info.basis[i][1])
        mmat = zeros(Float64, lb, lb)
        for j = 1:lb, k = j:lb
            bi = sadd(info.basis[i][1][j], info.basis[i][1][k])
            Locb = bfind(info.tsupp, bi)
            if Locb !== nothing
                mmat[j,k] = moment[Locb]
            end
        end
        MomMat[i] = Symmetric(mmat, :U)
    end
    return MomMat
end
