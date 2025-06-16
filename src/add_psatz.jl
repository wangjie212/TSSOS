mutable struct struct_data
    cql # number of cliques
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
    info = add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order; CS=false, cliques=[], TS="block", 
    SO=1, GroebnerBasis=false, QUIET=false, constrs=nothing)

Add a Putinar's style SOS representation of the polynomial `nonneg` to the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `nonneg`: a nonnegative polynomial constrained to be a Putinar's style SOS on a semialgebraic set
- `vars`: the set of POP variables
- `ineq_cons`: inequality constraints
- `eq_cons`: equality constraints
- `order`: relaxation order
- `CS`: method of chordal extension for correlative sparsity (`"MF"`, `"MD"`, `"NC"`, `false`)
- `cliques`: the set of cliques used in correlative sparsity
- `TS`: type of term sparsity (`"block"`, `"signsymmetry"`, `"MD"`, `"MF"`, `false`)
- `SO`: sparse order
- `GroebnerBasis`: exploit the quotient ring structure or not (`true`, `false`)
- `QUIET`: run in the quiet mode (`true`, `false`)
- `constrs`: the constraint name used in the JuMP model

# Output arguments
- `info`: auxiliary data
"""
function add_psatz!(model, nonneg::Poly{T}, vars, ineq_cons, eq_cons, order; CS=false, cliques=[], blocks=[], TS="block", SO=1, 
    GroebnerBasis=false, QUIET=false, constrs=nothing) where {T<:Union{Number,AffExpr}}
    n = length(vars)
    m = length(ineq_cons)
    if ineq_cons != []
        gsupp,gcoe = npolys_info(ineq_cons, vars)
        glt = length.(gcoe)
        dg = maxdegree.(ineq_cons)
    else
        gsupp = Matrix{UInt8}[]
        glt = dg = Int[]
    end
    if eq_cons != []
        hsupp,hcoe = npolys_info(eq_cons, vars)
        hlt = length.(hcoe)
        dh = maxdegree.(eq_cons)
    else
        hsupp = Matrix{UInt8}[]
        hlt = dh = Int[]
    end
    if GroebnerBasis == true && eq_cons != []
        l = 0
        gb = convert.(Poly{Float64}, eq_cons)
        SemialgebraicSets.gröbner_basis!(gb)
        nonneg = rem(nonneg, gb)
        leadm = SemialgebraicSets.leading_monomial.(gb)
        llead = length(leadm)
        lead = zeros(UInt8, n, llead)
        for i = 1:llead, j = 1:n
            @inbounds lead[j,i] = MP.degree(leadm[i], vars[j])
        end
    else
        gb = []
        l = length(eq_cons)
    end
    fsupp,fcoe = poly_info(nonneg, vars)
    if CS != false
        if cliques == []
            CS = CS == true ? "MF" : CS
            cliques,cql,cliquesize = clique_decomp(n, m, length(eq_cons), fsupp, gsupp, hsupp, alg=CS, QUIET=QUIET)
        else
            cql = length(cliques)
            cliquesize = length.(cliques)
        end
    else
        cliques,cql,cliquesize = [Vector(1:n)],1,[n]
    end
    ss = nothing
    if TS == "signsymmetry"
        temp = fsupp
        if ineq_cons != []
            temp = [temp reduce(hcat, gsupp)]
        end
        if eq_cons != []
            temp = [temp reduce(hcat, hsupp)]
        end
        ss = get_signsymmetry(permutedims(temp, [2,1]))
    end
    dmin = ceil(Int, maximum([maxdegree(nonneg); dg; dh])/2)
    order = order < dmin ? dmin : order
    I,J = assign_constraint(m, l, gsupp, hsupp, cliques, cql)
    basis = Vector{Vector{Matrix{UInt8}}}(undef, cql)
    for t = 1:cql
        basis[t] = Vector{Matrix{UInt8}}(undef, length(I[t])+length(J[t])+1)
        basis[t][1] = get_basis(n, order, var=cliques[t])
        for s = 1:length(I[t])
            basis[t][s+1] = get_basis(n, order-ceil(Int, dg[I[t][s]]/2), var=cliques[t])
        end
        for s = 1:length(J[t])
            basis[t][s+length(I[t])+1] = get_basis(n, 2*order-dh[J[t][s]], var=cliques[t])
        end
    end
    if isempty(blocks)
        blocks,cl,blocksize,eblocks = get_blocks(n, I, J, fsupp, gsupp, hsupp, basis, cliques, cql, TS=TS, SO=SO, signsymmetry=ss)
    else
        eblocks = nothing
        blocksize = [[length.(blocks[1][i]) for i = 1:length(blocks[1])]]
        cl = [length.(blocksize[1])]
    end
    ne = 0
    for t = 1:cql
        ne += sum(numele(blocksize[t][1]))
        if I[t] != []
            ne += sum(glt[I[t][k]]*numele(blocksize[t][k+1]) for k=1:length(I[t]))
        end
        if J[t] != []
            ne += sum(hlt[J[t][k]]*length(eblocks[t][k]) for k=1:length(J[t]))
        end
    end
    tsupp = zeros(UInt8, n, ne)
    q = 1
    for i = 1:cql
        for j = 1:cl[i][1], k = 1:blocksize[i][1][j], r = k:blocksize[i][1][j]
            @inbounds bi = basis[i][1][:, blocks[i][1][j][k]] + basis[i][1][:, blocks[i][1][j][r]]
            tsupp[:, q] = bi
            q += 1
        end
        if TS != false && TS != "signsymmetry"
            for (j, w) in enumerate(I[i]), p = 1:cl[i][j+1], t = 1:blocksize[i][j+1][p], r = t:blocksize[i][j+1][p], s = 1:glt[w]
                ind1 = blocks[i][j+1][p][t]
                ind2 = blocks[i][j+1][p][r]
                @inbounds bi = basis[i][j+1][:, ind1] + basis[i][j+1][:, ind2] + gsupp[w][:, s]
                tsupp[:, q] = bi
                q += 1
            end
            for (j, w) in enumerate(J[i]), t in eblocks[i][j], s = 1:hlt[w]
                @inbounds bi = basis[i][j+length(I[i])+1][:, t] + hsupp[w][:, s]
                tsupp[:, q] = bi
                q += 1
            end
        end
    end
    if !isempty(gb)
        tsupp = unique(tsupp, dims=2)
        nsupp = zeros(UInt8, n)
        for col in eachcol(tsupp)
            if divide(col, lead, n, llead)
                temp = reminder(col, vars, gb, n)[2]
                nsupp = [nsupp temp]
            else
                nsupp = [nsupp col]
            end
        end
        tsupp = nsupp
    end
    tsupp = sortslices(tsupp, dims=2)
    tsupp = unique(tsupp, dims=2)
    ltsupp = size(tsupp, 2)
    cons = [AffExpr(0) for i=1:ltsupp]
    pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, cql)
    mul = nothing
    if l > 0
        mul = Vector{Vector{Vector{VariableRef}}}(undef, cql)
    end
    for t = 1:cql
        pos[t] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, 1+length(I[t]))
        pos[t][1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[t][1])
        for i = 1:cl[t][1]
            bs = blocksize[t][1][i]
            if bs == 1
               pos[t][1][i] = @variable(model, lower_bound=0)
               bi = 2*basis[t][1][:, blocks[t][1][i][1]]
               if !isempty(gb) && divide(bi, lead, n, llead)
                    bi_lm,bi_supp,bi_coe = reminder(bi, vars, gb, n)
                    for z = 1:bi_lm
                        Locb = bfind(tsupp, ltsupp, bi_supp[:,z])
                        @inbounds add_to_expression!(cons[Locb], bi_coe[z], pos[t][1][i])
                    end
               else
                    Locb = bfind(tsupp, ltsupp, bi)
                    @inbounds add_to_expression!(cons[Locb], pos[t][1][i])
               end
            else
               pos[t][1][i] = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:bs, r = j:bs
                   bi = basis[t][1][:, blocks[t][1][i][j]] + basis[t][1][:, blocks[t][1][i][r]]
                   if !isempty(gb) && divide(bi, lead, n, llead)
                        bi_lm,bi_supp,bi_coe = reminder(bi, vars, gb, n)
                        for z = 1:bi_lm
                            Locb = bfind(tsupp, ltsupp, bi_supp[:,z])
                            if j == r
                                @inbounds add_to_expression!(cons[Locb], bi_coe[z], pos[t][1][i][j,r])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*bi_coe[z], pos[t][1][i][j,r])
                            end
                        end
                   else
                        Locb = bfind(tsupp, ltsupp, bi)
                        if j == r
                            @inbounds add_to_expression!(cons[Locb], pos[t][1][i][j,r])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2, pos[t][1][i][j,r])
                        end
                    end
               end
            end
        end
        for k = 1:length(I[t])
            pos[t][k+1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[t][k+1])
            for i = 1:cl[t][k+1]
                bs = blocksize[t][k+1][i]
                if bs == 1
                    pos[t][k+1][i] = @variable(model, lower_bound=0)
                    for s = 1:glt[I[t][k]]
                        bi = 2*basis[t][k+1][:, blocks[t][k+1][i][1]] + gsupp[I[t][k]][:,s]
                        if !isempty(gb) && divide(bi, lead, n, llead)
                            bi_lm,bi_supp,bi_coe = reminder(bi, vars, gb, n)
                            for z = 1:bi_lm
                                Locb = bfind(tsupp, ltsupp, bi_supp[:,z])
                                @inbounds add_to_expression!(cons[Locb], gcoe[I[t][k]][s]*bi_coe[z], pos[t][k+1][i])
                            end
                        else
                            Locb = bfind(tsupp, ltsupp, bi)
                            @inbounds add_to_expression!(cons[Locb], gcoe[I[t][k]][s], pos[t][k+1][i])
                        end
                    end
                else
                    pos[t][k+1][i] = @variable(model, [1:bs, 1:bs], PSD)
                    for j = 1:bs, r = j:bs, s = 1:glt[I[t][k]]
                        bi = basis[t][k+1][:, blocks[t][k+1][i][j]] + basis[t][k+1][:, blocks[t][k+1][i][r]] + gsupp[I[t][k]][:,s]
                        if !isempty(gb) && divide(bi, lead, n, llead)
                            bi_lm,bi_supp,bi_coe = reminder(bi, vars, gb, n)
                            for z = 1:bi_lm
                                Locb = bfind(tsupp, ltsupp, bi_supp[:,z])
                                if j == r
                                    @inbounds add_to_expression!(cons[Locb], gcoe[I[t][k]][s]*bi_coe[z], pos[t][k+1][i][j,r])
                                else
                                    @inbounds add_to_expression!(cons[Locb], 2*gcoe[I[t][k]][s]*bi_coe[z], pos[t][k+1][i][j,r])
                                end
                            end
                        else
                            Locb = bfind(tsupp, ltsupp, bi)
                            if j == r
                                @inbounds add_to_expression!(cons[Locb], gcoe[I[t][k]][s], pos[t][k+1][i][j,r])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*gcoe[I[t][k]][s], pos[t][k+1][i][j,r])
                            end
                        end
                    end
                end
            end
        end
        if l > 0
            mul[t] = Vector{Vector{VariableRef}}(undef, length(J[t]))
            for k = 1:length(J[t])
                bs = length(eblocks[t][k])
                mul[t][k] = @variable(model, [1:bs])
                for i = 1:bs, s = 1:hlt[J[t][k]]
                    bi = basis[t][k+length(I[t])+1][:, eblocks[t][k][i]] + hsupp[J[t][k]][:, s]
                    Locb = bfind(tsupp, ltsupp, bi)
                    @inbounds add_to_expression!(cons[Locb], hcoe[J[t][k]][s], mul[t][k][i])
                end
            end
        end
    end
    for i = 1:size(fsupp, 2)
        Locb = bfind(tsupp, ltsupp, fsupp[:, i])
        if Locb === nothing
            @error "The monomial basis is not enough!"
        else
            cons[Locb] -= fcoe[i]
        end
    end
    if constrs !== nothing
        @constraint(model, cons==zeros(ltsupp), base_name=constrs)
    else
        @constraint(model, cons==zeros(ltsupp))
    end
    info = struct_data(cql, cliquesize, cliques, basis, nothing, blocksize, blocks, eblocks, tsupp, I, J, pos, mul, constrs)
    return info
end

"""
    info = add_complex_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order; ipart=true, CS=false, cliques=[], TS="block", 
    SO=1, ConjugateBasis=false, normality=!ConjugateBasis, QUIET=false)

Add a complex Putinar's style Hermitian SOS representation of the complex polynomial `nonneg` to the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `nonneg`: a nonnegative complex polynomial constrained to be a complex Putinar's style Hermitian SOS on a semialgebraic set
- `vars`: the set of POP variables
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
function add_complex_psatz!(model, nonneg::Poly{T}, vars, ineq_cons, eq_cons, order; ipart=true, CS=false, cliques=[], TS="block", SO=1, 
    ConjugateBasis=false, normality=!ConjugateBasis, QUIET=false) where {T<:Union{Number,AffExpr,GenericAffExpr}}
    n = length(vars)
    m = length(ineq_cons)
    l = length(eq_cons)
    coe_type = ipart == true ? ComplexF64 : Float64
    if ineq_cons != []
        gsupp,gcoe = cpolys_info(ineq_cons, vars, ctype=coe_type)
        glt = length.(gcoe)
    else
        gsupp = Vector{Vector{Vector{UInt16}}}[]
        glt = Int[]
    end
    if eq_cons != []
        hsupp,hcoe = cpolys_info(eq_cons, vars, ctype=coe_type)
        hlt = length.(hcoe)
    else
        hsupp = Vector{Vector{Vector{UInt16}}}[]
        hlt = Int[]
    end
    fsupp,fcoe = cpolys_info([nonneg], vars, ctype=Union{coe_type, AffExpr, GenericAffExpr})
    fsupp,fcoe = fsupp[1],fcoe[1]
    dg = zeros(Int, m)
    for i = 1:m
        if ConjugateBasis == false
            dg[i] = maximum([max(length(item[1]), length(item[2])) for item in gsupp[i]])
        else
            dg[i] = ceil(Int, maximum([length(item[1]) + length(item[2]) for item in gsupp[i]])/2)
        end
    end
    dh = zeros(Int, l)
    for i = 1:l
        if ConjugateBasis == false
            dh[i] = maximum([max(length(item[1]), length(item[2])) for item in hsupp[i]])
        else
            dh[i] = maximum([length(item[1]) + length(item[2]) for item in hsupp[i]])
        end
    end
    if cliques != []
        cql = length(cliques)
        cliquesize = length.(cliques)
    else
        time = @elapsed begin
        CS = CS == true ? "MF" : CS
        cliques,cql,cliquesize = clique_decomp(n, m+l, [dg;dh], [[fsupp];gsupp;hsupp], order=order, alg=CS, ReducedCS=false)
        end
        if CS != false && QUIET == false
            mc = maximum(cliquesize)
            println("Obtained the variable cliques in $time seconds.\nThe maximal size of cliques is $mc.")
        end
    end
    I = assign_constraint(gsupp, cliques, cql)
    J = assign_constraint(hsupp, cliques, cql)
    ebasis = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, cql)
    if ConjugateBasis == false
        if normality == 0
            basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        else
            basis = Vector{Vector{Vector{Union{Vector{UInt16}, Vector{Vector{UInt16}}}}}}(undef, cql)
        end
    else
        basis = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, cql)
    end
    for i = 1:cql
        ebasis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(J[i]))
        if ConjugateBasis == false
            if normality == 0
                basis[i] = Vector{Vector{Vector{UInt16}}}(undef, length(I[i])+1)
                basis[i][1] = get_basis(cliques[i], order)
                for s = 1:length(I[i])
                    basis[i][s+1] = get_basis(cliques[i], order-dg[I[i][s]])
                end
            else
                basis[i] = Vector{Vector{Union{Vector{UInt16}, Vector{Vector{UInt16}}}}}(undef, length(I[i])+1+cliquesize[i])
                basis[i][1] = get_basis(cliques[i], order)
                for s = 1:cliquesize[i]
                    temp = get_basis(cliques[i], Int(normality))
                    basis[i][s+1] = [[[item, UInt16[]] for item in temp]; [[item, UInt16[cliques[i][s]]] for item in temp]]
                end
                for s = 1:length(I[i])
                    basis[i][s+1+cliquesize[i]] = get_basis(cliques[i], order-dg[I[i][s]])
                end
            end
            for s = 1:length(J[i])
                if order < dh[J[i][s]]
                    @error "The relaxation order is too small!"
                end
                temp = get_basis(cliques[i], order-dh[J[i][s]])
                ebasis[i][s] = vec([[item1, item2] for item1 in temp, item2 in temp])
                sort!(ebasis[i][s])
            end
        else
            if normality < order
                basis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(I[i])+1)
                basis[i][1] = get_conjugate_basis(cliques[i], order)
                for s = 1:length(I[i])
                    basis[i][s+1] = get_conjugate_basis(cliques[i], order-dg[I[i][s]])
                end
            else
                basis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(I[i])+1+cliquesize[i])
                basis[i][1] = get_conjugate_basis(cliques[i], order)
                for s = 1:cliquesize[i]
                    temp = get_basis(cliques[i], normality)
                    basis[i][s+1] = [[[item, UInt16[]] for item in temp]; [[item, UInt16[cliques[i][s]]] for item in temp]]
                end
                for s = 1:length(I[i])
                    basis[i][s+1+cliquesize[i]] = get_conjugate_basis(cliques[i], order-dg[I[i][s]])
                end
            end
            for s = 1:length(J[i])
                ebasis[i][s] = get_conjugate_basis(cliques[i], 2*order-dh[J[i][s]])
            end
        end
    end
    blocks,cl,blocksize,eblocks = get_blocks(order, I, J, fsupp, gsupp, hsupp, basis, ebasis, cliques, cliquesize, cql, TS=TS, SO=SO, ConjugateBasis=ConjugateBasis, normality=normality)
    tsupp = Vector{Vector{UInt16}}[]
    for i = 1:cql
        if ConjugateBasis == false
            a = normality > 0 ? 1 + cliquesize[i] : 1
        else
            a = normality >= order ? 1 + cliquesize[i] : 1
        end
        for s = 1:a, j = 1:cl[i][s], k = 1:blocksize[i][s][j], r = k:blocksize[i][s][j]
            if ConjugateBasis == false && s == 1
                @inbounds bi = [basis[i][s][blocks[i][s][j][k]], basis[i][s][blocks[i][s][j][r]]]
            else
                @inbounds bi = [sadd(basis[i][s][blocks[i][s][j][k]][1], basis[i][s][blocks[i][s][j][r]][2]), sadd(basis[i][s][blocks[i][s][j][k]][2], basis[i][s][blocks[i][s][j][r]][1])]
            end
            bi[1] <= bi[2] ? push!(tsupp, bi) : push!(tsupp, bi[2:-1:1])
        end
    end
    if TS != false
        csupp = get_csupp(order, basis, ebasis, gsupp, hsupp, cql, I, J, blocks, eblocks, cl, blocksize, cliquesize, ConjugateBasis=ConjugateBasis, normality=normality)
        append!(tsupp, csupp)
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp = length(tsupp)
    rcons = [AffExpr(0) for i=1:ltsupp]
    if ipart == true
        icons = [AffExpr(0) for i=1:ltsupp]
    end
    pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, cql)
    for i = 1:cql
        if ConjugateBasis == false
            a = normality > 0 ? 1 + cliquesize[i] : 1
        else
            a = normality >= order ? 1 + cliquesize[i] : 1
        end
        pos[i] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, a+length(I[i]))
        for j = 1:a
            pos[i][j] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i][j])
            for l = 1:cl[i][j]
                @inbounds bs = blocksize[i][j][l]
                if bs == 1
                    pos[i][j][l] = @variable(model, lower_bound=0)
                    if ConjugateBasis == false && j == 1
                        @inbounds bi = [basis[i][j][blocks[i][j][l][1]], basis[i][j][blocks[i][j][l][1]]]
                    else
                        @inbounds bi = [sadd(basis[i][j][blocks[i][j][l][1]][1], basis[i][j][blocks[i][j][l][1]][2]), sadd(basis[i][j][blocks[i][j][l][1]][1], basis[i][j][blocks[i][j][l][1]][2])]
                    end
                    Locb = bfind(tsupp, ltsupp, bi)
                    @inbounds add_to_expression!(rcons[Locb], pos[i][j][l])
                else
                    if ipart == true
                        pos[i][j][l] = @variable(model, [1:2bs, 1:2bs], PSD)
                    else
                        pos[i][j][l] = @variable(model, [1:bs, 1:bs], PSD)
                    end
                    for t = 1:bs, r = t:bs
                        @inbounds ind1 = blocks[i][j][l][t]
                        @inbounds ind2 = blocks[i][j][l][r]
                        if ConjugateBasis == false && j == 1
                            @inbounds bi = [basis[i][j][ind1], basis[i][j][ind2]]
                        else
                            @inbounds bi = [sadd(basis[i][j][ind1][1], basis[i][j][ind2][2]), sadd(basis[i][j][ind1][2], basis[i][j][ind2][1])]
                        end
                        if bi[1] <= bi[2]
                            Locb = bfind(tsupp, ltsupp, bi)
                            if ipart == true && bi[1] != bi[2]
                                @inbounds add_to_expression!(icons[Locb], pos[i][j][l][t,r+bs]-pos[i][j][l][r,t+bs])
                            end
                        else
                            Locb = bfind(tsupp, ltsupp, bi[2:-1:1])
                            if ipart == true
                                @inbounds add_to_expression!(icons[Locb], -1, pos[i][j][l][t,r+bs]-pos[i][j][l][r,t+bs])
                            end
                        end
                        if ipart == true
                            if bi[1] == bi[2] && t != r
                                @inbounds add_to_expression!(rcons[Locb], 2*pos[i][j][l][t,r]+2*pos[i][j][l][t+bs,r+bs])
                            else
                                @inbounds add_to_expression!(rcons[Locb], pos[i][j][l][t,r]+pos[i][j][l][t+bs,r+bs])
                            end
                        else
                            if bi[1] == bi[2] && t != r
                                @inbounds add_to_expression!(rcons[Locb], 2*pos[i][j][l][t,r])
                            else
                                @inbounds add_to_expression!(rcons[Locb], pos[i][j][l][t,r])
                            end
                        end
                    end
                end
            end
        end
    end
    for i = 1:cql, (j, w) in enumerate(I[i])
        if ConjugateBasis == false
            a = normality > 0 ? 1 + cliquesize[i] : 1
        else
            a = normality >= order ? 1 + cliquesize[i] : 1
        end
        pos[i][j+a] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i][j+a])
        for l = 1:cl[i][j+a]
            bs = blocksize[i][j+a][l]
            if bs == 1
                pos[i][j+a][l] = @variable(model, lower_bound=0)
                ind = blocks[i][j+a][l][1]
                for s = 1:length(gsupp[w])
                    if ConjugateBasis == false
                        @inbounds bi = [sadd(basis[i][j+a][ind], gsupp[w][s][1]), sadd(basis[i][j+a][ind], gsupp[w][s][2])]
                    else
                        @inbounds bi = [sadd(sadd(basis[i][j+a][ind][1], gsupp[w][s][1]), basis[i][j+a][ind][2]), sadd(sadd(basis[i][j+a][ind][2], gsupp[w][s][2]), basis[i][j+a][ind][1])]
                    end
                    if bi[1] <= bi[2]
                        Locb = bfind(tsupp, ltsupp, bi)
                        if ipart == true
                            @inbounds add_to_expression!(icons[Locb], imag(gcoe[w][s]), pos[i][j+a][l])
                        end
                        @inbounds add_to_expression!(rcons[Locb], real(gcoe[w][s]), pos[i][j+a][l])
                    end
                end
            else
                if ipart == true
                    pos[i][j+a][l] = @variable(model, [1:2bs, 1:2bs], PSD)
                else
                    pos[i][j+a][l] = @variable(model, [1:bs, 1:bs], PSD)
                end
                for t = 1:bs, r = 1:bs
                    ind1 = blocks[i][j+a][l][t]
                    ind2 = blocks[i][j+a][l][r]
                    for s = 1:length(gsupp[w])
                        if ConjugateBasis == false
                            @inbounds bi = [sadd(basis[i][j+a][ind1], gsupp[w][s][1]), sadd(basis[i][j+a][ind2], gsupp[w][s][2])]
                        else
                            @inbounds bi = [sadd(sadd(basis[i][j+a][ind1][1], gsupp[w][s][1]), basis[i][j+a][ind2][2]), sadd(sadd(basis[i][j+a][ind1][2], gsupp[w][s][2]), basis[i][j+a][ind2][1])]
                        end
                        if bi[1] <= bi[2]
                            Locb = bfind(tsupp, ltsupp, bi)
                            if ipart == true
                                @inbounds add_to_expression!(rcons[Locb], real(gcoe[w][s]), pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs])
                                @inbounds add_to_expression!(rcons[Locb], -imag(gcoe[w][s]), pos[i][j+a][l][t,r+bs]-pos[i][j+a][l][r,t+bs])
                                if bi[1] != bi[2]
                                    @inbounds add_to_expression!(icons[Locb], imag(gcoe[w][s]), pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs])
                                    @inbounds add_to_expression!(icons[Locb], real(gcoe[w][s]), pos[i][j+a][l][t,r+bs]-pos[i][j+a][l][r,t+bs])
                                end
                            else
                                @inbounds add_to_expression!(rcons[Locb], real(gcoe[w][s]), pos[i][j+a][l][t,r])
                            end
                        end
                    end
                end
            end
        end
    end
    free = Vector{Vector{Vector{VariableRef}}}(undef, cql)
    for i = 1:cql
        if !isempty(J[i])
            free[i] = Vector{Vector{VariableRef}}(undef, length(J[i]))
            for (j, w) in enumerate(J[i])
                mons = ebasis[i][j][eblocks[i][j]]
                temp = mons[[item[1] <= item[2] for item in mons]]
                lb = length(temp)
                if ipart == true
                    free[i][j] = @variable(model, [1:2*lb])
                else
                    free[i][j] = @variable(model, [1:lb])
                end
                for k in eblocks[i][j], s = 1:length(hsupp[w])
                    @inbounds bi = [sadd(ebasis[i][j][k][1], hsupp[w][s][1]), sadd(ebasis[i][j][k][2], hsupp[w][s][2])]
                    if bi[1] <= bi[2]
                        Locb = bfind(tsupp, ltsupp, bi)
                        if ebasis[i][j][k][1] <= ebasis[i][j][k][2]
                            loc = bfind(temp, lb, ebasis[i][j][k])
                            tag = 1
                            if ebasis[i][j][k][1] == ebasis[i][j][k][2]
                                tag = 0
                            end
                        else
                            loc = bfind(temp, lb, ebasis[i][j][k][2:-1:1])
                            tag = -1
                        end
                        if ipart == true
                            @inbounds add_to_expression!(rcons[Locb], real(hcoe[w][s])*free[i][j][loc]-tag*imag(hcoe[w][s])*free[i][j][loc+lb])
                            if bi[1] != bi[2]
                                @inbounds add_to_expression!(icons[Locb], tag*real(hcoe[w][s])*free[i][j][loc+lb]+imag(hcoe[w][s])*free[i][j][loc])
                            end
                        else
                            @inbounds add_to_expression!(rcons[Locb], real(hcoe[w][s]), free[i][j][loc])
                        end
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
    ind = [item[1] <= item[2] for item in fsupp]
    nsupp,ncoe = fsupp[ind],fcoe[ind]
    for (i,item) in enumerate(nsupp)
        Locb = bfind(tsupp, ltsupp, item)
        if Locb === nothing
            @error "The monomial basis is not enough!"
        else
            rcons[Locb] -= real(ncoe[i])
            if ipart == true && item[1] != item[2]
                Locb = bfind(itsupp, length(itsupp), item)
                icons[Locb] -= imag(ncoe[i])
            end
        end
    end
    @constraint(model, rcons .== 0)
    if ipart == true
        @constraint(model, icons .== 0)
    end
    info = struct_data(cql, cliquesize, cliques, basis, ebasis, blocksize, blocks, eblocks, tsupp, I, J, pos, free, nothing)
    return info
end

function assign_constraint(supp::Vector{Vector{Vector{Vector{UInt16}}}}, cliques, cql)
    I = [Int[] for i=1:cql]
    for (i, g) in enumerate(supp)
        ind = findall(k -> issubset(unique(reduce(vcat, [item[1] for item in g])), cliques[k]), 1:cql)
        push!.(I[ind], i)
    end
    return I
end

function get_csupp(order::Int, basis, ebasis, gsupp, hsupp, cql, I, J, blocks, eblocks, cl, blocksize, cliquesize; ConjugateBasis=false, normality=1)
    tsupp = Vector{Vector{UInt16}}[]
    for i = 1:cql
        if ConjugateBasis == false
            a = normality > 0 ? 1 + cliquesize[i] : 1
        else
            a = normality >= order ? 1 + cliquesize[i] : 1
        end
        for (j, w) in enumerate(I[i]), l = 1:cl[i][j+a], t = 1:blocksize[i][j+a][l], r = t:blocksize[i][j+a][l], item in gsupp[w]
            ind1 = blocks[i][j+a][l][t]
            ind2 = blocks[i][j+a][l][r]
            if ConjugateBasis == false
                @inbounds bi = [sadd(basis[i][j+a][ind1], item[1]), sadd(basis[i][j+a][ind2], item[2])]
            else
                @inbounds bi = [sadd(sadd(basis[i][j+a][ind1][1], item[1]), basis[i][j+a][ind2][2]), sadd(sadd(basis[i][j+a][ind1][2], item[2]), basis[i][j+a][ind2][1])]
            end
            bi[1] <= bi[2] ? push!(tsupp, bi) : push!(tsupp, bi[2:-1:1])
        end
        for (j, w) in enumerate(J[i]), k in eblocks[i][j], item in hsupp[w]
            @inbounds bi = [sadd(ebasis[i][j][k][1], item[1]), sadd(ebasis[i][j][k][2], item[2])]
            bi[1] <= bi[2] ? push!(tsupp, bi) : push!(tsupp, bi[2:-1:1])
        end
    end
    return tsupp
end

function clique_decomp(n, m, l, fsupp::Matrix{UInt8}, gsupp::Vector{Matrix{UInt8}}, hsupp::Vector{Matrix{UInt8}}; alg="MF", QUIET=false)
    G = SimpleGraph(n)
    for item in eachcol(fsupp)
        add_clique!(G, findall(item .!= 0))
    end
    for i = 1:m
        temp = findall(gsupp[i][:,1] .!= 0)
        for j = 2:size(gsupp[i], 2)
            append!(temp, findall(gsupp[i][:,j] .!= 0))
        end
        add_clique!(G, unique(temp))
    end
    for i = 1:l
        temp = findall(hsupp[i][:,1] .!= 0)
        for j = 2:size(hsupp[i], 2)
            append!(temp, findall(hsupp[i][:,j] .!= 0))
        end
        add_clique!(G, unique(temp))
    end
    if alg == "NC"
        cliques,cql,cliquesize = max_cliques(G)
    else
        cliques,cql,cliquesize = chordal_cliques!(G, method=alg, minimize=true)
    end
    uc = unique(cliquesize)
    sizes = [sum(cliquesize.== i) for i in uc]
    if QUIET == false
        println("-----------------------------------------------------------------------------")
        println("The clique sizes of varibles:\n$uc\n$sizes")
        println("-----------------------------------------------------------------------------")
    end
    return cliques,cql,cliquesize
end

function assign_constraint(m, l, gsupp::Vector{Matrix{UInt8}}, hsupp::Vector{Matrix{UInt8}}, cliques, cql)
    I = [Int[] for i=1:cql]
    J = [Int[] for i=1:cql]
    for i = 1:m
        rind = findall(gsupp[i][:,1] .!= 0)
        for j = 2:size(gsupp[i], 2)
            append!(rind, findall(gsupp[i][:,j] .!= 0))
        end
        unique!(rind)
        ind = findall(k->issubset(rind, cliques[k]), 1:cql)
        push!.(I[ind], i)
    end
    for i = 1:l
        rind = findall(hsupp[i][:,1] .!= 0)
        for j = 2:size(hsupp[i], 2)
            append!(rind, findall(hsupp[i][:,j] .!= 0))
        end
        unique!(rind)
        ind = findall(k->issubset(rind, cliques[k]), 1:cql)
        push!.(J[ind], i)
    end
    return I,J
end

function get_blocks(n, I, J, fsupp::Matrix{UInt8}, gsupp::Vector{Matrix{UInt8}}, hsupp::Vector{Matrix{UInt8}}, basis, cliques, cql; TS="block", SO=1, signsymmetry=nothing)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
    status = ones(Int, cql)
    tsupp = hcat(fsupp, gsupp..., hsupp...)
    tsupp = sortslices(tsupp, dims=2)
    tsupp = unique(tsupp, dims=2)
    for i = 1:cql
        supp = nothing
        if TS != false
            ind = [issubset(findall(item .!= 0), cliques[i]) for item in eachcol(tsupp)]
            supp = [tsupp[:, ind] UInt8(2)*basis[i][1]]
            supp = sortslices(supp, dims=2)
            supp = unique(supp, dims=2)
        end
        blocks[i],cl[i],blocksize[i],eblocks[i],status[i] = get_blocks(n, length(I[i]), length(J[i]), supp, [gsupp[I[i]]; hsupp[J[i]]], basis[i], TS=TS, SO=SO, signsymmetry=signsymmetry)
    end
    if minimum(status) == 1
        println("No higher TS step of the CS-TSSOS hierarchy!")
    end
    return blocks,cl,blocksize,eblocks
end

function get_blocks(order, I, J, fsupp::Vector{Vector{Vector{UInt16}}}, gsupp::Vector{Vector{Vector{Vector{UInt16}}}}, hsupp::Vector{Vector{Vector{Vector{UInt16}}}}, basis, ebasis, cliques, 
    cliquesize, cql; TS="block", SO=1, ConjugateBasis=false, normality=1)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
    status = ones(Int, cql)
    tsupp = filter(item -> item[1] <= item[2], vcat([[fsupp]; gsupp; hsupp]...))
    sort!(tsupp)
    unique!(tsupp)
    for i = 1:cql
        ksupp = TS == false ? nothing : tsupp[[issubset(union(item[1], item[2]), cliques[i]) for item in tsupp]]
        blocks[i],cl[i],blocksize[i],eblocks[i],status[i] = get_blocks(order, length(I[i]), length(J[i]), ksupp, gsupp[I[i]], hsupp[J[i]], basis[i], ebasis[i], 
        TS=TS, SO=SO, ConjugateBasis=ConjugateBasis, normality=normality, nvar=cliquesize[i])
    end
    if minimum(status) == 1
        println("No higher TS step of the CS-TSSOS hierarchy!")
    end
    return blocks,cl,blocksize,eblocks
end

function get_blocks(n::Int, m::Int, l::Int, tsupp, supp::Vector{Matrix{UInt8}}, basis::Vector{Matrix{UInt8}}; TS="block", SO=1, merge=false, md=3, signsymmetry=nothing)
    blocks = Vector{Vector{Vector{Int}}}(undef, m+1)
    blocksize = Vector{Vector{Int}}(undef, m+1)
    cl = Vector{Int}(undef, m+1)
    eblocks = Vector{Vector{Int}}(undef, l)
    status = 0
    if TS == false
        for k = 1:m+1
            blocks[k],blocksize[k],cl[k] = [Vector(1:size(basis[k],2))],[size(basis[k],2)],1       
        end
        for k = 1:l
            eblocks[k] = Vector(1:size(basis[k+m+1],2))
        end
    else       
        for i = 1:SO
            if i > 1
                oblocksize = deepcopy(blocksize)
                oeblocks = deepcopy(eblocks)
            end
            for k = 1:m+1
                if k == 1
                    G = get_graph(tsupp, basis[1], signsymmetry=signsymmetry)
                else
                    G = get_graph(tsupp, supp[k-1], basis[k], signsymmetry=signsymmetry)
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
                eblocks[k] = get_eblock(tsupp, supp[k+m], basis[k+m+1], signsymmetry=signsymmetry)
            end
            if i > 1 && blocksize == oblocksize && eblocks == oeblocks
                status = 1
                break
            end
            if i < SO
                tsupp = zeros(UInt8, n, numele(blocksize[1]))
                k = 1
                for t = 1:length(blocks[1]), j = 1:blocksize[1][t], r = j:blocksize[1][t]
                    tsupp[:,k] = basis[1][:,blocks[1][t][j]] + basis[1][:,blocks[1][t][r]]
                    k += 1
                end
                tsupp = sortslices(tsupp, dims=2)
                tsupp = unique(tsupp, dims=2)
            end
        end
    end
    return blocks,cl,blocksize,eblocks,status
end

function get_blocks(order, m::Int, l::Int, tsupp, gsupp::Vector{Vector{Vector{Vector{UInt16}}}}, hsupp::Vector{Vector{Vector{Vector{UInt16}}}}, basis, ebasis; TS="block", SO=1, ConjugateBasis=false, normality=1, nvar=0)
    if (ConjugateBasis == false && normality > 0) || (ConjugateBasis == true && normality >= order)
        uk = m + 1 + nvar
    else
        uk = m + 1
    end
    blocks = Vector{Vector{Vector{Int}}}(undef, uk)
    blocksize = Vector{Vector{Int}}(undef, uk)
    cl = Vector{Int}(undef, uk)
    eblocks = Vector{Vector{Int}}(undef, l)
    status = 0
    if TS == false
        for k = 1:uk
            lb = length(basis[k])
            blocks[k],blocksize[k],cl[k] = [Vector(1:lb)],[lb],1
        end
        for k = 1:l
            eblocks[k] = Vector(1:length(ebasis[k]))
        end
    else  
        for i = 1:SO
            if i > 1
                oblocksize = deepcopy(blocksize)
                oeblocks = deepcopy(eblocks)
            end
            for k = 1:uk
                if k == 1
                    G = get_graph(tsupp, basis[1], ConjugateBasis=ConjugateBasis)
                elseif ((ConjugateBasis == false && normality > 0) || (ConjugateBasis == true && normality >= order)) && k <= 1 + nvar
                    G = get_graph(tsupp, basis[k], ConjugateBasis=true)
                elseif ((ConjugateBasis == false && normality > 0) || (ConjugateBasis == true && normality >= order)) && k > 1 + nvar
                    G = get_graph(tsupp, gsupp[k-1-nvar], basis[k], ConjugateBasis=ConjugateBasis)
                else
                    G = get_graph(tsupp, gsupp[k-1], basis[k], ConjugateBasis=ConjugateBasis)
                end
                if TS == "block"
                    blocks[k] = connected_components(G)
                    blocksize[k] = length.(blocks[k])
                    cl[k] = length(blocksize[k])
                else
                    blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS)
                end
            end
            for k = 1:l
                 eblocks[k] = get_eblock(tsupp, hsupp[k], ebasis[k])
            end
            if i > 1 && blocksize == oblocksize && eblocks == oeblocks
                status = 1
                break
            end
            if i < SO
                tsupp = Vector{Vector{UInt16}}[]
                if ConjugateBasis == false
                    a = normality > 0 ? 1 + nvar : 1
                else
                    a = normality >= order ? 1 + nvar : 1
                end
                for s = 1:a, j = 1:cl[s], k = 1:blocksize[s][j], r = k:blocksize[s][j]
                    if ConjugateBasis == false && s == 1
                        @inbounds bi = [basis[s][blocks[s][j][k]], basis[s][blocks[s][j][r]]]
                    else
                        @inbounds bi = [sadd(basis[s][blocks[s][j][k]][1], basis[s][blocks[s][j][r]][2]), sadd(basis[s][blocks[s][j][k]][2], basis[s][blocks[s][j][r]][1])]
                    end
                    bi[1] <= bi[2] ? push!(tsupp, bi) : push!(tsupp, bi[2:-1:1])
                end
            end
            csupp = get_csupp(order, [basis], [ebasis], gsupp, hsupp, 1, Vector(1:m), Vector(1:l), [blocks], [eblocks], [cl], [blocksize], [nvar], ConjugateBasis=ConjugateBasis, normality=normality)
            append!(tsupp, csupp)
            sort!(tsupp)
            unique!(tsupp)
        end
    end
    return blocks,cl,blocksize,eblocks,status
end

function get_moment(supp::Array{UInt8, 2}, lb, ub)
    return [prod([(ub[j]^(item[j]+1)-lb[j]^(item[j]+1))/(item[j]+1) for j=1:length(lb)]) for item in eachcol(supp)]
end

function get_moment(mons::Vector{DP.Monomial{V, M}}, lb, ub) where {V, M}
    supp = exponents.(mons)
    return [prod([(ub[j]^(item[j]+1)-lb[j]^(item[j]+1))/(item[j]+1) for j=1:length(lb)]) for item in supp]
end

function get_moment_matrix(moment, info)
    MomMat = Vector{Union{Float64, Symmetric{Float64}, Array{Float64,2}}}(undef, info.cql)
    ltsupp = size(info.tsupp, 2)
    for i = 1:info.cql
        lb = size(info.basis[i][1], 2)
        MomMat[i] = zeros(Float64, lb, lb)
        for j = 1:lb, k = j:lb
            bi = info.basis[i][1][:, j] + info.basis[i][1][:, k]
            Locb = bfind(info.tsupp, ltsupp, bi)
            if Locb !== nothing
                MomMat[i][j,k] = moment[Locb]
            end
        end
        MomMat[i] = Symmetric(MomMat[i], :U)
    end
    return MomMat
end
