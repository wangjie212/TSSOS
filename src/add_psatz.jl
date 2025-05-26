mutable struct struct_data
    cql # number of cliques
    cliquesize # size of cliques
    cliques # clique structrue
    basis # monomial basis
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
function add_psatz!(model, nonneg::DP.Polynomial{V, M, T}, vars, ineq_cons, eq_cons, order; CS=false, cliques=[], blocks=[], TS="block", SO=1, GroebnerBasis=false, QUIET=false, constrs=nothing) where {V, M, T<:Union{Number,AffExpr}}
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
        dh = maxdegree.(ineq_cons)
    else
        hsupp = Matrix{UInt8}[]
        hlt = dh = Int[]
    end
    if GroebnerBasis == true && eq_cons != []
        l = 0
        gb = convert.(DP.Polynomial{V, M, Float64}, eq_cons)
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
        blocks,cl,blocksize,eblocks = get_blocks(n, I, J, m, l, fsupp, gsupp, hsupp, basis, cliques, cql, tsupp=[], TS=TS, SO=SO, signsymmetry=ss)
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
    bc = [AffExpr(0) for i=1:ltsupp]
    for i = 1:size(fsupp, 2)
        Locb = bfind(tsupp, ltsupp, fsupp[:, i])
        if Locb === nothing
            @error "The monomial basis is not enough!"
        else
            bc[Locb] = fcoe[i]
        end
    end
    if constrs !== nothing
        @constraint(model, cons==bc, base_name=constrs)
    else
        @constraint(model, cons==bc)
    end
    info = struct_data(cql, cliquesize, cliques, basis, blocksize, blocks, eblocks, tsupp, I, J, pos, mul, constrs)
    return info
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

function get_blocks(n, I, J, m, l, fsupp::Matrix{UInt8}, gsupp::Vector{Matrix{UInt8}}, hsupp::Vector{Matrix{UInt8}}, basis, cliques, cql; tsupp=[], TS="block", SO=1, signsymmetry=nothing)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    status = ones(Int, cql)
    if isempty(tsupp)
        tsupp = copy(fsupp)
        for i = 1:m
            tsupp = [tsupp gsupp[i]]
        end
        for i = 1:l
            tsupp = [tsupp hsupp[i]]
        end
        tsupp = sortslices(tsupp, dims=2)
        tsupp = unique(tsupp, dims=2)
    end
    for i = 1:cql
        lc = length(I[i])
        ind = [issubset(findall(item .!= 0), cliques[i]) for item in eachcol(tsupp)]
        supp = [tsupp[:, ind] UInt8(2)*basis[i][1]]
        supp = sortslices(supp, dims=2)
        supp = unique(supp, dims=2)
        blocks[i] = Vector{Vector{Vector{Int}}}(undef, lc+1)
        eblocks[i] = Vector{Vector{Int}}(undef, length(J[i]))
        cl[i] = Vector{Int}(undef, lc+1)
        blocksize[i] = Vector{Vector{Int}}(undef, lc+1)
        blocks[i],cl[i],blocksize[i],eblocks[i],status[i] = get_blocks(n, lc, length(J[i]), supp, [gsupp[I[i]]; hsupp[J[i]]], basis[i], TS=TS, SO=SO, signsymmetry=signsymmetry)
    end
    if minimum(status) == 1
        println("No higher TS step of the CS-TSSOS hierarchy!")
    end
    return blocks,cl,blocksize,eblocks
end

function get_blocks(n::Int, m::Int, l::Int, tsupp, supp::Vector{Array{UInt8, 2}}, basis::Vector{Array{UInt8, 2}}; TS="block", SO=1, merge=false, md=3, signsymmetry=nothing)
    blocks = Vector{Vector{Vector{Int}}}(undef, m+1)
    eblocks = Vector{Vector{Int}}(undef, l)
    blocksize = Vector{Vector{Int}}(undef, m+1)
    cl = Vector{Int}(undef, m+1)
    status = 0
    if TS == false
        for k = 1:m+1
            blocks[k],blocksize[k],cl[k] = [[i for i=1:size(basis[k],2)]],[size(basis[k],2)],1       
        end
        for k = 1:l
            eblocks[k] = [i for i=1:size(basis[k+m+1],2)]
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
