mutable struct struct_data
    cql # number of cliques
    cliquesize # size of cliques
    cliques # clique structrue
    basis # monomial basis
    cl # number of blocks
    blocksize # size of blocks
    blocks # block structrue for inequality constraints
    eblocks # block structrue for equality constraints
    tsupp # total support
    I # index sets of inequality constraints
    J # index sets of equality constraints
    gram # Gram variables
    constrs # constraint name
end

"""
model,info = add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order; CS=false, cliques=[], TS="block", 
SO=1, Groebnerbasis=false, QUIET=false, constrs=nothing)

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
- `TS`: type of term sparsity (`"block"`, `"MD"`, `"MF"`, `false`)
- `SO`: sparse order
- `Groebnerbasis`: exploit the quotient ring structure or not (`true`, `false`)
- `QUIET`: run in the quiet mode (`true`, `false`)
- `constrs`: the constraint name used in the JuMP model

# Output arguments
- `model`: the modified JuMP model
- `info`: other auxiliary data
"""
function add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order; CS=false, cliques=[], TS="block", SO=1, Groebnerbasis=false, QUIET=false, constrs=nothing)
    n = length(vars)
    m = length(ineq_cons)
    if ineq_cons != []
        gsupp,gcoe,glt,dg = npolys_info(ineq_cons, vars)
    else
        gsupp = Matrix{UInt8}[]
        glt = dg = Int[]
    end
    if eq_cons != []
        hsupp,hcoe,hlt,dh = npolys_info(eq_cons, vars)
    else
        hsupp = Matrix{UInt8}[]
        hlt = dh = Int[]
    end
    if Groebnerbasis == true && eq_cons != []
        l = 0
        gb = convert.(Polynomial{true,Float64}, eq_cons)
        SemialgebraicSets.gröbnerbasis!(gb)
        nonneg = rem(nonneg, gb)
        leadm = leadingmonomial.(gb)
        llead = length(leadm)
        lead = zeros(UInt8, n, llead)
        for i = 1:llead, j = 1:n
            @inbounds lead[j,i] = MultivariatePolynomials.degree(leadm[i], vars[j])
        end
    else
        gb = []
        l = length(eq_cons)
    end
    fsupp,fcoe = poly_info(nonneg, vars)
    if CS != false
        if cliques == []
            cliques,cql,cliquesize = clique_decomp(n, m, length(eq_cons), fsupp, gsupp, hsupp, alg=CS, QUIET=false)
        else
            cql = length(cliques)
            cliquesize = length.(cliques)
        end
    else
        cliques,cql,cliquesize = [Vector(1:n)],1,[n]
    end
    dmin = ceil(Int, maximum([maxdegree(nonneg); dg; dh])/2)
    order = order < dmin ? dmin : order
    I,J = assign_constraint(m, l, gsupp, hsupp, cliques, cql)
    basis = Vector{Vector{Matrix{UInt8}}}(undef, cql)
    for t = 1:cql
        basis[t] = Vector{Matrix{UInt8}}(undef, length(I[t])+length(J[t])+1)
        basis[t][1] = get_nbasis(n, order, var=cliques[t])
        for s = 1:length(I[t])
            basis[t][s+1] = get_nbasis(n, order-ceil(Int, dg[I[t][s]]/2), var=cliques[t])
        end
        for s = 1:length(J[t])
            basis[t][s+length(I[t])+1] = get_nbasis(n, 2*order-dh[J[t][s]], var=cliques[t])
        end
    end
    blocks,cl,blocksize,eblocks,_,_,_ = get_cblocks_mix(n, I, J, m, l, fsupp, gsupp, glt, hsupp, hlt, basis, cliques, cql, tsupp=[], TS=TS, SO=SO, QUIET=QUIET)
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
            return model,info
        else
            bc[Locb] = fcoe[i]
        end
    end
    if constrs !== nothing
        @constraint(model, [i=1:ltsupp], cons[i]==bc[i], base_name=constrs)
    else
        @constraint(model, cons.==bc)
    end
    info = struct_data(cql,cliquesize,cliques,basis,cl,blocksize,blocks,eblocks,tsupp,I,J,pos,constrs)
    return model,info
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
    I = [UInt32[] for i=1:cql]
    J = [UInt32[] for i=1:cql]
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

function get_cblocks_mix(n, I, J, m, l, fsupp::Matrix{UInt8}, gsupp::Vector{Matrix{UInt8}}, glt, hsupp::Vector{Matrix{UInt8}}, hlt, basis, cliques, cql; tsupp=[], TS="block", SO=1, QUIET=false)
    blocks = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
    eblocks = Vector{Vector{Vector{UInt16}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    sb = Vector{Vector{Int}}(undef, cql)
    numb = Vector{Vector{Int}}(undef, cql)
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
    status = ones(Int, cql)
    for i = 1:cql
        lc = length(I[i])
        ind = [issubset(findall(item .!= 0), cliques[i]) for item in eachcol(tsupp)]
        supp = [tsupp[:, ind] UInt8(2)*basis[i][1]]
        supp = sortslices(supp, dims=2)
        supp = unique(supp, dims=2)
        blocks[i] = Vector{Vector{Vector{UInt16}}}(undef, lc+1)
        eblocks[i] = Vector{Vector{UInt16}}(undef, length(J[i]))
        cl[i] = Vector{Int}(undef, lc+1)
        blocksize[i] = Vector{Vector{Int}}(undef, lc+1)
        sb[i] = Vector{Int}(undef, lc+1)
        numb[i] = Vector{Int}(undef, lc+1)
        blocks[i],cl[i],blocksize[i],eblocks[i],sb[i],numb[i],status[i] = get_blocks(n, lc, length(J[i]), supp, [gsupp[I[i]]; hsupp[J[i]]], [glt[I[i]]; hlt[J[i]]], basis[i], TS=TS, SO=SO, QUIET=QUIET)
    end
    return blocks,cl,blocksize,eblocks,sb,numb,maximum(status)
end

function get_cgraph(tsupp::Array{UInt8, 2}, gsupp::Array{UInt8, 2}, glt, basis::Array{UInt8, 2})
    lb = size(basis, 2)
    G = SimpleGraph(lb)
    ltsupp = size(tsupp, 2)
    for i = 1:lb, j = i+1:lb
        ind = findfirst(x -> bfind(tsupp, ltsupp, basis[:,i] + basis[:,j] + gsupp[:,x]) !== nothing, 1:glt)
        if ind !== nothing
           add_edge!(G, i, j)
        end
    end
    return G
end

function get_eblock(tsupp::Array{UInt8, 2}, hsupp::Array{UInt8, 2}, hlt, basis::Array{UInt8, 2})
    ltsupp = size(tsupp, 2)
    eblock = UInt16[]
    for i = 1:size(basis, 2)
        ind = findfirst(x -> bfind(tsupp, ltsupp, basis[:,i] + hsupp[:,x]) !== nothing, 1:hlt)
        if ind !== nothing
           push!(eblock, i)
        end
    end
    return eblock
end

function get_blocks(n::Int, m::Int, l::Int, tsupp, supp::Vector{Array{UInt8, 2}}, lt, basis::Vector{Array{UInt8, 2}}; TS="block", SO=1, merge=false, md=3, QUIET=false)
    blocks = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    eblocks = Vector{Vector{UInt16}}(undef, l)
    blocksize = Vector{Vector{Int}}(undef, m+1)
    cl = Vector{Int}(undef, m+1)
    if TS == false
        for k = 1:m+1
            blocks[k] = [[i for i=1:size(basis[k],2)]]
            blocksize[k] = [size(basis[k],2)]
            cl[k] = 1          
        end
        for k = 1:l
            eblocks[k] = [i for i=1:size(basis[k+m+1],2)]
        end
        sb = blocksize[1]
        numb = [1]
        status = 1
    else
        status = 1
        blocks[1] = Vector{UInt16}[]
        for i = 1:SO
            G = get_graph(tsupp, basis[1])
            if TS == "block"
                nblock = connected_components(G)
            else
                nblock = chordal_cliques!(G, method=TS)[1]
                if merge == true
                    nblock = clique_merge!(nblock, QUIET=true, d=md)[1]
                end
            end
            if nblock != blocks[1]
                blocks[1] = nblock
                if i < SO
                    blocksize[1] = length.(blocks[1])
                    tsupp = zeros(UInt8, n, numele(blocksize[1]))
                    k = 1
                    for i = 1:length(blocks[1]), j = 1:blocksize[1][i], r = j:blocksize[1][i]
                        tsupp[:,k] = basis[1][:,blocks[1][i][j]] + basis[1][:,blocks[1][i][r]]
                        k += 1
                    end
                    tsupp = sortslices(tsupp, dims=2)
                    tsupp = unique(tsupp, dims=2)
                end
            else
                println("No higher TS step of the TSSOS hierarchy!")
                status = 0
                sb = numb = nothing
                break
            end
        end
        if status == 1
            blocksize[1] = length.(blocks[1])
            cl[1] = length(blocksize[1])
            bz = sort(blocksize[1], rev=true)
            sb = unique(bz)
            numb = [sum(bz.== i) for i in sb]
            if QUIET == false
                println("------------------------------------------------------")
                println("The sizes of PSD blocks:\n$sb\n$numb")
                println("------------------------------------------------------")
            end
            for k = 1:m
                G = get_cgraph(tsupp, supp[k], lt[k], basis[k+1])
                blocks[k+1] = connected_components(G)
                blocksize[k+1] = length.(blocks[k+1])
                cl[k+1] = length(blocksize[k+1])
            end
            for k = 1:l
                eblocks[k] = get_eblock(tsupp, supp[k+m], lt[k+m], basis[k+m+1])
            end
        end
    end
    return blocks,cl,blocksize,eblocks,sb,numb,status
end

function get_moment(n, tsupp, lb, ub)
    ltsupp = size(tsupp, 2)
    moment = zeros(ltsupp)
    for i = 1:ltsupp
        moment[i] = prod([(ub[j]^(tsupp[j,i]+1)-lb[j]^(tsupp[j,i]+1))/(tsupp[j,i]+1) for j=1:n])
    end
    return moment
end

function get_moment_matrix(moment, tsupp, cql, basis)
    MomMat = Vector{Union{Float64, Symmetric{Float64}, Array{Float64,2}}}(undef, cql)
    ltsupp = size(tsupp, 2)
    for i = 1:cql
        lb = size(basis[i][1], 2)
        MomMat[i] = zeros(Float64, lb, lb)
        for j = 1:lb, k = j:lb
            bi = basis[i][1][:, j] + basis[i][1][:, k]
            Locb = bfind(tsupp, ltsupp, bi)
            if Locb !== nothing
                MomMat[i][j,k] = moment[Locb]
            end
        end
        MomMat[i] = Symmetric(MomMat[i], :U)
    end
    return MomMat
end
