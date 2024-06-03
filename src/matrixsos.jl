mutable struct poly_data
    n::Int
    supp::Vector{Vector{UInt16}}
    coe::Vector{Number}
end

mutable struct poly_matrix
    m::Int
    poly::Vector{poly_data} # store the upper triangular part by colomn
end

mutable struct mpop_data
    b
    obj_matrix
    cons_matrix
    basis # monomial basis
    gbasis
    ksupp # extended support at the k-th step
    cl
    blocksize
    blocks # block structrue
    cql
    cliquesize
    cliques
    I
    sb # sizes of different blocks
    numb # numbers of different blocks
    SDP_status
end

function tssos_first(F::Matrix{Polynomial{true, T}}, G, x, d; TS="block", QUIET=false, solve=true) where {T<:Number}
    return cs_tssos_first(F, G, x, d, CS=false, TS=TS, QUIET=QUIET, solve=solve)
end

function tssos_higher!(data::mpop_data; TS="block", QUIET=false, solve=true)
    return cs_tssos_higher!(data, TS=TS, QUIET=QUIET, solve=solve)
end

function cs_tssos_first(F::Matrix{Polynomial{true, T}}, G, x, d; CS="MF", TS="block", QUIET=false, solve=true) where {T<:Number}
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    n = length(x)
    m = length(G)
    dG = [maximum(maxdegree.(vec(G[i]))) for i=1:m]
    obj_matrix = poly_matrix(size(F,1), Vector{poly_data}(undef, Int((size(F,1)+1)*size(F,1)/2)))
    for i = 1:obj_matrix.m, j = i:obj_matrix.m
        _,supp,coe = polys_info([F[i,j]], x)
        obj_matrix.poly[i+Int(j*(j-1)/2)] = poly_data(n, supp[1], coe[1])
    end
    cons_matrix = Vector{poly_matrix}(undef, m)
    for k = 1:m
        cons_matrix[k] = poly_matrix(size(G[k],1), Vector{poly_data}(undef, Int((size(G[k],1)+1)*size(G[k],1)/2)))
        for i = 1:cons_matrix[k].m, j = i:cons_matrix[k].m
            _,supp,coe = polys_info([G[k][i,j]], x)
            cons_matrix[k].poly[i+Int(j*(j-1)/2)] = poly_data(n, supp[1], coe[1])
        end
    end
    cliques,cql,cliquesize = clique_decomp(n, m, d, dG, obj_matrix, cons_matrix, alg=CS, minimize=true)
    I,ncc = assign_constraint(m, cons_matrix, cliques, cql)
    basis = Vector{Vector{Vector{UInt16}}}(undef, cql)
    gbasis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
    for i = 1:cql
        basis[i] = get_sbasis(cliques[i], d)
        gbasis[i] = Vector{Vector{Vector{UInt16}}}(undef, length(I[i]))
        for (s,k) in enumerate(I[i])
            gbasis[i][s] = get_sbasis(cliques[i], d-Int(ceil(dG[k]/2)))
        end
    end
    ksupp = Vector{Vector{Vector{UInt16}}}(undef, Int((obj_matrix.m+1)*obj_matrix.m/2))
    if TS != false
        for i = 1:obj_matrix.m, j = i:obj_matrix.m
            ind = i + Int(j*(j-1)/2)
            ksupp[ind] = copy(obj_matrix.poly[ind].supp)
            if i == j
                for k = 1:cql, item in basis[k]
                    push!(ksupp[ind], sadd(item, item))
                end
            end
        end
        unique!.(ksupp)
        sort!.(ksupp)
        if QUIET == false
            println("Starting to compute the block structure...")
        end
    end
    time = @elapsed begin
    blocks,cl,blocksize,sb,numb,_ = get_cmblocks_mix(I, obj_matrix.m, cons_matrix, cliques, cql, ksupp, basis, gbasis, blocks=[], cl=[], blocksize=[], sb=[], numb=[], TS=TS)
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,SDP_status = pmo_sdp(obj_matrix, cons_matrix, basis, gbasis, blocks, cl, blocksize, cql, I, QUIET=QUIET, solve=solve)
    data = mpop_data(nothing, obj_matrix, cons_matrix, basis, gbasis, ksupp, cl, blocksize, blocks, cql, cliquesize, cliques, I, sb, numb, SDP_status)
    return opt,data
end

function cs_tssos_higher!(data::mpop_data; TS="block", QUIET=false, solve=true)
    basis = data.basis
    gbasis = data.gbasis
    ksupp = data.ksupp
    obj_matrix = data.obj_matrix
    cons_matrix = data.cons_matrix
    blocks = data.blocks
    cl = data.cl
    blocksize = data.blocksize
    sb = data.sb
    numb = data.numb
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,data.sb,data.numb,status = get_cmblocks_mix(data.I, obj_matrix.m, cons_matrix, data.cliques, data.cql, ksupp, basis, gbasis, blocks=blocks, 
    cl=cl, blocksize=blocksize, sb=sb, numb=numb, TS=TS, QUIET=QUIET)
    end
    opt = nothing
    if status == 1
        if TS != false && QUIET == false
            mb = maximum(maximum.(data.sb))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,SDP_status = pmo_sdp(obj_matrix, cons_matrix, basis, gbasis, blocks, cl, blocksize, data.cql, data.I, QUIET=QUIET, solve=solve)
        data.ksupp = ksupp
        data.blocks = blocks
        data.cl = cl
        data.blocksize = blocksize
        data.SDP_status = SDP_status
    end
    return opt,data
end

function clique_decomp(n, m, d, dG, obj_matrix, cons_matrix; alg="MF", minimize=false)
    if alg == false
        cliques,cql,cliquesize = [UInt16[i for i=1:n]],1,[n]
    else
        G = SimpleGraph(n)
        for i = 1:Int((obj_matrix.m + 1)*obj_matrix.m/2)
            foreach(x -> add_clique!(G, unique(x)), obj_matrix.poly[i].supp)
        end
        for k = 1:m
            if d == ceil(Int, dG[k]/2)
                for i = 1:Int((cons_matrix[k].m + 1)*cons_matrix[k].m/2)
                    foreach(x -> add_clique!(G, unique(x)), cons_matrix[k].poly[i].supp)
                end
            else
                for i = 1:Int((cons_matrix[k].m + 1)*cons_matrix[k].m/2)
                    add_clique!(G, unique(reduce(vcat, cons_matrix[k].poly[i].supp)))
                end
            end
        end
        if alg == "NC"
            cliques,cql,cliquesize = max_cliques(G)
        else
            cliques,cql,cliquesize = chordal_cliques!(G, method=alg, minimize=minimize)
        end
    end
    uc = unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("-----------------------------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("-----------------------------------------------------------------------------")
    return cliques,cql,cliquesize
end

function assign_constraint(m, cons_matrix, cliques, cql)
    I = [UInt16[] for i=1:cql]
    ncc = UInt16[]
    for i = 1:m
        ind = findall(k->issubset(unique(vcat([reduce(vcat, cons_matrix[i].poly[s].supp) for s=1:Int((cons_matrix[i].m + 1)*cons_matrix[i].m/2)]...)), cliques[k]), 1:cql)
        isempty(ind) ? push!(ncc, i) : push!.(I[ind], i)
    end
    return I,ncc
end

function get_mgraph(tsupp, basis, om)
    lb = length(basis)
    G = SimpleGraph(lb*om)
    for i = 1:om, j = i:om
        lt = length(tsupp[i+Int(j*(j-1)/2)])
        for k = 1:lb, l = k:lb
            bi = sadd(basis[k], basis[l])
            if bfind(tsupp[i+Int(j*(j-1)/2)], lt, bi) !== nothing
               add_edge!(G, i+(k-1)*om, j+(l-1)*om)
            end
        end
    end
    return G
end

function get_cmgraph(tsupp, cons_matrix, gbasis, om)
    com = cons_matrix.m*om
    lb = length(gbasis)*com
    G = SimpleGraph(lb)
    for j = 1:lb, k = j+1:lb
        p = mod(j, com)
        p = p != 0 ? p : com
        q = mod(k, com)
        q = q != 0 ? q : com
        p1 = ceil(Int, p/cons_matrix.m)
        q1 = ceil(Int, q/cons_matrix.m)
        ind = p1 <= q1 ? p1 + Int(q1*(q1-1)/2) : q1 + Int(p1*(p1-1)/2)
        t = mod(j, cons_matrix.m)
        t = t != 0 ? t : cons_matrix.m
        r = mod(k, cons_matrix.m)
        r = r != 0 ? r : cons_matrix.m
        loc = t <= r ? t + Int(r*(r-1)/2) : r + Int(t*(t-1)/2)
        flag = 0
        for w = 1:length(cons_matrix.poly[loc].supp)
            if bfind(tsupp[ind], length(tsupp[ind]), sadd(sadd(gbasis[ceil(Int, j/com)], gbasis[ceil(Int, k/com)]), cons_matrix.poly[loc].supp[w])) !== nothing
                flag = 1
                break
            end
        end
        if flag == 1
           add_edge!(G, j, k)
        end
    end
    return G
end

function get_mblocks(om, cons_matrix, tsupp, basis, gbasis; TS="block", blocks=[], cl=[], blocksize=[], sb=[], numb=[], QUIET=false)
    if isempty(blocks)
        blocks = Vector{Vector{Vector{UInt16}}}(undef, length(cons_matrix)+1)
        blocksize = Vector{Vector{UInt16}}(undef, length(cons_matrix)+1)
        cl = Vector{UInt16}(undef, length(cons_matrix)+1)
    end
    if TS == false
        for i = 1:length(cons_matrix) + 1
            lb = i == 1 ? om*length(basis) : om*cons_matrix[i-1].m*length(gbasis[i-1])
            blocks[i] = [Vector(1:lb)]
            blocksize[i] = [lb]
            cl[i] = 1
        end
        status = 1
        nsb = [om*length(basis)]
        nnumb = [1]
    else
        G = get_mgraph(tsupp, basis, om)
        if TS == "block"
            blocks[1] = connected_components(G)
            blocksize[1] = length.(blocks[1])
            cl[1] = length(blocksize[1])
        else
            blocks[1],cl[1],blocksize[1] = chordal_cliques!(G, method=TS, minimize=false)
        end
        nsb = sort(Int.(unique(blocksize[1])), rev=true)
        nnumb = [sum(blocksize[1].== i) for i in nsb]
        if isempty(sb) || nsb!=sb || nnumb!=numb
            status = 1
            for k = 1:length(cons_matrix)
                G = get_cmgraph(tsupp, cons_matrix[k], gbasis[k],om)
                if TS == "block"
                    blocks[k+1] = connected_components(G)
                    blocksize[k+1] = length.(blocks[k+1])
                    cl[k+1] = length(blocksize[k+1])
                else
                    blocks[k+1],cl[k+1],blocksize[k+1] = chordal_cliques!(G, method=TS, minimize=false)
                end
            end
        else
            status = 0
            if QUIET == false
                println("No higher TS step of the TSSOS hierarchy!")
            end
        end
    end
    if status == 1 && QUIET == false
        println("-----------------------------------------------------------------------------")
        println("The sizes of PSD blocks:\n$nsb\n$nnumb")
        println("-----------------------------------------------------------------------------")
    end
    return blocks,cl,blocksize,nsb,nnumb,status
end

function get_cmblocks_mix(I, om, cons_matrix, cliques, cql, tsupp, basis, gbasis; blocks=[], cl=[], blocksize=[], sb=[], numb=[], TS="block", QUIET=true)
    status = ones(Int, cql)
    if isempty(blocks)
        blocks = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        cl = Vector{Vector{UInt16}}(undef, cql)
        blocksize = Vector{Vector{Vector{UInt16}}}(undef, cql)
        sb = Vector{Vector{UInt16}}(undef, cql)
        numb = Vector{Vector{UInt16}}(undef, cql)
        for i = 1:cql
            blocks[i] = Vector{Vector{Vector{UInt16}}}(undef, length(I[i])+1)
            cl[i] = Vector{UInt16}(undef, length(I[i])+1)
            blocksize[i] = Vector{Vector{UInt16}}(undef, length(I[i])+1)
            sb[i] = Vector{UInt16}(undef, length(I[i])+1)
            numb[i] = Vector{UInt16}(undef, length(I[i])+1)
            nsupp = nothing
            if TS != false
                ind = [[issubset(item[j], cliques[i]) for j in eachindex(item)] for item in tsupp]
                nsupp = [tsupp[k][ind[k]] for k = 1:length(tsupp)]
            end
            blocks[i],cl[i],blocksize[i],sb[i],numb[i],status[i] = get_mblocks(om, cons_matrix[I[i]], nsupp, basis[i], gbasis[i], TS=TS, QUIET=true)
        end
    else
        for i = 1:cql
            ind = [[issubset(item[j], cliques[i]) for j in eachindex(item)] for item in tsupp]
            nsupp = [tsupp[k][ind[k]] for k = 1:length(tsupp)]
            blocks[i],cl[i],blocksize[i],sb[i],numb[i],status[i] = get_mblocks(om, cons_matrix[I[i]], nsupp, basis[i], gbasis[i], blocks=blocks[i], 
            cl=cl[i], blocksize=blocksize[i], sb=sb[i], numb=numb[i], TS=TS, QUIET=true)
        end
    end
    return blocks,cl,blocksize,sb,numb,maximum(status)
end

function pmo_sdp(obj_matrix, cons_matrix, basis, gbasis, blocks, cl, blocksize, cql, I; solve=true, QUIET=false)
    om = obj_matrix.m
    ksupp = [Vector{UInt16}[] for i = 1:length(obj_matrix.poly)]
    for u = 1:cql, i = 1:cl[u][1], j = 1:blocksize[u][1][i], k = j:blocksize[u][1][i]
        bi = sadd(basis[u][ceil(Int, blocks[u][1][i][j]/om)], basis[u][ceil(Int, blocks[u][1][i][k]/om)])
        ind1 = mod(blocks[u][1][i][j], om)
        ind1 = ind1 != 0 ? ind1 : om
        ind2 = mod(blocks[u][1][i][k], om)
        ind2 = ind2 != 0 ? ind2 : om
        push!(ksupp[ind1+Int(ind2*(ind2-1)/2)], bi)
    end
    for u = 1:cql, (s,v) = enumerate(I[u])
        com = cons_matrix[v].m*om
        for i = 1:cl[u][s+1], j = 1:blocksize[u][s+1][i], k = j:blocksize[u][s+1][i]
            p = mod(blocks[u][s+1][i][j], com)
            p = p != 0 ? p : com
            q = mod(blocks[u][s+1][i][k], com)
            q = q != 0 ? q : com
            p1 = ceil(Int, p/cons_matrix[v].m)
            q1 = ceil(Int, q/cons_matrix[v].m)
            ind = p1 <= q1 ? p1 + Int(q1*(q1-1)/2) : q1 + Int(p1*(p1-1)/2)
            t = mod(blocks[u][s+1][i][j], cons_matrix[v].m)
            t = t != 0 ? t : cons_matrix[v].m
            r = mod(blocks[u][s+1][i][k], cons_matrix[v].m)
            r = r != 0 ? r : cons_matrix[v].m
            loc = t <= r ? t + Int(r*(r-1)/2) : r + Int(t*(t-1)/2)
            for w = 1:length(cons_matrix[v].poly[loc].supp)
                bi = sadd(sadd(gbasis[u][s][ceil(Int, blocks[u][s+1][i][j]/com)], gbasis[u][s][ceil(Int, blocks[u][s+1][i][k]/com)]), cons_matrix[v].poly[loc].supp[w])
                push!(ksupp[ind], bi)
            end
        end
    end
    sort!.(ksupp)
    unique!.(ksupp)
    objv = SDP_status = nothing
    if solve == true
        if QUIET == false
            ncons = sum(length.(ksupp))
            println("Assembling the SDP...")
            println("There are $ncons affine constraints.")
        end
        model = Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = Vector{Vector{AffExpr}}(undef, length(obj_matrix.poly))
        for i = 1:length(obj_matrix.poly)
            cons[i] = [AffExpr(0) for j=1:length(ksupp[i])]
        end
        for u = 1:cql
            pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[u][1])
            for i = 1:cl[u][1]
                pos[i] = @variable(model, [1:blocksize[u][1][i], 1:blocksize[u][1][i]], PSD)
                for j = 1:blocksize[u][1][i], k = j:blocksize[u][1][i]
                    p = mod(blocks[u][1][i][j], om)
                    p = p != 0 ? p : om
                    q = mod(blocks[u][1][i][k], om)
                    q = q != 0 ? q : om
                    ind = p <= q ? p + Int(q*(q-1)/2) : q + Int(p*(p-1)/2)
                    Locb = bfind(ksupp[ind], length(ksupp[ind]), sadd(basis[u][ceil(Int, blocks[u][1][i][j]/om)], basis[u][ceil(Int, blocks[u][1][i][k]/om)]))
                    if p != q || j == k
                        @inbounds add_to_expression!(cons[ind][Locb], pos[i][j,k])
                    else
                        @inbounds add_to_expression!(cons[ind][Locb], 2, pos[i][j,k])
                    end
                end
            end
            gpos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(I[u]))
            for (s,v) = enumerate(I[u])
                gpos[s] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[u][s+1])
                com = cons_matrix[v].m*om
                for i = 1:cl[u][s+1]
                    gpos[s][i] = @variable(model, [1:blocksize[u][s+1][i], 1:blocksize[u][s+1][i]], PSD)
                    for j = 1:blocksize[u][s+1][i], k = j:blocksize[u][s+1][i]
                        p = mod(blocks[u][s+1][i][j], com)
                        p = p != 0 ? p : com
                        q = mod(blocks[u][s+1][i][k], com)
                        q = q != 0 ? q : com
                        p1 = ceil(Int, p/cons_matrix[v].m)
                        q1 = ceil(Int, q/cons_matrix[v].m)
                        ind = p1 <= q1 ? p1 + Int(q1*(q1-1)/2) : q1 + Int(p1*(p1-1)/2)
                        t = mod(blocks[u][s+1][i][j], cons_matrix[v].m)
                        t = t != 0 ? t : cons_matrix[v].m
                        r = mod(blocks[u][s+1][i][k], cons_matrix[v].m)
                        r = r != 0 ? r : cons_matrix[v].m
                        loc = t <= r ? t + Int(r*(r-1)/2) : r + Int(t*(t-1)/2)
                        for w = 1:length(cons_matrix[v].poly[loc].supp)
                            Locb = bfind(ksupp[ind], length(ksupp[ind]), sadd(sadd(gbasis[u][s][ceil(Int, blocks[u][s+1][i][j]/com)], 
                            gbasis[u][s][ceil(Int, blocks[u][s+1][i][k]/com)]), cons_matrix[v].poly[loc].supp[w]))
                            if p != q || j == k
                                @inbounds add_to_expression!(cons[ind][Locb], cons_matrix[v].poly[loc].coe[w], gpos[s][i][j,k])
                            else
                                @inbounds add_to_expression!(cons[ind][Locb], 2*cons_matrix[s].poly[loc].coe[w], gpos[s][i][j,k])
                            end
                        end
                    end
                end
            end
        end
        @variable(model, lower)
        for i = 1:om, j = i:om
            ind = i + Int(j*(j-1)/2)
            bc = zeros(length(ksupp[ind]))
            for k = 1:length(obj_matrix.poly[ind].supp)
                Locb = bfind(ksupp[ind], length(ksupp[ind]), obj_matrix.poly[ind].supp[k])
                if Locb === nothing
                   @error "The monomial basis is not enough!"
                   return nothing,nothing,nothing
                else
                   bc[Locb] = obj_matrix.poly[ind].coe[k]
                end
            end
            if i == j
                cons[ind][1] += lower
            end
            @constraint(model, cons[ind].==bc)
        end
        @objective(model, Max, lower)
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
            println("Solving the SDP...")
        end
        time = @elapsed begin
        optimize!(model)
        end
        if QUIET == false
            println("SDP solving time: $time seconds.")
        end
        SDP_status = termination_status(model)
        objv = objective_value(model)
        if SDP_status != MOI.OPTIMAL
            println("termination status: $SDP_status")
            status = primal_status(model)
            println("solution status: $status")
        end
        println("optimum = $objv")
    end
    return objv,ksupp,SDP_status
end

function LinearPMI_first(b, F::Vector{Matrix{Polynomial{true, T}}}, G, x, d; TS="block", QUIET=false, solve=true) where {T<:Number}
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    n = length(x)
    s = length(F)
    m = length(G)
    dG = [maximum(maxdegree.(vec(G[i]))) for i=1:m]
    obj_matrix = Vector{poly_matrix}(undef, s)
    for k = 1:s
        obj_matrix[k] = poly_matrix(size(F[k],1), Vector{poly_data}(undef, Int((size(F[k],1)+1)*size(F[k],1)/2)))
        for i = 1:obj_matrix[k].m, j = i:obj_matrix[k].m
            _,supp,coe = polys_info([F[k][i,j]], x)
            obj_matrix[k].poly[i+Int(j*(j-1)/2)] = poly_data(n, supp[1], coe[1])
        end
    end
    basis = get_sbasis(Vector(1:n), d)
    cons_matrix = Vector{poly_matrix}(undef, m)
    gbasis = Vector{Vector{Vector{UInt16}}}(undef, m)
    for k = 1:m
        gbasis[k] = get_sbasis(Vector(1:n), d-Int(ceil(dG[k]/2)))
        cons_matrix[k] = poly_matrix(size(G[k],1), Vector{poly_data}(undef, Int((size(G[k],1)+1)*size(G[k],1)/2)))
        for i = 1:cons_matrix[k].m, j = i:cons_matrix[k].m
            _,supp,coe = polys_info([G[k][i,j]], x)
            cons_matrix[k].poly[i+Int(j*(j-1)/2)] = poly_data(n, supp[1], coe[1])
        end
    end
    ksupp = Vector{Vector{Vector{UInt16}}}(undef, Int((obj_matrix[1].m+1)*obj_matrix[1].m/2))
    if TS != false
        for i = 1:obj_matrix[1].m, j = i:obj_matrix[1].m
            ind = i + Int(j*(j-1)/2)
            ksupp[ind] = reduce(vcat, [obj_matrix[k].poly[ind].supp for k=1:s])
            if i == j
                for item in basis
                    push!(ksupp[ind], sadd(item, item))
                end
            end
        end
        unique!.(ksupp)
        sort!.(ksupp)
        if QUIET == false
            println("Starting to compute the block structure...")
        end
    end
    time = @elapsed begin
    blocks,cl,blocksize,sb,numb,_ = get_mblocks(obj_matrix[1].m, cons_matrix, ksupp, basis, gbasis, TS=TS, QUIET=QUIET)
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,SDP_status = LinearPMI_sdp(b, obj_matrix, cons_matrix, basis, gbasis, blocks, cl, blocksize, QUIET=QUIET, solve=solve)
    data = mpop_data(b, obj_matrix, cons_matrix, basis, gbasis, ksupp, cl, blocksize, blocks, nothing, nothing, nothing, nothing, sb, numb, SDP_status)
    return opt,data
end

function LinearPMI_higher!(data::mpop_data; TS="block", QUIET=false, solve=true)
    basis = data.basis
    gbasis = data.gbasis
    ksupp = data.ksupp
    obj_matrix = data.obj_matrix
    cons_matrix = data.cons_matrix
    blocks = data.blocks
    cl = data.cl
    blocksize = data.blocksize
    sb = data.sb
    numb = data.numb
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,data.sb,data.numb,status = get_mblocks(obj_matrix[1].m, cons_matrix, ksupp, basis, gbasis, blocks=blocks, cl=cl, blocksize=blocksize, sb=sb, numb=numb, TS=TS, QUIET=QUIET)
    end
    opt = nothing
    if status == 1
        if TS != false && QUIET == false
            mb = maximum(maximum.(data.sb))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,SDP_status = LinearPMI_sdp(data.b, obj_matrix, cons_matrix, basis, gbasis, blocks, cl, blocksize, QUIET=QUIET, solve=solve)
        data.ksupp = ksupp
        data.blocks = blocks
        data.cl = cl
        data.blocksize = blocksize
        data.SDP_status = SDP_status
    end
    return opt,data
end

function LinearPMI_sdp(b, obj_matrix, cons_matrix, basis, gbasis, blocks, cl, blocksize; solve=true, QUIET=false)
    om = obj_matrix[1].m
    ksupp = [Vector{UInt16}[] for i = 1:length(obj_matrix[1].poly)]
    for i = 1:cl[1], j = 1:blocksize[1][i], k = j:blocksize[1][i]
        bi = sadd(basis[ceil(Int, blocks[1][i][j]/om)], basis[ceil(Int, blocks[1][i][k]/om)])
        ind1 = mod(blocks[1][i][j], om)
        ind1 = ind1 != 0 ? ind1 : om
        ind2 = mod(blocks[1][i][k], om)
        ind2 = ind2 != 0 ? ind2 : om
        push!(ksupp[ind1+Int(ind2*(ind2-1)/2)], bi)
    end
    for s = 1:length(cons_matrix)
        com = cons_matrix[s].m*om
        for i = 1:cl[s+1], j = 1:blocksize[s+1][i], k = j:blocksize[s+1][i]
            p = mod(blocks[s+1][i][j], com)
            p = p != 0 ? p : com
            q = mod(blocks[s+1][i][k], com)
            q = q != 0 ? q : com
            p1 = ceil(Int, p/cons_matrix[s].m)
            q1 = ceil(Int, q/cons_matrix[s].m)
            ind = p1 <= q1 ? p1 + Int(q1*(q1-1)/2) : q1 + Int(p1*(p1-1)/2)
            t = mod(blocks[s+1][i][j], cons_matrix[s].m)
            t = t != 0 ? t : cons_matrix[s].m
            r = mod(blocks[s+1][i][k], cons_matrix[s].m)
            r = r != 0 ? r : cons_matrix[s].m
            loc = t <= r ? t + Int(r*(r-1)/2) : r + Int(t*(t-1)/2)
            for w = 1:length(cons_matrix[s].poly[loc].supp)
                bi = sadd(sadd(gbasis[s][ceil(Int, blocks[s+1][i][j]/com)], gbasis[s][ceil(Int, blocks[s+1][i][k]/com)]), cons_matrix[s].poly[loc].supp[w])
                push!(ksupp[ind], bi)
            end
        end
    end
    unique!.(ksupp)
    sort!.(ksupp)
    objv = SDP_status = nothing
    if solve == true
        if QUIET == false
            ncons = sum(length.(ksupp))
            println("Assembling the SDP...")
            println("There are $ncons affine constraints.")
        end
        model = Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = Vector{Vector{AffExpr}}(undef, length(obj_matrix[1].poly))
        for i = 1:length(obj_matrix[1].poly)
            cons[i] = [AffExpr(0) for j=1:length(ksupp[i])]
        end
        pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1])
        for i = 1:cl[1]
            pos[i] = @variable(model, [1:blocksize[1][i], 1:blocksize[1][i]], PSD)
            for j = 1:blocksize[1][i], k = j:blocksize[1][i]
                p = mod(blocks[1][i][j], om)
                p = p != 0 ? p : om
                q = mod(blocks[1][i][k], om)
                q = q != 0 ? q : om
                ind = p <= q ? p + Int(q*(q-1)/2) : q + Int(p*(p-1)/2)
                Locb = bfind(ksupp[ind], length(ksupp[ind]), sadd(basis[ceil(Int, blocks[1][i][j]/om)], basis[ceil(Int, blocks[1][i][k]/om)]))
                if p != q || j == k
                    @inbounds add_to_expression!(cons[ind][Locb], pos[i][j,k])
                else
                    @inbounds add_to_expression!(cons[ind][Locb], 2, pos[i][j,k])
                end
            end
        end
        gpos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(cons_matrix))
        for s = 1:length(cons_matrix)
            gpos[s] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[s+1])
            com = cons_matrix[s].m*om
            for i = 1:cl[s+1]
                gpos[s][i] = @variable(model, [1:blocksize[s+1][i], 1:blocksize[s+1][i]], base_name="P", PSD)
                for j = 1:blocksize[s+1][i], k = j:blocksize[s+1][i]
                    p = mod(blocks[s+1][i][j], com)
                    p = p != 0 ? p : com
                    q = mod(blocks[s+1][i][k], com)
                    q = q != 0 ? q : com
                    p1 = ceil(Int, p/cons_matrix[s].m)
                    q1 = ceil(Int, q/cons_matrix[s].m)
                    ind = p1 <= q1 ? p1 + Int(q1*(q1-1)/2) : q1 + Int(p1*(p1-1)/2)
                    t = mod(blocks[s+1][i][j], cons_matrix[s].m)
                    t = t != 0 ? t : cons_matrix[s].m
                    r = mod(blocks[s+1][i][k], cons_matrix[s].m)
                    r = r != 0 ? r : cons_matrix[s].m
                    loc = t <= r ? t + Int(r*(r-1)/2) : r + Int(t*(t-1)/2)
                    for w = 1:length(cons_matrix[s].poly[loc].supp)
                        Locb = bfind(ksupp[ind], length(ksupp[ind]), sadd(sadd(gbasis[s][ceil(Int, blocks[s+1][i][j]/com)], 
                        gbasis[s][ceil(Int, blocks[s+1][i][k]/com)]), cons_matrix[s].poly[loc].supp[w]))
                        if p != q || j == k
                            @inbounds add_to_expression!(cons[ind][Locb], cons_matrix[s].poly[loc].coe[w], gpos[s][i][j,k])
                        else
                            @inbounds add_to_expression!(cons[ind][Locb], 2*cons_matrix[s].poly[loc].coe[w], gpos[s][i][j,k])
                        end
                    end
                end
            end
        end
        λ = @variable(model, [1:length(b)])
        for i = 1:om, j = i:om
            ind = i + Int(j*(j-1)/2)
            bc = [AffExpr(0) for k = 1:length(ksupp[ind])]
            for k = 1:length(obj_matrix[1].poly[ind].supp)
                Locb = bfind(ksupp[ind], length(ksupp[ind]), obj_matrix[1].poly[ind].supp[k])
                if Locb === nothing
                    @error "The monomial basis is not enough!"
                    return nothing,nothing,nothing
                else
                    @inbounds add_to_expression!(bc[Locb], obj_matrix[1].poly[ind].coe[k])
                end
            end
            for t = 2:length(obj_matrix), k = 1:length(obj_matrix[t].poly[ind].supp)
                Locb = bfind(ksupp[ind], length(ksupp[ind]), obj_matrix[t].poly[ind].supp[k])
                if Locb === nothing
                    @error "The monomial basis is not enough!"
                    return nothing,nothing,nothing
                else
                    @inbounds add_to_expression!(bc[Locb], λ[t-1], obj_matrix[t].poly[ind].coe[k])
                end
            end
            @constraint(model, cons[ind].==bc)
        end
        @objective(model, Min, b'*λ)
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
            println("Solving the SDP...")
        end
        time = @elapsed begin
        optimize!(model)
        end
        if QUIET == false
            println("SDP solving time: $time seconds.")
        end
        SDP_status = termination_status(model)
        objv = objective_value(model)
        if SDP_status != MOI.OPTIMAL
            println("termination status: $SDP_status")
            status = primal_status(model)
            println("solution status: $status")
        end
        println("optimum = $objv")
    end
    return objv,ksupp,SDP_status
end