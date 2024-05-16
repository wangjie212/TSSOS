mutable struct poly_data
    n::Int
    supp::Vector{Vector{UInt16}}
    coe::Vector{Number}
end

mutable struct poly_matrix
    m::Int
    poly::Vector{poly_data}
end

mutable struct mpop_data
    b
    obj_matrix
    cons_matrix
    basis # monomial basis
    gbasis
    ksupp # extended support at the k-th step
    blocks # block structrue
    cl
    blocksize
    sb # sizes of different blocks
    numb # numbers of different blocks
    SDP_status
end

function tssos_first(F::Matrix{Polynomial{true, T}}, G, x, d; TS="block", QUIET=false, solve=true) where {T<:Number}
    println("*********************************** TSSOS ***********************************")
    println("Version 1.1.2, developed by Jie Wang, 2020--2024")
    println("TSSOS is launching...")
    n = length(x)
    m = length(G)
    dG = [maximum(maxdegree.(vec(G[i]))) for i=1:m]
    obj_matrix = poly_matrix(size(F,1), Vector{poly_data}(undef, Int((size(F,1)+1)*size(F,1)/2)))
    basis = get_sbasis(Vector(1:n), d)
    for i = 1:obj_matrix.m, j = i:obj_matrix.m
        _,supp,coe = polys_info([F[i,j]], x)
        obj_matrix.poly[i+Int(j*(j-1)/2)] = poly_data(n, supp[1], coe[1])
    end
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
    ksupp = Vector{Vector{Vector{UInt16}}}(undef, Int((obj_matrix.m+1)*obj_matrix.m/2))
    if TS != false
        for i = 1:obj_matrix.m, j = i:obj_matrix.m
            ind = i + Int(j*(j-1)/2)
            ksupp[ind] = copy(obj_matrix.poly[ind].supp)
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
    blocks,cl,blocksize,sb,numb,_ = get_mblocks(obj_matrix.m, cons_matrix, ksupp, basis, gbasis, TS=TS, QUIET=QUIET)
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,SDP_status = pmo_sdp(obj_matrix, cons_matrix, basis, gbasis, blocks, cl, blocksize, QUIET=QUIET, solve=solve)
    data = mpop_data(nothing, obj_matrix, cons_matrix, basis, gbasis, ksupp, blocks, cl, blocksize, sb, numb, SDP_status)
    return opt,data
end

function tssos_higher!(data::mpop_data; TS="block", QUIET=false, solve=true)
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
    blocks,cl,blocksize,data.sb,data.numb,status = get_mblocks(obj_matrix.m, cons_matrix, ksupp, basis, gbasis, blocks=blocks, cl=cl, blocksize=blocksize, sb=sb, numb=numb, TS=TS, QUIET=QUIET)
    end
    opt = nothing
    if status == 1
        if TS != false && QUIET == false
            mb = maximum(maximum.(data.sb))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,SDP_status = pmo_sdp(obj_matrix, cons_matrix, basis, gbasis, blocks, cl, blocksize, QUIET=QUIET, solve=solve)
        data.ksupp = ksupp
        data.blocks = blocks
        data.cl = cl
        data.blocksize = blocksize
        data.SDP_status = SDP_status
    end
    return opt,data
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
            blocks[i] = [[j for j=1:lb]]
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

function pmo_sdp(obj_matrix, cons_matrix, basis, gbasis, blocks, cl, blocksize; solve=true, QUIET=false)
    om = obj_matrix.m
    ksupp = [Vector{UInt16}[] for i = 1:length(obj_matrix.poly)]
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
        cons = Vector{Vector{AffExpr}}(undef, length(obj_matrix.poly))
        for i = 1:length(obj_matrix.poly)
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
    println("Version 1.1.2, developed by Jie Wang, 2020--2024")
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
    data = mpop_data(b, obj_matrix, cons_matrix, basis, gbasis, ksupp, blocks, cl, blocksize, sb, numb, SDP_status)
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