mutable struct ncmpop_data
    n::Int # number of all variables
    m::Int # number of all constraints
    numeq::Int # number of equality constraints
    d::Int # relaxation order
    supp # support data
    coe # coefficient data
    obj # "eigen" or "trace"
    ksupp # extending support at the k-th step
    basis # monomial bses
    cql # number of cliques
    cliques # cliques of variables
    cliquesize # numbers of cliques
    J # constraints associated to each clique
    ncc # constraints associated to no clique
    blocks # block structure
    cl # numbers of blocks
    blocksize # sizes of blocks
    sb # sizes of different blocks
    numb # numbers of different blocks
end

function cs_nctssos_first(f, x; d=0, CS="MF", minimize=false, TS="block", merge=false, md=3,
    QUIET=false, obj="eigen", solve=true)
    println("***************************NCTSSOS***************************")
    println("NCTSSOS is launching...")
    n,supp,coe = poly_info(f, x)
    if d == 0
        d = ceil(Int, maxdegree(f)/2)
    end
    opt,data = cs_nctssos_first(supp, coe, n, d=d, CS=CS, minimize=minimize, TS=TS, merge=merge, md=md, QUIET=QUIET, obj=obj, solve=solve)
    return opt,data
end

function cs_nctssos_first(supp::Vector{Vector{UInt16}}, coe, n::Int; d=0, CS="MF",
    minimize=false, TS="block", merge=false, md=3, QUIET=false, obj="eigen", solve=true)
    if obj == "trace"
        supp,coe = cyclic_canon(supp, coe)
    else
        supp,coe = sym_canon(supp, coe)
    end
    time = @elapsed begin
    cliques,cql,cliquesize = clique_decomp(n, supp)
    end
    if CS != false && QUIET == false
        mc = maximum(cliquesize)
        println("Obtained the variable cliques in $time seconds. The maximal size of cliques is $mc.")
    end
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,sb,numb,basis,_ = get_blocks_mix(d, supp, cliques, cql, cliquesize, TS=TS, merge=merge, md=md, obj=obj)
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
    end
    opt,ksupp = blockupop_mix(n, supp, coe, basis, cliques, cql, cliquesize, blocks, cl, blocksize, obj=obj, solve=solve, QUIET=QUIET)
    data = ncmpop_data(n, 0, 0, d, supp, coe, obj, ksupp, basis, cql, cliques, cliquesize, [], [], blocks, cl, blocksize, sb, numb)
    return opt,data
end

"""
    opt,data = cs_nctssos_first(pop, x, d; numeq=0, CS="MF", TS="block", merge=false, md=3,
    QUIET=false, obj="eigen", solve=true)

Compute the first step of the CS-NCTSSOS hierarchy for constrained noncommutative polynomial
optimization with relaxation order `d`. Return the optimum and other auxiliary data.

# Arguments
- `pop`: the vector of the objective function, inequality constraints, and equality constraints.
- `x`: the set of noncommuting variables.
- `d`: the relaxation order of the moment-SOHS hierarchy.
- `numeq`: the number of equality constraints.
"""
function cs_nctssos_first(pop, x, d; numeq=0, CS="MF", minimize=false, assign="first", TS="block", merge=false, md=3, QUIET=false, obj="eigen", solve=true)
    n,supp,coe = polys_info(pop, x)
    opt,data = cs_nctssos_first(supp, coe, n, d, numeq=numeq, CS=CS, minimize=minimize, assign=assign, TS=TS, QUIET=QUIET, obj=obj, solve=solve)
    return opt,data
end

"""
    opt,data = cs_nctssos_first(supp::Vector{Vector{Vector{UInt16}}}, coe, n::Int, d::Int; numeq=0,
    CS="MF", TS="block", merge=false, md=3, QUIET=false, obj="eigen", solve=true)

Compute the first step of the CS-NCTSSOS hierarchy for constrained noncommutative polynomial optimization
with relaxation order `d`. Here the polynomial optimization problem is defined by `supp` and `coe`,
corresponding to the supports and coeffients of `pop` respectively. Return the optimum and other auxiliary data.

# Arguments
- `supp`: the supports of the polynomial optimization problem.
- `coe`: the coeffients of the polynomial optimization problem.
- `d`: the relaxation order of the moment-SOHS hierarchy.
- `numeq`: the number of equality constraints.
"""
function cs_nctssos_first(supp::Vector{Vector{Vector{UInt16}}}, coe, n::Int, d::Int; numeq=0, CS="MF",
    minimize=false, assign="first", TS="block", merge=false, md=3, QUIET=false, obj="eigen", solve=true)
    println("***************************NCTSSOS***************************")
    println("NCTSSOS is launching...")
    m = length(supp)-1
    dg = [maximum(length.(supp[i])) for i=2:m+1]
    if obj == "trace"
        supp[1],coe[1] = cyclic_canon(supp[1], coe[1])
    else
        supp[1],coe[1] = sym_canon(supp[1], coe[1])
    end
    time = @elapsed begin
    cliques,cql,cliquesize = clique_decomp(n, m, d, dg, supp, alg=CS, minimize=minimize)
    end
    if CS != false && QUIET == false
        mc = maximum(cliquesize)
        println("Obtained the variable cliques in $time seconds. The maximal size of cliques is $mc.")
    end
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    J,ncc = assign_constraint(m, supp, cliques, cql, cliquesize, assign=assign)
    time = @elapsed begin
    blocks,cl,blocksize,sb,numb,basis,status = get_cblocks_mix(d, dg, J, m, supp, cliques, cql, cliquesize, TS=TS, obj=obj)
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
    end
    opt,ksupp = blockcpop_mix(n, m, supp, coe, basis, cliques, cql, cliquesize, J, ncc, blocks, cl, blocksize, numeq=numeq, QUIET=QUIET, obj=obj, solve=solve)
    data = ncmpop_data(n, m, numeq, d, supp, coe, obj, ksupp, basis, cql, cliques, cliquesize, J, ncc, blocks, cl, blocksize, sb, numb)
    return opt,data
end

"""
    opt,data = cs_nctssos_higher!(data; TS="block", QUIET=false, merge=false, md=3, solve=true)

Compute higher steps of the CS-NCTSSOS hierarchy.
Return the optimum and other auxiliary data.
"""
function cs_nctssos_higher!(data::ncmpop_data; TS="block", QUIET=false, merge=false, md=3, solve=true)
    n = data.n
    m = data.m
    numeq = data.numeq
    d = data.d
    supp = data.supp
    coe = data.coe
    obj = data.obj
    ksupp = data.ksupp
    basis = data.basis
    cql = data.cql
    cliques = data.cliques
    cliquesize = data.cliquesize
    J = data.J
    ncc = data.ncc
    blocks = data.blocks
    cl = data.cl
    blocksize = data.blocksize
    sb = data.sb
    numb = data.numb
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    if m == 0
        time = @elapsed begin
        blocks,cl,blocksize,sb,numb,basis,status = get_blocks_mix(d, supp, cliques, cql, cliquesize, basis=basis, sb=sb, numb=numb, TS=TS, merge=merge, md=md, obj=obj)
        end
        if status == 1
            if QUIET == false
                mb = maximum(maximum.(sb))
                println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
            end
            opt,ksupp = blockupop_mix(n, supp, coe, basis, cliques, cql, cliquesize, blocks, cl, blocksize, obj=obj, solve=solve, QUIET=QUIET)
        end
    else
        time = @elapsed begin
        blocks,cl,blocksize,sb,numb,basis,status = get_cblocks_mix(d, [], J, m, supp, cliques, cql, cliquesize, ksupp=ksupp, basis=basis, blocks=blocks, cl=cl, blocksize=blocksize, sb=sb, numb=numb, TS=TS, merge=merge, md=md, obj=obj)
        end
        if status==1
            if QUIET == false
                mb = maximum(maximum.(sb))
                println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
            end
            opt,ksupp = blockcpop_mix(n, m, supp, coe, basis, cliques, cql, cliquesize, J, ncc, blocks, cl, blocksize, numeq=numeq, QUIET=QUIET, obj=obj, solve=solve)
        end
    end
    if status == 0
        opt = nothing
        println("No higher CS-NCTSSOS hierarchy!")
    end
    data.ksupp = ksupp
    data.blocks = blocks
    data.cl = cl
    data.blocksize = blocksize
    data.sb = sb
    data.numb = numb
    return opt,data
end

function blockupop_mix(n, supp, coe, basis, cliques, cql, cliquesize, blocks, cl, blocksize; QUIET=false, obj="eigen", solve=true)
    ksupp = Vector{UInt16}[]
    for i = 1:cql, j = 1:cl[i], k = 1:blocksize[i][j], r = k:blocksize[i][j]
        @inbounds bi = [basis[i][blocks[i][j][k]][end:-1:1]; basis[i][blocks[i][j][r]]]
        push!(ksupp, bi)
    end
    ksupp = reduce!.(ksupp, obj=obj)
    sort!(ksupp)
    unique!(ksupp)
    lksupp = length(ksupp)
    if QUIET == false
        println("There are $lksupp affine constraints.")
    end
    objv = nothing
    if solve == true
        if QUIET == false
            println("Assembling the SDP...")
        end
        model = Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = [AffExpr(0) for i=1:lksupp]
        pos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, cql)
        for i = 1:cql
            pos[i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i])
            for k = 1:cl[i]
                if blocksize[i][k] == 1
                   pos[i][k] = @variable(model, lower_bound=0)
                   @inbounds bi = [basis[i][blocks[i][k][1]][end:-1:1]; basis[i][blocks[i][k][1]]]
                   bi = reduce!(bi, obj=obj)
                   Locb = ncbfind(ksupp, lksupp, bi)
                   @inbounds add_to_expression!(cons[Locb], pos[i][k])
                else
                   pos[i][k] = @variable(model, [1:blocksize[i][k], 1:blocksize[i][k]], PSD)
                   for j = 1:blocksize[i][k], r = j:blocksize[i][k]
                       @inbounds ind1 = blocks[i][k][j]
                       @inbounds ind2 = blocks[i][k][r]
                       @inbounds bi = [basis[i][ind1][end:-1:1]; basis[i][ind2]]
                       bi = reduce!(bi, obj=obj)
                       Locb = ncbfind(ksupp, lksupp, bi)
                       if j == r
                           @inbounds add_to_expression!(cons[Locb], pos[i][k][j,r])
                       else
                           @inbounds add_to_expression!(cons[Locb], 2, pos[i][k][j,r])
                       end
                   end
                end
            end
        end
        bc = zeros(lksupp)
        for i = 1:length(supp)
            Locb = ncbfind(ksupp, lksupp, supp[i])
            if Locb == 0
               @error "The monomial basis is not enough!"
               return nothing,nothing
            else
               bc[Locb] = coe[i]
            end
        end
        @variable(model, lower)
        cons[1] += lower
        @constraint(model, cons.==bc)
        @objective(model, Max, lower)
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
            println("Solving the SDP...")
        end
        time=@elapsed begin
        optimize!(model)
        end
        if QUIET == false
            println("SDP solving time: $time seconds.")
        end
        status = termination_status(model)
        objv = objective_value(model)
        if status != MOI.OPTIMAL
           println("termination status: $status")
           status = primal_status(model)
           println("solution status: $status")
        end
        println("optimum = $objv")
    end
    return objv,ksupp
end

function blockcpop_mix(n, m, supp, coe, basis, cliques, cql, cliquesize, J, ncc, blocks, cl, blocksize; numeq=0, QUIET=false, obj="eigen", solve=true)
    ksupp = Vector{UInt16}[]
    for i = 1:cql
        for j = 1:cl[i][1], k = 1:blocksize[i][1][j], r = k:blocksize[i][1][j]
            @inbounds bi = [basis[i][1][blocks[i][1][j][k]][end:-1:1]; basis[i][1][blocks[i][1][j][r]]]
            push!(ksupp, bi)
        end
        for (j, w) in enumerate(J[i])
            for l = 1:cl[i][j+1], t = 1:blocksize[i][j+1][l], r = t:blocksize[i][j+1][l], s = 1:length(supp[w+1])
                ind1 = blocks[i][j+1][l][t]
                ind2 = blocks[i][j+1][l][r]
                @inbounds bi = [basis[i][j+1][ind1][end:-1:1]; supp[w+1][s]; basis[i][j+1][ind2]]
                push!(ksupp, bi)
            end
        end
    end
    for i ∈ ncc
        append!(ksupp, supp[i+1])
    end
    ksupp = reduce!.(ksupp, obj=obj)
    sort!(ksupp)
    unique!(ksupp)
    lksupp = length(ksupp)
    if QUIET == false
        println("There are $lksupp affine constraints.")
    end
    objv = nothing
    if solve == true
        if QUIET == false
            println("Assembling the SDP...")
        end
        model = Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = [AffExpr(0) for i=1:lksupp]
        for i = 1:cql, l = 1:cl[i][1]
            if blocksize[i][1][l] == 1
               @inbounds pos = @variable(model, lower_bound=0)
               @inbounds bi = [basis[i][1][blocks[i][1][l][1]][end:-1:1]; basis[i][1][blocks[i][1][l][1]]]
               bi = reduce!(bi, obj=obj)
               Locb = ncbfind(ksupp, lksupp,bi)
               @inbounds add_to_expression!(cons[Locb], pos)
            else
               @inbounds bs = blocksize[i][1][l]
               @inbounds pos = @variable(model, [1:bs, 1:bs], PSD)
               for t = 1:bs, r = t:bs
                   @inbounds ind1 = blocks[i][1][l][t]
                   @inbounds ind2 = blocks[i][1][l][r]
                   @inbounds bi = [basis[i][1][ind1][end:-1:1]; basis[i][1][ind2]]
                   bi = reduce!(bi, obj=obj)
                   Locb = ncbfind(ksupp, lksupp, bi)
                   if t == r
                      @inbounds add_to_expression!(cons[Locb], pos[t,r])
                   else
                      @inbounds add_to_expression!(cons[Locb], 2, pos[t,r])
                   end
               end
            end
        end
        for k = 1:length(ncc)
            i = ncc[k]
            if i <= m-numeq
                pos = @variable(model, lower_bound=0)
            else
                pos = @variable(model)
            end
            for j = 1:length(supp[i+1])
                bi = reduce!(supp[i+1][j], obj=obj)
                Locb = ncbfind(ksupp, lksupp, bi)
                @inbounds add_to_expression!(cons[Locb], coe[i+1][j], pos)
            end
        end
        for i = 1:cql, (j, w) in enumerate(J[i])
            for l = 1:cl[i][j+1]
                bs = blocksize[i][j+1][l]
                if bs == 1
                    if j <= m-numeq
                        pos = @variable(model, lower_bound=0)
                    else
                        pos = @variable(model)
                    end
                    ind1 = blocks[i][j+1][l][1]
                    for s = 1:length(supp[w+1])
                        @inbounds bi = [basis[i][j+1][ind1][end:-1:1]; supp[w+1][s]; basis[i][j+1][ind1]]
                        bi = reduce!(bi, obj=obj)
                        Locb = ncbfind(ksupp, lksupp, bi)
                        @inbounds add_to_expression!(cons[Locb], coe[w+1][s], pos)
                    end
                else
                    if j <= m-numeq
                        pos = @variable(model, [1:bs, 1:bs], PSD)
                    else
                        pos = @variable(model, [1:bs, 1:bs], Symmetric)
                    end
                    for t = 1:bs, r = t:bs
                        ind1 = blocks[i][j+1][l][t]
                        ind2 = blocks[i][j+1][l][r]
                        for s = 1:length(supp[w+1])
                            @inbounds bi = [basis[i][j+1][ind1][end:-1:1]; supp[w+1][s]; basis[i][j+1][ind2]]
                            bi = reduce!(bi, obj=obj)
                            Locb = ncbfind(ksupp, lksupp, bi)
                            if t == r
                                @inbounds add_to_expression!(cons[Locb], coe[w+1][s], pos[t,r])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*coe[w+1][s], pos[t,r])
                            end
                        end
                    end
                end
            end
        end
        bc = zeros(lksupp)
        for i = 1:length(supp[1])
            Locb = ncbfind(ksupp, lksupp, supp[1][i])
            if Locb == 0
               @error "The monomial basis is not enough!"
               return nothing,nothing
            else
               bc[Locb] = coe[1][i]
            end
        end
        @variable(model, lower)
        cons[1] += lower
        @constraint(model, cons.==bc)
        @objective(model, Max, lower)
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
            println("Solving the SDP...")
        end
        time=@elapsed begin
        optimize!(model)
        end
        if QUIET == false
            println("SDP solving time: $time seconds.")
        end
        status = termination_status(model)
        objv = objective_value(model)
        if status != MOI.OPTIMAL
           println("termination status: $status")
           status = primal_status(model)
           println("solution status: $status")
        end
        println("optimum = $objv")
    end
    return objv,ksupp
end

function get_blocks_mix(d, supp, cliques, cql, cliquesize; basis=[], sb=[], numb=[], TS="block", merge=false, md=3, obj="eigen")
    blocks = Vector{Vector{Vector{UInt16}}}(undef, cql)
    cl = Vector{UInt16}(undef, cql)
    blocksize = Vector{Vector{UInt16}}(undef, cql)
    status = ones(UInt8, cql)
    if isempty(basis)
        sb = Vector{Vector{UInt16}}(undef, cql)
        numb = Vector{Vector{UInt16}}(undef, cql)
        basis = Vector{Vector{Vector{UInt16}}}(undef, cql)
        flag = 1
    else
        flag = 0
    end
    for i = 1:cql
        nvar = cliquesize[i]
        ind = [issubset(supp[j], cliques[i]) for j=1:length(supp)]
        ksupp = copy(supp[ind])
        if flag == 1
            basis[i] = get_ncbasis(nvar, d, ind=cliques[i])
            if obj == "trace"
                append!(ksupp, [_cyclic_canon([basis[i][k][end:-1:1]; basis[i][k]]) for k=1:length(basis[i])])
            else
                append!(ksupp, [[basis[i][k][end:-1:1]; basis[i][k]] for k=1:length(basis[i])])
            end
            sort!(ksupp)
            unique!(ksupp)
            blocks[i],cl[i],blocksize[i],sb[i],numb[i],status[i] = get_ncblocks(ksupp, basis[i], TS=TS, obj=obj, QUIET=true, merge=merge, md=md)
        else
            blocks[i],cl[i],blocksize[i],sb[i],numb[i],status[i] = get_ncblocks(ksupp, basis[i], sb=sb[i], numb=numb[i], TS=TS, obj=obj, QUIET=true, merge=merge, md=md)
        end
    end
    return blocks,cl,blocksize,sb,numb,basis,maximum(status)
end

function get_cblocks_mix(d, dg, J, m, supp, cliques, cql, cliquesize; ksupp=[], basis=[], blocks=[], cl=[], blocksize=[], sb=[], numb=[], TS="block", merge=false, md=3, obj="eigen")
    if isempty(basis)
        blocks = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        cl = Vector{Vector{UInt16}}(undef, cql)
        blocksize = Vector{Vector{Vector{UInt16}}}(undef, cql)
        sb = Vector{Vector{UInt16}}(undef, cql)
        numb = Vector{Vector{UInt16}}(undef, cql)
        basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        ksupp = copy(supp[1])
        for i = 2:m+1
            append!(ksupp,  _sym_canon.(supp[i]))
        end
        sort!(ksupp)
        unique!(ksupp)
        flag = 1
    else
        flag = 0
    end
    status = ones(UInt8, cql)
    for i = 1:cql
        lc = length(J[i])
        nvar = cliquesize[i]
        ind = [issubset(ksupp[j], cliques[i]) for j=1:length(ksupp)]
        fsupp = copy(ksupp[ind])
        if flag == 1
            basis[i] = Vector{Vector{Vector{UInt16}}}(undef, lc+1)
            basis[i][1] = get_ncbasis(cliquesize[i], d, ind=cliques[i])
            for s = 1:lc
                basis[i][s+1] = get_ncbasis(nvar, d-ceil(Int, dg[J[i][s]]/2), ind=cliques[i])
            end
            blocks[i] = Vector{Vector{Vector{UInt16}}}(undef, lc+1)
            cl[i] = Vector{UInt16}(undef, lc+1)
            blocksize[i] = Vector{Vector{UInt16}}(undef, lc+1)
            sb[i] = Vector{UInt16}(undef, lc+1)
            numb[i] = Vector{UInt16}(undef, lc+1)
            if obj == "trace"
                append!(fsupp, [_cyclic_canon([basis[i][1][k][end:-1:1]; basis[i][1][k]]) for k=1:length(basis[i][1])])
            else
                append!(fsupp, [[basis[i][1][k][end:-1:1]; basis[i][1][k]] for k=1:length(basis[i][1])])
            end
            sort!(fsupp)
            unique!(fsupp)
            blocks[i],cl[i],blocksize[i],sb[i],numb[i],status[i] = get_nccblocks(lc, fsupp, supp[J[i].+1], basis[i], TS=TS, QUIET=true, merge=merge, md=md, obj=obj)
        else
            blocks[i],cl[i],blocksize[i],sb[i],numb[i],status[i] = get_nccblocks(lc, fsupp, supp[J[i].+1], basis[i], blocks=blocks[i], cl=cl[i], blocksize=blocksize[i], sb=sb[i], numb=numb[i], TS=TS, QUIET=true, merge=merge, md=md, obj=obj)
        end
    end
    return blocks,cl,blocksize,sb,numb,basis,maximum(status)
end

function assign_constraint(m, supp, cliques, cql, cliquesize; assign="first")
    J = [UInt16[] for i=1:cql]
    ncc = UInt16[]
    for i = 2:m+1
        rind = copy(supp[i][1])
        for j = 2:length(supp[i])
            append!(rind, supp[i][j])
        end
        rind = unique(rind)
        if assign == "first"
            ind = findfirst(k->issubset(rind, cliques[k]), 1:cql)
            if ind != nothing
                push!(J[ind], i-1)
            else
                push!(ncc, i-1)
            end
        else
            temp = UInt16[]
            for j = 1:cql
                if issubset(rind, cliques[j])
                    push!(temp,j)
                end
            end
            if !isempty(temp)
                if assign == "min"
                    push!(J[temp[argmin(cliquesize[temp])]], i-1)
                else
                    push!(J[temp[argmax(cliquesize[temp])]], i-1)
                end
            else
                push!(ncc, i-1)
            end
        end
    end
    return J,ncc
end

function clique_decomp(n::Int, supp::Vector{Vector{UInt16}}; alg="MF", minimize=false)
    if alg == false
        cliques = [UInt16[i for i=1:n]]
        cql = 1
        cliquesize = [n]
    else
        G = SimpleGraph(n)
        for j = 1:length(supp)
            add_clique!(G, unique(supp[j]))
        end
        if alg == "NC"
            cliques,cql,cliquesize = max_cliques(G)
        else
            cliques,cql,cliquesize = chordal_cliques!(G, method=alg, minimize=minimize)
        end
    end
    uc = unique(cliquesize)
    sizes = [sum(cliquesize.== i) for i in uc]
    println("------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("------------------------------------------------------")
    return cliques,cql,cliquesize
end

function clique_decomp(n::Int, m::Int, d::Int, dg::Vector{Int}, supp::Vector{Vector{Vector{UInt16}}}; alg="MF", minimize=false)
    if alg == false
        cliques = [UInt16[i for i=1:n]]
        cql = 1
        cliquesize=[n]
    else
        G = SimpleGraph(n)
        for i = 1:m+1
            if i == 1 || d == ceil(Int, dg[i-1]/2)
                for j = 1:length(supp[i])
                    add_clique!(G, unique(supp[i][j]))
                end
            else
                temp = copy(supp[i][1])
                for j = 2:length(supp[i])
                    append!(temp, supp[i][j])
                end
                add_clique!(G, unique(temp))
            end
        end
        if alg == "NC"
            cliques,cql,cliquesize = max_cliques(G)
        else
            cliques,cql,cliquesize = chordal_cliques!(G, method=alg, minimize=minimize)
        end
    end
    uc = unique(cliquesize)
    sizes = [sum(cliquesize.== i) for i in uc]
    println("------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("------------------------------------------------------")
    return cliques,cql,cliquesize
end
