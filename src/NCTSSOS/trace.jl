mutable struct traceopt_type
    supp # support data
    coe # coefficient data
    constraint # "projection" or "unipotent"
    ptsupp # pure trace support
    wbasis # word basis
    tbasis # trace basis
    basis # non-trace basis
    ksupp # extending support at the k-th step
    sb # sizes of different blocks
    numb # numbers of different blocks
end

function get_sbasis(n, d, ptsupp)
    basis = Vector{UInt16}[[]]
    i = 0
    temp = UInt16[]
    while i < d+1
        if sum(temp) == n*i
            temp = ones(UInt16, i+1)
            if i < d && sum(length.(ptsupp[temp])) <= d
                push!(basis, temp)
            end
            i += 1
        else
            temp2 = copy(temp)
            j = temp[1]
            ind = findfirst(x->temp[x]!=j, 1:length(temp))
            if ind == nothing
                ind = length(temp)+1
            end
            if j != 1
                temp2[1:ind-2] = ones(UInt16, ind-2)
            end
            temp2[ind-1] = j+1
            temp = temp2
            if sum(length.(ptsupp[temp])) <= d
                push!(basis, temp)
            end
        end
    end
    return basis
end

function trace_basis(n, d, ptsupp, bsupp)
    ind = [length(bsupp[i]) <= d for i = 1:length(bsupp)]
    basis = bsupp[ind]
    tbasis = get_sbasis(length(ptsupp), d, ptsupp)
    wbasis = Vector{UInt16}[]
    for i = 1:length(tbasis), j = 1:length(basis)
        if sum(length.(ptsupp[tbasis[i]])) + length(basis[j]) <= d
            push!(wbasis, [i,j])
        end
    end
    return wbasis,tbasis,basis
end

function sym_cyclic(word)
    return min(_cyclic_canon(word), _cyclic_canon(reverse(word)))
end

function ptraceopt_first(tr_supp, coe, n, d; TS="block", monosquare=false, QUIET=false, constraint="unipotent", solve=true)
    println("********************************** NCTSSOS **********************************")
    println("Version 0.2.0, developed by Jie Wang, 2020--2022")
    println("NCTSSOS is launching...")
    bsupp = get_ncbasis(n, 2d)
    ind = [length(bsupp[i]) <= 1 || (bsupp[i][1] != bsupp[i][end] && findfirst(j -> bsupp[i][j]
    == bsupp[i][j+1], 1:length(bsupp[i])-1) == nothing) for i = 1:length(bsupp)]
    bsupp = bsupp[ind]
    ind = [sym_cyclic(bsupp[i])==bsupp[i] for i=1:length(bsupp)]
    ptsupp = bsupp[ind]
    ptsupp = ptsupp[2:end]
    sort!(ptsupp)
    supp = Vector{Vector{UInt16}}(undef, length(tr_supp))
    for i = 1:length(tr_supp)
        supp[i] = sort([ncbfind(ptsupp, length(ptsupp), tr_supp[i][j]) for j=1:length(tr_supp[i])])
    end
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    wbasis,tbasis,basis = trace_basis(n, d, ptsupp, bsupp)
    ksupp = copy(supp)
    if monosquare == true
        for i = 1:length(wbasis)
            bi1 = sort([tbasis[wbasis[i][1]]; tbasis[wbasis[i][1]]])
            bi2 = [reverse(basis[wbasis[i][2]]); basis[wbasis[i][2]]]
            constraint_reduce!(bi2, constraint=constraint)
            bi = trace_reduce(bi1, bi2, ptsupp)
            push!(ksupp, bi)
        end
    end
    sort!(ksupp)
    unique!(ksupp)
    blocks,cl,blocksize,sb,numb,_ = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, TS=TS, QUIET=QUIET, constraint=constraint)
    end
    if QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
    end
    opt,ksupp = ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, QUIET=QUIET, constraint=constraint, solve=solve)
    data = traceopt_type(supp, coe, constraint, ptsupp, wbasis, tbasis, basis, ksupp, sb, numb)
    return opt,data
end

function ptraceopt_higher!(data; TS="block", QUIET=false, solve=true)
    supp = data.supp
    coe = data.coe
    constraint = data.constraint
    ptsupp = data.ptsupp
    wbasis = data.wbasis
    tbasis = data.tbasis
    basis = data.basis
    ksupp = data.ksupp
    sb = data.sb
    numb = data.numb
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    blocks,cl,blocksize,sb,numb,status = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, sb=sb, numb=numb, TS=TS, QUIET=QUIET, constraint=constraint)
    opt = nothing
    if status == 1
        if QUIET == false
            mb = maximum(maximum.(sb))
            println("Obtained the block structure. The maximal size of blocks is $mb.")
        end
        opt,ksupp = ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, QUIET=QUIET, constraint=constraint, solve=solve)
    end
    data.ksupp = ksupp
    data.sb = sb
    data.numb = numb
    return opt,data
end

function ptrace_SDP(supp, coe, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize; QUIET=false, constraint="unipotent", solve=true)
    ksupp = Vector{Vector{UInt16}}(undef, Int(sum(Int.(blocksize).^2+blocksize)/2))
    k = 1
    for i = 1:cl, j = 1:blocksize[i], r = j:blocksize[i]
        @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][r]][1]]])
        @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][r]][2]]]
        constraint_reduce!(bi2, constraint=constraint)
        @inbounds ksupp[k] = trace_reduce(bi1, bi2, ptsupp)
        k += 1
    end
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
        cons = [AffExpr(0) for i=1:lksupp]
        for i = 1:cl
            bs = blocksize[i]
            if bs == 1
               @inbounds pos = @variable(model, lower_bound=0)
               @inbounds bi1 = sort([tbasis[wbasis[blocks[i][1]][1]]; tbasis[wbasis[blocks[i][1]][1]]])
               @inbounds bi2 = [reverse(basis[wbasis[blocks[i][1]][2]]); basis[wbasis[blocks[i][1]][2]]]
               constraint_reduce!(bi2, constraint=constraint)
               bi = trace_reduce(bi1, bi2, ptsupp)
               Locb = ncbfind(ksupp, lksupp, bi)
               @inbounds add_to_expression!(cons[Locb], pos)
            else
               @inbounds pos = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:blocksize[i], r = j:blocksize[i]
                   @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][r]][1]]])
                   @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][r]][2]]]
                   constraint_reduce!(bi2, constraint=constraint)
                   bi = trace_reduce(bi1, bi2, ptsupp)
                   Locb = ncbfind(ksupp, lksupp, bi)
                   if Locb == 0
                       @error "The word does not exist!"
                       return nothing,nothing
                   end
                   if j == r
                       @inbounds add_to_expression!(cons[Locb], pos[j,r])
                   else
                       @inbounds add_to_expression!(cons[Locb], 2, pos[j,r])
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
        if QUIET == false
            println("Solving the SDP...")
        end
        time = @elapsed begin
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

function trace_reduce(word1, word2, ptsupp)
    if isempty(word2)
        ind = UInt16[]
    else
        ind = UInt16(ncbfind(ptsupp, length(ptsupp), sym_cyclic(word2)))
    end
    return sort([word1; ind])
end

function constraint_reduce!(word; constraint="unipotent")
    while length(word) > 1 && word[1] == word[end]
        deleteat!(word, 1)
        if constraint == "unipotent"
            deleteat!(word, length(word))
        end
    end
    i = 1
    while i < length(word)
        if word[i] == word[i+1]
            deleteat!(word, i)
            if constraint == "unipotent"
                deleteat!(word, i)
            end
        else
            i += 1
        end
    end
    return word
end

function get_ncgraph(ksupp, ptsupp, wbasis, tbasis, basis; constraint="unipotent")
    lb = length(wbasis)
    G = SimpleGraph(lb)
    lksupp = length(ksupp)
    for i = 1:lb, j = i+1:lb
        @inbounds bi1 = sort([tbasis[wbasis[i][1]]; tbasis[wbasis[j][1]]])
        @inbounds bi2 = [reverse(basis[wbasis[i][2]]); basis[wbasis[j][2]]]
        constraint_reduce!(bi2, constraint=constraint)
        bi = trace_reduce(bi1, bi2, ptsupp)
        if ncbfind(ksupp, lksupp, bi) != 0
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis; sb=[], numb=[], TS="block", QUIET=false, constraint="unipotent")
    if TS == false
        blocksize = [length(wbasis)]
        blocks = [[i for i=1:length(wbasis)]]
        cl = 1
    else
        G = get_ncgraph(ksupp, ptsupp, wbasis, tbasis, basis, constraint=constraint)
        if TS == "block"
            blocks = connected_components(G)
            blocksize = length.(blocks)
            cl = length(blocksize)
        else
            blocks,cl,blocksize = chordal_cliques!(G, method=TS)
        end
    end
    nsb = sort(unique(blocksize), rev=true)
    nnumb = [sum(blocksize.== i) for i in nsb]
    if isempty(sb) || nsb!=sb || nnumb!=numb
        status = 1
        if QUIET == false
            println("-----------------------------------------------------------------------------")
            println("The sizes of PSD blocks:\n$nsb\n$nnumb")
            println("-----------------------------------------------------------------------------")
        end
    else
        status = 0
        println("No higher TS step of the NCTSSOS hierarchy!")
    end
    return blocks,cl,blocksize,nsb,nnumb,status
end

function Werner_witness_first(dY, sigma, n, d; TS="block", monosquare=false, QUIET=false, solve=true)
    println("********************************** NCTSSOS **********************************")
    println("Version 0.2.0, developed by Jie Wang, 2020--2022")
    println("NCTSSOS is launching...")
    bsupp = get_ncbasis(n, 2d)
    ind = [length(bsupp[i]) <= 1 || (bsupp[i][1] != bsupp[i][end] && findfirst(j -> bsupp[i][j]
    == bsupp[i][j+1], 1:length(bsupp[i])-1) == nothing) for i = 1:length(bsupp)]
    bsupp = bsupp[ind]
    ind = [sym_cyclic(bsupp[i])==bsupp[i] for i=1:length(bsupp)]
    ptsupp = bsupp[ind]
    ptsupp = ptsupp[2:end]
    sort!(ptsupp)
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    wbasis,tbasis,basis = trace_basis(n, d, ptsupp, bsupp)
    htrace = generate_htrace(n)
    supp = Vector{Vector{UInt16}}(undef, length(htrace))
    for i = 1:length(htrace)
        supp[i] = sort([ncbfind(ptsupp, length(ptsupp), sym_cyclic(htrace[i][j])) for j=1:length(htrace[i])])
    end
    if monosquare == true
        for i = 1:length(wbasis)
            bi1 = sort([tbasis[wbasis[i][1]]; tbasis[wbasis[i][1]]])
            bi2 = [reverse(basis[wbasis[i][2]]); basis[wbasis[i][2]]]
            constraint_reduce!(bi2, constraint="projection")
            bi = trace_reduce(bi1, bi2, ptsupp)
            push!(supp, bi)
        end
    end
    sort!(supp)
    unique!(supp)
    blocks,cl,blocksize,sb,numb,_ = get_ncblocks(supp, ptsupp, wbasis, tbasis, basis, TS=TS, QUIET=QUIET, constraint="projection")
    end
    if QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds. The maximal size of blocks is $mb.")
    end
    opt,ksupp = Werner_SDP(dY, sigma, htrace, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, QUIET=QUIET, solve=solve)
    data = traceopt_type(htrace, dY, sigma, ptsupp, wbasis, tbasis, basis, ksupp, sb, numb)
    return opt,data
end

function Werner_witness_higher!(data; TS="block", QUIET=false, solve=true)
    htrace = data.supp
    dY = data.coe
    sigma = data.constraint
    ptsupp = data.ptsupp
    wbasis = data.wbasis
    tbasis = data.tbasis
    basis = data.basis
    ksupp = data.ksupp
    sb = data.sb
    numb = data.numb
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    blocks,cl,blocksize,sb,numb,status = get_ncblocks(ksupp, ptsupp, wbasis, tbasis, basis, sb=sb, numb=numb, TS=TS, QUIET=QUIET, constraint="projection")
    opt = nothing
    if status == 1
        if QUIET == false
            mb = maximum(maximum.(sb))
            println("Obtained the block structure. The maximal size of blocks is $mb.")
        end
        opt,ksupp = Werner_SDP(dY, sigma, htrace, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize, QUIET=QUIET, solve=solve)
    end
    data.ksupp = ksupp
    data.sb = sb
    data.numb = numb
    return opt,data
end

function Werner_SDP(dY, sigma, htrace, ptsupp, wbasis, tbasis, basis, blocks, cl, blocksize; QUIET=false, solve=true)
    ksupp = Vector{Vector{UInt16}}(undef, Int(sum(Int.(blocksize).^2+blocksize)/2))
    k = 1
    for i = 1:cl, j = 1:blocksize[i], r = j:blocksize[i]
        @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][r]][1]]])
        @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][r]][2]]]
        constraint_reduce!(bi2, constraint="projection")
        @inbounds ksupp[k] = trace_reduce(bi1, bi2, ptsupp)
        k += 1
    end
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
        cons = [AffExpr(0) for i=1:lksupp]
        for i = 1:cl
            bs = blocksize[i]
            if bs == 1
               @inbounds pos = @variable(model, lower_bound=0)
               @inbounds bi1 = sort([tbasis[wbasis[blocks[i][1]][1]]; tbasis[wbasis[blocks[i][1]][1]]])
               @inbounds bi2 = [reverse(basis[wbasis[blocks[i][1]][2]]); basis[wbasis[blocks[i][1]][2]]]
               constraint_reduce!(bi2, constraint="projection")
               bi = trace_reduce(bi1, bi2, ptsupp)
               Locb = ncbfind(ksupp, lksupp, bi)
               @inbounds add_to_expression!(cons[Locb], pos)
            else
               @inbounds pos = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:blocksize[i], r = j:blocksize[i]
                   @inbounds bi1 = sort([tbasis[wbasis[blocks[i][j]][1]]; tbasis[wbasis[blocks[i][r]][1]]])
                   @inbounds bi2 = [reverse(basis[wbasis[blocks[i][j]][2]]); basis[wbasis[blocks[i][r]][2]]]
                   constraint_reduce!(bi2, constraint="projection")
                   bi = trace_reduce(bi1, bi2, ptsupp)
                   Locb = ncbfind(ksupp, lksupp, bi)
                   if Locb == 0
                       @error "The word does not exist!"
                       return nothing,nothing
                   end
                   if j == r
                       @inbounds add_to_expression!(cons[Locb], pos[j,r])
                   else
                       @inbounds add_to_expression!(cons[Locb], 2, pos[j,r])
                   end
               end
            end
        end
        bc = [AffExpr(0) for i=1:lksupp]
        coe = @variable(model, [1:length(htrace)])
        tcons = AffExpr(0)
        for i = 1:length(htrace)
            for k = 1:length(dY)
                tcons += sigma[k]*prod(tr(prod(dY[k][htrace[i][j]])) for j=1:length(htrace[i]))*coe[i]
            end
            temp = sort([ncbfind(ptsupp, length(ptsupp), sym_cyclic(htrace[i][j])) for j=1:length(htrace[i])])
            Locb = ncbfind(ksupp, lksupp, temp)
            if Locb == 0
               @error "The monomial basis is not enough!"
               return nothing,nothing
            else
               bc[Locb] += coe[i]
            end
        end
        @constraint(model, tcons==-16)
        @variable(model, epsilon)
        cons[1] -= epsilon
        @constraint(model, cons.==bc)
        @objective(model, Min, epsilon)
        if QUIET == false
            println("Solving the SDP...")
        end
        time = @elapsed begin
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

function generate_htrace(n)
    return _generate_htrace(UInt16[i for i=1:n])
end

function _generate_htrace(var)
    if isempty(var)
         htrace = [Vector{UInt16}[]]
    else
        htrace = Vector{Vector{UInt16}}[]
        for i = 1:length(var)
            sset = _selete(var[2:end], i-1)
            pushfirst!.(sset, var[1])
            for mem in sset
                cbasis = _cyclic_basis(mem)
                for emem in cbasis
                    sub_htrace = _generate_htrace(UInt16.(setdiff(var, mem)))
                    for j = 1:length(sub_htrace)
                        pushfirst!(sub_htrace[j], emem)
                    end
                    append!(htrace, sub_htrace)
                end
            end
        end
    end
    return htrace
end

function _selete(var, d)
    if d > 0
        set = Vector{UInt16}[]
        for i = 1:length(var)-d+1
            iset = _selete(var[i+1:end], d-1)
            pushfirst!.(iset, var[i])
            append!(set, iset)
        end
    else
        set = [UInt16[]]
    end
    return set
end

function _cyclic_basis(var)
    basis = _permutation(var, ones(length(var)))
    return unique(_cyclic_canon.(basis))
end
