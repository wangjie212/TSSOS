mutable struct ncupop_type
    supp # support data
    coe # coefficient data
    partition # the first 'partition' variables commutes with the remaining variables
    constraint # "projection" or "unipotent"
    basis # basis
    obj # "eigen" or "trace"
    ksupp # extending support at the k-th step
    sb # sizes of different blocks
    numb # numbers of different blocks
end

"""
    opt,data = nctssos_first(f::Polynomial{false, T} where T<:Number, x::Vector{PolyVar{false}};
        newton=true, reducebasis=true, TS="block", obj="eigen", merge=false, md=3, solve=true, QUIET=false)

Compute the first step of the NCTSSOS hierarchy for unconstrained noncommutative polynomial optimization.
If `newton=true`, then compute a monomial basis by the Newton chip method.
If `reducebasis=true`, then remove monomials from the monomial basis by diagonal inconsistency.
If `TS="block"`, use maximal chordal extensions; if `TS="MD"`, use approximately smallest chordal
extensions. If obj="eigen", minimize the eigenvalue; if obj="trace", then minimize the trace.
If `merge=true`, perform the PSD block merging. Return the optimum and other auxiliary data.

# Arguments
- `f`: the objective function for unconstrained noncommutative polynomial optimization.
- `x`: the set of noncommuting variables.
- `md`: the tunable parameter for merging blocks.
"""
function nctssos_first(f::Polynomial{false, T} where T<:Number, x::Vector{PolyVar{false}}; d=0, newton=true, reducebasis=true,
    monosquare=true, TS="block", obj="eigen", partition=0, constraint=nothing, merge=false, md=3, solve=true, QUIET=false)
    n,supp,coe = poly_info(f, x)
    opt,data = nctssos_first(supp, coe, n, d=d, newton=newton, reducebasis=reducebasis, monosquare=monosquare, TS=TS, obj=obj,
    partition=partition, constraint=constraint, merge=merge, md=md, solve=solve, QUIET=QUIET)
    return opt,data
end

function nctssos_first(supp::Vector{Vector{UInt16}}, coe, n::Int; d=0, newton=true, reducebasis=true, monosquare=true,
    TS="block", obj="eigen", partition=0, constraint=nothing, merge=false, md=3, solve=true, QUIET=false)
    println("************************NCTSSOS************************")
    println("NCTSSOS is launching...")
    if d == 0
        d = Int(maximum(length.(supp))/2)
    end
    if obj == "trace"
        supp,coe = cyclic_canon(supp, coe)
    else
        supp,coe = sym_canon(supp, coe)
    end
    if newton == true && constraint == nothing
        if obj == "trace"
            basis = newton_cyclic(supp, n, d)
        else
            basis = newton_ncbasis(supp)
        end
    else
        basis = get_ncbasis(n, d)
    end
    if partition > 0
        ind = [_comm(basis[i], partition) == basis[i] for i=1:length(basis)]
        basis = basis[ind]
    end
    if constraint != nothing
        ind = [findfirst(j -> basis[i][j]== basis[i][j+1], 1:length(basis[i])-1) == nothing for i=1:length(basis)]
        basis = basis[ind]
    end
    ksupp = copy(supp)
    if monosquare == true && TS != false
        if obj == "trace"
            append!(ksupp, [_cyclic_canon([basis[i][end:-1:1]; basis[i]]) for i=1:length(basis)])
        else
            append!(ksupp, [[basis[i][end:-1:1]; basis[i]] for i=1:length(basis)])
        end
        if partition > 0
            ksupp = _comm.(ksupp, partition)
        end
        if constraint != nothing
            reduce_cons!.(ksupp, constraint = constraint)
        end
        sort!(ksupp)
        unique!(ksupp)
    end
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    blocks,cl,blocksize,sb,numb,_ = get_ncblocks(ksupp, basis, TS=TS, QUIET=QUIET, merge=merge, md=md, obj=obj, partition=partition, constraint=constraint)
    if reducebasis == true && obj == "eigen" && constraint == nothing
        psupp = copy(supp)
        psupp = psupp[is_sym.(psupp)]
        push!(psupp, UInt16[])
        basis,flag = reducebasis!(psupp, basis, blocks, cl, blocksize)
        if flag == 1
            ksupp = copy(supp)
            if monosquare == true
                if obj == "trace"
                    append!(ksupp, [_cyclic_canon([basis[i][end:-1:1]; basis[i]]) for i=1:length(basis)])
                else
                    append!(ksupp, [[basis[i][end:-1:1]; basis[i]] for i=1:length(basis)])
                end
                sort!(ksupp)
                unique!(ksupp)
            end
            blocks,cl,blocksize,sb,numb,_ = get_ncblocks(ksupp, basis, TS=TS, QUIET=QUIET, merge=merge, md=md, obj=obj)
        end
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure. The maximal size of blocks is $mb.")
    end
    opt,ksupp = ncblockupop(supp, coe, basis, blocks, cl, blocksize, QUIET=QUIET, obj=obj, partition=partition, constraint=constraint, solve=solve)
    data = ncupop_type(supp, coe, partition, constraint, basis, obj, ksupp, sb, numb)
    return opt,data
end

function poly_info(f, x)
    n = length(x)
    mon = monomials(f)
    coe = coefficients(f)
    lm = length(mon)
    supp = [UInt16[] for i=1:lm]
    for i = 1:lm
        ind = mon[i].z .> 0
        vars = mon[i].vars[ind]
        exp = mon[i].z[ind]
        for j = 1:length(vars)
            k = ncbfind(x, n, vars[j], rev=true)
            append!(supp[i], k*ones(UInt16, exp[j]))
        end
    end
    return n,supp,coe
end

"""
    opt,data = nctssos_higher!(data, TS="block", merge=false, md=3, QUIET=false)

Compute higher steps of the NCTSSOS hierarchy.
Return the optimum and other auxiliary data.
"""
function nctssos_higher!(data::ncupop_type; TS="block", merge=false, md=3, solve=true, QUIET=false)
    supp = data.supp
    basis = data.basis
    coe = data.coe
    partition = data.partition
    constraint = data.constraint
    obj = data.obj
    ksupp = data.ksupp
    sb = data.sb
    numb = data.numb
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    blocks,cl,blocksize,sb,numb,status = get_ncblocks(ksupp, basis, sb=sb, numb=numb, TS=TS, QUIET=QUIET, merge=merge, md=md, obj=obj, partition=partition, constraint=constraint)
    opt = nothing
    if status == 1
        if QUIET == false
            mb = maximum(maximum.(sb))
            println("Obtained the block structure. The maximal size of blocks is $mb.")
        end
        opt,ksupp = ncblockupop(supp, coe, basis, blocks, cl, blocksize, QUIET=QUIET, obj=obj, partition=partition, constraint=constraint, solve=solve)
    end
    data.ksupp = ksupp
    data.sb = sb
    data.numb = numb
    return opt,data
end

function cc(a::Vector{UInt16}, n::Int)
    ua = unique(a)
    ca = zeros(UInt8, n)
    for i = 1:length(ua)
        ca[ua[i]] = count(x->isequal(ua[i], x), a)
    end
    return ca
end

function remove(csupp, dw, n)
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    t = @variable(model)
    alpha = @variable(model, [1:n])
    @constraint(model, [i=1:length(csupp)], alpha'*(csupp[i].-2*dw)<=t)
    @objective(model, Min, t)
    optimize!(model)
    if objective_value(model) >= 0
        return true
    else
        return false
    end
end

function permutation(a)
    b = sparse(a)
    ua = convert(Vector{UInt16}, b.nzind)
    na = convert(Vector{UInt16}, b.nzval)
    return _permutation(ua, na)
end

function _permutation(ua, na)
    if !isempty(ua)
        perm = Vector{UInt16}[]
        for i = 1:length(ua)
            nua = copy(ua)
            nna = copy(na)
            if na[i] == 1
                deleteat!(nua, i)
                deleteat!(nna, i)
            else
                nna[i] -= 1
            end
            temp = _permutation(nua, nna)
            push!.(temp, ua[i])
            append!(perm, temp)
        end
        return perm
    else
        return [UInt16[]]
    end
end

function newton_cyclic(supp, n, d)
    pbasis = get_basis(n,d)
    basis = [UInt16[]]
    csupp = cc.(supp, n)
    pushfirst!(csupp, zeros(UInt8, n))
    sort!(csupp)
    unique!(csupp)
    for i = 2:size(pbasis,2)
        if remove(csupp, pbasis[:,i], n)
            append!(basis, permutation(pbasis[:,i]))
        end
    end
    sort!(basis)
    return basis
end

function get_basis(n, d)
    lb = binomial(n+d, d)
    basis = zeros(UInt8, n, lb)
    i = 0
    t = 1
    while i < d+1
        t += 1
        if basis[n,t-1] == i
           if i < d
              basis[1,t] = i+1
           end
           i += 1
        else
            j = findfirst(x->basis[x,t-1]!=0, 1:n)
            basis[:,t] = basis[:,t-1]
            if j == 1
               basis[1,t] -= 1
               basis[2,t] += 1
            else
               basis[1,t] = basis[j,t] - 1
               basis[j,t] = 0
               basis[j+1,t] += 1
            end
        end
    end
    return basis
end

function _cyclic_canon(a::Vector{UInt16})
    if isempty(a)
        return a
    else
        return minimum([[a[i+1:length(a)]; a[1:i]] for i=0:length(a)-1])
    end
end

function cyclic_canon(supp, coe)
    nsupp = [min(_cyclic_canon(word), _cyclic_canon(reverse(word))) for word in supp]
    sort!(nsupp)
    unique!(nsupp)
    l = length(nsupp)
    ncoe = zeros(l)
    for i = 1:length(supp)
        bi = min(_cyclic_canon(supp[i]), _cyclic_canon(reverse(supp[i])))
        Locb = ncbfind(nsupp, l, bi)
        ncoe[Locb] += coe[i]
    end
    return nsupp,ncoe
end

function _sym_canon(a::Vector{UInt16})
    i = 1
    while i <= Int(ceil((length(a)-1)/2))
        if a[i] < a[end+1-i]
            return a
        elseif a[i] > a[end+1-i]
            return reverse(a)
        else
            i += 1
        end
    end
    return a
end

function is_sym(a::Vector{UInt16})
    l = Int(ceil((length(a)-1)/2))
    return isequal(a[1:l], a[end:-1:end-l+1])
end

function sym_canon(supp, coe)
    nsupp = copy(supp)
    nsupp = _sym_canon.(nsupp)
    sort!(nsupp)
    unique!(nsupp)
    l = length(nsupp)
    ncoe = zeros(l)
    for i = 1:length(supp)
        Locb = ncbfind(nsupp, l, _sym_canon(supp[i]))
        ncoe[Locb] += coe[i]
    end
    return nsupp,ncoe
end

function get_ncbasis(n, d; ind=UInt16[i for i=1:n])
    basis = [UInt16[]]
    for i = 1:d
        append!(basis, _get_ncbasis_deg(n, i, ind=ind))
    end
    return basis
end

function _get_ncbasis_deg(n, d; ind=UInt16[i for i=1:n])
    if d > 0
        basis = Vector{UInt16}[]
        for i = 1:n
            temp = _get_ncbasis_deg(n, d-1, ind=ind)
            push!.(temp, ind[i])
            append!(basis, temp)
        end
        return basis
    else
        return [UInt16[]]
    end
end

function newton_ncbasis(supp)
    nbasis = [UInt16[]]
    for bi in supp
        if iseven(length(bi))
            k = Int(length(bi)/2)
            w = bi[end-k+1:end]
            if isequal(w, bi[k:-1:1])
                for j = 1:k
                    push!(nbasis, w[end-j+1:end])
                end
            end
        end
    end
    sort!(nbasis)
    unique!(nbasis)
    return nbasis
end

function ncbfind(A, l, a; rev=false)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if isequal(A[mid], a)
           return mid
        elseif isless(A[mid], a)
            if rev == false
                low = mid+1
            else
                high = mid-1
            end
        else
            if rev == false
                high = mid-1
            else
                low = mid+1
            end
        end
    end
    return 0
end

function get_ncgraph(ksupp, basis; obj="eigen", partition=0, constraint=nothing)
    lb = length(basis)
    G = SimpleGraph(lb)
    lksupp = length(ksupp)
    for i = 1:lb, j = i+1:lb
        bi = [basis[i][end:-1:1]; basis[j]]
        bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint)
        if ncbfind(ksupp, lksupp, bi) != 0
           add_edge!(G, i, j)
        end
    end
    return G
end

function get_ncblocks(ksupp, basis; sb=[], numb=[], TS="block", obj="eigen", partition=0, constraint=nothing, QUIET=true, merge=false, md=3)
    if TS == false
        blocksize = [length(basis)]
        blocks = [[i for i=1:length(basis)]]
        cl = 1
    else
        G = get_ncgraph(ksupp, basis, obj=obj, partition=partition, constraint=constraint)
        if TS == "block"
            blocks = connected_components(G)
            blocksize = length.(blocks)
            cl = length(blocksize)
        else
            blocks,cl,blocksize = chordal_cliques!(G, method=TS)
            if merge == true
                blocks,cl,blocksize = clique_merge!(blocks, d=md, QUIET=true)
            end
        end
    end
    nsb = sort(unique(blocksize), rev=true)
    nnumb = [sum(blocksize.== i) for i in nsb]
    if isempty(sb) || nsb!=sb || nnumb!=numb
        status = 1
        if QUIET == false
            println("------------------------------------------------------")
            println("The sizes of PSD blocks:\n$nsb\n$nnumb")
            println("------------------------------------------------------")
        end
    else
        status = 0
        println("No higher NCTSSOS hierarchy!")
    end
    return blocks,cl,blocksize,nsb,nnumb,status
end

function ncblockupop(supp, coe, basis, blocks, cl, blocksize; QUIET=true, obj="eigen", partition=0, constraint=nothing, solve=true)
    ksupp = Vector{Vector{UInt16}}(undef, Int(sum(Int.(blocksize).^2+blocksize)/2))
    k = 1
    for i = 1:cl, j = 1:blocksize[i], r = j:blocksize[i]
        @inbounds bi = [basis[blocks[i][j]][end:-1:1]; basis[blocks[i][r]]]
        @inbounds ksupp[k] = bi
        k += 1
    end
    ksupp = reduce!.(ksupp, obj=obj, partition=partition, constraint=constraint)
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
               @inbounds bi = [basis[blocks[i][1]][end:-1:1]; basis[blocks[i][1]]]
               bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint)
               Locb = ncbfind(ksupp,lksupp,bi)
               @inbounds add_to_expression!(cons[Locb], pos)
            else
               @inbounds pos = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:blocksize[i], r = j:blocksize[i]
                   @inbounds bi = [basis[blocks[i][j]][end:-1:1]; basis[blocks[i][r]]]
                   bi = reduce!(bi, obj=obj, partition=partition, constraint=constraint)
                   Locb = ncbfind(ksupp, lksupp, bi)
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

function reduce_cons!(word; constraint="unipotent")
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

function reduce!(word; obj="eigen", partition=0, constraint=nothing)
    if obj == "trace"
        word = min(_cyclic_canon(word), _cyclic_canon(reverse(word)))
    else
        if partition > 0 && constraint == nothing
            word = min(_comm(word, partition), _comm(reverse(word), partition))
        elseif partition == 0 && constraint != nothing
            cword = copy(word)
            word = min(reduce_cons!(word, constraint = constraint), reduce_cons!(reverse(cword), constraint = constraint))
        elseif partition > 0 && constraint != nothing
            word = min(reduce_cons!(_comm(word, partition), constraint = constraint), reduce_cons!(_comm(reverse(word), partition), constraint = constraint))
        else
            word = _sym_canon(word)
        end
    end
    return word
end

function _comm(word::Vector{UInt16}, partition)
    ind1 = word .<= partition
    ind2 = word .> partition
    return [word[ind1]; word[ind2]]
end
