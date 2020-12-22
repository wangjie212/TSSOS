mutable struct mdata_type
    n # the number of all variables
    m # the number of all constraints
    dg # the degree vector
    supp # the support data
    coe # the coefficient data
    numeq # the number of equality constraints
    rlorder # the relaxation order
    supp0 # the extending support at each sparse order
    fbasis # the whole basis
    gbasis # the basis of constraint multiplier
    cql # the number of cliques
    cliques # cliques of variables
    cliquesize # the sizes of cliques
    J # constraints associated to each clique
    ncc # constraints associated to no clique
    blocks # the block structure
    cl # the number of blocks
    blocksize # the sizes of blocks
    ub # the unique sizes of blocks
    sizes # the number of different blocks
end

function cs_tssos_first(supp::Vector{Vector{Vector{UInt16}}},coe,n,d,dg;numeq=0,foc=100,CS="MF",minimize=true,assign="first",TS="block",QUIET=false,solve=true,solution=false,MomentOne=true)
    m=length(supp)-1
    cliques,cql,cliquesize=clique_decomp(n,m,dg,supp,order=d,alg=CS,minimize=minimize)
    J,ncc=assign_constraint(m,supp,cliques,cql,cliquesize,assign=assign)
    rlorder=init_order(dg,J,cliquesize,cql,foc=foc,order=d)
    blocks,cl,blocksize,ub,sizes,fbasis,gbasis,status=get_cblocks_mix!(dg,J,rlorder,m,supp,cliques,cql,cliquesize,TS=TS)
    opt,supp0=blockcpop_mix(n,m,dg,supp,coe,fbasis,gbasis,cliques,cql,cliquesize,J,ncc,blocks,cl,blocksize,numeq=numeq,QUIET=QUIET,TS=TS,solve=solve,solution=solution,MomentOne=MomentOne)
    data=mdata_type(n,m,dg,supp,coe,numeq,rlorder,supp0,fbasis,gbasis,cql,cliques,cliquesize,J,ncc,blocks,cl,blocksize,ub,sizes)
    return opt,data
end

function cs_tssos_higher!(data;TS="block",QUIET=false,solve=true,solution=false,MomentOne=false)
    n=data.n
    m=data.m
    dg=data.dg
    supp=data.supp
    coe=data.coe
    numeq=data.numeq
    rlorder=data.rlorder
    supp0=data.supp0
    fbasis=data.fbasis
    gbasis=data.gbasis
    cql=data.cql
    cliques=data.cliques
    cliquesize=data.cliquesize
    J=data.J
    ncc=data.ncc
    blocks=data.blocks
    cl=data.cl
    blocksize=data.blocksize
    ub=data.ub
    sizes=data.sizes
    blocks,cl,blocksize,ub,sizes,fbasis,gbasis,status=get_cblocks_mix!(dg,J,rlorder,m,supp,cliques,cql,cliquesize,tsupp=supp0,fbasis=fbasis,gbasis=gbasis,blocks=blocks,cl=cl,blocksize=blocksize,ub=ub,sizes=sizes,TS=TS)
    if status==1
        opt,supp0=blockcpop_mix(n,m,dg,supp,coe,fbasis,gbasis,cliques,cql,cliquesize,J,ncc,blocks,cl,blocksize,numeq=numeq,QUIET=QUIET,solve=solve,solution=solution,MomentOne=MomentOne)
        data.supp0=supp0
        data.blocks=blocks
        data.cl=cl
        data.blocksize=blocksize
        data.ub=ub
        data.sizes=sizes
    else
        opt=nothing
        println("No higher CS-TSSOS hierarchy!")
    end
    return opt,data
end

function blockcpop_mix(n,m,dg,supp,coe,fbasis,gbasis,cliques,cql,cliquesize,J,ncc,blocks,cl,blocksize;numeq=0,QUIET=false,TS="block",solve=true,solution=false,MomentOne=false)
    tsupp=Vector{UInt16}[]
    for i=1:cql, j=1:cl[i][1], k=1:blocksize[i][1][j], r=k:blocksize[i][1][j]
        @inbounds bi=sadd(fbasis[i][blocks[i][1][j][k]], fbasis[i][blocks[i][1][j][r]])
        push!(tsupp, bi)
    end
    if (MomentOne==true||solution==true)&&TS!=false
        supp0=copy(tsupp)
    end
    for i=1:cql, (j, w) in enumerate(J[i])
        for l=1:cl[i][j+1], t=1:blocksize[i][j+1][l], r=t:blocksize[i][j+1][l], s=1:length(supp[w+1])
            ind1=blocks[i][j+1][l][t]
            ind2=blocks[i][j+1][l][r]
            @inbounds bi=sadd(sadd(gbasis[i][j][ind1], supp[w+1][s]), gbasis[i][j][ind2])
            push!(tsupp, bi)
        end
    end
    for i ∈ ncc, j=1:length(supp[i+1])
        push!(tsupp, supp[i+1][j])
    end
    if (MomentOne==true||solution==true)&&TS!=false
        for i=1:cql, j=1:cliquesize[i]
            push!(tsupp, UInt16[cliques[i][j]])
            for k=j+1:cliquesize[i]
                push!(tsupp, UInt16[cliques[i][j];cliques[i][k]])
            end
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    if (MomentOne==true||solution==true)&&TS!=false
        sort!(supp0)
        unique!(supp0)
    else
        supp0=tsupp
    end
    objv=nothing
    if solve==true
        ltsupp=length(tsupp)
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:ltsupp]
        for i=1:cql
            if (MomentOne==true||solution==true)&&TS!=false
                bs=cliquesize[i]+1
                pos=@variable(model, [1:bs, 1:bs], PSD)
                for t=1:bs, r=t:bs
                    if t==1&&r==1
                        bi=UInt16[]
                    elseif t==1&&r>1
                        bi=UInt16[cliques[i][r-1]]
                    else
                        bi=UInt16[cliques[i][t-1];cliques[i][r-1]]
                    end
                    Locb=bfind(tsupp, ltsupp, bi)
                    if t==r
                        @inbounds cons[Locb]+=pos[t,r]
                    else
                        @inbounds cons[Locb]+=2*pos[t,r]
                    end
                end
            end
            for l=1:cl[i][1]
                @inbounds bs=blocksize[i][1][l]
                if bs==1
                    pos=@variable(model, lower_bound=0)
                    @inbounds bi=sadd(fbasis[i][blocks[i][1][l][1]], fbasis[i][blocks[i][1][l][1]])
                    Locb=bfind(tsupp,ltsupp,bi)
                    @inbounds cons[Locb]+=pos
                else
                    pos=@variable(model, [1:bs, 1:bs], PSD)
                    for t=1:bs, r=t:bs
                        @inbounds ind1=blocks[i][1][l][t]
                        @inbounds ind2=blocks[i][1][l][r]
                        @inbounds bi=sadd(fbasis[i][ind1], fbasis[i][ind2])
                        Locb=bfind(tsupp,ltsupp,bi)
                        if t==r
                            @inbounds cons[Locb]+=pos[t,r]
                        else
                            @inbounds cons[Locb]+=2*pos[t,r]
                        end
                    end
                end
            end
        end
        for i in ncc
            if i<=m-numeq
                pos=@variable(model, lower_bound=0)
            else
                pos=@variable(model)
            end
            for j=1:length(supp[i+1])
                Locb=bfind(tsupp,ltsupp,supp[i+1][j])
                @inbounds cons[Locb]+=coe[i+1][j]*pos
            end
        end
        for i=1:cql, (j, w) in enumerate(J[i]), l=1:cl[i][j+1]
            bs=blocksize[i][j+1][l]
            if bs==1
                if w<=m-numeq
                    pos=@variable(model, lower_bound=0)
                else
                    pos=@variable(model)
                end
                ind=blocks[i][j+1][l][1]
                for s=1:length(supp[w+1])
                    @inbounds bi=sadd(sadd(gbasis[i][j][ind], supp[w+1][s]), gbasis[i][j][ind])
                    Locb=bfind(tsupp,ltsupp,bi)
                    @inbounds cons[Locb]+=coe[w+1][s]*pos
                end
            else
                if w<=m-numeq
                    pos=@variable(model, [1:bs, 1:bs], PSD)
                else
                    pos=@variable(model, [1:bs, 1:bs], Symmetric)
                end
                for t=1:bs, r=t:bs
                    ind1=blocks[i][j+1][l][t]
                    ind2=blocks[i][j+1][l][r]
                    for s=1:length(supp[w+1])
                        @inbounds bi=sadd(sadd(gbasis[i][j][ind1], supp[w+1][s]), gbasis[i][j][ind2])
                        Locb=bfind(tsupp,ltsupp,bi)
                        if t==r
                            @inbounds cons[Locb]+=coe[w+1][s]*pos[t,r]
                        else
                            @inbounds cons[Locb]+=2*coe[w+1][s]*pos[t,r]
                        end
                    end
                end
            end
        end
        bc=zeros(ltsupp)
        for i=1:length(supp[1])
            Locb=bfind(tsupp,ltsupp,supp[1][i])
            if Locb==0
               @error "The monomial basis is not enough!"
               return nothing,nothing
            else
               bc[Locb]=coe[1][i]
            end
        end
        @variable(model, lower)
        @constraint(model, cons[2:end].==bc[2:end])
        @constraint(model, cons[1]+lower==bc[1])
        @objective(model, Max, lower)
        optimize!(model)
        status=termination_status(model)
        objv = objective_value(model)
        if status!=MOI.OPTIMAL
           println("termination status: $status")
           status=primal_status(model)
           println("solution status: $status")
        end
        println("optimum = $objv")
    end
    return objv,supp0
end

function init_order(dg,J,cliquesize,cql;foc=100,order="min")
    rlorder=ones(Int,cql)
    if order=="min"
        for i=1:cql
            if !isempty(J[i])
                rlorder[i]=ceil(Int, maximum(dg[J[i]])/2)
            end
        end
    else
        for i=1:cql
            if cliquesize[i]<=foc
                rlorder[i]=order
            end
        end
    end
    return rlorder
end

function get_cblocks_mix!(dg,J,rlorder,m,supp,cliques,cql,cliquesize;tsupp=[],fbasis=[],gbasis=[],blocks=[],cl=[],blocksize=[],ub=[],sizes=[],TS="block",merge=false)
    if isempty(tsupp)
        blocks=Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        cl=Vector{Vector{UInt16}}(undef, cql)
        blocksize=Vector{Vector{Vector{Int}}}(undef, cql)
        ub=Vector{Vector{UInt16}}(undef, cql)
        sizes=Vector{Vector{UInt16}}(undef, cql)
        fbasis=Vector{Vector{Vector{UInt16}}}(undef, cql)
        gbasis=Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        tsupp=copy(supp[1])
        for i=2:m+1, j=1:length(supp[i])
            push!(tsupp, supp[i][j])
        end
        sort!(tsupp)
        unique!(tsupp)
        flag=1
    else
        flag=0
    end
    status=ones(Int, cql)
    for i=1:cql
        lc=length(J[i])
        lt=length.(supp[J[i].+1])
        nvar=cliquesize[i]
        ind=[issubset(tsupp[j], cliques[i]) for j=1:length(tsupp)]
        fsupp=copy(tsupp[ind])
        if flag==1
            fbasis[i]=get_sbasis(cliques[i], rlorder[i])
            for j=1:length(fbasis[i])
                push!(fsupp, sadd(fbasis[i][j], fbasis[i][j]))
            end
            sort!(fsupp)
            unique!(fsupp)
            gbasis[i]=Vector{Vector{Vector{UInt16}}}(undef, lc)
            for s=1:lc
                gbasis[i][s]=get_sbasis(cliques[i], rlorder[i]-ceil(Int, dg[J[i][s]]/2))
            end
            blocks[i]=Vector{Vector{Vector{UInt16}}}(undef, lc+1)
            cl[i]=Vector{UInt16}(undef, lc+1)
            blocksize[i]=Vector{Vector{Int}}(undef, lc+1)
            ub[i]=Vector{UInt16}(undef, lc+1)
            sizes[i]=Vector{UInt16}(undef, lc+1)
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],ub[i],sizes[i],status[i]=get_cblocks!(lc,fsupp,supp[J[i].+1],lt,fbasis[i],gbasis[i],TS=TS,QUIET=true,merge=merge)
        else
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],ub[i],sizes[i],status[i]=get_cblocks!(lc,fsupp,supp[J[i].+1],lt,fbasis[i],gbasis[i],gblocks=blocks[i][2:end],gcl=cl[i][2:end],gblocksize=blocksize[i][2:end],ub=ub[i],sizes=sizes[i],TS=TS,QUIET=true,merge=merge)
        end
    end
    return blocks,cl,blocksize,ub,sizes,fbasis,gbasis,maximum(status)
end

function assign_constraint(m,supp,cliques,cql,cliquesize;assign="first")
    J=[UInt32[] for i=1:cql]
    ncc=UInt32[]
    for i=2:m+1
        rind=copy(supp[i][1])
        for j=2:length(supp[i])
            append!(rind, supp[i][j])
        end
        unique!(rind)
        if assign=="first"
            ind=findfirst(k->issubset(rind, cliques[k]), 1:cql)
            if ind!=nothing
                push!(J[ind], i-1)
            else
                push!(ncc, i-1)
            end
        else
            temp=UInt32[]
            for j=1:cql
                if issubset(rind, cliques[j])
                    push!(temp, j)
                end
            end
            if !isempty(temp)
                if assign=="min"
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

function get_sbasis(var,d)
    n=length(var)
    lb=binomial(n+d,d)
    basis=Vector{Vector{UInt16}}(undef, lb)
    basis[1]=UInt16[]
    i=0
    t=1
    while i<d+1
        t+=1
        if length(basis[t-1])>=i&&basis[t-1][end-i+1:end]==var[n]*ones(UInt16, i)
           if i<d
               basis[t]=var[1]*ones(UInt16, i+1)
           end
           i+=1
        else
            j=bfind(var, n, basis[t-1][1])
            basis[t]=copy(basis[t-1])
            ind=findfirst(x->basis[t][x]!=var[j], 1:length(basis[t]))
            if ind==nothing
                ind=length(basis[t])+1
            end
            if j!=1
                basis[t][1:ind-2]=var[1]*ones(UInt16, ind-2)
            end
            basis[t][ind-1]=var[j+1]
        end
    end
    return basis
end

function get_graph(tsupp::Vector{Vector{UInt16}},basis::Vector{Vector{UInt16}})
    lb=length(basis)
    G=SimpleGraph(lb)
    ltsupp=length(tsupp)
    for i = 1:lb, j = i+1:lb
        bi=sadd(basis[i], basis[j])
        if bfind(tsupp, ltsupp, bi)!=0
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_cgraph(tsupp::Vector{Vector{UInt16}},supp::Vector{Vector{UInt16}},lt::Int,basis::Vector{Vector{UInt16}})
    lb=length(basis)
    ltsupp=length(tsupp)
    G=SimpleGraph(lb)
    for i = 1:lb, j = i+1:lb
        r=1
        while r<=lt
            bi=sadd(sadd(basis[i], supp[r]), basis[j])
            if bfind(tsupp,ltsupp,bi)!=0
               break
            else
                r+=1
            end
        end
        if r<=lt
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_cblocks!(m::Int,tsupp::Vector{Vector{UInt16}},gsupp::Vector{Vector{Vector{UInt16}}},lt::Vector{Int},fbasis::Vector{Vector{UInt16}},gbasis::Vector{Vector{Vector{UInt16}}};gblocks=[],gcl=[],gblocksize=[],ub=[],sizes=[],TS="block",minimize=false,QUIET=true,merge=false)
    if isempty(gblocks)
        gblocks=Vector{Vector{Vector{UInt16}}}(undef, m)
        gblocksize=Vector{Vector{UInt16}}(undef, m)
        gcl=Vector{UInt16}(undef, m)
    end
    if TS==false
        fblocksize=[length(fbasis)]
        fblocks=[[i for i=1:length(fbasis)]]
        fcl=1
        for k=1:m
            gblocks[k]=[[i for i=1:length(gbasis[k])]]
            gblocksize[k]=[length(gbasis[k])]
            gcl[k]=1
        end
        status=0
        nub=fblocksize
        nsizes=[1]
    else
        G=get_graph(tsupp,fbasis)
        if TS=="block"
            fblocks=connected_components(G)
            fblocksize=length.(fblocks)
            fcl=length(fblocksize)
        else
            fblocks,fcl,fblocksize=chordal_cliques!(G, method=TS, minimize=minimize)
            if merge==true
                fblocks,fcl,fblocksize=clique_merge!(fblocks,fcl,QUIET=true)
            end
        end
        nub=unique(fblocksize)
        nsizes=[sum(fblocksize.== i) for i in nub]
        if isempty(ub)||nub!=ub||nsizes!=sizes
            status=1
            for k=1:m
                G=get_cgraph(tsupp,gsupp[k],lt[k],gbasis[k])
                if TS=="block"
                    gblocks[k]=connected_components(G)
                    gblocksize[k]=length.(gblocks[k])
                    gcl[k]=length(gblocksize[k])
                else
                    gblocks[k],gcl[k],gblocksize[k]=chordal_cliques!(G, method=TS, minimize=minimize)
                    if merge==true
                        gblocks[k],gcl[k],gblocksize[k]=clique_merge!(gblocks[k],gcl[k],QUIET=true)
                    end
                end
            end
        else
            status=0
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,status
end

function clique_decomp(n,m,dg,supp;order="min",alg="MF",minimize=false)
    if alg==false
        cliques=[UInt16[i for i=1:n]]
        cql=1
        cliquesize=[n]
    else
        G=SimpleGraph(n)
        for i=1:m+1
            if order=="min"||i==1||order==ceil(Int, dg[i-1]/2)
                for j = 1:length(supp[i])
                    add_clique!(G, unique(supp[i][j]))
                end
            else
                temp=copy(supp[i][1])
                for j=2:length(supp[i])
                    append!(temp, supp[i][j])
                end
                add_clique!(G, unique(temp))
            end
        end
        if alg=="NC"
            cliques,cql,cliquesize=max_cliques(G)
        else
            cliques,cql,cliquesize=chordal_cliques!(G, method=alg, minimize=minimize)
        end
    end
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("------------------------------------------------------")
    return cliques,cql,cliquesize
end

function sadd(a, b)
    c=[a;b]
    return sort!(c)
end

function bfind(A, l, a)
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        if A[mid]==a
           return mid
        elseif A[mid]<a
            low=mid+1
        else
            high=mid-1
        end
    end
    return 0
end
