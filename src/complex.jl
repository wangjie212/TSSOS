function cs_tssos_first(supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe, n, d, dg; numeq=0, foc=100,
    CS="MF", minimize=true, assign="first", TS="block", QUIET=false, solve=true, solution=false,
    MomentOne=false)
    m=length(supp)-1
    ind=[supp[1][i][1]<=supp[1][i][2] for i=1:length(supp[1])]
    supp[1]=supp[1][ind]
    coe[1]=coe[1][ind]
    cliques,cql,cliquesize=clique_decomp(n,m,dg,supp,order=d,alg=CS,minimize=minimize)
    J,ncc=assign_constraint(m,supp,cliques,cql,cliquesize,assign=assign)
    rlorder=init_order(dg,J,cliquesize,cql,foc=foc,order=d)
    blocks,cl,blocksize,sb,numb,basis,status=get_cblocks_mix!(dg,J,rlorder,m,supp,cliques,cql,cliquesize,TS=TS)
    opt,ksupp=blockcpop_mix(n,m,supp,coe,basis,cliques,cql,cliquesize,J,ncc,blocks,cl,blocksize,numeq=numeq,QUIET=QUIET,TS=TS,solve=solve,solution=solution,MomentOne=MomentOne)
    data=mcpop_data(n,0,m,numeq,supp,coe,dg,basis,rlorder,ksupp,cql,cliques,cliquesize,J,ncc,sb,numb,blocks,cl,blocksize,"Mosek",1e-4,1)
    return opt,data
end

# function cs_tssos_higher!(data;TS="block",QUIET=false,solve=true,solution=false,MomentOne=true)
#     n=data.n
#     m=data.m
#     numeq=data.numeq
#     dg=data.dg
#     supp=data.supp
#     coe=data.coe
#     basis=data.basis
#     rlorder=data.rlorder
#     ksupp=data.ksupp
#     cql=data.cql
#     cliques=data.cliques
#     cliquesize=data.cliquesize
#     J=data.J
#     ncc=data.ncc
#     blocks=data.blocks
#     cl=data.cl
#     blocksize=data.blocksize
#     sb=data.sb
#     numb=data.numb
#     blocks,cl,blocksize,sb,numb,basis,status=get_cblocks_mix!(dg,J,rlorder,m,supp,cliques,cql,cliquesize,tsupp=ksupp,basis=basis,blocks=blocks,cl=cl,blocksize=blocksize,sb=sb,numb=numb,TS=TS)
#     if status==1
#         opt,ksupp=blockcpop_mix(n,m,supp,coe,basis,cliques,cql,cliquesize,J,ncc,blocks,cl,blocksize,numeq=numeq,QUIET=QUIET,solve=solve,solution=solution,MomentOne=MomentOne)
#         data.ksupp=ksupp
#         data.blocks=blocks
#         data.cl=cl
#         data.blocksize=blocksize
#         data.sb=sb
#         data.numb=numb
#     else
#         opt=nothing
#         println("No higher CS-TSSOS hierarchy!")
#     end
#     return opt,data
# end

function blockcpop_mix(n, m, supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe, basis, cliques,
    cql, cliquesize, J, ncc, blocks, cl, blocksize; numeq=0, QUIET=false, TS="block", solve=true,
    solution=false, MomentOne=false)
    tsupp=Vector{Vector{UInt16}}[]
    for i=1:cql
        for j=1:cl[i][1], k=1:blocksize[i][1][j], r=k:blocksize[i][1][j]
            @inbounds bi=[basis[i][1][blocks[i][1][j][k]], basis[i][1][blocks[i][1][j][r]]]
            if bi[1]<=bi[2]
                push!(tsupp, bi)
            else
                push!(tsupp, bi[2:-1:1])
            end
        end
        for (j, w) in enumerate(J[i])
            for l=1:cl[i][j+1], t=1:blocksize[i][j+1][l], r=t:blocksize[i][j+1][l], s=1:length(supp[w+1])
                ind1=blocks[i][j+1][l][t]
                ind2=blocks[i][j+1][l][r]
                @inbounds bi=[sadd(basis[i][j+1][ind1], supp[w+1][s][1]), sadd(basis[i][j+1][ind2], supp[w+1][s][2])]
                if bi[1]<=bi[2]
                    push!(tsupp, bi)
                else
                    push!(tsupp, bi[2:-1:1])
                end
            end
        end
    end
    for i ∈ ncc, j=1:length(supp[i+1])
        if supp[i+1][j][1]<=supp[i+1][j][2]
            push!(tsupp, supp[i+1][j])
        end
    end
    if (MomentOne==true||solution==true)&&TS!=false
        ksupp=copy(tsupp)
        for i=1:cql, j=1:cliquesize[i]
            push!(tsupp, [UInt16[], UInt16[cliques[i][j]]])
            for k=j+1:cliquesize[i]
                push!(tsupp, [UInt16[cliques[i][j]], UInt16[cliques[i][k]]])
            end
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    if (MomentOne==true||solution==true)&&TS!=false
        sort!(ksupp)
        unique!(ksupp)
    else
        ksupp=tsupp
    end
    objv=nothing
    if solve==true
        ltsupp=length(tsupp)
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        rcons=[AffExpr(0) for i=1:ltsupp]
        icons=[AffExpr(0) for i=1:ltsupp]
        for i=1:cql
            if (MomentOne==true||solution==true)&&TS!=false
                bs=cliquesize[i]+1
                pos=@variable(model, [1:2*bs, 1:2*bs], PSD)
                for t=1:bs, r=t:bs
                    @constraint(model, pos[t,r]==pos[t+bs,r+bs])
                    @constraint(model, pos[r,t+bs]+pos[t,r+bs]==0)
                    if t==1&&r==1
                        bi=[UInt16[], UInt16[]]
                    elseif t==1&&r>1
                        bi=[UInt16[], UInt16[cliques[i][r-1]]]
                    else
                        bi=[UInt16[cliques[i][t-1]], UInt16[cliques[i][r-1]]]
                    end
                    Locb=bfind(tsupp, ltsupp, bi)
                    @inbounds add_to_expression!(icons[Locb], pos[t+bs,r])
                    @inbounds add_to_expression!(rcons[Locb], pos[t,r])
                end
            end
            for l=1:cl[i][1]
                @inbounds bs=blocksize[i][1][l]
                if bs==1
                    pos=@variable(model, lower_bound=0)
                    @inbounds bi=[basis[i][1][blocks[i][1][l][1]], basis[i][1][blocks[i][1][l][1]]]
                    Locb=bfind(tsupp,ltsupp,bi)
                    @inbounds add_to_expression!(rcons[Locb], pos)
                else
                    pos=@variable(model, [1:2*bs, 1:2*bs], PSD)
                    for t=1:bs, r=t:bs
                        @constraint(model, pos[t,r]==pos[t+bs,r+bs])
                        @constraint(model, pos[r,t+bs]+pos[t,r+bs]==0)
                        @inbounds ind1=blocks[i][1][l][t]
                        @inbounds ind2=blocks[i][1][l][r]
                        @inbounds bi=[basis[i][1][ind1], basis[i][1][ind2]]
                        if bi[1]<=bi[2]
                            Locb=bfind(tsupp,ltsupp,bi)
                            @inbounds add_to_expression!(icons[Locb], pos[t+bs,r])
                        else
                            Locb=bfind(tsupp,ltsupp,bi[2:-1:1])
                            @inbounds add_to_expression!(icons[Locb], -1, pos[t+bs,r])
                        end
                        @inbounds add_to_expression!(rcons[Locb], pos[t,r])
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
                if supp[i+1][j][1]<=supp[i+1][j][2]
                    Locb=bfind(tsupp,ltsupp,supp[i+1][j])
                    @inbounds add_to_expression!(icons[Locb], imag(coe[i+1][j]), pos)
                    @inbounds add_to_expression!(rcons[Locb], real(coe[i+1][j]), pos)
                end
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
                    @inbounds bi=[sadd(basis[i][j+1][ind], supp[w+1][s][1]), sadd(basis[i][j+1][ind], supp[w+1][s][2])]
                    if bi[1]<=bi[2]
                        Locb=bfind(tsupp,ltsupp,bi)
                        @inbounds add_to_expression!(icons[Locb], imag(coe[w+1][s]), pos)
                        @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s]), pos)
                    end
                end
            else
                if w<=m-numeq
                    pos=@variable(model, [1:2*bs, 1:2*bs], PSD)
                else
                    pos=@variable(model, [1:2*bs, 1:2*bs], Symmetric)
                end
                for t=1:bs, r=t:bs
                    @constraint(model, pos[t,r]==pos[t+bs,r+bs])
                    @constraint(model, pos[r,t+bs]+pos[t,r+bs]==0)
                end
                for t=1:bs, r=1:bs
                    ind1=blocks[i][j+1][l][t]
                    ind2=blocks[i][j+1][l][r]
                    for s=1:length(supp[w+1])
                        @inbounds bi=[sadd(basis[i][j+1][ind1], supp[w+1][s][1]), sadd(basis[i][j+1][ind2], supp[w+1][s][2])]
                        if bi[1]<=bi[2]
                            Locb=bfind(tsupp,ltsupp,bi)
                            @inbounds add_to_expression!(icons[Locb], imag(coe[w+1][s]), pos[t,r])
                            @inbounds add_to_expression!(icons[Locb], real(coe[w+1][s]), pos[t+bs,r])
                            @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s]), pos[t,r])
                            @inbounds add_to_expression!(rcons[Locb], -imag(coe[w+1][s]), pos[t+bs,r])
                        end
                    end
                end
            end
        end
        rbc=zeros(ltsupp)
        ibc=zeros(ltsupp)
        for i=1:length(supp[1])
            Locb=bfind(tsupp,ltsupp,supp[1][i])
            if Locb==0
               @error "The monomial basis is not enough!"
               return nothing,nothing
            else
               rbc[Locb]=real(coe[1][i])
               ibc[Locb]=imag(coe[1][i])
            end
        end
        @variable(model, lower)
        @constraint(model, rcons[2:end].==rbc[2:end])
        @constraint(model, rcons[1]+lower==rbc[1])
        @constraint(model, icons.==ibc)
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
    return objv,ksupp
end

function get_cblocks_mix!(dg, J, rlorder, m, supp::Vector{Vector{Vector{Vector{UInt16}}}}, cliques, cql,
    cliquesize; tsupp=[], basis=[], blocks=[], cl=[], blocksize=[], sb=[], numb=[],
    TS="block", merge=false)
    if isempty(tsupp)
        blocks=Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        cl=Vector{Vector{UInt16}}(undef, cql)
        blocksize=Vector{Vector{Vector{Int}}}(undef, cql)
        sb=Vector{Vector{UInt16}}(undef, cql)
        numb=Vector{Vector{UInt16}}(undef, cql)
        basis=Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        tsupp=copy(supp[1])
        for i=2:m+1, j=1:length(supp[i])
            if supp[i][j][1]<=supp[i][j][2]
                push!(tsupp, supp[i][j])
            end
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
        ind=[issubset(union(tsupp[j][1], tsupp[j][2]), cliques[i]) for j=1:length(tsupp)]
        fsupp=tsupp[ind]
        if flag==1
            basis[i]=Vector{Vector{Vector{UInt16}}}(undef, lc+1)
            basis[i][1]=get_basis(cliques[i], rlorder[i])
            for s=1:lc
                basis[i][s+1]=get_basis(cliques[i], rlorder[i]-ceil(Int, dg[J[i][s]]/2))
            end
            blocks[i]=Vector{Vector{Vector{UInt16}}}(undef, lc+1)
            cl[i]=Vector{UInt16}(undef, lc+1)
            blocksize[i]=Vector{Vector{Int}}(undef, lc+1)
            sb[i]=Vector{UInt16}(undef, lc+1)
            numb[i]=Vector{UInt16}(undef, lc+1)
            blocks[i],cl[i],blocksize[i],sb[i],numb[i],status[i]=get_cblocks!(lc,fsupp,supp[J[i].+1],lt,basis[i],TS=TS,QUIET=true,merge=merge)
        else
            blocks[i],cl[i],blocksize[i],sb[i],numb[i],status[i]=get_cblocks!(lc,fsupp,supp[J[i].+1],lt,basis[i],blocks=blocks[i],cl=cl[i],blocksize=blocksize[i],sb=sb[i],numb=numb[i],TS=TS,QUIET=true,merge=merge)
        end
    end
    return blocks,cl,blocksize,sb,numb,basis,maximum(status)
end

function assign_constraint(m, supp::Vector{Vector{Vector{Vector{UInt16}}}}, cliques, cql, cliquesize;
    assign="first")
    J=[UInt32[] for i=1:cql]
    ncc=UInt32[]
    for i=2:m+1
        rind=copy(supp[i][1][1])
        for j=2:length(supp[i])
            append!(rind, supp[i][j][1])
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

function get_basis(var::Vector{UInt16}, d)
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

function get_graph(tsupp::Vector{Vector{Vector{UInt16}}}, basis; nb=0)
    lb=length(basis)
    G=SimpleGraph(lb)
    ltsupp=length(tsupp)
    for i = 1:lb, j = i+1:lb
        bi=[basis[i], basis[j]]
        if bi[1]>bi[2]
            bi=bi[2:-1:1]
        end
        if bfind(tsupp, ltsupp, bi)!=0
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_cgraph(tsupp::Vector{Vector{Vector{UInt16}}}, supp, lt, basis; nb=0)
    lb=length(basis)
    ltsupp=length(tsupp)
    G=SimpleGraph(lb)
    for i = 1:lb, j = i+1:lb
        r=1
        while r<=lt
            bi=[sadd(basis[i], supp[r][1]), sadd(basis[j], supp[r][2])]
            if bi[1]>bi[2]
                bi=bi[2:-1:1]
            end
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

function clique_decomp(n, m, dg, supp::Vector{Vector{Vector{Vector{UInt16}}}}; order="min", alg="MF",
    minimize=true)
    if alg==false
        cliques=[UInt16[i for i=1:n]]
        cql=1
        cliquesize=[n]
    else
        G=SimpleGraph(n)
        for i=1:m+1
            if order=="min"||i==1||order==ceil(Int, dg[i-1]/2)
                for j = 1:length(supp[i])
                    add_clique!(G, unique([supp[i][j][1];supp[i][j][2]]))
                end
            else
                temp=copy(supp[i][1][1])
                for j=2:length(supp[i])
                    append!(temp, supp[i][j][1])
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
