mutable struct cpop_data
    n # the number of variables
    nb # the number of binary variables
    m # the number of constraints
    numeq # the number of equality constraints
    x # the set of variables
    pop # the polynomial optimization problem
    gb # the Grobner basis
    leadsupp # the leader terms of the Grobner basis
    supp # the support data
    coe # the coefficient data
    lt # the numbers of terms
    dg # the degrees of constraints
    basis # monomial bases
    ksupp # the extending support at the k-th step
    sb # the sizes of different blocks
    numb # the numbers of different blocks
    blocks # the block structure
    cl # the numbers of blocks
    blocksize # the sizes of blocks
    solver # the SDP solver
    tol # the tolerance to certify global optimality
    flag # 0 if global optimality is certified; 1 otherwise
end

"""
    opt,sol,data = tssos_first(pop, x, d; nb=0, numeq=0, quotient=true, basis=nothing,
    reducebasis=false, TS="block", merge=false, solver="Mosek", QUIET=false, solve=true,
    MomentOne=false, solution=false, tol=1e-4)

Compute the first step of the TSSOS hierarchy for constrained polynomial optimization with
relaxation order `d`.
If `quotient=true`, then exploit the quotient ring structure defined by the equality constraints.
Return the optimum, the (near) optimal solution (if `solution=true`) and other data.

# Arguments
- `pop`: the vector of the objective function, inequality constraints, and equality constraints.
- `x`: the set of variables.
- `d`: the relaxation order of the moment-SOS hierarchy.
- `nb`: the number of binary variables in `x`.
- `numeq`: the number of equality constraints.
"""
function tssos_first(pop, x, d; nb=0, numeq=0, quotient=true, basis=nothing, reducebasis=false,
    TS="block", merge=false, solver="Mosek", QUIET=false, solve=true, MomentOne=false,
    solution=false, tol=1e-4)
    n=length(x)
    if quotient==true
        cpop=copy(pop)
        gb=cpop[end-numeq+1:end]
        cpop=cpop[1:end-numeq]
        SemialgebraicSets.gröbnerbasis!(gb)
        cpop[1]=rem(cpop[1], gb)
        lead=leadingmonomial.(gb)
        llead=length(lead)
        leadsupp=zeros(UInt8,n,llead)
        for i=1:llead, j=1:n
            @inbounds leadsupp[j,i]=MultivariatePolynomials.degree(lead[i],x[j])
        end
    else
        cpop=pop
        gb=[]
        leadsupp=[]
    end
    m=length(cpop)-1
    dg=zeros(Int,m)
    coe=Vector{Vector{Float64}}(undef, m+1)
    supp=Vector{Array{UInt8,2}}(undef, m+1)
    lt=Vector{UInt32}(undef, m+1)
    for k=1:m+1
        mon=monomials(cpop[k])
        coe[k]=coefficients(cpop[k])
        lt[k]=length(mon)
        supp[k]=zeros(UInt8,n,lt[k])
        for i=1:lt[k], j=1:n
            @inbounds supp[k][j,i]=MultivariatePolynomials.degree(mon[i],x[j])
        end
    end
    isupp=supp[1]
    for i=2:m+1
        dg[i-1]=maxdegree(pop[i])
        isupp=[isupp supp[i]]
    end
    if basis==nothing
        basis=Vector{Array{UInt8,2}}(undef, m+1)
        basis[1]=get_basis(n,d,nb=nb,lead=leadsupp)
        for k=1:m
            basis[k+1]=get_basis(n,d-Int(ceil(dg[k]/2)),nb=nb,lead=leadsupp)
        end
    end
    tsupp=[isupp bin_add(basis[1],basis[1],nb)]
    tsupp=sortslices(tsupp,dims=2)
    tsupp=unique(tsupp,dims=2)
    blocks,cl,blocksize,sb,numb,_=get_cblocks!(m,tsupp,supp[2:end],lt[2:end],basis,nb=nb,TS=TS,QUIET=QUIET,merge=merge)
    if reducebasis==true
        gsupp=get_gsupp(n,m,lt,supp,basis[2:end],blocks[2:end],cl[2:end],blocksize[2:end],nb=nb)
        psupp=[supp[1] zeros(UInt8,n)]
        psupp=[psupp gsupp]
        basis[1],flag=reducebasis!(psupp,basis[1],blocks[1],cl[1],blocksize[1],nb=nb)
        if flag==1
            tsupp=[isupp bin_add(basis[1],basis[1],nb)]
            tsupp=sortslices(tsupp,dims=2)
            tsupp=unique(tsupp,dims=2)
            blocks,cl,blocksize,sb,numb,_=get_cblocks!(m,tsupp,supp[2:end],lt[2:end],basis,nb=nb,TS=TS,QUIET=QUIET,merge=merge)
        end
    end
    opt,ksupp,moment=blockcpop(n,m,supp,coe,lt,basis,blocks,cl,blocksize,nb=nb,numeq=numeq,gb=gb,x=x,lead=leadsupp,solver=solver,QUIET=QUIET,solve=solve,solution=solution,MomentOne=MomentOne)
    if solution==true
        sol,flag=extract_solutions(moment,opt,pop,x,numeq=numeq,tol=tol)
    else
        sol=nothing
        flag=1
    end
    data=cpop_data(n,nb,m,numeq,x,pop,gb,leadsupp,supp,coe,lt,dg,basis,ksupp,sb,numb,blocks,cl,blocksize,solver,tol,flag)
    return opt,sol,data
end

function tssos_higher!(data::cpop_data; TS="block", merge=false, QUIET=false, solve=true,
    MomentOne=false, solution=false)
    n=data.n
    nb=data.nb
    m=data.m
    numeq=data.numeq
    x=data.x
    pop=data.pop
    gb=data.gb
    leadsupp=data.leadsupp
    supp=data.supp
    coe=data.coe
    lt=data.lt
    dg=data.dg
    basis=data.basis
    ksupp=data.ksupp
    sb=data.sb
    numb=data.numb
    blocks=data.blocks
    cl=data.cl
    blocksize=data.blocksize
    solver=data.solver
    tol=data.tol
    ksupp=sortslices(ksupp,dims=2)
    ksupp=unique(ksupp,dims=2)
    blocks,cl,blocksize,sb,numb,status=get_cblocks!(m,ksupp,supp[2:end],lt[2:end],basis,blocks=blocks,cl=cl,blocksize=blocksize,sb=sb,numb=numb,nb=nb,TS=TS,QUIET=QUIET,merge=merge)
    opt=nothing
    sol=nothing
    if status==1
        opt,ksupp,moment=blockcpop(n,m,supp,coe,lt,basis,blocks,cl,blocksize,nb=nb,numeq=numeq,gb=gb,x=x,lead=leadsupp,solver=solver,QUIET=QUIET,solve=solve,solution=solution,MomentOne=MomentOne)
        if solution==true
            sol,data.flag=extract_solutions(moment,opt,pop,x,numeq=numeq,tol=tol)
        end
    end
    data.ksupp=ksupp
    data.sb=sb
    data.numb=numb
    data.blocks=blocks
    data.cl=cl
    data.blocksize=blocksize
    return opt,sol,data
end

function get_gsupp(n, m, lt, supp, gbasis, blocks, cl, blocksize; nb=0)
    gsupp=zeros(UInt8,n,sum(lt[k+1]*Int(sum(Int.(blocksize[k]).^2+blocksize[k])/2) for k=1:m))
    l=1
    for k=1:m, i=1:cl[k], j=1:blocksize[k][i], r=j:blocksize[k][i], s=1:lt[k+1]
        @inbounds bi=bin_add(gbasis[k][:,blocks[k][i][j]],gbasis[k][:,blocks[k][i][r]],nb)
        @inbounds bi=bin_add(bi,supp[k+1][:,s],nb)
        @inbounds gsupp[:,l]=bi
        l+=1
    end
    return gsupp
end

function reducebasis!(supp, basis, blocks, cl, blocksize; nb=0)
    esupp=supp[:,all.(iseven, eachcol(supp))]
    init,flag,check=0,0,0
    while init==0||check>0
        init,check=1,0
        tsupp=esupp
        for i=1:cl
            if blocksize[i]>1
                for j=1:blocksize[i], r=j+1:blocksize[i]
                    @inbounds bi=bin_add(basis[:,blocks[i][j]],basis[:,blocks[i][r]],nb)
                    tsupp=[tsupp bi]
                end
            end
        end
        tsupp=unique(tsupp,dims=2)
        tsupp=sortslices(tsupp,dims=2)
        ltsupp=size(tsupp,2)
        for i=1:cl
            lo=blocksize[i]
            indexb=[k for k=1:lo]
            j=1
            while lo>=j
                bi=bin_add(basis[:,blocks[i][indexb[j]]],basis[:,blocks[i][indexb[j]]],nb)
                Locb=bfind(tsupp,ltsupp,bi)
                if Locb==0
                   check,flag=1,1
                   deleteat!(indexb,j)
                   lo=lo-1
                else
                   j=j+1
                end
            end
            blocks[i]=blocks[i][indexb]
            blocksize[i]=lo
        end
    end
    if flag==1
       indexb=blocks[1]
       for i=2:cl
           indexb=append!(indexb, blocks[i])
       end
       sort!(indexb)
       unique!(indexb)
       return basis[:,indexb],flag
    else
       return basis,flag
    end
end

function get_cgraph(tsupp::Array{UInt8, 2}, supp::Array{UInt8, 2}, lt, basis::Array{UInt8, 2}; nb=0)
    lb=size(basis,2)
    G=SimpleGraph(lb)
    ltsupp=size(tsupp,2)
    for i = 1:lb, j = i+1:lb
        r=1
        while r<=lt
            bi=bin_add(basis[:,i],basis[:,j],nb)
            bi=bin_add(bi,supp[:,r],nb)
            if bfind(tsupp,ltsupp,bi)!=0
               break
            else
                r+=1
            end
        end
        if r<=lt
           add_edge!(G,i,j)
        end
    end
    return G
end

function get_cblocks!(m, tsupp, supp, lt, basis; blocks=[], cl=[], blocksize=[], sb=[], numb=[], nb=0,
    TS="block", QUIET=true, merge=false)
    if isempty(blocks)
        blocks=Vector{Vector{Vector{UInt16}}}(undef, m+1)
        blocksize=Vector{Vector{UInt16}}(undef, m+1)
        cl=Vector{UInt16}(undef, m+1)
    end
    if TS==false
        blocksize[1]=[size(basis[1],2)]
        blocks[1]=[[i for i=1:size(basis[1],2)]]
        cl[1]=1
        for k=1:m
            blocks[k+1]=[[i for i=1:size(basis[k+1],2)]]
            blocksize[k+1]=[size(basis[k+1],2)]
            cl[k+1]=1
        end
        status=1
        nsb=Int.(blocksize[1])
        nnumb=[1]
        if QUIET==false
            println("------------------------------------------------------")
            println("The sizes of PSD blocks:\n$nsb\n$nnumb")
            println("------------------------------------------------------")
        end
    else
        G=get_graph(tsupp, basis[1], nb=nb)
        if TS=="block"
            blocks[1]=connected_components(G)
            blocksize[1]=length.(blocks[1])
            cl[1]=length(blocksize[1])
        else
            blocks[1],cl[1],blocksize[1]=chordal_cliques!(G, method=TS, minimize=false)
            if merge==true
                blocks[1],cl[1],blocksize[1]=clique_merge!(blocks[1], cl[1], QUIET=true)
            end
        end
        nsb=Int.(unique(blocksize[1]))
        nnumb=[sum(blocksize[1].== i) for i in nsb]
        if isempty(sb)||nsb!=sb||nnumb!=numb
            status=1
            if QUIET==false
                println("------------------------------------------------------")
                println("The sizes of PSD blocks:\n$nsb\n$nnumb")
                println("------------------------------------------------------")
            end
            for k=1:m
                G=get_cgraph(tsupp,supp[k],lt[k],basis[k+1],nb=nb)
                if TS=="block"
                    blocks[k+1]=connected_components(G)
                    blocksize[k+1]=length.(blocks[k+1])
                    cl[k+1]=length(blocksize[k+1])
                else
                    blocks[k+1],cl[k+1],blocksize[k+1]=chordal_cliques!(G, method=TS, minimize=false)
                    if merge==true
                        blocks[k+1],cl[k+1],blocksize[k+1]=clique_merge!(blocks[k+1], cl[k+1], QUIET=true)
                    end
                end
            end
        else
            status=0
            if QUIET==false
               println("No higher TSSOS hierarchy!")
            end
        end
    end
    return blocks,cl,blocksize,nsb,nnumb,status
end

function blockcpop(n, m, supp, coe, lt, basis, blocks, cl, blocksize; nb=0, numeq=0, gb=[],
    x=[], lead=[], solver="Mosek", QUIET=true, solve=true, solution=false, MomentOne=false)
    ksupp=zeros(UInt8, n, Int(sum(Int.(blocksize[1]).^2+blocksize[1])/2))
    k=1
    for i=1:cl[1], j=1:blocksize[1][i], r=j:blocksize[1][i]
        @inbounds bi=bin_add(basis[1][:,blocks[1][i][j]], basis[1][:,blocks[1][i][r]], nb)
        @inbounds ksupp[:,k]=bi
        k+=1
    end
    objv=nothing
    moment=nothing
    if solve==true
        tsupp=ksupp
        if m>0
            gsupp=get_gsupp(n, m, lt, supp, basis[2:end], blocks[2:end], cl[2:end], blocksize[2:end], nb=nb)
            tsupp=[tsupp gsupp]
        end
        if MomentOne==true||solution==true
            tsupp=[tsupp get_basis(n,2,nb=nb)]
        end
        tsupp=unique(tsupp,dims=2)
        if !isempty(gb)
            nsupp=zeros(UInt8, n)
            llead=size(lead, 2)
            for i=1:size(tsupp, 2)
                if divide(tsupp[:,i], lead, n, llead)
                    _,temp,_=reminder(tsupp[:,i],x,gb,n)
                    nsupp=[nsupp temp]
                else
                    nsupp=[nsupp tsupp[:,i]]
                end
            end
            tsupp=nsupp
        end
        tsupp=sortslices(tsupp,dims=2)
        tsupp=unique(tsupp,dims=2)
        ltsupp=size(tsupp,2)
        if solver=="Mosek"
            model=Model(optimizer_with_attributes(Mosek.Optimizer))
        elseif solver=="SDPT3"
            model=Model(optimizer_with_attributes(SDPT3.Optimizer))
        else
            @error "The solver is currently not supported!"
            return nothing,nothing,nothing
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:ltsupp]
        pos=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1])
        for i=1:cl[1]
            if MomentOne==true||solution==true
                pos0=@variable(model, [1:n+1, 1:n+1], PSD)
                for j=1:n+1, k=j:n+1
                    @inbounds bi=bin_add(basis[1][:,j],basis[1][:,k],nb)
                    if !isempty(gb)&&divide(bi, lead, n, llead)
                        bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                        for l=1:bi_lm
                            Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                            if j==k
                               @inbounds add_to_expression!(cons[Locb], bi_coe[l], pos0[j,k])
                            else
                               @inbounds add_to_expression!(cons[Locb], 2*bi_coe[l], pos0[j,k])
                            end
                        end
                    else
                        Locb=bfind(tsupp,ltsupp,bi)
                        if j==k
                           @inbounds add_to_expression!(cons[Locb], pos0[j,k])
                        else
                           @inbounds add_to_expression!(cons[Locb], 2, pos0[j,k])
                        end
                    end
                end
            end
            bs=blocksize[1][i]
            if bs==1
               @inbounds pos[i]=@variable(model, lower_bound=0)
               @inbounds bi=bin_add(basis[1][:,blocks[1][i][1]],basis[1][:,blocks[1][i][1]],nb)
               if !isempty(gb)&&divide(bi, lead, n, llead)
                   bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                   for l=1:bi_lm
                       Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                       @inbounds add_to_expression!(cons[Locb], bi_coe[l], pos[i])
                   end
               else
                   Locb=bfind(tsupp,ltsupp,bi)
                   @inbounds add_to_expression!(cons[Locb], pos[i])
               end
            else
               @inbounds pos[i]=@variable(model, [1:bs, 1:bs], PSD)
               for j=1:bs, r=j:bs
                   @inbounds bi=bin_add(basis[1][:,blocks[1][i][j]],basis[1][:,blocks[1][i][r]],nb)
                   if !isempty(gb)&&divide(bi, lead, n, llead)
                       bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                       for l=1:bi_lm
                           Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                           if j==r
                              @inbounds add_to_expression!(cons[Locb], bi_coe[l], pos[i][j,r])
                           else
                              @inbounds add_to_expression!(cons[Locb], 2*bi_coe[l], pos[i][j,r])
                           end
                       end
                   else
                       Locb=bfind(tsupp,ltsupp,bi)
                       if j==r
                          @inbounds add_to_expression!(cons[Locb], pos[i][j,r])
                       else
                          @inbounds add_to_expression!(cons[Locb], 2, pos[i][j,r])
                       end
                   end
               end
            end
        end
        gpos=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m)
        for k=1:m
            gpos[k]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k])
            for i=1:cl[k+1]
                bs=blocksize[k+1][i]
                if bs==1
                    if k<=m-numeq
                        gpos[k][i]=@variable(model, lower_bound=0)
                    else
                        gpos[k][i]=@variable(model)
                    end
                    for s=1:lt[k+1]
                        @inbounds bi=bin_add(basis[k+1][:,blocks[k+1][i][1]],basis[k+1][:,blocks[k+1][i][1]],nb)
                        @inbounds bi=bin_add(bi,supp[k+1][:,s],nb)
                        if !isempty(gb)&&divide(bi, lead, n, llead)
                            bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                            for l=1:bi_lm
                                Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                                @inbounds add_to_expression!(cons[Locb], coe[k+1][s]*bi_coe[l], gpos[k][i])
                            end
                        else
                            Locb=bfind(tsupp,ltsupp,bi)
                            @inbounds add_to_expression!(cons[Locb], coe[k+1][s], gpos[k][i])
                        end
                    end
                else
                    if k<=m-numeq
                       gpos[k][i]=@variable(model, [1:bs, 1:bs], PSD)
                    else
                       gpos[k][i]=@variable(model, [1:bs, 1:bs], Symmetric)
                    end
                    for j=1:bs, r=j:bs, s=1:lt[k+1]
                        @inbounds bi=bin_add(basis[k+1][:,blocks[k+1][i][j]],basis[k+1][:,blocks[k+1][i][r]],nb)
                        @inbounds bi=bin_add(bi,supp[k+1][:,s],nb)
                        if !isempty(gb)&&divide(bi, lead, n, llead)
                            bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                            for l=1:bi_lm
                                Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                                if j==r
                                   @inbounds add_to_expression!(cons[Locb], coe[k+1][s]*bi_coe[l], gpos[k][i][j,r])
                                else
                                   @inbounds add_to_expression!(cons[Locb], 2*coe[k+1][s]*bi_coe[l], gpos[k][i][j,r])
                                end
                            end
                        else
                            Locb=bfind(tsupp,ltsupp,bi)
                            if j==r
                               @inbounds add_to_expression!(cons[Locb], coe[k+1][s], gpos[k][i][j,r])
                            else
                               @inbounds add_to_expression!(cons[Locb], 2*coe[k+1][s], gpos[k][i][j,r])
                            end
                        end
                    end
                end
            end
        end
        bc=zeros(ltsupp)
        for i=1:lt[1]
            Locb=bfind(tsupp,ltsupp,supp[1][:,i])
            if Locb==0
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing
            else
               bc[Locb]=coe[1][i]
           end
        end
        @variable(model, lower)
        cons[1]+=lower
        @constraint(model, con[i=1:ltsupp], cons[i]==bc[i])
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
        if solution==true
            dual_var=-dual.(con)
            moment=zeros(Float64,n+1,n+1)
            for j=1:n+1, k=j:n+1
                bi=bin_add(basis[1][:,j],basis[1][:,k],nb)
                if !isempty(gb)&&divide(bi, lead, n, llead)
                    bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                    moment[j,k]=0
                    for l=1:bi_lm
                        Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                        moment[j,k]+=bi_coe[l]*dual_var[Locb]
                    end
                else
                    Locb=bfind(tsupp,ltsupp,bi)
                    moment[j,k]=dual_var[Locb]
                end
            end
            moment=Symmetric(moment,:U)
        end
    end
    return objv,ksupp,moment
end
