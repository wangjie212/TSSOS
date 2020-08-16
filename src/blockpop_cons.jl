mutable struct cdata_type
    n
    nb
    m
    numeq
    x
    pop
    gb
    leadsupp
    ssupp
    coe
    lt
    d
    dg
    fbasis
    gbasis
    fsupp
    ub
    sizes
    gblocks
    gcl
    gblocksize
end

function tssos_first(pop,x,d;nb=0,numeq=0,quotient=true,basis=nothing,reducebasis=false,TS="block",minimize=false,merge=false,QUIET=false,solve=true,extra_sos=false,solution=false,tol=1e-5)
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
        m=length(cpop)-1
    else
        m=length(pop)-1
        gb=[]
        leadsupp=[]
    end
    dg=zeros(Int,m)
    coe=Vector{Vector{Float64}}(undef, m+1)
    ssupp=Vector{Array{UInt8,2}}(undef, m+1)
    lt=Vector{UInt32}(undef,m+1)
    for k=1:m+1
        mon=monomials(cpop[k])
        coe[k]=coefficients(cpop[k])
        lt[k]=length(mon)
        ssupp[k]=zeros(UInt8,n,lt[k])
        for i=1:lt[k], j=1:n
            @inbounds ssupp[k][j,i]=MultivariatePolynomials.degree(mon[i],x[j])
        end
    end
    supp=ssupp[1]
    for i=2:m+1
        dg[i-1]=maxdegree(pop[i])
        supp=[supp ssupp[i]]
    end
    if basis==nothing
        fbasis=get_basis(n,d,nb=nb,lead=leadsupp)
        gbasis=Vector{Array{UInt8,2}}(undef,m)
        for k=1:m
            gbasis[k]=get_basis(n,d-Int(ceil(dg[k]/2)),nb=nb,lead=leadsupp)
        end
    else
        fbasis=basis[1]
        gbasis=basis[2:end]
    end
    tsupp=[supp bin_add(fbasis,fbasis,nb)]
    tsupp=sortslices(tsupp,dims=2)
    tsupp=unique(tsupp,dims=2)
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,_=get_cblocks!(m,tsupp,ssupp,lt,fbasis,gbasis,nb=nb,TS=TS,minimize=minimize,QUIET=QUIET,merge=merge)
    if reducebasis==true
        gsupp=get_gsupp(n,m,lt,ssupp,gbasis,gblocks,gcl,gblocksize,nb=nb)
        psupp=[ssupp[1] zeros(UInt8,n)]
        psupp=[psupp gsupp]
        fbasis,flag=reducebasis!(psupp,fbasis,fblocks,fcl,fblocksize,nb=nb)
        if flag==1
            tsupp=[supp bin_add(fbasis,fbasis,nb)]
            tsupp=sortslices(tsupp,dims=2)
            tsupp=unique(tsupp,dims=2)
            fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,_=get_cblocks!(m,tsupp,ssupp,lt,fbasis,gbasis,nb=nb,TS=TS,minimize=minimize,QUIET=QUIET,merge=merge)
        end
    end
    opt,fsupp,moment=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nb=nb,numeq=numeq,gb=gb,x=x,lead=leadsupp,QUIET=QUIET,solve=solve,solution=solution,extra_sos=extra_sos)
    if solution==true
        sol=extract_solutions(moment,opt,pop,x,numeq=numeq,tol=tol)
    else
        sol=nothing
    end
    data=cdata_type(n,nb,m,numeq,x,cpop,gb,leadsupp,ssupp,coe,lt,d,dg,fbasis,gbasis,fsupp,ub,sizes,gblocks,gcl,gblocksize)
    return opt,sol,data
end

function tssos_higher!(data::cdata_type;TS="block",minimize=false,merge=false,QUIET=false,solve=true,extra_sos=false,solution=false,tol=1e-5)
    n=data.n
    nb=data.nb
    m=data.m
    numeq=data.numeq
    x=data.x
    pop=data.pop
    gb=data.gb
    leadsupp=data.leadsupp
    ssupp=data.ssupp
    coe=data.coe
    lt=data.lt
    d=data.d
    dg=data.dg
    fbasis=data.fbasis
    gbasis=data.gbasis
    fsupp=data.fsupp
    ub=data.ub
    sizes=data.sizes
    gblocks=data.gblocks
    gcl=data.gcl
    gblocksize=data.gblocksize
    fsupp=sortslices(fsupp,dims=2)
    fsupp=unique(fsupp,dims=2)
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_cblocks!(m,fsupp,ssupp,lt,fbasis,gbasis,gblocks=gblocks,gcl=gcl,gblocksize=gblocksize,ub=ub,sizes=sizes,nb=nb,TS=TS,minimize=minimize,QUIET=QUIET,merge=merge)
    opt=nothing
    sol=nothing
    if status==1
        opt,fsupp,moment=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nb=nb,numeq=numeq,gb=gb,x=x,lead=leadsupp,QUIET=QUIET,solve=solve,solution=solution,extra_sos=extra_sos)
        if solution==true
            sol=extract_solutions(moment,opt,pop,x,numeq=numeq,tol=tol)
        end
    end
    data.fsupp=fsupp
    data.ub=ub
    data.sizes=sizes
    data.gblocks=gblocks
    data.gcl=gcl
    data.gblocksize=gblocksize
    return opt,sol,data
end

function get_gsupp(n,m,lt,ssupp,gbasis,gblocks,gcl,gblocksize;nb=0)
    gsupp=zeros(UInt8,n,sum(lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2) for k=1:m))
    l=1
    for k=1:m, i=1:gcl[k], j=1:gblocksize[k][i], r=j:gblocksize[k][i], s=1:lt[k+1]
        @inbounds bi=bin_add(gbasis[k][:,gblocks[k][i][j]],gbasis[k][:,gblocks[k][i][r]],nb)
        @inbounds bi=bin_add(bi,ssupp[k+1][:,s],nb)
        @inbounds gsupp[:,l]=bi
        l+=1
    end
    return gsupp
end

function reducebasis!(supp,basis,blocks,cl,blocksize;nb=0)
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

function get_cgraph(tsupp,ssupp,lt,basis;nb=0)
    lb=size(basis,2)
    G=SimpleGraph(lb)
    ltsupp=size(tsupp,2)
    for i = 1:lb, j = i+1:lb
        r=1
        while r<=lt
            bi=bin_add(basis[:,i],basis[:,j],nb)
            bi=bin_add(bi,ssupp[:,r],nb)
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

function get_cblocks!(m,tsupp,ssupp,lt,fbasis,gbasis;gblocks=[],gcl=[],gblocksize=[],ub=[],sizes=[],nb=0,TS="block",minimize=false,QUIET=true,merge=false)
    if isempty(gblocks)
        gblocks=Vector{Vector{Vector{UInt16}}}(undef,m)
        gblocksize=Vector{Vector{UInt16}}(undef, m)
        gcl=Vector{UInt16}(undef,m)
    end
    if TS==false
        fblocksize=[size(fbasis,2)]
        fblocks=[[i for i=1:size(fbasis,2)]]
        fcl=1
        for k=1:m
            gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
            gblocksize[k]=[size(gbasis[k],2)]
            gcl[k]=1
        end
        status=1
        nub=fblocksize
        nsizes=[1]
        if QUIET==false
            println("------------------------------------------------------")
            println("The sizes of blocks:\n$nub\n$nsizes")
            println("------------------------------------------------------")
        end
    elseif TS=="block"
        G=get_graph(tsupp,fbasis,nb=nb)
        fblocks=connected_components(G)
        fblocksize=length.(fblocks)
        fcl=length(fblocksize)
        nub=unique(fblocksize)
        nsizes=[sum(fblocksize.== i) for i in nub]
        if isempty(ub)||nub!=ub||nsizes!=sizes
            status=1
            if QUIET==false
                println("------------------------------------------------------")
                println("The sizes of blocks:\n$nub\n$nsizes")
                println("------------------------------------------------------")
            end
            for k=1:m
                G=get_cgraph(tsupp,ssupp[k+1],lt[k+1],gbasis[k],nb=nb)
                gblocks[k]=connected_components(G)
                gblocksize[k]=length.(gblocks[k])
                gcl[k]=length(gblocksize[k])
            end
        else
            status=0
            if QUIET==false
               println("No higher TSSOS hierarchy!")
            end
        end
    else
        G=get_graph(tsupp,fbasis,nb=nb)
        fblocks,fcl,fblocksize=chordal_cliques!(G, method=TS, minimize=minimize)
        if merge==true
            fblocks,fcl,fblocksize=clique_merge!(fblocks,fcl,QUIET=true)
        end
        fblocksize=length.(fblocks)
        fcl=length(fblocksize)
        nub=unique(fblocksize)
        nsizes=[sum(fblocksize.== i) for i in nub]
        if isempty(ub)||nub!=ub||nsizes!=sizes
            status=1
            if QUIET==false
                println("------------------------------------------------------")
                println("The sizes of blocks:\n$nub\n$nsizes")
                println("------------------------------------------------------")
            end
            for k=1:m
                G=get_cgraph(tsupp,ssupp[k+1],lt[k+1],gbasis[k],nb=nb)
                gblocks[k],gcl[k],gblocksize[k]=chordal_cliques!(G, method=TS, minimize=minimize)
                if merge==true
                    gblocks[k],gcl[k],gblocksize[k]=clique_merge!(gblocks[k],gcl[k],QUIET=true)
                end
            end
        else
            status=0
            if QUIET==false
               println("No higher TSSOS hierarchy!")
            end
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,status
end

function blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize;nb=0,numeq=0,gb=[],x=[],lead=[],QUIET=true,solve=true,solution=false,extra_sos=false)
    fsupp=zeros(UInt8,n,Int(sum(fblocksize.^2+fblocksize)/2))
    k=1
    for i=1:fcl, j=1:fblocksize[i], r=j:fblocksize[i]
        @inbounds bi=bin_add(fbasis[:,fblocks[i][j]],fbasis[:,fblocks[i][r]],nb)
        @inbounds fsupp[:,k]=bi
        k+=1
    end
    objv=nothing
    moment=nothing
    if solve==true
        gsupp=get_gsupp(n,m,lt,ssupp,gbasis,gblocks,gcl,gblocksize,nb=nb)
        tsupp=[fsupp gsupp]
        if extra_sos==true||solution==true
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
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:ltsupp]
        pos=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, fcl)
        for i=1:fcl
            if extra_sos==true||solution==true
                pos0=@variable(model, [1:n+1, 1:n+1], PSD)
                for j=1:n+1, k=j:n+1
                    @inbounds bi=bin_add(fbasis[:,j],fbasis[:,k],nb)
                    if !isempty(gb)&&divide(bi, lead, n, llead)
                        bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                        for l=1:bi_lm
                            Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                            if j==k
                               @inbounds cons[Locb]+=bi_coe[l]*pos0[j,k]
                            else
                               @inbounds cons[Locb]+=2*bi_coe[l]*pos0[j,k]
                            end
                        end
                    else
                        Locb=bfind(tsupp,ltsupp,bi)
                        if j==k
                           @inbounds cons[Locb]+=pos0[j,k]
                        else
                           @inbounds cons[Locb]+=2*pos0[j,k]
                        end
                    end
                end
            end
            bs=fblocksize[i]
            if bs==1
               @inbounds pos[i]=@variable(model, lower_bound=0)
               @inbounds bi=bin_add(fbasis[:,fblocks[i][1]],fbasis[:,fblocks[i][1]],nb)
               if !isempty(gb)&&divide(bi, lead, n, llead)
                   bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                   for l=1:bi_lm
                       Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                       @inbounds cons[Locb]+=bi_coe[l]*pos[i]
                   end
               else
                   Locb=bfind(tsupp,ltsupp,bi)
                   @inbounds cons[Locb]+=pos[i]
               end
            else
               @inbounds pos[i]=@variable(model, [1:bs, 1:bs], PSD)
               for j=1:bs, r=j:bs
                   @inbounds bi=bin_add(fbasis[:,fblocks[i][j]],fbasis[:,fblocks[i][r]],nb)
                   if !isempty(gb)&&divide(bi, lead, n, llead)
                       bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                       for l=1:bi_lm
                           Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                           if j==r
                              @inbounds cons[Locb]+=bi_coe[l]*pos[i][j,r]
                           else
                              @inbounds cons[Locb]+=2*bi_coe[l]*pos[i][j,r]
                           end
                       end
                   else
                       Locb=bfind(tsupp,ltsupp,bi)
                       if j==r
                          @inbounds cons[Locb]+=pos[i][j,r]
                       else
                          @inbounds cons[Locb]+=2*pos[i][j,r]
                       end
                   end
               end
            end
        end
        gpos=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m)
        for k=1:m
            gpos[k]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, gcl[k])
            for i=1:gcl[k]
                bs=gblocksize[k][i]
                if bs==1
                    if k<=m-numeq
                        gpos[k][i]=@variable(model, lower_bound=0)
                    else
                        gpos[k][i]=@variable(model)
                    end
                    for s=1:lt[k+1]
                        @inbounds bi=bin_add(gbasis[k][:,gblocks[k][i][1]],gbasis[k][:,gblocks[k][i][1]],nb)
                        @inbounds bi=bin_add(bi,ssupp[k+1][:,s],nb)
                        if !isempty(gb)&&divide(bi, lead, n, llead)
                            bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                            for l=1:bi_lm
                                Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                                @inbounds cons[Locb]+=coe[k+1][s]*bi_coe[l]*gpos[k][i]
                            end
                        else
                            Locb=bfind(tsupp,ltsupp,bi)
                            @inbounds cons[Locb]+=coe[k+1][s]*gpos[k][i]
                        end
                    end
                else
                    if k<=m-numeq
                       gpos[k][i]=@variable(model, [1:bs, 1:bs], PSD)
                    else
                       gpos[k][i]=@variable(model, [1:bs, 1:bs], Symmetric)
                    end
                    for j=1:bs, r=j:bs, s=1:lt[k+1]
                        @inbounds bi=bin_add(gbasis[k][:,gblocks[k][i][j]],gbasis[k][:,gblocks[k][i][r]],nb)
                        @inbounds bi=bin_add(bi,ssupp[k+1][:,s],nb)
                        if !isempty(gb)&&divide(bi, lead, n, llead)
                            bi_lm,bi_supp,bi_coe=reminder(bi,x,gb,n)
                            for l=1:bi_lm
                                Locb=bfind(tsupp,ltsupp,bi_supp[:,l])
                                if j==r
                                   @inbounds cons[Locb]+=coe[k+1][s]*bi_coe[l]*gpos[k][i][j,r]
                                else
                                   @inbounds cons[Locb]+=2*coe[k+1][s]*bi_coe[l]*gpos[k][i][j,r]
                                end
                            end
                        else
                            Locb=bfind(tsupp,ltsupp,bi)
                            if j==r
                               @inbounds cons[Locb]+=coe[k+1][s]*gpos[k][i][j,r]
                            else
                               @inbounds cons[Locb]+=2*coe[k+1][s]*gpos[k][i][j,r]
                            end
                        end
                    end
                end
            end
        end
        bc=zeros(ltsupp)
        for i=1:lt[1]
            Locb=bfind(tsupp,ltsupp,ssupp[1][:,i])
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
        if  status==MOI.OPTIMAL
            objv=objective_value(model)
            println("optimum = $objv")
        else
            objv=objective_value(model)
            println("termination status: $status")
            sstatus=primal_status(model)
            println("solution status: $sstatus")
            println("optimum = $objv")
        end
        if solution==true
            dual_var=-dual.(con)
            moment=zeros(Float64,n+1,n+1)
            for j=1:n+1, k=j:n+1
                bi=bin_add(fbasis[:,j],fbasis[:,k],nb)
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
    return objv,fsupp,moment
end
