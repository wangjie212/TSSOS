mutable struct data_type
    n
    nb
    x
    f
    supp
    basis
    coe
    supp0
    ub
    sizes
end

function tssos_first(f,x;nb=0,newton=true,reducebasis=false,TS="block",merge=false,QUIET=false,solve=true,extra_sos=false,solution=false,tol=1e-5)
    n=length(x)
    mon=monomials(f)
    coe=coefficients(f)
    lm=length(mon)
    supp=zeros(UInt8,n,lm)
    for i=1:lm, j=1:n
        @inbounds supp[j,i]=MultivariatePolynomials.degree(mon[i],x[j])
    end
    d=Int(maxdegree(f)/2)
    if newton==true
       if sum(supp[:,end])!=0
          supp=[supp zeros(UInt8,n)]
          coe=[coe;0]
       end
       basis=newton_basis(n,d,supp)
    else
       basis=get_basis(n,d,nb=nb)
    end
    tsupp=[supp bin_add(basis,basis,nb)]
    tsupp=sortslices(tsupp,dims=2)
    tsupp=unique(tsupp,dims=2)
    blocks,cl,blocksize,ub,sizes,_=get_blocks(tsupp,basis,nb=nb,TS=TS,QUIET=QUIET,merge=merge)
    if reducebasis==true
        psupp=[supp zeros(UInt8,n)]
        basis,flag=reducebasis!(psupp,basis,blocks,cl,blocksize,nb=nb)
        if flag==1
            tsupp=[supp bin_add(basis,basis,nb)]
            tsupp=sortslices(tsupp,dims=2)
            tsupp=unique(tsupp,dims=2)
            blocks,cl,blocksize,ub,sizes,_=get_blocks(tsupp,basis,nb=nb,TS=TS,QUIET=QUIET,merge=merge)
        end
    end
    opt,supp0,moment=blockupop(n,supp,coe,basis,blocks,cl,blocksize,nb=nb,QUIET=QUIET,solve=solve,solution=solution,extra_sos=extra_sos)
    if solution==true
        sol=extract_solutions(moment,opt,[f],x,tol=tol)
    else
        sol=nothing
    end
    data=data_type(n,nb,x,f,supp,basis,coe,supp0,ub,sizes)
    return opt,sol,data
end

function tssos_higher!(data::data_type;TS="block",merge=false,QUIET=false,solve=true,extra_sos=false,solution=false,tol=1e-5)
    n=data.n
    nb=data.nb
    x=data.x
    f=data.f
    supp=data.supp
    basis=data.basis
    coe=data.coe
    supp0=data.supp0
    ub=data.ub
    sizes=data.sizes
    supp0=sortslices(supp0,dims=2)
    supp0=unique(supp0,dims=2)
    blocks,cl,blocksize,ub,sizes,status=get_blocks(supp0,basis,ub=ub,sizes=sizes,nb=nb,TS=TS,QUIET=QUIET,merge=merge)
    opt=nothing
    sol=nothing
    if status==1
        opt,supp0,moment=blockupop(n,supp,coe,basis,blocks,cl,blocksize,nb=nb,QUIET=QUIET,solve=solve,solution=solution,extra_sos=extra_sos)
        if solution==true
            sol=extract_solutions(moment,opt,[f],x,tol=tol)
        else
            sol=nothing
        end
    end
    data.supp0=supp0
    data.ub=ub
    data.sizes=sizes
    return opt,sol,data
end

function bin_add(bi,bj,nb)
    bs=bi+bj
    if nb>0
        bs[1:nb,:]=isodd.(bs[1:nb,:])
    end
    return bs
end

function get_basis(n,d;nb=0,lead=[])
    lb=binomial(n+d,d)
    basis=zeros(UInt8,n,lb)
    i=0
    t=1
    while i<d+1
        t+=1
        if basis[n,t-1]==i
           if i<d
              basis[1,t]=i+1
           end
           i+=1
        else
            j=findfirst(x->basis[x,t-1]!=0,1:n)
            basis[:,t]=basis[:,t-1]
            if j==1
               basis[1,t]-=1
               basis[2,t]+=1
            else
               basis[1,t]=basis[j,t]-1
               basis[j,t]=0
               basis[j+1,t]+=1
            end
        end
    end
    if nb>0
        basis_bin=basis[1:nb,:]
        basis_valid=all.(x->x<=1, eachcol(basis_bin))
        basis=basis[:, basis_valid]
    end
    if !isempty(lead)
        basis_valid=map(a->!divide(a, lead, n, size(lead, 2)), eachcol(basis))
        basis=basis[:, basis_valid]
    end
    return basis
end

function divide(a, lead, n, llead)
    return any(j->all(i->lead[i,j]<=a[i], 1:n), 1:llead)
end

function reminder(a,x,gb,n)
    remind=rem(prod(x.^a), gb)
    mon=monomials(remind)
    coe=coefficients(remind)
    lm=length(mon)
    supp=zeros(UInt8,n,lm)
    for i=1:lm, j=1:n
        @inbounds supp[j,i]=MultivariatePolynomials.degree(mon[i],x[j])
    end
    return lm,supp,coe
end

function newton_basis(n,d,supp;e=1e-5)
    lsupp=size(supp,2)
    basis=get_basis(n,d)
    lb=size(basis,2)
    A0=[-1/2*supp' ones(lsupp,1)]
    t=1
    indexb=[i for i=1:lb]
    temp=sortslices(supp,dims=2)
    while t<=lb
          i=indexb[t]
          if bfind(temp,lsupp,UInt8(2)*basis[:,i])!=0
             t=t+1
          else
             model=Model(optimizer_with_attributes(Mosek.Optimizer))
             set_optimizer_attribute(model, MOI.Silent(), true)
             @variable(model, x[1:n+1], lower_bound=-10, upper_bound=10)
             @constraint(model, [A0; [basis[:,i]' -1]]*x .<= zeros(lsupp+1))
             @objective(model, Min, [basis[:,i]' -1]*x)
             optimize!(model)
             vx=value.(x)
             if abs(objective_value(model))<=e&&sum(abs.(vx))<=e
                t=t+1
             else
                if abs(objective_value(model))<=e&&sum(abs.(vx))>e
                   t=t+1
                else
                   lb=lb-1
                   indexb=deleteat!(indexb,t)
                end
                r=t
                while lb>=r
                      j=indexb[r]
                      if [basis[:,j]' -1]*vx<=-e
                         lb=lb-1
                         indexb=deleteat!(indexb,r)
                      else
                         r=r+1
                      end
                end
             end
          end
    end
    return basis[:,indexb]
end

function generate_basis!(supp,basis)
    supp=sortslices(supp,dims=2)
    supp=unique(supp,dims=2)
    lsupp=size(supp,2)
    lb=size(basis,2)
    indexb=UInt32[]
    for i = 1:lb, j = i:lb
        bi=basis[:,i]+basis[:,j]
        if bfind(supp,lsupp,bi)!=0
             push!(indexb,i,j)
        end
    end
    sort!(indexb)
    unique!(indexb)
    return basis[:,indexb]
end

function comp(a, b)
    i=1
    while i<=length(a)
          if a[i]<b[i]
             return -1
          elseif a[i]>b[i]
             return 1
          else
             i+=1
          end
    end
    return 0
end

function bfind(A::Array{T, 2}, l::Int, a::Array{T, 1}) where {T<:Number}
    if l==0
        return 0
    end
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        order=comp(A[:, mid], a)
        if order==0
           return mid
        elseif order<0
           low=mid+1
        else
           high=mid-1
        end
    end
    return 0
end

function bfind(A::Vector{T}, l::Int, a) where {T<:Number}
    if l==0
        return 0
    end
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        order=comp(A[mid], a)
        if order==0
           return mid
        elseif order<0
           low=mid+1
        else
           high=mid-1
        end
    end
    return 0
end

function get_graph(tsupp,basis;nb=0)
    lb=size(basis,2)
    G=SimpleGraph(lb)
    ltsupp=size(tsupp,2)
    for i = 1:lb, j = i+1:lb
        bi=bin_add(basis[:,i],basis[:,j],nb)
        if bfind(tsupp,ltsupp,bi)!=0
           add_edge!(G,i,j)
        end
    end
    return G
end

function get_blocks(tsupp,basis;ub=[],sizes=[],nb=0,TS="block",minimize=false,QUIET=true,merge=false)
    if TS==false
        blocksize=[size(basis,2)]
        blocks=[[i for i=1:size(basis,2)]]
        cl=1
    else
        G=get_graph(tsupp,basis,nb=nb)
        if TS=="block"
            blocks=connected_components(G)
            blocksize=length.(blocks)
            cl=length(blocksize)
        else
            blocks,cl,blocksize=chordal_cliques!(G, method=TS, minimize=minimize)
            if merge==true
                blocks,cl,blocksize=clique_merge!(blocks,cl,QUIET=true)
            end
        end
    end
    nub=unique(blocksize)
    nsizes=[sum(blocksize.== i) for i in nub]
    if isempty(ub)||nub!=ub||nsizes!=sizes
        status=1
        if QUIET==false
            println("------------------------------------------------------")
            println("The sizes of blocks:\n$nub\n$nsizes")
            println("------------------------------------------------------")
        end
    else
        status=0
        if QUIET==false
            println("No higher TSSOS hierarchy!")
        end
    end
    return blocks,cl,blocksize,nub,nsizes,status
end

function blockupop(n,supp,coe,basis,blocks,cl,blocksize;nb=0,QUIET=true,solve=true,solution=false,extra_sos=false)
    tsupp=zeros(UInt8,n,Int(sum(blocksize.^2+blocksize)/2))
    k=1
    for i=1:cl, j=1:blocksize[i], r=j:blocksize[i]
        @inbounds bi=bin_add(basis[:,blocks[i][j]],basis[:,blocks[i][r]],nb)
        @inbounds tsupp[:,k]=bi
        k+=1
    end
    supp0=tsupp
    if extra_sos==true||solution==true
        tsupp=[tsupp get_basis(n,2,nb=nb)]
    end
    tsupp=unique(tsupp,dims=2)
    tsupp=sortslices(tsupp,dims=2)
    objv=nothing
    moment=nothing
    if solve==true
        ltsupp=size(tsupp,2)
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:ltsupp]
        if extra_sos==true||solution==true
            pos0=@variable(model, [1:n+1, 1:n+1], PSD)
            for j=1:n+1, k=j:n+1
                @inbounds bi=bin_add(basis[:,j],basis[:,k],nb)
                Locb=bfind(tsupp,ltsupp,bi)
                if j==k
                   @inbounds cons[Locb]+=pos0[j,k]
                else
                   @inbounds cons[Locb]+=2*pos0[j,k]
                end
            end
        end
        pos=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl)
        for i=1:cl
            bs=blocksize[i]
            if bs==1
               @inbounds pos[i]=@variable(model, lower_bound=0)
               @inbounds bi=bin_add(basis[:,blocks[i][1]],basis[:,blocks[i][1]],nb)
               Locb=bfind(tsupp,ltsupp,bi)
               @inbounds cons[Locb]+=pos[i]
            else
               @inbounds pos[i]=@variable(model, [1:bs, 1:bs], PSD)
               for j=1:blocksize[i], r=j:blocksize[i]
                   @inbounds bi=bin_add(basis[:,blocks[i][j]],basis[:,blocks[i][r]],nb)
                   Locb=bfind(tsupp,ltsupp,bi)
                   if j==r
                       @inbounds cons[Locb]+=pos[i][j,r]
                   else
                       @inbounds cons[Locb]+=2*pos[i][j,r]
                   end
               end
            end
        end
        bc=zeros(ltsupp)
        for i=1:size(supp,2)
            Locb=bfind(tsupp,ltsupp,supp[:,i])
            if Locb==0
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing
            else
               bc[Locb]=coe[i]
            end
        end
        @variable(model, lower)
        cons[1]+=lower
        @constraint(model, con[i=1:ltsupp], cons[i]==bc[i])
        @objective(model, Max, lower)
        optimize!(model)
        status=termination_status(model)
        if status == MOI.OPTIMAL
           objv = objective_value(model)
           println("optimum = $objv")
        else
           objv = objective_value(model)
           println("termination status: $status")
           sstatus=primal_status(model)
           println("solution status: $sstatus")
           println("optimum = $objv")
        end
        if solution==true
            dual_var=-dual.(con)
            moment=zeros(Float64,n+1,n+1)
            for j=1:n+1, k=j:n+1
                bi=bin_add(basis[:,j],basis[:,k],nb)
                Locb=bfind(tsupp,ltsupp,bi)
                moment[j,k]=dual_var[Locb]
            end
            moment=Symmetric(moment,:U)
        end
    end
    return objv,supp0,moment
end

function extract_solutions(moment,opt,pop,x;numeq=0,tol=1e-5)
    n=length(x)
    m=length(pop)-1
    F=eigen(moment, n+1:n+1)
    sol=sqrt(F.values[1])*F.vectors[:,1]
    sol=sol[2:end]/sol[1]
    flag=0
    if abs(opt-polynomial(pop[1])(x => sol))>=tol
        flag=1
    end
    for i=1:m-numeq
        if polynomial(pop[i+1])(x => sol)<=-tol
            flag=1
        end
    end
    for i=m-numeq+1:m
        if abs(polynomial(pop[i+1])(x => sol))>=tol
            flag=1
        end
    end
    if flag==0
        println("Global optimality certified!")
    end
    return sol
end
