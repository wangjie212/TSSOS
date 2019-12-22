mutable struct data_type
    supp
    basis
    coe
    supp1
    ub
    sizes
end

function blockupop_first(n,pd,supp,coe;newton=1,method="block",reducebasis=0)
d=Int(pd/2)
if newton==1
   if sum(supp[:,end])!=0
      supp=[supp zeros(UInt8,n,1)]
      coe=[coe 0]
   end
   basis=newton_basis(n,d,supp)
else
   basis=get_basis(n,d)
end
if method=="block"&&reducebasis==0
blocks,cl,blocksize,ub,sizes,status=get_blocks(n,supp,basis)
elseif method=="block"&&reducebasis==1
    flag=1
    while flag==1
          blocks,cl,blocksize,ub,sizes,status=get_blocks(n,supp,basis,reduce=1)
          tsupp=[supp zeros(UInt8,n,1)]
          basis,flag=reducebasis!(n,tsupp,basis,blocks,cl,blocksize)
    end
elseif method=="clique"&&reducebasis==0
blocks,cl,blocksize,ub,sizes,status=get_cliques(n,supp,basis)
else
    flag=1
    while flag==1
          blocks,cl,blocksize,ub,sizes,status=get_cliques(n,supp,basis,reduce=1)
          tsupp=[supp zeros(UInt8,n,1)]
          basis,flag=reducebasis!(n,tsupp,basis,blocks,cl,blocksize)
    end
end
if status==1
   opt,supp1=blockupop(n,supp,coe,basis,blocks,cl,blocksize)
else
   opt=0
   supp1=0
end
data=data_type(supp,basis,coe,supp1,ub,sizes)
return opt,data,status
end

function blockupop_higher!(n,data;method="block",reducebasis=0)
supp=data.supp
basis=data.basis
coe=data.coe
supp1=data.supp1
ub=data.ub
sizes=data.sizes
opt=0
if method=="block"&&reducebasis==0
blocks,cl,blocksize,ub,sizes,status=get_hblocks!(n,supp1,basis,ub,sizes)
elseif method=="block"&&reducebasis==1
    flag=1
    while flag==1
          blocks,cl,blocksize,ub,sizes,status=get_hblocks!(n,supp1,basis,ub,sizes,redcue=1)
          tsupp=[supp zeros(UInt8,n,1)]
          basis,flag=reducebasis!(n,tsupp,basis,blocks,cl,blocksize)
    end
elseif method=="clique"&&reducebasis==0
blocks,cl,blocksize,ub,sizes,status=get_hcliques!(n,supp1,basis,ub,sizes)
else
    flag=1
    while flag==1
          blocks,cl,blocksize,ub,sizes,status=get_hcliques!(n,supp1,basis,ub,sizes,reduce=1)
          tsupp=[supp zeros(UInt8,n,1)]
          basis,flag=reducebasis!(n,tsupp,basis,blocks,cl,blocksize)
    end
end
if status==1
   opt,supp1=blockupop(n,supp,coe,basis,blocks,cl,blocksize)
end
data.supp1=supp1
data.ub=ub
data.sizes=sizes
return opt,data,status
end

function get_basis(n,d)
lb=binomial(n+d,d)
basis=zeros(UInt8,n,lb)
i=UInt16(0)
t=UInt16(1)
while i<d+1
    if basis[n,t]==i
       if i<d
          t=t+1
          basis[1,t]=i+1
          i=i+1
       else i=i+1
       end
    else j=1
         while basis[j,t]==0
               j=j+1
         end
         if j==1
            t=t+1
            basis[:,t]=basis[:,t-1]
            basis[1,t]=basis[1,t]-1
            basis[2,t]=basis[2,t]+1
         else t=t+1
              basis[:,t]=basis[:,t-1]
              basis[1,t]=basis[j,t]-1
              basis[j,t]=0
              basis[j+1,t]=basis[j+1,t]+1
         end
    end
end
return basis
end

function newton_basis(n,d,supp)
lsupp=size(supp,2)
lb=binomial(n+d,d)
basis=get_basis(n,d)
A0=[-1/2*supp' ones(lsupp,1)]
t=1
indexb=[i for i=1:lb]
temp=sortslices(supp,dims=2)
while t<=lb
      i=indexb[t]
      if bfind(temp,lsupp,2*basis[:,i],n)!=0
         t=t+1
      else
         model=Model(with_optimizer(Mosek.Optimizer, QUIET=true))
         @variable(model, x[1:n+1])
         @constraint(model, con, [A0; [basis[:,i]' -1]]*x.<=zeros(lsupp+1,1))
         @objective(model, Min, [basis[:,i]' -1]*x)
         optimize!(model)
         if termination_status(model) == MOI.OPTIMAL
            t=t+1
         else
            vx = value.(x)
            lb=lb-1
            indexb=deleteat!(indexb,t)
            r=t
            while lb>=r
                  j=indexb[r]
                  if [basis[:,j]' -1]*vx<=-0.001&&bfind(temp,lsupp,2*basis[:,i],n)==0
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

function generate_basis!(n,supp,basis)
supp=sortslices(supp,dims=2)
lsupp=size(supp,2)
lb=size(basis,2)
indexb=[0]
for i = 1:lb
    for j = i:lb
        bi=basis[:,i]+basis[:,j]
         if bfind(supp,lsupp,bi,n)!=0
            indexb=[indexb i j]
         end
    end
end
indexb=indexb[2:end]
indexb=sort(indexb)
indexb=unique(indexb)
return basis[:,indexb]
end

function odd_supp(n,supp)
lo=size(supp,2)
indexb=[i for i=1:lo]
i=Int(1)
while lo>=i
      bi=supp[:,indexb[i]]
      if sum(Int[iseven(bi[j]) for j=1:n])==n
         deleteat!(indexb,i)
         lo=lo-1
      else
         i=i+1
      end
end
return supp[:,indexb]
end

function even_supp(n,supp)
lo=size(supp,2)
indexb=[i for i=1:lo]
i=Int(1)
while lo>=i
      bi=supp[:,indexb[i]]
      if sum(Int[iseven(bi[j]) for j=1:n])<n
         deleteat!(indexb,i)
         lo=lo-1
      else
         i=i+1
      end
end
return supp[:,indexb]
end

function comp(a,b,n)
    i=Int(1)
    while i<=n
          if a[i]<b[i]
             return -1
          elseif a[i]>b[i]
             return 1
          else
             i+=1
          end
    end
    if i==n+1
       return 0
    end
end

function bfind(A,l,a,n)
    if l==0
        return 0
    end
    low=Int(1)
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        order=comp(A[:,mid],a,n)
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

function cliquesFromSpMatD(A)
ms=MSession()
mat"lb = size($A,1);
A = spones($A) + (2*lb+1)*speye(lb);
I = amd(A);
R = chol(A(I,I));
Cliques = spones(R);
[value,orig_idx] = sort(I);
remainIdx = 1;
for i=2:lb
    idx = i:lb;
    one = find(Cliques(i,idx));
    noOfone = length(one);
    cliqueResult = sum(Cliques(remainIdx,idx(one)),2);
    if isempty(find(cliqueResult == noOfone,1))
       remainIdx = [remainIdx;i];
    end
end
cSet = Cliques(remainIdx,orig_idx);
cl = length(remainIdx);
[Elem,~] = find(cSet');
NoElem = full(sum(cSet,2));"
cl=jscalar(get_mvariable(ms,:cl))
Elem=jarray(get_mvariable(ms,:Elem))
blocksize=jarray(get_mvariable(ms,:NoElem))
cl=convert(UInt16,cl)
Elem=convert(Array{UInt16},Elem)
blocksize=convert(Array{Int},blocksize)
blocks=Array{Any}(undef,cl)
blocks[1]=Elem[1:blocksize[1]]
for i=2:cl
    idx=sum(blocksize[1:i-1])
    blocks[i]=Elem[idx+1:idx+blocksize[i]]
end
return blocks,cl,blocksize
end

function get_cliques(n,supp,basis;reduce=0)
if reduce==1
supp1=[supp 2*basis]
supp1=sortslices(supp1,dims=2)
supp1=unique(supp1,dims=2)
lsupp1=size(supp1,2)
lb=size(basis,2)
A=zeros(UInt8,lb,lb)
for i = 1:lb
    for j = i:lb
        bi=basis[:,i]+basis[:,j]
         if bfind(supp1,lsupp1,bi,n)!=0
           A[i,j]=1
           A[j,i]=1
        end
    end
end
else
    osupp=odd_supp(n,supp)
    osupp=sortslices(osupp,dims=2)
    lo=size(osupp,2)
    lb=size(basis,2)
    A=zeros(UInt8,lb,lb)
    for i = 1:lb
        for j = i:lb
            bi=basis[:,i]+basis[:,j]
            if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(osupp,lo,bi,n)!=0
               A[i,j]=1
               A[j,i]=1
            end
        end
    end
end
blocks,cl,blocksize=cliquesFromSpMatD(A)
if cl==1
   println("Unblockable")
   return blocks,cl,blocksize,blocksize,[1],0
else
   ub=unique(blocksize)
   sizes=[sum(blocksize.== i) for i in ub]
   println("blocksizes:\n$ub\n$sizes")
   return blocks,cl,blocksize,ub,sizes,1
end
end

function get_hcliques!(n,supp,basis,ub,sizes;reduce=0)
if reduce==1
supp1=[supp 2*basis]
supp1=sortslices(supp1,dims=2)
supp1=unique(supp1,dims=2)
lsupp1=size(supp1,2)
lb=size(basis,2)
A=zeros(UInt8,lb,lb)
for i = 1:lb
    for j = i:lb
        bi=basis[:,i]+basis[:,j]
         if bfind(supp1,lsupp1,bi,n)!=0
           A[i,j]=1
           A[j,i]=1
        end
    end
end
else
    osupp=odd_supp(n,supp)
    osupp=sortslices(osupp,dims=2)
    lo=size(osupp,2)
    lb=size(basis,2)
    A=zeros(UInt8,lb,lb)
    for i = 1:lb
        for j = i:lb
            bi=basis[:,i]+basis[:,j]
            if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(osupp,lo,bi,n)!=0
               A[i,j]=1
               A[j,i]=1
            end
        end
    end
end
blocks,cl,blocksize=cliquesFromSpMatD(A)
if cl==1
   println("$blocksize\n[1]")
   return blocks,cl,blocksize,blocksize,[1],1
else
   nub=unique(blocksize)
   nsizes=[sum(blocksize.== i) for i in nub]
   if nub!=ub||nsizes!=sizes
      ub=nub
      sizes=nsizes
      println("$ub\n$sizes")
      return blocks,cl,blocksize,ub,sizes,1
   else
      println("No higher blocking hierarchy")
      return 0,0,0,0,0,0
   end
end
end

function get_blocks(n,supp,basis;reduce=0)
if reduce==1
supp1=[supp 2*basis]
supp1=sortslices(supp1,dims=2)
supp1=unique(supp1,dims=2)
lsupp1=size(supp1,2)
lb=size(basis,2)
G=SimpleGraph(lb)
for i = 1:lb
    for j = i:lb
         bi=basis[:,i]+basis[:,j]
         if bfind(supp1,lsupp1,bi,n)!=0
            add_edge!(G,i,j)
         end
    end
end
else
    osupp=odd_supp(n,supp)
    osupp=sortslices(osupp,dims=2)
    lo=size(osupp,2)
    lb=size(basis,2)
    G=SimpleGraph(lb)
    for i = 1:lb
        for j = i:lb
            bi=basis[:,i]+basis[:,j]
            if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(osupp,lo,bi,n)!=0
               add_edge!(G,i,j)
            end
        end
    end
end
blocks=connected_components(G)
cl=size(blocks,1)
blocksize=zeros(Int,1,cl)
for i=1:cl
    blocksize[i]=length(blocks[i])
end
if cl==1
   println("Unblockable")
   return blocks,cl,blocksize,blocksize,[1],0
else
   ub=unique(blocksize)
   sizes=[sum(blocksize.== i) for i in ub]
   println("blocksizes:\n$ub\n$sizes")
   return blocks,cl,blocksize,ub,sizes,1
end
end

function blockupop(n,supp,coe,basis,blocks,cl,blocksize)
    lsupp=size(supp,2)
    supp1=zeros(UInt8,n,1)
    for i=1:cl
        for j=1:blocksize[i]
            for r=j:blocksize[i]
                bi=basis[:,blocks[i][j]]+basis[:,blocks[i][r]]
                supp1=[supp1 bi]
            end
        end
    end
    supp1=sortslices(supp1,dims=2)
    supp1=unique(supp1,dims=2)
    lsupp1=size(supp1,2)
    model=Model(with_optimizer(Mosek.Optimizer, QUIET=true))
    cons=Array{Any}(undef, lsupp1)
    cons.=AffExpr(0)
    pos=Array{Any}(undef, cl)
    for i=1:cl
        if blocksize[i]==1
           pos[i]=@variable(model, lower_bound=0)
           bi=UInt8(2)*basis[:,blocks[i]]
           Locb=bfind(supp1,lsupp1,bi,n)
           cons[Locb]+=pos[i]
        else
           pos[i]=@variable(model, [1:blocksize[i], 1:blocksize[i]], PSD)
           for j=1:blocksize[i]
               for r=j:blocksize[i]
                   bi=basis[:,blocks[i][j]]+basis[:,blocks[i][r]]
                   Locb=bfind(supp1,lsupp1,bi,n)
                   if j==r
                      cons[Locb]+=pos[i][j,r]
                   else
                      cons[Locb]+=2*pos[i][j,r]
                   end
               end
           end
        end
    end
    bc=zeros(1,lsupp1)
    for i=1:lsupp
        Locb=bfind(supp1,lsupp1,supp[:,i],n)
        bc[Locb]=coe[i]
    end
    @constraint(model, cons[2:end].==bc[2:end])
    @variable(model, lower)
    @constraint(model, cons[1]+lower==bc[1])
    @objective(model, Max, lower)
    optimize!(model)
    status=termination_status(model)
    if status == MOI.OPTIMAL
       objv = objective_value(model)
       println("optimum = $objv")
    else
       objv = objective_value(model)
       println("$status")
       println("optimum = $objv")
    end
    return objv,supp1
end

function blockupopm(n,supp,coe,basis,blocks,cl,blocksize)
    lsupp=size(supp,2)
    supp1=zeros(UInt8,n,1)
    for i=1:cl
        for j=1:blocksize[i]
            for r=j:blocksize[i]
                bi=basis[:,blocks[i][j]]+basis[:,blocks[i][r]]
                supp1=[supp1 bi]
            end
        end
    end
    supp1=sortslices(supp1,dims=2)
    supp1=unique(supp1,dims=2)
    lsupp1=size(supp1,2)
    # printstream(msg::String)=print(msg)
    indexb=[i for i=1:cl]
    oneb=indexb[[blocksize[k]==1 for k=1:cl]]
    semb=indexb[[blocksize[k]!=1 for k=1:cl]]
    sblocksize=blocksize[semb]
    oblocks=blocks[oneb]
    sblocks=blocks[semb]
    lone=length(oneb)
    scl=cl-lone
    bkc=[MSK_BK_FX for i=1:lsupp1]
    bc=zeros(1,lsupp1)[1:end]
    for i=1:lsupp
        Locb=bfind(supp1,lsupp1,supp[:,i],n)
        bc[Locb]=coe[i]
    end
    consi=Array{Any}(undef,lsupp1,scl)
    consj=Array{Any}(undef,lsupp1,scl)
    consk=Array{Any}(undef,lsupp1,scl)
    dims=zeros(UInt16,lsupp1,scl)
    for i=1:lsupp1
        for j=1:scl
            dims[i,j]=sblocksize[j]
            consi[i,j]=[UInt32(0)]
            consj[i,j]=[UInt32(0)]
            consk[i,j]=[0]
        end
    end
    for i=1:scl
        for j=1:sblocksize[i]
            for r=j:sblocksize[i]
                bi=basis[:,sblocks[i][j]]+basis[:,sblocks[i][r]]
                Locb=bfind(supp1,lsupp1,bi,n)
                consi[Locb,i]=[consi[Locb,i];r]
                consj[Locb,i]=[consj[Locb,i];j]
                consk[Locb,i]=[consk[Locb,i];1.0]
           end
        end
    end
    for i=1:lsupp1
        for j=1:scl
            consi[i,j]=consi[i,j][2:end]
            consj[i,j]=consj[i,j][2:end]
            consk[i,j]=consk[i,j][2:end]
        end
    end
    oLocb=zeros(UInt32,lone,1)[1:end]
    for i=1:lone
        bi=2*basis[:,oblocks[i][1]]
        oLocb[i]=bfind(supp1,lsupp1,bi,n)
    end
    A=sparse(oLocb,[i for i=1:lone],[1.0 for i=1:lone])
    maketask() do task
    # putstreamfunc(task,MSK_STREAM_LOG,printstream)
    appendvars(task,1+lone)
    appendcons(task,lsupp1)
    appendbarvars(task,sblocksize[1:end])
    putcj(task,1,1.0)
    putvarboundslice(task,1,2,[MSK_BK_FR],[-Inf],[+Inf])
    putvarboundslice(task,2,2+lone,[MSK_BK_LO for i=1:lone],[0 for i=1:lone],[+Inf for i=1:lone])
    putconboundslice(task,1,lsupp1+1,bkc,bc,bc)
    putacolslice(task,1,2,[1],[2],[1],[1.0])
    putacolslice(task,2,2+lone,A.colptr[1:lone],A.colptr[2:lone+1],A.rowval,A.nzval)
    for k=1:lsupp1
        for j=1:scl
            cons=appendsparsesymmat(task,dims[k,j],consi[k,j],consj[k,j],consk[k,j])
            putbaraij(task,k,j,[cons],[1.0])
        end
    end
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)
    solsta=getsolsta(task,MSK_SOL_ITR)
    if solsta==MSK_SOL_STA_OPTIMAL
       opt=getprimalobj(task,MSK_SOL_ITR)
       barx=getbarxj(task,MSK_SOL_ITR,1)
    #   @printf("Optimal solution: \n%s\n", opt')
    elseif solsta==MSK_SOL_STA_DUAL_INFEAS_CER
       println("Primal or dual infeasibility.\n")
    elseif solsta==MSK_SOL_STA_PRIM_INFEAS_CER
       println("Primal or dual infeasibility.\n")
    elseif solsta==MSK_SOL_STA_UNKNOWN
       println("Unknown solution status")
    else
       println("Other solution status")
    end
    end
    return supp1
end

function get_hblocks!(n,supp,basis,ub,sizes;reduce=0)
if reduce==1
supp1=[supp 2*basis]
supp1=sortslices(supp1,dims=2)
supp1=unique(supp1,dims=2)
lsupp1=size(supp1,2)
lb=size(basis,2)
G=SimpleGraph(lb)
for i = 1:lb
    for j = i:lb
        bi=basis[:,i]+basis[:,j]
        if bfind(supp1,lsupp1,bi,n)!=0
           add_edge!(G,i,j)
        end
    end
end
else
    osupp=odd_supp(n,supp)
    lo=size(osupp,2)
    lb=size(basis,2)
    G=SimpleGraph(lb)
    for i = 1:lb
        for j = i:lb
            bi=basis[:,i]+basis[:,j]
            if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(osupp,lo,bi,n)!=0
               add_edge!(G,i,j)
            end
        end
    end
end
blocks=connected_components(G)
cl=size(blocks,1)
blocksize=zeros(Int,1,cl)
for i=1:cl
    blocksize[i]=length(blocks[i])
end
if cl==1
   println("blocksizes:\n$blocksize\n[1]")
   return blocks,cl,blocksize,blocksize,[1],1
else
   nub=unique(blocksize)
   nsizes=[sum(blocksize.== i) for i in nub]
   if nub!=ub||nsizes!=sizes
      ub=nub
      sizes=nsizes
      println("blocksizes:\n$ub\n$sizes")
      return blocks,cl,blocksize,ub,sizes,1
   else
      println("No higher block hierarchy")
      return 0,0,0,0,0,0
   end
end
end
