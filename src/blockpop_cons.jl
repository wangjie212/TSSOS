mutable struct cdata_type
    ssupp
    coe
    lt
    fbasis
    gbasis
    fsupp
    gsupp
    ub
    sizes
end

function blockcpop_first(n,m,d,dg,supp,ssupp,coe,lt)
fbasis=get_basis(n,d)
gbasis=Array{Any}(undef,m)
opt=0
for k=1:m
    gbasis[k]=get_basis(n,d-Int(ceil(dg[k]/2)))
end
fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis)
if status==1
   opt,fsupp,gsupp=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize)
end
data=cdata_type(ssupp,coe,lt,fbasis,gbasis,fsupp,gsupp,ub,sizes)
return opt,data,status
end

function blockcpop_higher!(n,m,data)
ssupp=data.ssupp
coe=data.coe
lt=data.lt
fbasis=data.fbasis
gbasis=data.gbasis
fsupp=data.fsupp
gsupp=data.gsupp
ub=data.ub
sizes=data.sizes
opt=0
fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chblocks(n,m,ssupp,lt,fbasis,gbasis,fsupp,gsupp,ub,sizes)
if status==1
   opt,fsupp,gsupp=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize)
end
data.fsupp=fsupp
data.gsupp=gsupp
data.ub=ub
data.sizes=sizes
return opt,data,status
end

function get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis)
osupp=odd_supp(n,supp)
osupp=sortslices(osupp,dims=2)
lo=size(osupp,2)
flb=binomial(n+d,d)
fG=SimpleGraph(flb)
for i = 1:flb
    for j = i:flb
        bi=fbasis[:,i]+fbasis[:,j]
        if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(osupp,lo,bi,n)!=0
           add_edge!(fG,i,j)
        end
    end
end
fblocks=connected_components(fG)
fcl=size(fblocks,1)
fblocksize=zeros(Int,1,fcl)
for i=1:fcl
    fblocksize[i]=length(fblocks[i])
end
if fcl==1
   println("Unblockable")
   return 0,0,0,0,0,0,0,0,0
else
   ub=unique(fblocksize)
   sizes=[sum(fblocksize.== i) for i in ub]
   println("$ub\n$sizes")
end

glb=Array{UInt16}(undef, m)
gcl=Array{UInt8}(undef, m)
gblocks=Array{Any}(undef, m)
gblocksize=Array{Any}(undef, m)
for k=1:m
    glb[k]=size(gbasis[k],2)
    gG=SimpleGraph(glb[k])
    for i = 1:glb[k]
        for j = i:glb[k]
            r=1
            while r<=lt[k+1]
                  bi=ssupp[k+1][:,r]+gbasis[k][:,i]+gbasis[k][:,j]
                  if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(osupp,lo,bi,n)!=0
                     break
                  else
                     r=r+1
                  end
            end
            if r<=lt[k+1]
               add_edge!(gG,i,j)
            end
        end
    end
    gblocks[k]=connected_components(gG)
    gcl[k]=size(gblocks[k],1)
    gblocksize[k]=zeros(Int,1,gcl[k])
    for i=1:gcl[k]
        gblocksize[k][i]=length(gblocks[k][i])
    end
    if gcl[k]==1
       gbk=gblocksize[k]
       println("$gbk\n[1]")
    else
       gub=unique(gblocksize[k])
       gsizes=[sum(gblocksize[k].== i) for i in gub]
       println("$gub\n$gsizes")
    end
end
return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,1
end

function blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize)
    fsupp=zeros(UInt8,n,1)
    for i=1:fcl
        for j=1:fblocksize[i]
            for r=j:fblocksize[i]
                bi=fbasis[:,fblocks[i][j]]+fbasis[:,fblocks[i][r]]
                fsupp=[fsupp bi]
            end
        end
    end
    gsupp=Array{Any}(undef, m)
    supp1=fsupp
    for k=1:m
        gsupp[k]=zeros(UInt8,n,1)
        for i=1:gcl[k]
            for j=1:gblocksize[k][i]
                for r=j:gblocksize[k][i]
                    for s=1:lt[k+1]
                        bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                        gsupp[k]=[gsupp[k] bi]
                    end
                end
            end
        end
        supp1=[supp1 gsupp[k]]
    end
    supp1=sortslices(supp1,dims=2)
    supp1=unique(supp1,dims=2)
    lsupp1=size(supp1,2)
    model=Model(with_optimizer(Mosek.Optimizer, QUIET=true))
    cons=Array{Any}(undef, lsupp1)
    cons.=AffExpr(0)
    pos=Array{Any}(undef, fcl)
    for i=1:fcl
        if fblocksize[i]==1
           pos[i]=@variable(model, lower_bound=0)
           bi=UInt8(2)*fbasis[:,fblocks[i]]
           Locb=bfind(supp1,lsupp1,bi,n)
           cons[Locb]+=pos[i]
        else
           pos[i]=@variable(model, [1:fblocksize[i], 1:fblocksize[i]], PSD)
           for j=1:fblocksize[i]
               for r=j:fblocksize[i]
                   bi=fbasis[:,fblocks[i][j]]+fbasis[:,fblocks[i][r]]
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
    gpos=Array{Any}(undef, m)
    for k=1:m
        gpos[k]=Array{Any}(undef, gcl[k])
        for i=1:gcl[k]
            if gblocksize[k][i]==1
               gpos[k][i]=@variable(model, lower_bound=0)
               for s=1:lt[k+1]
                   bi=ssupp[k+1][:,[s]]+UInt8(2)*gbasis[k][:,gblocks[k][i]]
                   Locb=bfind(supp1,lsupp1,bi,n)
                   cons[Locb]+=coe[k+1][s]*gpos[k][i]
               end
            else
               gpos[k][i]=@variable(model, base_name="gx$k$i", [1:gblocksize[k][i], 1:gblocksize[k][i]], PSD)
               for j=1:gblocksize[k][i]
                   for r=j:gblocksize[k][i]
                       for s=1:lt[k+1]
                           bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                           Locb=bfind(supp1,lsupp1,bi,n)
                           if j==r
                              cons[Locb]+=coe[k+1][s]*gpos[k][i][j,r]
                           else
                              cons[Locb]+=2*coe[k+1][s]*gpos[k][i][j,r]
                           end
                       end
                   end
               end
           end
       end
    end
    bc=zeros(1,lsupp1)
    for i=1:lt[1]
        Locb=bfind(supp1,lsupp1,ssupp[1][:,i],n)
        bc[Locb]=coe[1][i]
    end
    @constraint(model, cons[2:end].==bc[2:end])
    @variable(model, lower)
    @constraint(model, cons[1]+lower==bc[1])
    @objective(model, Max, lower)
    optimize!(model)
    status=termination_status(model)
    if  status==MOI.OPTIMAL
        objv=objective_value(model)
        println("$objv")
    else
        println("$status")
    end
    return objv,fsupp,gsupp
end

function get_chblocks(n,m,ssupp,lt,fbasis,gbasis,fsupp,gsupp,ub,sizes)
ofsupp=odd_supp(n,fsupp)
ofsupp=sortslices(ofsupp,dims=2)
lfo=size(ofsupp,2)
flb=size(fbasis,2)
fG=SimpleGraph(flb)
for i = 1:flb
    for j = i:flb
        bi=fbasis[:,i]+fbasis[:,j]
        if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(ofsupp,lfo,bi,n)!=0
           add_edge!(fG,i,j)
        end
    end
end
fblocks=connected_components(fG)
fcl=size(fblocks,1)
fblocksize=zeros(Int,1,fcl)
for i=1:fcl
    fblocksize[i]=length(fblocks[i])
end
if fcl==1
   println("$fblocksize\n[1]")
   gblocksize=Array{Any}(undef, m)
   gcl=Array{UInt8}(undef, m)
   for k=1:m
       gblocksize[k]=[i for i=1:size(gbasis[k],1)]
       gcl[k]=size(gblocks[k],1)
   end
   return fblocks,fcl,fblocksize,gblocks,gcl,[1],fblocksize,1
else
   nub=unique(fblocksize)
   nsizes=[sum(fblocksize.== i) for i in nub]
   if nub!=ub||nsizes!=sizes
      ub=nub
      sizes=nsizes
      println("$ub\n$sizes")
   else
      println("No higher blocking hierarchy")
      return 0,0,0,0,0,0,0,0,0,0,0
   end
end

glb=Array{UInt16}(undef, m)
gcl=Array{UInt8}(undef, m)
gblocks=Array{Any}(undef, m)
gblocksize=Array{Any}(undef, m)
for k=1:m
    glb[k]=size(gbasis[k],2)
    gG=SimpleGraph(glb[k])
    ggsupp=[fsupp gsupp[k]]
    ogsupp=odd_supp(n,ggsupp)
    ogsupp=sortslices(ogsupp,dims=2)
    lgo=size(ogsupp,2)
    for i = 1:glb[k]
        for j = i:glb[k]
            r=1
            while r<=lt[k+1]
                  bi=ssupp[k+1][:,r]+gbasis[k][:,i]+gbasis[k][:,j]
                  if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(ogsupp,lgo,bi,n)!=0
                     break
                  else
                     r=r+1
                  end
            end
            if r<=lt[k+1]
               add_edge!(gG,i,j)
            end
        end
    end
    gblocks[k]=connected_components(gG)
    gcl[k]=size(gblocks[k],1)
    gblocksize[k]=zeros(Int,1,gcl[k])
    for i=1:gcl[k]
        gblocksize[k][i]=length(gblocks[k][i])
    end
    if gcl[k]==1
       gbk=gblocksize[k]
       println("$gbk\n[1]")
    else
       gub=unique(gblocksize[k])
       gsizes=[sum(gblocksize[k].== i) for i in gub]
       println("$gub\n$gsizes")
    end
end
return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,1
end
