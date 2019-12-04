mutable struct cdata_type
    ssupp
    coe
    lt
    d
    dg
    fbasis
    gbasis
    fsupp
    gsupp
    ub
    sizes
end

function blockcpop_first(n,m,d,dg,supp,ssupp,coe,lt;method="block",reducebasis=0)
fbasis=get_basis(n,d)
gbasis=Array{Any}(undef,m)
for k=1:m
    gbasis[k]=get_basis(n,d-Int(ceil(dg[k]/2)))
end
if method=="block"&&reducebasis==0
   fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis)
elseif method=="block"&&reducebasis==1
    flag=1
    while flag==1
          fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis,reduce=1)
          tsupp=[ssupp[1] zeros(UInt8,n,1)]
          for k=1:m
              for i=1:gcl[k]
                  for j=1:gblocksize[k][i]
                      for r=j:gblocksize[k][i]
                          for s=1:lt[k+1]
                              bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                              tsupp=[tsupp bi]
                          end
                      end
                  end
              end
          end
          fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize)
    end
elseif method=="clique"&&reducebasis==0
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis)
else
    flag=1
    while flag==1
          fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis,reduce=1)
          tsupp=[ssupp[1] zeros(UInt8,n,1)]
          for k=1:m
              for i=1:gcl[k]
                  for j=1:gblocksize[k][i]
                      for r=j:gblocksize[k][i]
                          for s=1:lt[k+1]
                              bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                              tsupp=[tsupp bi]
                          end
                      end
                  end
              end
          end
          fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize)
    end
end
if status==1
   opt,fsupp,gsupp=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize)
else
   opt=0
   fsupp=0
   gsupp=0
end
data=cdata_type(ssupp,coe,lt,d,dg,fbasis,gbasis,fsupp,gsupp,ub,sizes)
return opt,data,status
end

function blockcpop_higher!(n,m,data;method="block",reducebasis=0)
ssupp=data.ssupp
coe=data.coe
lt=data.lt
d=data.d
dg=data.dg
gbasis=data.gbasis
fsupp=data.fsupp
gsupp=data.gsupp
ub=data.ub
sizes=data.sizes
opt=0
if method=="block"&&reducebasis==0
   fbasis=data.fbasis
   fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gsupp,ub,sizes)
elseif method=="block"&&reducebasis==1
    fbasis=get_basis(n,d)
    flag=1
    while flag==1
          fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gsupp,ub,sizes,reduce=1)
          tsupp=[ssupp[1] zeros(UInt8,n,1)]
          for k=1:m
              for i=1:gcl[k]
                  for j=1:gblocksize[k][i]
                      for r=j:gblocksize[k][i]
                          for s=1:lt[k+1]
                              bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                              tsupp=[tsupp bi]
                          end
                      end
                  end
              end
          end
          fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize)
    end
elseif method=="clique"&&reducebasis==0
    fbasis=data.fbasis
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gsupp,ub,sizes)
else
    fbasis=get_basis(n,d)
    flag=1
    while flag==1
          fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gsupp,ub,sizes,reduce=1)
          tsupp=[ssupp[1] zeros(UInt8,n,1)]
          for k=1:m
              for i=1:gcl[k]
                  for j=1:gblocksize[k][i]
                      for r=j:gblocksize[k][i]
                          for s=1:lt[k+1]
                              bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                              tsupp=[tsupp bi]
                          end
                      end
                  end
              end
          end
          fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize)
    end
end
if status==1
   opt,fsupp,gsupp=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize)
end
data.fsupp=fsupp
data.gsupp=gsupp
data.fbasis=fbasis
data.ub=ub
data.sizes=sizes
return opt,data,status
end

function reducebasis!(n,supp,basis,blocks,cl,blocksize)
esupp=even_supp(n,supp)
init=0
flag=0
check=0
while init==0||check>0
init=1
check=0
supp1=esupp
for i=1:cl
    if blocksize[i]>1
    for j=1:blocksize[i]
        for r=j+1:blocksize[i]
            bi=basis[:,blocks[i][j]]+basis[:,blocks[i][r]]
            supp1=[supp1 bi]
        end
    end
    end
end
supp1=sortslices(supp1,dims=2)
supp1=unique(supp1,dims=2)
lsupp1=size(supp1,2)
for i=1:cl
    lo=blocksize[i]
    indexb=[k for k=1:lo]
    j=Int(1)
    while lo>=j
        bi=2*basis[:,blocks[i][indexb[j]]]
        Locb=bfind(supp1,lsupp1,bi,n)
        if Locb==0
           check=1
           flag=1
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
       indexb=[indexb;blocks[i]]
   end
   indexb=sort(indexb)
   indexb=unique(indexb)
   return basis[:,indexb],flag
else
   return basis,flag
end
end

function get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis;reduce=0)
if reduce==1
supp1=[supp 2*fbasis]
supp1=sortslices(supp1,dims=2)
supp1=unique(supp1,dims=2)
lsupp1=size(supp1,2)
flb=size(fbasis,2)
fG=SimpleGraph(flb)
for i = 1:flb
    for j = i:flb
        bi=fbasis[:,i]+fbasis[:,j]
        if bfind(supp1,lsupp1,bi,n)!=0
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
   println("fblocksizes:\n$ub\n$sizes")
end
println("gblocksizes:")
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
                  if bfind(supp1,lsupp1,bi,n)!=0
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
else
    osupp=odd_supp(n,supp)
    osupp=sortslices(osupp,dims=2)
    lo=size(osupp,2)
    flb=size(fbasis,2)
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
       println("fblocksizes:\n$ub\n$sizes")
    end
    println("gblocksizes:")
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
end
return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,1
end

function get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis;reduce=0)
if reduce==1
supp1=[supp 2*fbasis]
supp1=sortslices(supp1,dims=2)
supp1=unique(supp1,dims=2)
lsupp1=size(supp1,2)
flb=size(fbasis,2)
A=zeros(UInt8,flb,flb)
for i = 1:flb
    for j = i:flb
        bi=fbasis[:,i]+fbasis[:,j]
        if bfind(supp1,lsupp1,bi,n)!=0
            A[i,j]=1
            A[j,i]=1
        end
    end
end
fblocks,fcl,fblocksize=cliquesFromSpMatD(A)
if fcl==1
   println("Unblockable")
   return 0,0,0,0,0,0,0,0,0
else
   ub=unique(fblocksize)
   sizes=[sum(fblocksize.== i) for i in ub]
   println("fblocksizes:\n$ub\n$sizes")
end
println("gblocksizes:")
glb=Array{UInt16}(undef, m)
gcl=Array{UInt8}(undef, m)
gblocks=Array{Any}(undef, m)
gblocksize=Array{Any}(undef, m)
for k=1:m
    glb[k]=size(gbasis[k],2)
    A=zeros(UInt8,glb[k],glb[k])
    for i = 1:glb[k]
        for j = i:glb[k]
            r=1
            while r<=lt[k+1]
                  bi=ssupp[k+1][:,r]+gbasis[k][:,i]+gbasis[k][:,j]
                  if bfind(supp1,lsupp1,bi,n)!=0
                     break
                  else
                     r=r+1
                  end
            end
            if r<=lt[k+1]
                A[i,j]=1
                A[j,i]=1
            end
        end
    end
    gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A)
    if gcl[k]==1
       gbk=gblocksize[k]
       println("$gbk\n[1]")
    else
       gub=unique(gblocksize[k])
       gsizes=[sum(gblocksize[k].== i) for i in gub]
       println("$gub\n$gsizes")
    end
end
else
    osupp=odd_supp(n,supp)
    osupp=sortslices(osupp,dims=2)
    lo=size(osupp,2)
    flb=size(fbasis,2)
    A=zeros(UInt8,flb,flb)
    for i = 1:flb
        for j = i:flb
            bi=fbasis[:,i]+fbasis[:,j]
            if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(osupp,lo,bi,n)!=0
                A[i,j]=1
                A[j,i]=1
            end
        end
    end
    fblocks,fcl,fblocksize=cliquesFromSpMatD(A)
    if fcl==1
       println("Unblockable")
       return 0,0,0,0,0,0,0,0,0
    else
       ub=unique(fblocksize)
       sizes=[sum(fblocksize.== i) for i in ub]
       println("fblocksizes:\n$ub\n$sizes")
    end
    println("gblocksizes:")
    glb=Array{UInt16}(undef, m)
    gcl=Array{UInt8}(undef, m)
    gblocks=Array{Any}(undef, m)
    gblocksize=Array{Any}(undef, m)
    for k=1:m
        glb[k]=size(gbasis[k],2)
        A=zeros(UInt8,glb[k],glb[k])
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
                    A[i,j]=1
                    A[j,i]=1
                end
            end
        end
        gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A)
        if gcl[k]==1
           gbk=gblocksize[k]
           println("$gbk\n[1]")
        else
           gub=unique(gblocksize[k])
           gsizes=[sum(gblocksize[k].== i) for i in gub]
           println("$gub\n$gsizes")
        end
    end
end
return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,1
end

function get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gsupp,ub,sizes;reduce=0)
if reduce==1
    supp1=[fsupp 2*fbasis]
    supp1=sortslices(supp1,dims=2)
    supp1=unique(supp1,dims=2)
    lsupp1=size(supp1,2)
    flb=size(fbasis,2)
    A=zeros(UInt8,flb,flb)
    for i = 1:flb
        for j = i:flb
            bi=fbasis[:,i]+fbasis[:,j]
            if bfind(supp1,lsupp1,bi,n)!=0
               A[i,j]=1
               A[j,i]=1
            end
        end
    end
    fblocks,fcl,fblocksize=cliquesFromSpMatD(A)
    if fcl==1
       println("fblocksizes:\n$fblocksize\n[1]")
       println("gblocksizes:")
       gblocks=Array{Any}(undef, m)
       gblocksize=Array{Any}(undef, m)
       gcl=Array{UInt8}(undef, m)
       for k=1:m
           gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
           gblocksize[k]=[size(gbasis[k],2)]
           gcl[k]=1
           gbk=gblocksize[k]
           println("$gbk\n[1]")
       end
       return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
    else
       nub=unique(fblocksize)
       nsizes=[sum(fblocksize.== i) for i in nub]
       if nub!=ub||nsizes!=sizes
          ub=nub
          sizes=nsizes
          println("fblocksizes:\n$ub\n$sizes")
       else
          println("No higher block hierarchy")
          return 0,0,0,0,0,0,0,0,0
       end
    end
    println("gblocksizes:")
    glb=Array{UInt16}(undef, m)
    gcl=Array{UInt8}(undef, m)
    gblocks=Array{Any}(undef, m)
    gblocksize=Array{Any}(undef, m)
    for k=1:m
        glb[k]=size(gbasis[k],2)
        A=zeros(UInt8,glb[k],glb[k])
        ggsupp=[fsupp gsupp[k]]
        supp2=[ggsupp 2*fbasis]
        supp2=sortslices(supp2,dims=2)
        supp2=unique(supp2,dims=2)
        lsupp2=size(supp2,2)
        for i = 1:glb[k]
            for j = i:glb[k]
                r=1
                while r<=lt[k+1]
                      bi=ssupp[k+1][:,r]+gbasis[k][:,i]+gbasis[k][:,j]
                      if bfind(supp2,lsupp2,bi,n)!=0
                         break
                      else
                         r=r+1
                      end
                end
                if r<=lt[k+1]
                    A[i,j]=1
                    A[j,i]=1
                end
            end
        end
        gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A)
        if gcl[k]==1
           gbk=gblocksize[k]
           println("$gbk\n[1]")
        else
           gub=unique(gblocksize[k])
           gsizes=[sum(gblocksize[k].== i) for i in gub]
           println("$gub\n$gsizes")
        end
    end
else
ofsupp=odd_supp(n,fsupp)
ofsupp=sortslices(ofsupp,dims=2)
lfo=size(ofsupp,2)
flb=size(fbasis,2)
A=zeros(UInt8,flb,flb)
for i = 1:flb
    for j = i:flb
        bi=fbasis[:,i]+fbasis[:,j]
        if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(ofsupp,lfo,bi,n)!=0
           A[i,j]=1
           A[j,i]=1
        end
    end
end
fblocks,fcl,fblocksize=cliquesFromSpMatD(A)
if fcl==1
   println("fblocksizes:\n$fblocksize\n[1]")
   println("gblocksizes:")
   gblocks=Array{Any}(undef, m)
   gblocksize=Array{Any}(undef, m)
   gcl=Array{UInt8}(undef, m)
   for k=1:m
       gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
       gblocksize[k]=[size(gbasis[k],2)]
       gcl[k]=1
       gbk=gblocksize[k]
       println("$gbk\n[1]")
   end
   return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
else
   nub=unique(fblocksize)
   nsizes=[sum(fblocksize.== i) for i in nub]
   if nub!=ub||nsizes!=sizes
      ub=nub
      sizes=nsizes
      println("fblocksizes:\n$ub\n$sizes")
   else
      println("No higher block hierarchy")
      return 0,0,0,0,0,0,0,0,0
   end
end
println("gblocksizes:")
glb=Array{UInt16}(undef, m)
gcl=Array{UInt8}(undef, m)
gblocks=Array{Any}(undef, m)
gblocksize=Array{Any}(undef, m)
for k=1:m
    glb[k]=size(gbasis[k],2)
    A=zeros(UInt8,glb[k],glb[k])
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
                A[i,j]=1
                A[j,i]=1
            end
        end
    end
    gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A)
    if gcl[k]==1
       gbk=gblocksize[k]
       println("$gbk\n[1]")
    else
       gub=unique(gblocksize[k])
       gsizes=[sum(gblocksize[k].== i) for i in gub]
       println("$gub\n$gsizes")
    end
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
        println("optimum=\n$objv")
    else
        objv = 0
        println("Failed, status=\n$status")
    end
    return objv,fsupp,gsupp
end

function get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gsupp,ub,sizes;reduce=0)
if reduce==1
    supp1=[fsupp 2*fbasis]
    supp1=sortslices(supp1,dims=2)
    supp1=unique(supp1,dims=2)
    lsupp1=size(supp1,2)
    flb=size(fbasis,2)
    fG=SimpleGraph(flb)
    for i = 1:flb
        for j = i:flb
            bi=fbasis[:,i]+fbasis[:,j]
            if bfind(supp1,lsupp1,bi,n)!=0
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
       println("fblocksizes:\n$fblocksize\n[1]")
       println("gblocksizes:")
       gblocks=Array{Any}(undef, m)
       gblocksize=Array{Any}(undef, m)
       gcl=Array{UInt8}(undef, m)
       for k=1:m
           gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
           gblocksize[k]=[size(gbasis[k],2)]
           gcl[k]=1
           gbk=gblocksize[k]
           println("$gbk\n[1]")
       end
       return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
    else
       nub=unique(fblocksize)
       nsizes=[sum(fblocksize.== i) for i in nub]
       if nub!=ub||nsizes!=sizes
          ub=nub
          sizes=nsizes
          println("fblocksizes:\n$ub\n$sizes")
       else
          println("No higher block hierarchy")
          return 0,0,0,0,0,0,0,0,0
       end
    end
    println("gblocksizes:")
    glb=Array{UInt16}(undef, m)
    gcl=Array{UInt8}(undef, m)
    gblocks=Array{Any}(undef, m)
    gblocksize=Array{Any}(undef, m)
    for k=1:m
        glb[k]=size(gbasis[k],2)
        gG=SimpleGraph(glb[k])
        ggsupp=[fsupp gsupp[k]]
        supp2=[ggsupp 2*fbasis]
        supp2=sortslices(supp2,dims=2)
        supp2=unique(supp2,dims=2)
        lsupp2=size(supp2,2)
        for i = 1:glb[k]
            for j = i:glb[k]
                r=1
                while r<=lt[k+1]
                      bi=ssupp[k+1][:,r]+gbasis[k][:,i]+gbasis[k][:,j]
                      if bfind(supp2,lsupp2,bi,n)!=0
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
else
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
   println("fblocksizes:\n$fblocksize\n[1]")
   println("gblocksizes:")
   gblocks=Array{Any}(undef, m)
   gblocksize=Array{Any}(undef, m)
   gcl=Array{UInt8}(undef, m)
   for k=1:m
       gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
       gblocksize[k]=[size(gbasis[k],2)]
       gcl[k]=1
       gbk=gblocksize[k]
       println("$gbk\n[1]")
   end
   return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
else
   nub=unique(fblocksize)
   nsizes=[sum(fblocksize.== i) for i in nub]
   if nub!=ub||nsizes!=sizes
      ub=nub
      sizes=nsizes
      println("fblocksizes:\n$ub\n$sizes")
   else
      println("No higher block hierarchy")
      return 0,0,0,0,0,0,0,0,0
   end
end
println("gblocksizes:")
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
end
return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,1
end
