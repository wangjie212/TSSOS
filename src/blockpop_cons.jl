mutable struct cdata_type
    n
    m
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
end

function blockcpop_first(pop,x,d;method="block",reducebasis=0,numeq=0,QUIET=true,dense=10)
n=length(x)
m=length(pop)-1
dg=zeros(UInt8,1,m)
coe=Array{Array{Float64,1}}(undef, m+1)
mon=Array{Any}(undef, m+1)
ssupp=Array{Array{UInt8,2}}(undef, m+1)
lt=zeros(Int,1,m+1)
for k=1:m+1
    mon[k]=monomials(pop[k])
    coe[k]=coefficients(pop[k])
    lt[k]=length(mon[k])
    ssupp[k]=zeros(UInt8,n,lt[k])
    for i=1:lt[k]
        for j=1:n
            ssupp[k][j,i]=MultivariatePolynomials.degree(mon[k][i],x[j])
        end
    end
end
supp=ssupp[1]
for i=2:m+1
    dg[i-1]=maxdegree(pop[i])
    supp=[supp ssupp[i]]
end
supp=unique(supp,dims=2)
fbasis=get_basis(n,d)
gbasis=Array{Array{UInt8,2}}(undef,m)
for k=1:m
    gbasis[k]=get_basis(n,d-Int(ceil(dg[k]/2)))
end
if method=="block"&&reducebasis==0
   fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis)
elseif method=="block"&&reducebasis==1
    flag=1
    while flag==1
          fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis,reduce=1)
          tsupp=[ssupp[1] zeros(UInt8,n,1)]
          for k=1:m
              gsupp=zeros(UInt8,n,lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2))
              l=1
              for i=1:gcl[k]
                  for j=1:gblocksize[k][i]
                      for r=j:gblocksize[k][i]
                          for s=1:lt[k+1]
                              bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                              gsupp[:,l]=bi
                              l+=1
                          end
                      end
                  end
              end
              tsupp=[tsupp gsupp]
          end
          fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize)
    end
elseif method=="clique"&&reducebasis==0
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis,dense=dense)
else
    flag=1
    while flag==1
          fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis,reduce=1,dense=dense)
          tsupp=[ssupp[1] zeros(UInt8,n,1)]
          for k=1:m
              gsupp=zeros(UInt8,n,lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2))
              l=1
              for i=1:gcl[k]
                  for j=1:gblocksize[k][i]
                      for r=j:gblocksize[k][i]
                          for s=1:lt[k+1]
                              bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                              gsupp[:,l]=bi
                              l+=1
                          end
                      end
                  end
              end
              tsupp=[tsupp gsupp]
          end
          fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize)
    end
end
opt,fsupp,Gram=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,numeq=numeq,QUIET=QUIET)
data=cdata_type(n,m,ssupp,coe,lt,d,dg,fbasis,gbasis,fsupp,ub,sizes)
return opt,data
end

function blockcpop_higher!(data;method="block",reducebasis=0,numeq=0,QUIET=true,dense=10)
n=data.n
m=data.m
ssupp=data.ssupp
coe=data.coe
lt=data.lt
d=data.d
dg=data.dg
gbasis=data.gbasis
fsupp=data.fsupp
ub=data.ub
sizes=data.sizes
opt=0
if method=="block"&&reducebasis==0
   fbasis=data.fbasis
   fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,ub,sizes)
elseif method=="block"&&reducebasis==1
    fbasis=get_basis(n,d)
    flag=1
    while flag==1
          fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,ub,sizes,reduce=1)
          tsupp=[ssupp[1] zeros(UInt8,n,1)]
          for k=1:m
              gsupp=zeros(UInt8,n,lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2))
              l=1
              for i=1:gcl[k]
                  for j=1:gblocksize[k][i]
                      for r=j:gblocksize[k][i]
                          for s=1:lt[k+1]
                              bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                              gsupp[:,l]=bi
                              l+=1
                          end
                      end
                  end
              end
              tsupp=[tsupp gsupp]
          end
          fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize)
    end
elseif method=="clique"&&reducebasis==0
    fbasis=data.fbasis
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,ub,sizes,dense=dense)
else
    fbasis=get_basis(n,d)
    flag=1
    while flag==1
          fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,ub,sizes,reduce=1,dense=dense)
          tsupp=[ssupp[1] zeros(UInt8,n,1)]
          for k=1:m
              gsupp=zeros(UInt8,n,lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2))
              l=1
              for i=1:gcl[k]
                  for j=1:gblocksize[k][i]
                      for r=j:gblocksize[k][i]
                          for s=1:lt[k+1]
                              bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                              gsupp[:,l]=bi
                              l+=1
                          end
                      end
                  end
              end
              tsupp=[tsupp gsupp]
          end
          fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize)
    end
end
if status==1
    opt,fsupp,Gram=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,numeq=numeq,QUIET=QUIET)
end
data.fsupp=fsupp
data.fbasis=fbasis
data.ub=ub
data.sizes=sizes
return opt,data
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
    j=1
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
    println("fblocksizes:\n$fblocksize\n[1]")
    println("-----------------------------------------------")
    println("gblocksizes:")
    gblocks=Array{Array{Array{Int,1},1}}(undef, m)
    gblocksize=Array{Array{Int,1}}(undef, m)
    gcl=zeros(Int,1,m)
    for k=1:m
        gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
        gblocksize[k]=[size(gbasis[k],2)]
        gcl[k]=1
        gbk=gblocksize[k]
        println("$gbk\n[1]")
        println("-----------------------------------------------")
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1]
else
ub=unique(fblocksize)
sizes=[sum(fblocksize.== i) for i in ub]
println("fblocksizes:\n$ub\n$sizes")
println("-----------------------------------------------")
println("gblocksizes:")
glb=zeros(Int,1,m)
gblocks=Array{Array{Array{Int,1},1}}(undef, m)
gblocksize=Array{Array{Int,2}}(undef, m)
gcl=zeros(Int,1,m)
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
       println("-----------------------------------------------")
    else
       gub=unique(gblocksize[k])
       gsizes=[sum(gblocksize[k].== i) for i in gub]
       println("$gub\n$gsizes")
       println("-----------------------------------------------")
    end
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
        println("fblocksizes:\n$fblocksize\n[1]")
        println("-----------------------------------------------")
        println("gblocksizes:")
        gblocks=Array{Array{Array{Int,1},1}}(undef, m)
        gblocksize=Array{Array{Int,1}}(undef, m)
        gcl=zeros(Int,1,m)
        for k=1:m
            gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
            gblocksize[k]=[size(gbasis[k],2)]
            gcl[k]=1
            gbk=gblocksize[k]
            println("$gbk\n[1]")
            println("-----------------------------------------------")
        end
        return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1]
    else
       ub=unique(fblocksize)
       sizes=[sum(fblocksize.== i) for i in ub]
       println("fblocksizes:\n$ub\n$sizes")
       println("-----------------------------------------------")
    end
    println("gblocksizes:")
    glb=zeros(Int,1,m)
    gblocks=Array{Array{Array{Int,1},1}}(undef, m)
    gblocksize=Array{Array{Int,2}}(undef, m)
    gcl=zeros(Int,1,m)
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
           println("-----------------------------------------------")
        else
           gub=unique(gblocksize[k])
           gsizes=[sum(gblocksize[k].== i) for i in gub]
           println("$gub\n$gsizes")
           println("-----------------------------------------------")
        end
    end
end
return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes
end

function get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis;reduce=0,dense=10)
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
fblocks,fcl,fblocksize=cliquesFromSpMatD(A,dense=dense)
if fcl==1
    println("fblocksizes:\n$fblocksize\n[1]")
    println("-----------------------------------------------")
    println("gblocksizes:")
    gblocks=Array{Array{Array{Int,1},1}}(undef, m)
    gblocksize=Array{Array{Int,1}}(undef, m)
    gcl=zeros(Int,1,m)
    for k=1:m
        gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
        gblocksize[k]=[size(gbasis[k],2)]
        gcl[k]=1
        gbk=gblocksize[k]
        println("$gbk\n[1]")
        println("-----------------------------------------------")
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1]
else
   ub=unique(fblocksize)
   sizes=[sum(fblocksize.== i) for i in ub]
   println("fblocksizes:\n$ub\n$sizes")
   println("-----------------------------------------------")
end
println("gblocksizes:")
glb=zeros(Int,1,m)
gblocks=Array{Array{Array{Int,1},1}}(undef, m)
gblocksize=Array{Array{Int,2}}(undef, m)
gcl=zeros(Int,1,m)
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
    gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A,dense=dense)
    if gcl[k]==1
       gbk=gblocksize[k]
       println("$gbk\n[1]")
       println("-----------------------------------------------")
    else
       gub=unique(gblocksize[k])
       gsizes=[sum(gblocksize[k].== i) for i in gub]
       println("$gub\n$gsizes")
       println("-----------------------------------------------")
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
    fblocks,fcl,fblocksize=cliquesFromSpMatD(A,dense=dense)
    if fcl==1
        println("fblocksizes:\n$fblocksize\n[1]")
        println("-----------------------------------------------")
        println("gblocksizes:")
        gblocks=Array{Array{Array{Int,1},1}}(undef, m)
        gblocksize=Array{Array{Int,1}}(undef, m)
        gcl=zeros(Int,1,m)
        for k=1:m
            gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
            gblocksize[k]=[size(gbasis[k],2)]
            gcl[k]=1
            gbk=gblocksize[k]
            println("$gbk\n[1]")
            println("-----------------------------------------------")
        end
        return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1]
    else
       ub=unique(fblocksize)
       sizes=[sum(fblocksize.== i) for i in ub]
       println("fblocksizes:\n$ub\n$sizes")
       println("-----------------------------------------------")
    end
    println("gblocksizes:")
    glb=zeros(Int,1,m)
    gblocks=Array{Array{Array{Int,1},1}}(undef, m)
    gblocksize=Array{Array{Int,2}}(undef, m)
    gcl=zeros(Int,1,m)
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
        gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A,dense=dense)
        if gcl[k]==1
           gbk=gblocksize[k]
           println("$gbk\n[1]")
           println("-----------------------------------------------")
        else
           gub=unique(gblocksize[k])
           gsizes=[sum(gblocksize[k].== i) for i in gub]
           println("$gub\n$gsizes")
           println("-----------------------------------------------")
        end
    end
end
return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes
end

function get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,ub,sizes;reduce=0,dense=10)
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
    fblocks,fcl,fblocksize=cliquesFromSpMatD(A,dense=dense)
    if fcl==1
       println("fblocksizes:\n$fblocksize\n[1]")
       println("-----------------------------------------------")
       println("gblocksizes:")
       gblocks=Array{Array{Array{Int,1},1}}(undef, m)
       gblocksize=Array{Array{Int,1}}(undef, m)
       gcl=zeros(Int,1,m)
       for k=1:m
           gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
           gblocksize[k]=[size(gbasis[k],2)]
           gcl[k]=1
           gbk=gblocksize[k]
           println("$gbk\n[1]")
           println("-----------------------------------------------")
       end
       return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
    else
       nub=unique(fblocksize)
       nsizes=[sum(fblocksize.== i) for i in nub]
       if nub!=ub||nsizes!=sizes
          ub=nub
          sizes=nsizes
          println("fblocksizes:\n$ub\n$sizes")
          println("-----------------------------------------------")
       else
          println("No higher clique hierarchy")
          return nothing,nothing,nothing,nothing,nothing,nothing,nothing,nothing,0
       end
    end
    println("gblocksizes:")
    glb=zeros(Int,1,m)
    gblocks=Array{Array{Array{Int,1},1}}(undef, m)
    gblocksize=Array{Array{Int,2}}(undef, m)
    gcl=zeros(Int,1,m)
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
        gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A,dense=dense)
        if gcl[k]==1
           gbk=gblocksize[k]
           println("$gbk\n[1]")
           println("-----------------------------------------------")
        else
           gub=unique(gblocksize[k])
           gsizes=[sum(gblocksize[k].== i) for i in gub]
           println("$gub\n$gsizes")
           println("-----------------------------------------------")
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
fblocks,fcl,fblocksize=cliquesFromSpMatD(A,dense=dense)
if fcl==1
   println("fblocksizes:\n$fblocksize\n[1]")
   println("-----------------------------------------------")
   println("gblocksizes:")
   gblocks=Array{Array{Array{Int,1},1}}(undef, m)
   gblocksize=Array{Array{Int,1}}(undef, m)
   gcl=zeros(Int,1,m)
   for k=1:m
       gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
       gblocksize[k]=[size(gbasis[k],2)]
       gcl[k]=1
       gbk=gblocksize[k]
       println("$gbk\n[1]")
       println("-----------------------------------------------")
   end
   return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
else
   nub=unique(fblocksize)
   nsizes=[sum(fblocksize.== i) for i in nub]
   if nub!=ub||nsizes!=sizes
      ub=nub
      sizes=nsizes
      println("fblocksizes:\n$ub\n$sizes")
      println("-----------------------------------------------")
   else
      println("No higher clique hierarchy")
      return nothing,nothing,nothing,nothing,nothing,nothing,nothing,nothing,0
   end
end
println("gblocksizes:")
glb=zeros(Int,1,m)
gblocks=Array{Array{Array{Int,1},1}}(undef, m)
gblocksize=Array{Array{Int,2}}(undef, m)
gcl=zeros(Int,1,m)
for k=1:m
    glb[k]=size(gbasis[k],2)
    A=zeros(UInt8,glb[k],glb[k])
    for i = 1:glb[k]
        for j = i:glb[k]
            r=1
            while r<=lt[k+1]
                  bi=ssupp[k+1][:,r]+gbasis[k][:,i]+gbasis[k][:,j]
                  if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(ofsupp,lfo,bi,n)!=0
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
    gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A,dense=dense)
    if gcl[k]==1
       gbk=gblocksize[k]
       println("$gbk\n[1]")
       println("-----------------------------------------------")
    else
       gub=unique(gblocksize[k])
       gsizes=[sum(gblocksize[k].== i) for i in gub]
       println("$gub\n$gsizes")
       println("-----------------------------------------------")
    end
end
end
return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,1
end

function get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,ub,sizes;reduce=0)
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
       println("-----------------------------------------------")
       println("gblocksizes:")
       gblocks=Array{Array{Array{Int,1},1}}(undef, m)
       gblocksize=Array{Array{Int,1}}(undef, m)
       gcl=zeros(Int,1,m)
       for k=1:m
           gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
           gblocksize[k]=[size(gbasis[k],2)]
           gcl[k]=1
           gbk=gblocksize[k]
           println("$gbk\n[1]")
           println("-----------------------------------------------")
       end
       return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
    else
       nub=unique(fblocksize)
       nsizes=[sum(fblocksize.== i) for i in nub]
       if nub!=ub||nsizes!=sizes
          ub=nub
          sizes=nsizes
          println("fblocksizes:\n$ub\n$sizes")
          println("-----------------------------------------------")
       else
          println("No higher block hierarchy")
          return nothing,nothing,nothing,nothing,nothing,nothing,nothing,nothing,0
       end
    end
    println("gblocksizes:")
    glb=zeros(Int,1,m)
    gblocks=Array{Array{Array{Int,1},1}}(undef, m)
    gblocksize=Array{Array{Int,2}}(undef, m)
    gcl=zeros(Int,1,m)
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
           println("-----------------------------------------------")
        else
           gub=unique(gblocksize[k])
           gsizes=[sum(gblocksize[k].== i) for i in gub]
           println("$gub\n$gsizes")
           println("-----------------------------------------------")
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
   println("-----------------------------------------------")
   println("gblocksizes:")
   gblocks=Array{Array{Array{Int,1},1}}(undef, m)
   gblocksize=Array{Array{Int,1}}(undef, m)
   gcl=zeros(Int,1,m)
   for k=1:m
       gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
       gblocksize[k]=[size(gbasis[k],2)]
       gcl[k]=1
       gbk=gblocksize[k]
       println("$gbk\n[1]")
       println("-----------------------------------------------")
   end
   return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
else
   nub=unique(fblocksize)
   nsizes=[sum(fblocksize.== i) for i in nub]
   if nub!=ub||nsizes!=sizes
      ub=nub
      sizes=nsizes
      println("fblocksizes:\n$ub\n$sizes")
      println("-----------------------------------------------")
   else
      println("No higher block hierarchy")
      return nothing,nothing,nothing,nothing,nothing,nothing,nothing,nothing,0
   end
end
println("gblocksizes:")
glb=zeros(Int,1,m)
gblocks=Array{Array{Array{Int,1},1}}(undef, m)
gblocksize=Array{Array{Int,2}}(undef, m)
gcl=zeros(Int,1,m)
for k=1:m
    glb[k]=size(gbasis[k],2)
    gG=SimpleGraph(glb[k])
    for i = 1:glb[k]
        for j = i:glb[k]
            r=1
            while r<=lt[k+1]
                  bi=ssupp[k+1][:,r]+gbasis[k][:,i]+gbasis[k][:,j]
                  if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(ofsupp,lfo,bi,n)!=0
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
       println("-----------------------------------------------")
    else
       gub=unique(gblocksize[k])
       gsizes=[sum(gblocksize[k].== i) for i in gub]
       println("$gub\n$gsizes")
       println("-----------------------------------------------")
    end
end
end
return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,1
end

function blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize;numeq=0,QUIET=true)
    fsupp=zeros(UInt8,n,1)
    fsupp=zeros(UInt8,n,Int(sum(fblocksize.^2+fblocksize)/2))
    k=1
    for i=1:fcl
        for j=1:fblocksize[i]
            for r=j:fblocksize[i]
                bi=fbasis[:,fblocks[i][j]]+fbasis[:,fblocks[i][r]]
                fsupp[:,k]=bi
                k+=1
            end
        end
    end
    supp1=fsupp
    for k=1:m
        gsupp=zeros(UInt8,n,lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2))
        l=1
        for i=1:gcl[k]
            for j=1:gblocksize[k][i]
                for r=j:gblocksize[k][i]
                    for s=1:lt[k+1]
                        bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                        gsupp[:,l]=bi
                        l+=1
                    end
                end
            end
        end
        supp1=[supp1 gsupp]
    end
    supp1=sortslices(supp1,dims=2)
    supp1=unique(supp1,dims=2)
    lsupp1=size(supp1,2)
    model=Model(with_optimizer(Mosek.Optimizer, QUIET=QUIET))
    cons=[AffExpr(0) for i=1:lsupp1]
    pos=Array{Any}(undef, fcl)
    gram=Array{Any}(undef, fcl)
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
    for k=1:m-numeq
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
               gpos[k][i]=@variable(model, [1:gblocksize[k][i], 1:gblocksize[k][i]], PSD)
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
    if numeq>0
    for k=m-numeq+1:m
        gpos[k]=Array{Any}(undef, gcl[k])
        for i=1:gcl[k]
            if gblocksize[k][i]==1
               gpos[k][i]=@variable(model)
               for s=1:lt[k+1]
                   bi=ssupp[k+1][:,[s]]+UInt8(2)*gbasis[k][:,gblocks[k][i]]
                   Locb=bfind(supp1,lsupp1,bi,n)
                   cons[Locb]+=coe[k+1][s]*gpos[k][i]
               end
            else
               gpos[k][i]=@variable(model, [1:gblocksize[k][i], 1:gblocksize[k][i]])
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
        println("optimum = $objv")
        for i=1:fcl
            gram[i]=value.(pos[i])
        end
    else
        objv=objective_value(model)
        println("$status")
        println("optimum = $objv")
        for i=1:fcl
            gram[i]=value.(pos[i])
        end
    end
    return objv,fsupp,gram
end
