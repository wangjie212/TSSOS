mutable struct cdata_type
    n
    m
    x
    pop
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
    numeq
    gblocks
    gcl
    gblocksize
end

function blockcpop_first(pop,x,d;method="block",reducebasis=0,numeq=0,QUIET=false,dense=10,chor_alg="amd",solve=true)
    n=length(x)
    m=length(pop)-1
    dg=zeros(UInt8,1,m)
    coe=Vector{Vector{Float64}}(undef, m+1)
    ssupp=Vector{Array{UInt8,2}}(undef, m+1)
    lt=Vector{UInt32}(undef,m+1)
    for k=1:m+1
        mon=monomials(pop[k])
        coe[k]=coefficients(pop[k])
        lt[k]=length(mon)
        ssupp[k]=zeros(UInt8,n,lt[k])
        for i=1:lt[k]
            for j=1:n
                @inbounds ssupp[k][j,i]=MultivariatePolynomials.degree(mon[i],x[j])
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
    gbasis=Vector{Array{UInt8,2}}(undef,m)
    for k=1:m
        gbasis[k]=get_basis(n,d-Int(ceil(dg[k]/2)))
    end
    if method=="block"&&reducebasis==0
       fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis,QUIET=QUIET)
    elseif method=="block"&&reducebasis==1
        flag=1
        while flag==1
              fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis,reduce=1,QUIET=QUIET)
              tsupp=[ssupp[1] zeros(UInt8,n,1)]
              for k=1:m
                  gsupp=zeros(UInt8,n,lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2))
                  l=1
                  for i=1:gcl[k]
                      for j=1:gblocksize[k][i]
                          for r=j:gblocksize[k][i]
                              for s=1:lt[k+1]
                                  @inbounds bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                                  @inbounds gsupp[:,l]=bi
                                  l+=1
                              end
                          end
                      end
                  end
                  tsupp=[tsupp gsupp]
              end
              fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize)
        end
    elseif method=="chordal"&&reducebasis==0
        fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis,dense=dense,QUIET=QUIET,alg=chor_alg)
    else
        flag=1
        while flag==1
              fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis,reduce=1,dense=dense,QUIET=QUIET,alg=chor_alg)
              tsupp=[ssupp[1] zeros(UInt8,n,1)]
              for k=1:m
                  gsupp=zeros(UInt8,n,lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2))
                  l=1
                  for i=1:gcl[k]
                      for j=1:gblocksize[k][i]
                          for r=j:gblocksize[k][i]
                              for s=1:lt[k+1]
                                  @inbounds bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                                  @inbounds gsupp[:,l]=bi
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
    opt,fsupp,Gram=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,numeq=numeq,QUIET=QUIET,solve=solve)
    sol=extract_solutions(n,m,x,d,pop,numeq,opt,fbasis,fblocks,fcl,fblocksize,Gram,method=method)
    data=cdata_type(n,m,x,pop,ssupp,coe,lt,d,dg,fbasis,gbasis,fsupp,ub,sizes,numeq,gblocks,gcl,gblocksize)
    return opt,sol,data
end

function blockcpop_higher!(data;method="block",reducebasis=0,QUIET=false,dense=10,chor_alg="amd",solve=true)
    n=data.n
    m=data.m
    x=data.x
    pop=data.pop
    ssupp=data.ssupp
    coe=data.coe
    lt=data.lt
    d=data.d
    dg=data.dg
    gbasis=data.gbasis
    fsupp=data.fsupp
    ub=data.ub
    sizes=data.sizes
    numeq=data.numeq
    gblocks=data.gblocks
    gcl=data.gcl
    gblocksize=data.gblocksize
    opt=nothing
    if method=="block"&&reducebasis==0
       fbasis=data.fbasis
       fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes,QUIET=QUIET)
    elseif method=="block"&&reducebasis==1
        fbasis=get_basis(n,d)
        flag=1
        while flag==1
              fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes,reduce=1,QUIET=QUIET)
              tsupp=[ssupp[1] zeros(UInt8,n,1)]
              for k=1:m
                  gsupp=zeros(UInt8,n,lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2))
                  l=1
                  for i=1:gcl[k]
                      for j=1:gblocksize[k][i]
                          for r=j:gblocksize[k][i]
                              for s=1:lt[k+1]
                                  @inbounds bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                                  @inbounds gsupp[:,l]=bi
                                  l+=1
                              end
                          end
                      end
                  end
                  tsupp=[tsupp gsupp]
              end
              fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize)
        end
    elseif method=="chordal"&&reducebasis==0
        fbasis=data.fbasis
        fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes,dense=dense,QUIET=QUIET,alg=chor_alg)
    else
        fbasis=get_basis(n,d)
        flag=1
        while flag==1
              fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes,reduce=1,dense=dense,QUIET=QUIET,alg=chor_alg)
              tsupp=[ssupp[1] zeros(UInt8,n,1)]
              for k=1:m
                  gsupp=zeros(UInt8,n,lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2))
                  l=1
                  for i=1:gcl[k]
                      for j=1:gblocksize[k][i]
                          for r=j:gblocksize[k][i]
                              for s=1:lt[k+1]
                                  @inbounds bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                                  @inbounds gsupp[:,l]=bi
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
    sol=nothing
    if status==1
        opt,fsupp,Gram=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,numeq=numeq,QUIET=QUIET,solve=solve)
        sol=extract_solutions(n,m,x,d,pop,numeq,opt,fbasis,fblocks,fcl,fblocksize,Gram,method=method)
    end
    data.fsupp=fsupp
    data.fbasis=fbasis
    data.ub=ub
    data.sizes=sizes
    data.gblocks=gblocks
    data.gcl=gcl
    data.gblocksize=gblocksize
    return opt,sol,data
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
                        @inbounds bi=basis[:,blocks[i][j]]+basis[:,blocks[i][r]]
                        supp1=[supp1 bi]
                    end
                end
            end
        end
        supp1=sortslices(supp1,dims=2)
        supp1=copy(unique(supp1,dims=2))
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

function get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis;reduce=0,QUIET=QUIET)
    gblocks=Vector{Vector{Vector{UInt16}}}(undef,m)
    gblocksize=Vector{Vector{Int}}(undef, m)
    gcl=Vector{UInt16}(undef,m)
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
        fblocksize=Vector{Int}(undef,fcl)
        for i=1:fcl
            fblocksize[i]=length(fblocks[i])
        end
        ub=unique(fblocksize)
        sizes=[sum(fblocksize.== i) for i in ub]
        if QUIET==false
            println("fblocksizes:\n$ub\n$sizes")
            println("-----------------------------------------------")
            println("gblocksizes:")
        end
        if fcl==1
            for k=1:m
                gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                gblocksize[k]=[size(gbasis[k],2)]
                gcl[k]=1
                if QUIET==false
                    gbk=gblocksize[k]
                    println("$gbk\n[1]")
                    println("-----------------------------------------------")
                end
            end
            return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1]
        else
            for k=1:m
                glb=size(gbasis[k],2)
                gG=SimpleGraph(glb)
                for i = 1:glb
                    for j = i:glb
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
                gblocksize[k]=Vector{Int}(undef,gcl[k])
                for i=1:gcl[k]
                    gblocksize[k][i]=length(gblocks[k][i])
                end
                if QUIET==false
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
        fblocksize=Vector{Int}(undef,fcl)
        for i=1:fcl
            fblocksize[i]=length(fblocks[i])
        end
        ub=unique(fblocksize)
        sizes=[sum(fblocksize.== i) for i in ub]
        if QUIET==false
            println("fblocksizes:\n$ub\n$sizes")
            println("-----------------------------------------------")
            println("gblocksizes:")
        end
        if fcl==1
            for k=1:m
                gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                gblocksize[k]=[size(gbasis[k],2)]
                gcl[k]=1
                if QUIET==false
                    gbk=gblocksize[k]
                    println("$gbk\n[1]")
                    println("-----------------------------------------------")
                end
            end
            return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1]
        else
            for k=1:m
                glb=size(gbasis[k],2)
                gG=SimpleGraph(glb)
                for i = 1:glb
                    for j = i:glb
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
                gblocksize[k]=Vector{Int}(undef,gcl[k])
                for i=1:gcl[k]
                    gblocksize[k][i]=length(gblocks[k][i])
                end
                gub=unique(gblocksize[k])
                gsizes=[sum(gblocksize[k].== i) for i in gub]
                if QUIET==false
                     println("$gub\n$gsizes")
                     println("-----------------------------------------------")
                end
            end
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes
end

function get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes;reduce=0,QUIET=QUIET)
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
        fblocksize=Vector{Int}(undef,fcl)
        for i=1:fcl
            fblocksize[i]=length(fblocks[i])
        end
        if fcl==1
            if length(sizes)==1&&sizes[1]==1
                if QUIET==false
                   println("No higher block hierarchy!")
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],0
            else
                if QUIET==false
                    println("fblocksizes:\n$fblocksize\n[1]")
                    println("-----------------------------------------------")
                    println("gblocksizes:")
                end
                for k=1:m
                    gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                    gblocksize[k]=[size(gbasis[k],2)]
                    gcl[k]=1
                    if QUIET==false
                        gbk=gblocksize[k]
                        println("$gbk\n[1]")
                        println("-----------------------------------------------")
                    end
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
            end
        else
           nub=unique(fblocksize)
           nsizes=[sum(fblocksize.== i) for i in nub]
           if nub!=ub||nsizes!=sizes
                if QUIET==false
                    println("fblocksizes:\n$nub\n$nsizes")
                    println("-----------------------------------------------")
                end
           else
               if QUIET==false
                  println("No higher block hierarchy!")
               end
              return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,0
           end
        end
        if QUIET==false
            println("gblocksizes:")
        end
        for k=1:m
            glb=size(gbasis[k],2)
            gG=SimpleGraph(glb)
            for i = 1:glb
                for j = i:glb
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
            gblocksize[k]=Vector{Int}(undef,gcl[k])
            for i=1:gcl[k]
                gblocksize[k][i]=length(gblocks[k][i])
            end
            if QUIET==false
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
        fblocksize=Vector{Int}(undef,fcl)
        for i=1:fcl
            fblocksize[i]=length(fblocks[i])
        end
        if fcl==1
            if length(sizes)==1&&sizes[1]==1
                if QUIET==false
                   println("No higher block hierarchy!")
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],0
            else
                if QUIET==false
                    println("fblocksizes:\n$fblocksize\n[1]")
                    println("-----------------------------------------------")
                    println("gblocksizes:")
                end
                for k=1:m
                    gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                    gblocksize[k]=[size(gbasis[k],2)]
                    gcl[k]=1
                    if QUIET==false
                        gbk=gblocksize[k]
                        println("$gbk\n[1]")
                        println("-----------------------------------------------")
                    end
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
            end
        else
           nub=unique(fblocksize)
           nsizes=[sum(fblocksize.== i) for i in nub]
           if nub!=ub||nsizes!=sizes
              if QUIET==false
                  println("fblocksizes:\n$nub\n$nsizes")
                  println("-----------------------------------------------")
              end
           else
               if QUIET==false
                  println("No higher block hierarchy!")
               end
              return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,0
           end
        end
        if QUIET==false
            println("gblocksizes:")
        end
        for k=1:m
            glb=size(gbasis[k],2)
            gG=SimpleGraph(glb)
            for i = 1:glb
                for j = i:glb
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
            gblocksize[k]=Vector{Int}(undef,gcl[k])
            for i=1:gcl[k]
                gblocksize[k][i]=length(gblocks[k][i])
            end
            gub=unique(gblocksize[k])
            gsizes=[sum(gblocksize[k].== i) for i in gub]
            if QUIET==false
                println("$gub\n$gsizes")
                println("-----------------------------------------------")
            end
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,1
end

function get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis;reduce=0,dense=10,QUIET=QUIET,alg="amd")
    gblocks=Vector{Vector{Vector{UInt16}}}(undef,m)
    gblocksize=Vector{Vector{Int}}(undef, m)
    gcl=Vector{UInt16}(undef,m)
    if reduce==1
        supp1=[supp 2*fbasis]
        supp1=unique(supp1,dims=2)
        supp1=sortslices(supp1,dims=2)
        lsupp1=size(supp1,2)
        flb=size(fbasis,2)
        if alg=="greedy"
            G=CGraph()
            for i=1:flb
                cadd_node!(G)
            end
        else
            A=zeros(UInt8,flb,flb)
        end
        for i = 1:flb
            for j = i+1:flb
                bi=fbasis[:,i]+fbasis[:,j]
                if bfind(supp1,lsupp1,bi,n)!=0
                    if alg=="greedy"
                        cadd_edge!(G,i,j)
                    else
                        A[i,j]=1
                        A[j,i]=1
                    end
                end
            end
        end
        if alg=="greedy"
            fblocks,fcl,fblocksize=chordal_extension(G, GreedyFillIn())
        else
            fblocks,fcl,fblocksize=cliquesFromSpMatD(A,dense=dense)
        end
        ub=unique(fblocksize)
        sizes=[sum(fblocksize.== i) for i in ub]
        if QUIET==false
            println("fblocksizes:\n$ub\n$sizes")
            println("-----------------------------------------------")
            println("gblocksizes:")
        end
        if fcl==1
            for k=1:m
                gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                gblocksize[k]=[size(gbasis[k],2)]
                gcl[k]=1
                if QUIET==false
                    gbk=gblocksize[k]
                    println("$gbk\n[1]")
                    println("-----------------------------------------------")
                end
            end
            return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1]
        else
            for k=1:m
                glb=size(gbasis[k],2)
                if alg=="greedy"
                    G=CGraph()
                    for i=1:glb
                        cadd_node!(G)
                    end
                else
                    A=zeros(UInt8,glb,glb)
                end
                for i = 1:glb
                    for j = i+1:glb
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
                            if alg=="greedy"
                                cadd_edge!(G,i,j)
                            else
                                A[i,j]=1
                                A[j,i]=1
                            end
                        end
                    end
                end
                if alg=="greedy"
                    gblocks[k],gcl[k],gblocksize[k]=chordal_extension(G, GreedyFillIn())
                else
                    gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A,dense=dense)
                end
                if QUIET==false
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
        if alg=="greedy"
            G=CGraph()
            for i=1:flb
                cadd_node!(G)
            end
        else
            A=zeros(UInt8,flb,flb)
        end
        for i = 1:flb
            for j = i+1:flb
                bi=fbasis[:,i]+fbasis[:,j]
                if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(osupp,lo,bi,n)!=0
                    if alg=="greedy"
                        cadd_edge!(G,i,j)
                    else
                        A[i,j]=1
                        A[j,i]=1
                    end
                end
            end
        end
        if alg=="greedy"
            fblocks,fcl,fblocksize=chordal_extension(G, GreedyFillIn())
        else
            fblocks,fcl,fblocksize=cliquesFromSpMatD(A,dense=dense)
        end
        ub=unique(fblocksize)
        sizes=[sum(fblocksize.== i) for i in ub]
        if QUIET==false
            println("fblocksizes:\n$ub\n$sizes")
            println("-----------------------------------------------")
            println("gblocksizes:")
        end
        if fcl==1
            for k=1:m
                gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                gblocksize[k]=[size(gbasis[k],2)]
                gcl[k]=1
                if QUIET==false
                    gbk=gblocksize[k]
                    println("$gbk\n[1]")
                    println("-----------------------------------------------")
                end
            end
            return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1]
        else
            for k=1:m
                glb=size(gbasis[k],2)
                if alg=="greedy"
                    G=CGraph()
                    for i=1:glb
                        cadd_node!(G)
                    end
                else
                    A=zeros(UInt8,glb,glb)
                end
                for i = 1:glb
                    for j = i+1:glb
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
                            if alg=="greedy"
                                cadd_edge!(G,i,j)
                            else
                                A[i,j]=1
                                A[j,i]=1
                            end
                        end
                    end
                end
                if alg=="greedy"
                    gblocks[k],gcl[k],gblocksize[k]=chordal_extension(G, GreedyFillIn())
                else
                    gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A,dense=dense)
                end
                if QUIET==false
                    gub=unique(gblocksize[k])
                    gsizes=[sum(gblocksize[k].== i) for i in gub]
                    println("$gub\n$gsizes")
                    println("-----------------------------------------------")
                end
            end
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes
end

function get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes;reduce=0,dense=10,QUIET=QUIET,alg="amd")
    if reduce==1
        supp1=[fsupp 2*fbasis]
        supp1=unique(supp1,dims=2)
        supp1=sortslices(supp1,dims=2)
        lsupp1=size(supp1,2)
        flb=size(fbasis,2)
        if alg=="greedy"
            G=CGraph()
            for i=1:flb
                cadd_node!(G)
            end
        else
            A=zeros(UInt8,flb,flb)
        end
        for i = 1:flb
            for j = i+1:flb
                bi=fbasis[:,i]+fbasis[:,j]
                if bfind(supp1,lsupp1,bi,n)!=0
                    if alg=="greedy"
                        cadd_edge!(G,i,j)
                    else
                        A[i,j]=1
                        A[j,i]=1
                    end
                end
            end
        end
        if alg=="greedy"
            fblocks,fcl,fblocksize=chordal_extension(G, GreedyFillIn())
        else
            fblocks,fcl,fblocksize=cliquesFromSpMatD(A,dense=dense)
        end
        if fcl==1
            if length(sizes)==1&&sizes[1]==1
                if QUIET==false
                   println("No higher chordal hierarchy!")
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],0
            else
                if QUIET==false
                    println("fblocksizes:\n$fblocksize\n[1]")
                    println("-----------------------------------------------")
                    println("gblocksizes:")
                end
                for k=1:m
                    gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                    gblocksize[k]=[size(gbasis[k],2)]
                    gcl[k]=1
                    if QUIET==false
                        gbk=gblocksize[k]
                        println("$gbk\n[1]")
                        println("-----------------------------------------------")
                    end
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
            end
        else
           nub=unique(fblocksize)
           nsizes=[sum(fblocksize.== i) for i in nub]
           if nub!=ub||nsizes!=sizes
               if QUIET==false
                   println("fblocksizes:\n$nub\n$nsizes")
                   println("-----------------------------------------------")
               end
           else
               if QUIET==false
                  println("No higher chordal hierarchy!")
               end
              return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,0
           end
        end
        if QUIET==false
            println("gblocksizes:")
        end
        for k=1:m
            glb=size(gbasis[k],2)
            if alg=="greedy"
                G=CGraph()
                for i=1:glb
                    cadd_node!(G)
                end
            else
                A=zeros(UInt8,glb,glb)
            end
            for i = 1:glb
                for j = i+1:glb
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
                        if alg=="greedy"
                            cadd_edge!(G,i,j)
                        else
                            A[i,j]=1
                            A[j,i]=1
                        end
                    end
                end
            end
            if alg=="greedy"
                gblocks[k],gcl[k],gblocksize[k]=chordal_extension(G, GreedyFillIn())
            else
                gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A,dense=dense)
            end
            if QUIET==false
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
        if alg=="greedy"
            G=CGraph()
            for i=1:flb
                cadd_node!(G)
            end
        else
            A=zeros(UInt8,flb,flb)
        end
        for i = 1:flb
            for j = i+1:flb
                bi=fbasis[:,i]+fbasis[:,j]
                if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(ofsupp,lfo,bi,n)!=0
                    if alg=="greedy"
                        cadd_edge!(G,i,j)
                    else
                        A[i,j]=1
                        A[j,i]=1
                    end
                end
            end
        end
        if alg=="greedy"
            fblocks,fcl,fblocksize=chordal_extension(G, GreedyFillIn())
        else
            fblocks,fcl,fblocksize=cliquesFromSpMatD(A,dense=dense)
        end
        if fcl==1
            if length(sizes)==1&&sizes[1]==1
                if QUIET==false
                   println("No higher chordal hierarchy!")
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],0
            else
                if QUIET==false
                    println("fblocksizes:\n$fblocksize\n[1]")
                    println("-----------------------------------------------")
                    println("gblocksizes:")
                end
                for k=1:m
                    gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                    gblocksize[k]=[size(gbasis[k],2)]
                    gcl[k]=1
                    if QUIET==false
                        gbk=gblocksize[k]
                        println("$gbk\n[1]")
                        println("-----------------------------------------------")
                    end
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
            end
        else
           nub=unique(fblocksize)
           nsizes=[sum(fblocksize.== i) for i in nub]
           if nub!=ub||nsizes!=sizes
               if QUIET==false
                   println("fblocksizes:\n$nub\n$nsizes")
                   println("-----------------------------------------------")
               end
           else
              if QUIET==false
                 println("No higher chordal hierarchy!")
              end
              return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,0
           end
        end
        if QUIET==false
            println("gblocksizes:")
        end
        for k=1:m
            glb=size(gbasis[k],2)
            if alg=="greedy"
                G=CGraph()
                for i=1:glb
                    cadd_node!(G)
                end
            else
                A=zeros(UInt8,glb,glb)
            end
            for i = 1:glb
                for j = i+1:glb
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
                        if alg=="greedy"
                            cadd_edge!(G,i,j)
                        else
                            A[i,j]=1
                            A[j,i]=1
                        end
                    end
                end
            end
            if alg=="greedy"
                gblocks[k],gcl[k],gblocksize[k]=chordal_extension(G, GreedyFillIn())
            else
                gblocks[k],gcl[k],gblocksize[k]=cliquesFromSpMatD(A,dense=dense)
            end
            if QUIET==false
               gub=unique(gblocksize[k])
               gsizes=[sum(gblocksize[k].== i) for i in gub]
               println("$gub\n$gsizes")
               println("-----------------------------------------------")
            end
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,1
end

function blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize;numeq=0,QUIET=true,solve=true)
    fsupp=zeros(UInt8,n,Int(sum(fblocksize.^2+fblocksize)/2))
    k=1
    for i=1:fcl
        for j=1:fblocksize[i]
            for r=j:fblocksize[i]
                @inbounds bi=fbasis[:,fblocks[i][j]]+fbasis[:,fblocks[i][r]]
                @inbounds fsupp[:,k]=bi
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
                        @inbounds bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                        @inbounds gsupp[:,l]=bi
                        l+=1
                    end
                end
            end
        end
        supp1=[supp1 gsupp]
    end
    supp1=unique(supp1,dims=2)
    supp1=sortslices(supp1,dims=2)
    objv=nothing
    gram=nothing
    if solve==true
        lsupp1=size(supp1,2)
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:lsupp1]
        pos=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, fcl)
        gram=Vector{Union{Float64,Array{Float64,2}}}(undef, fcl)
        for i=1:fcl
            bs=fblocksize[i]
            if bs==1
               @inbounds pos[i]=@variable(model, lower_bound=0)
               @inbounds bi=2*fbasis[:,fblocks[i]]
               Locb=bfind(supp1,lsupp1,bi,n)
               @inbounds cons[Locb]+=pos[i]
            else
               @inbounds pos[i]=@variable(model, [1:bs, 1:bs], PSD)
               for j=1:bs
                   for r=j:bs
                       @inbounds bi=fbasis[:,fblocks[i][j]]+fbasis[:,fblocks[i][r]]
                       Locb=bfind(supp1,lsupp1,bi,n)
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
                        @inbounds bi=ssupp[k+1][:,[s]]+2*gbasis[k][:,gblocks[k][i]]
                        Locb=bfind(supp1,lsupp1,bi,n)
                        @inbounds cons[Locb]+=coe[k+1][s]*gpos[k][i]
                    end
                else
                    if k<=m-numeq
                       gpos[k][i]=@variable(model, [1:bs, 1:bs], PSD)
                    else
                       gpos[k][i]=@variable(model, [1:bs, 1:bs], Symmetric)
                    end
                    for j=1:bs
                        for r=j:bs
                            for s=1:lt[k+1]
                                @inbounds bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                                Locb=bfind(supp1,lsupp1,bi,n)
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
        end
        bc=zeros(1,lsupp1)
        for i=1:lt[1]
            Locb=bfind(supp1,lsupp1,ssupp[1][:,i],n)
            if Locb==0
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing
            else
               bc[Locb]=coe[1][i]
           end
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
        else
            objv=objective_value(model)
            println("termination status: $status")
            sstatus=primal_status(model)
            println("solution status: $sstatus")
            println("optimum = $objv")
        end
        for i=1:fcl
            gram[i]=value.(pos[i])
        end
    end
    return objv,fsupp,gram
end
