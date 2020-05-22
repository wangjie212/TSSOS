mutable struct cdata_type
    n
    nb
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

function blockcpop_first(pop,x,d;nb=0,method="block",reducebasis=false,numeq=0,QUIET=false,dense=10,chor_alg="amd",solve=true,extra_sos=false,solution=false,tol=1e-5,merge=false)
    n=length(x)
    m=length(pop)-1
    dg=zeros(Int,m)
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
    fbasis=get_basis(n,d,nb=nb)
    gbasis=Vector{Array{UInt8,2}}(undef,m)
    for k=1:m
        gbasis[k]=get_basis(n,d-Int(ceil(dg[k]/2)),nb=nb)
    end
    if method=="block"&&reducebasis==false
       fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis,nb=nb,QUIET=QUIET)
    elseif method=="block"&&reducebasis==true
        flag=1
        while flag==1
              fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis,nb=nb,reduce=true,QUIET=QUIET)
              gsupp=get_gsupp(n,m,lt,ssupp,gbasis,gblocks,gcl,gblocksize,nb=nb)
              tsupp=[ssupp[1] zeros(UInt8,n,1)]
              tsupp=[tsupp gsupp]
              fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize,nb=nb)
        end
    elseif method=="chordal"&&reducebasis==false
        fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis,nb=nb,dense=dense,QUIET=QUIET,alg=chor_alg,merge=merge)
    else
        flag=1
        while flag==1
              fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes=get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis,nb=nb,reduce=true,dense=dense,QUIET=QUIET,alg=chor_alg,merge=merge)
              gsupp=get_gsupp(n,m,lt,ssupp,gbasis,gblocks,gcl,gblocksize,nb=nb)
              tsupp=[ssupp[1] zeros(UInt8,n,1)]
              tsupp=[tsupp gsupp]
              fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize,nb=nb)
        end
    end
    opt,fsupp,moment=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nb=nb,numeq=numeq,QUIET=QUIET,solve=solve,solution=solution,extra_sos=extra_sos)
    sol=nothing
    if solution==true
        sol=extract_solutions(moment,opt,n,m,pop,x,numeq=numeq,tol=tol)
    end
    data=cdata_type(n,nb,m,x,pop,ssupp,coe,lt,d,dg,fbasis,gbasis,fsupp,ub,sizes,numeq,gblocks,gcl,gblocksize)
    return opt,sol,data
end

function blockcpop_higher!(data;nb=0,method="block",reducebasis=false,QUIET=false,dense=10,chor_alg="amd",solve=true,extra_sos=false,solution=false,tol=1e-5,merge=false)
    n=data.n
    nb=data.nb
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
    if method=="block"&&reducebasis==false
       fbasis=data.fbasis
       fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes,nb=nb,QUIET=QUIET)
   elseif method=="block"&&reducebasis==true
        fbasis=get_basis(n,d,nb=nb)
        flag=1
        while flag==1
              fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes,nb=nb,reduce=true,QUIET=QUIET)
              gsupp=get_gsupp(n,m,lt,ssupp,gbasis,gblocks,gcl,gblocksize,nb=nb)
              tsupp=[ssupp[1] zeros(UInt8,n,1)]
              tsupp=[tsupp gsupp]
              fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize,nb=nb)
        end
    elseif method=="chordal"&&reducebasis==false
        fbasis=data.fbasis
        fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes,nb=nb,dense=dense,QUIET=QUIET,alg=chor_alg,merge=merge)
    else
        fbasis=get_basis(n,d,nb=nb)
        flag=1
        while flag==1
              fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes,nb=nb,reduce=true,dense=dense,QUIET=QUIET,alg=chor_alg,merge=merge)
              gsupp=get_gsupp(n,m,lt,ssupp,gbasis,gblocks,gcl,gblocksize,nb=nb)
              tsupp=[ssupp[1] zeros(UInt8,n,1)]
              tsupp=[tsupp gsupp]
              fbasis,flag=reducebasis!(n,tsupp,fbasis,fblocks,fcl,fblocksize,nb=nb)
        end
    end
    sol=nothing
    if status==1
        opt,fsupp,moment=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nb=nb,numeq=numeq,QUIET=QUIET,solve=solve,solution=solution,extra_sos=extra_sos)
        if solution==true
            sol=extract_solutions(moment,opt,n,m,pop,x,numeq=numeq,tol=tol)
        end
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

function get_gsupp(n,m,lt,ssupp,gbasis,gblocks,gcl,gblocksize;nb=0)
    gsupp=zeros(UInt8,n,sum(lt[k+1]*Int(sum(gblocksize[k].^2+gblocksize[k])/2) for k=1:m))
    l=1
    for k=1:m
        for i=1:gcl[k]
            for j=1:gblocksize[k][i]
                for r=j:gblocksize[k][i]
                    for s=1:lt[k+1]
                        # @inbounds bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                        @inbounds bi=bin_add(gbasis[k][:,gblocks[k][i][j]],gbasis[k][:,gblocks[k][i][r]],nb)
                        @inbounds bi=bin_add(bi,ssupp[k+1][:,s],nb)
                        @inbounds gsupp[:,l]=bi
                        l+=1
                    end
                end
            end
        end
    end
    return gsupp
end

function reducebasis!(n,supp,basis,blocks,cl,blocksize;nb=0)
    esupp=even_supp(supp)
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
                        # @inbounds bi=basis[:,blocks[i][j]]+basis[:,blocks[i][r]]
                        @inbounds bi=bin_add(basis[:,blocks[i][j]],basis[:,blocks[i][r]],nb)
                        supp1=[supp1 bi]
                    end
                end
            end
        end
        supp1=unique(supp1,dims=2)
        supp1=sortslices(supp1,dims=2)
        lsupp1=size(supp1,2)
        for i=1:cl
            lo=blocksize[i]
            indexb=[k for k=1:lo]
            j=1
            while lo>=j
                # bi=2*basis[:,blocks[i][indexb[j]]]
                bi=bin_add(basis[:,blocks[i][indexb[j]]],basis[:,blocks[i][indexb[j]]],nb)
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

function get_cblocks(n,m,supp,ssupp,lt,fbasis,gbasis;nb=0,reduce=false,QUIET=true)
    gblocks=Vector{Vector{Vector{UInt16}}}(undef,m)
    gblocksize=Vector{Vector{Int}}(undef, m)
    gcl=Vector{UInt16}(undef,m)
    if reduce==true||nb>0
        supp1=[supp bin_add(fbasis,fbasis,nb)]
        supp1=unique(supp1,dims=2)
        supp1=sortslices(supp1,dims=2)
        lsupp1=size(supp1,2)
        flb=size(fbasis,2)
        fG=SimpleGraph(flb)
        for i = 1:flb
            for j = i:flb
                bi=bin_add(fbasis[:,i],fbasis[:,j],nb)
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
            println("------------------------------------------------------")
            println("The sizes of blocks:\n$ub\n$sizes")
            println("------------------------------------------------------")
        end
        if fcl==1
            for k=1:m
                gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                gblocksize[k]=[size(gbasis[k],2)]
                gcl[k]=1
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
                              bi=bin_add(gbasis[k][:,i],gbasis[k][:,j],nb)
                              bi=bin_add(bi,ssupp[k+1][:,r],nb)
                              if bfind(supp1,lsupp1,bi,n)!=0
                                 break
                              else
                                  r+=1
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
            end
        end
    else
        osupp=odd_supp(supp)
        osupp=unique(osupp,dims=2)
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
            println("------------------------------------------------------")
            println("The sizes of blocks:\n$ub\n$sizes")
            println("------------------------------------------------------")
        end
        if fcl==1
            for k=1:m
                gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                gblocksize[k]=[size(gbasis[k],2)]
                gcl[k]=1
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
                                  r+=1
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
            end
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes
end

function get_chblocks!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes;nb=0,reduce=false,QUIET=true)
    if reduce==true||nb>0
        supp1=[fsupp bin_add(fbasis, fbasis, nb)]
        supp1=unique(supp1,dims=2)
        supp1=sortslices(supp1,dims=2)
        lsupp1=size(supp1,2)
        flb=size(fbasis,2)
        fG=SimpleGraph(flb)
        for i = 1:flb
            for j = i:flb
                # bi=fbasis[:,i]+fbasis[:,j]
                bi=bin_add(fbasis[:,i],fbasis[:,j],nb)
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
                   println("No higher TSSOS hierarchy!")
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],0
            else
                if QUIET==false
                    println("------------------------------------------------------")
                    println("The sizes of blocks:\n$fblocksize\n[1]")
                    println("------------------------------------------------------")
                end
                for k=1:m
                    gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                    gblocksize[k]=[size(gbasis[k],2)]
                    gcl[k]=1
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
            end
        else
           nub=unique(fblocksize)
           nsizes=[sum(fblocksize.== i) for i in nub]
           if nub!=ub||nsizes!=sizes
                if QUIET==false
                    println("------------------------------------------------------")
                    println("The sizes of blocks:\n$nub\n$nsizes")
                    println("------------------------------------------------------")
                end
           else
               if QUIET==false
                  println("No higher TSSOS hierarchy!")
               end
              return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,0
           end
        end
        for k=1:m
            glb=size(gbasis[k],2)
            gG=SimpleGraph(glb)
            for i = 1:glb
                for j = i:glb
                    r=1
                    while r<=lt[k+1]
                          # bi=ssupp[k+1][:,r]+gbasis[k][:,i]+gbasis[k][:,j]
                          bi=bin_add(gbasis[k][:,i],gbasis[k][:,j],nb)
                          bi=bin_add(bi,ssupp[k+1][:,r],nb)
                          if bfind(supp1,lsupp1,bi,n)!=0
                             break
                          else
                              r+=1
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
        end
    else
        ofsupp=odd_supp(fsupp)
        ofsupp=unique(ofsupp,dims=2)
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
                   println("No higher TSSOS hierarchy!")
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],0
            else
                if QUIET==false
                    println("------------------------------------------------------")
                    println("The sizes of blocks:\n$fblocksize\n[1]")
                    println("------------------------------------------------------")
                end
                for k=1:m
                    gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                    gblocksize[k]=[size(gbasis[k],2)]
                    gcl[k]=1
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
            end
        else
           nub=unique(fblocksize)
           nsizes=[sum(fblocksize.== i) for i in nub]
           if nub!=ub||nsizes!=sizes
              if QUIET==false
                  println("------------------------------------------------------")
                  println("The sizes of blocks:\n$nub\n$nsizes")
                  println("------------------------------------------------------")
              end
           else
               if QUIET==false
                  println("No higher TSSOS hierarchy!")
               end
              return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,0
           end
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
                             r+=1
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
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,1
end

function get_ccliques(n,m,supp,ssupp,lt,fbasis,gbasis;nb=0,reduce=false,dense=10,QUIET=true,alg="amd",merge=false)
    gblocks=Vector{Vector{Vector{UInt16}}}(undef,m)
    gblocksize=Vector{Vector{Int}}(undef, m)
    gcl=Vector{UInt16}(undef,m)
    if reduce==true||nb>0
        # supp1=[supp 2*fbasis]
        supp1=[supp bin_add(fbasis, fbasis, nb)]
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
                # bi=fbasis[:,i]+fbasis[:,j]
                bi=bin_add(fbasis[:,i],fbasis[:,j],nb)
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
        if merge==true
            fblocks,fcl,fblocksize=clique_merge!(fblocks,fcl,QUIET=true)
        end
        ub=unique(fblocksize)
        sizes=[sum(fblocksize.== i) for i in ub]
        if QUIET==false
            println("------------------------------------------------------")
            println("The sizes of blocks:\n$ub\n$sizes")
            println("------------------------------------------------------")
        end
        if fcl==1
            for k=1:m
                gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                gblocksize[k]=[size(gbasis[k],2)]
                gcl[k]=1
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
                              # bi=ssupp[k+1][:,r]+gbasis[k][:,i]+gbasis[k][:,j]
                              bi=bin_add(gbasis[k][:,i],gbasis[k][:,j],nb)
                              bi=bin_add(bi,ssupp[k+1][:,r],nb)
                              if bfind(supp1,lsupp1,bi,n)!=0
                                 break
                              else
                                  r+=1
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
                if merge==true
                    gblocks[k],gcl[k],gblocksize[k]=clique_merge!(gblocks[k],gcl[k],QUIET=true)
                end
            end
        end
    else
        osupp=odd_supp(supp)
        osupp=unique(osupp,dims=2)
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
        if merge==true
            fblocks,fcl,fblocksize=clique_merge!(fblocks,fcl,QUIET=true)
        end
        ub=unique(fblocksize)
        sizes=[sum(fblocksize.== i) for i in ub]
        if QUIET==false
            println("------------------------------------------------------")
            println("The sizes of blocks:\n$ub\n$sizes")
            println("------------------------------------------------------")
        end
        if fcl==1
            for k=1:m
                gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                gblocksize[k]=[size(gbasis[k],2)]
                gcl[k]=1
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
                                 r+=1
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
                if merge==true
                    gblocks[k],gcl[k],gblocksize[k]=clique_merge!(gblocks[k],gcl[k],QUIET=true)
                end
            end
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes
end

function get_chcliques!(n,m,ssupp,lt,fbasis,gbasis,fsupp,gblocks,gcl,gblocksize,ub,sizes;nb=0,reduce=false,dense=10,QUIET=true,alg="amd",merge=false)
    if reduce==true||nb>0
        # supp1=[fsupp 2*fbasis]
        supp1=[fsupp bin_add(fbasis, fbasis, nb)]
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
                # bi=fbasis[:,i]+fbasis[:,j]
                bi=bin_add(fbasis[:,i],fbasis[:,j],nb)
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
        if merge==true
            fblocks,fcl,fblocksize=clique_merge!(fblocks,fcl,QUIET=true)
        end
        if fcl==1
            if length(sizes)==1&&sizes[1]==1
                if QUIET==false
                   println("No higher TSSOS hierarchy!")
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],0
            else
                if QUIET==false
                    println("------------------------------------------------------")
                    println("The sizes of blocks:\n$fblocksize\n[1]")
                    println("------------------------------------------------------")
                end
                for k=1:m
                    gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                    gblocksize[k]=[size(gbasis[k],2)]
                    gcl[k]=1
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
            end
        else
           nub=unique(fblocksize)
           nsizes=[sum(fblocksize.== i) for i in nub]
           if nub!=ub||nsizes!=sizes
               if QUIET==false
                   println("------------------------------------------------------")
                   println("The sizes of blocks:\n$nub\n$nsizes")
                   println("------------------------------------------------------")
               end
           else
               if QUIET==false
                  println("No higher TSSOS hierarchy!")
               end
              return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,0
           end
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
                          # bi=ssupp[k+1][:,r]+gbasis[k][:,i]+gbasis[k][:,j]
                          bi=bin_add(gbasis[k][:,i],gbasis[k][:,j],nb)
                          bi=bin_add(bi,ssupp[k+1][:,r],nb)
                          if bfind(supp1,lsupp1,bi,n)!=0
                             break
                          else
                              r+=1
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
            if merge==true
                gblocks[k],gcl[k],gblocksize[k]=clique_merge!(gblocks[k],gcl[k],QUIET=true)
            end
        end
    else
        ofsupp=odd_supp(fsupp)
        ofsupp=unique(ofsupp,dims=2)
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
        if merge==true
            fblocks,fcl,fblocksize=clique_merge!(fblocks,fcl,QUIET=true)
        end
        if fcl==1
            if length(sizes)==1&&sizes[1]==1
                if QUIET==false
                   println("No higher TSSOS hierarchy!")
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],0
            else
                if QUIET==false
                    println("------------------------------------------------------")
                    println("The sizes of blocks:\n$fblocksize\n[1]")
                    println("------------------------------------------------------")
                end
                for k=1:m
                    gblocks[k]=[[i for i=1:size(gbasis[k],2)]]
                    gblocksize[k]=[size(gbasis[k],2)]
                    gcl[k]=1
                end
                return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,fblocksize,[1],1
            end
        else
           nub=unique(fblocksize)
           nsizes=[sum(fblocksize.== i) for i in nub]
           if nub!=ub||nsizes!=sizes
               if QUIET==false
                   println("------------------------------------------------------")
                   println("The sizes of blocks:\n$nub\n$nsizes")
                   println("------------------------------------------------------")
               end
           else
              if QUIET==false
                 println("No higher TSSOS hierarchy!")
              end
              return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,0
           end
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
                              r+=1
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
            if merge==true
                gblocks[k],gcl[k],gblocksize[k]=clique_merge!(gblocks[k],gcl[k],QUIET=true)
            end
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,1
end

function blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize;nb=0,numeq=0,QUIET=true,solve=true,solution=false,extra_sos=false)
    fsupp=zeros(UInt8,n,Int(sum(fblocksize.^2+fblocksize)/2))
    k=1
    for i=1:fcl
        for j=1:fblocksize[i]
            for r=j:fblocksize[i]
                # @inbounds bi=fbasis[:,fblocks[i][j]]+fbasis[:,fblocks[i][r]]
                @inbounds bi=bin_add(fbasis[:,fblocks[i][j]],fbasis[:,fblocks[i][r]],nb)
                @inbounds fsupp[:,k]=bi
                k+=1
            end
        end
    end
    objv=nothing
    moment=nothing
    if solve==true
        gsupp=get_gsupp(n,m,lt,ssupp,gbasis,gblocks,gcl,gblocksize,nb=nb)
        supp1=[fsupp gsupp]
        if extra_sos==true||solution==true
            supp1=[supp1 get_basis(n,2,nb=nb)]
        end
        supp1=unique(supp1,dims=2)
        supp1=sortslices(supp1,dims=2)
        lsupp1=size(supp1,2)
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:lsupp1]
        pos=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, fcl)
        for i=1:fcl
            if extra_sos==true||solution==true
                pos0=@variable(model, [1:n+1, 1:n+1], PSD)
                for j=1:n+1
                    for k=j:n+1
                        # @inbounds bi=fbasis[:,j]+fbasis[:,k]
                        @inbounds bi=bin_add(fbasis[:,j],fbasis[:,k],nb)
                        Locb=bfind(supp1,lsupp1,bi,n)
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
               # @inbounds bi=2*fbasis[:,fblocks[i]]
               @inbounds bi=bin_add(fbasis[:,fblocks[i]],fbasis[:,fblocks[i]],nb)
               Locb=bfind(supp1,lsupp1,bi,n)
               @inbounds cons[Locb]+=pos[i]
            else
               @inbounds pos[i]=@variable(model, [1:bs, 1:bs], PSD)
               for j=1:bs
                   for r=j:bs
                       # @inbounds bi=fbasis[:,fblocks[i][j]]+fbasis[:,fblocks[i][r]]
                       @inbounds bi=bin_add(fbasis[:,fblocks[i][j]],fbasis[:,fblocks[i][r]],nb)
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
                        # @inbounds bi=ssupp[k+1][:,[s]]+2*gbasis[k][:,gblocks[k][i]]
                        @inbounds bi=bin_add(gbasis[k][:,gblocks[k][i]],gbasis[k][:,gblocks[k][i]],nb)
                        @inbounds bi=bin_add(bi,ssupp[k+1][:,[s]],nb)
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
                                # @inbounds bi=ssupp[k+1][:,s]+gbasis[k][:,gblocks[k][i][j]]+gbasis[k][:,gblocks[k][i][r]]
                                @inbounds bi=bin_add(gbasis[k][:,gblocks[k][i][j]],gbasis[k][:,gblocks[k][i][r]],nb)
                                @inbounds bi=bin_add(bi,ssupp[k+1][:,s],nb)
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
        @variable(model, lower)
        cons[1]+=lower
        @constraint(model, con[i=1:lsupp1], cons[i]==bc[i])
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
            for j=1:n+1
                for k=j:n+1
                    # bi=fbasis[:,j]+fbasis[:,k]
                    bi=bin_add(fbasis[:,j],fbasis[:,k],nb)
                    Locb=bfind(supp1,lsupp1,bi,n)
                    moment[j,k]=dual_var[Locb]
                end
            end
            moment=Symmetric(moment,:U)
        end
    end
    return objv,fsupp,moment
end
