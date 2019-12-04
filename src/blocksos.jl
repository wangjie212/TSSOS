function get_acliques(n,supp,basis)
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
blocks=maximal_cliques(G)
cl=size(blocks,1)
blocksize=zeros(Int,1,cl)
for i=1:cl
    blocksize[i]=length(blocks[i])
end
if cl==1
   println("Unblockable")
   return 0,0,0,0,0,0,0
else
   ub=unique(blocksize)
   sizes=[sum(blocksize.== i) for i in ub]
   println("$ub\n$sizes")
   return blocks,cl,blocksize,ub,sizes,1
end
end

function blockcsos(n,m,ssupp,coe,lt,basis,blocks,cl,blocksize)
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
    model=Model(with_optimizer(Mosek.Optimizer, QUIET=false))
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
    @variable(model, gpos[1:m])
    for i=1:m
        for j=1:lt[i+1]
            Locb=bfind(supp1,lsupp1,ssupp[i+1][:,[j]],n)
            cons[Locb]+=coe[i+1][j]*gpos[i]
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
    if status == MOI.OPTIMAL
       objv = objective_value(model)
       println("optimum=\n$objv")
    else
       objv = 0
       println("Failed, status=\n$status")
    end
    return objv,supp1
end
