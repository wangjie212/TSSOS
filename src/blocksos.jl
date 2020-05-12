# function get_acliques(n,supp,basis;QUIET=false)
#     osupp=odd_supp(n,supp)
#     osupp=sortslices(osupp,dims=2)
#     lo=size(osupp,2)
#     lb=size(basis,2)
#     G=SimpleGraph(lb)
#     for i = 1:lb
#         for j = i:lb
#             bi=basis[:,i]+basis[:,j]
#             if sum(Int[iseven(bi[k]) for k=1:n])==n||bfind(osupp,lo,bi,n)!=0
#                add_edge!(G,i,j)
#             end
#         end
#     end
#     blocks=maximal_cliques(G)
#     cl=size(blocks,1)
#     blocksize=zeros(Int,1,cl)
#     for i=1:cl
#         blocksize[i]=length(blocks[i])
#     end
#     ub=unique(blocksize)
#     sizes=[sum(blocksize.== i) for i in ub]
#     if QUIET==false
#         println("The size of blocks:\n$nub\n$nsizes")
#     end
#     return blocks,cl,blocksize,ub,sizes
# end

function blockcsos(n,m,ssupp,coe,lt,basis,blocks,cl,blocksize;QUIET=false,solve=true)
    supp1=zeros(UInt8,n,Int(sum(blocksize.^2+blocksize)/2))
    k=1
    for i=1:cl
        for j=1:blocksize[i]
            for r=j:blocksize[i]
                @inbounds bi=basis[:,blocks[i][j]]+basis[:,blocks[i][r]]
                @inbounds supp1[:,k]=bi
                k+=1
            end
        end
    end
    supp1=sortslices(supp1,dims=2)
    supp1=unique(supp1,dims=2)
    objv=nothing
    if solve==true
        lsupp1=size(supp1,2)
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:lsupp1]
        pos=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl)
        for i=1:cl
            bs=blocksize[i]
            if bs==1
                @inbounds pos[i]=@variable(model, lower_bound=0)
                @inbounds bi=2*basis[:,blocks[i]]
                Locb=bfind(supp1,lsupp1,bi,n)
                @inbounds cons[Locb]+=pos[i]
            else
                @inbounds pos[i]=@variable(model, [1:bs, 1:bs], PSD)
                for j=1:blocksize[i]
                    for r=j:blocksize[i]
                        @inbounds bi=basis[:,blocks[i][j]]+basis[:,blocks[i][r]]
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
           println("optimum = $objv")
        else
           objv = objective_value(model)
           println("termination status: $status")
           sstatus=primal_status(model)
           println("solution status: $sstatus")
           println("optimum = $objv")
        end
    end
    return objv,supp1
end
