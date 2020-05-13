mutable struct mdata_type
    n
    m
    dg
    supp
    coe
    numeq
    rlorder
    supp0
    ssupp
    lt
    fbasis
    gbasis
    cql
    cliques
    cliquesize
    I
    ncc
    blocks
    cl
    blocksize
    ub
    sizes
end

function cs_tssos_first(n,m,dg,supp,coe,order;numeq=0,CS="amd",assign="min",TS="block",QUIET=false,solve=true,solution=false,extra_sos=true)
    cliques,cql,cliquesize=clique_cdecomp(n,m,dg,supp,order=order,alg=CS)
    I,ncc=assign_constraint(m,supp,cliques,cql,cliquesize,assign=assign)
    rlorder=init_order(dg,I,cql,order=order)
    if TS==false
        opt,supp0,supp1,measure,moment=blockcpop_mix(n,m,dg,rlorder,supp,coe,cliques,cql,cliquesize,I,ncc,nothing,nothing,nothing,numeq=numeq,mix=false,QUIET=QUIET,solve=solve,solution=solution,extra_sos=false)
        blocks=nothing
        cl=nothing
        blocksize=nothing
        fbasis=nothing
        gbasis=nothing
        ssupp=nothing
        lt=nothing
        ub=nothing
        sizes=nothing
    else
        if TS=="block"
            blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis=get_cblocks_mix(dg,I,rlorder,m,supp,cliques,cql,cliquesize,method="block",chor_alg=nothing)
        else
            blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis=get_cblocks_mix(dg,I,rlorder,m,supp,cliques,cql,cliquesize,method="chordal",chor_alg=TS)
        end
        opt,supp0,supp1,measure,moment=blockcpop_mix(n,m,dg,rlorder,supp,coe,cliques,cql,cliquesize,I,ncc,blocks,cl,blocksize,numeq=numeq,mix=true,QUIET=QUIET,solve=solve,solution=solution,extra_sos=extra_sos)
    end
    if solution==true
        sol=approx_sol(moment,n,cliques,cql,cliquesize)
    else
        sol=nothing
    end
    data=mdata_type(n,m,dg,supp,coe,numeq,rlorder,supp0,ssupp,lt,fbasis,gbasis,cql,cliques,cliquesize,I,ncc,blocks,cl,blocksize,ub,sizes)
    return opt,sol,data
end

function cs_tssos_higher!(data;TS="block",QUIET=false,solve=true,solution=false,extra_sos=true)
    n=data.n
    m=data.m
    dg=data.dg
    supp=data.supp
    coe=data.coe
    numeq=data.numeq
    rlorder=data.rlorder
    supp0=data.supp0
    ssupp=data.ssupp
    lt=data.lt
    fbasis=data.fbasis
    gbasis=data.gbasis
    cql=data.cql
    cliques=data.cliques
    cliquesize=data.cliquesize
    I=data.I
    ncc=data.ncc
    blocks=data.blocks
    cl=data.cl
    blocksize=data.blocksize
    ub=data.ub
    sizes=data.sizes
    if TS=="block"
        blocks,cl,blocksize,ub,sizes,status=get_chblocks_mix!(m,I,supp0,ssupp,lt,fbasis,gbasis,cql,cliques,cliquesize,blocks,cl,blocksize,ub,sizes,method="block",chor_alg=nothing)
    else
        blocks,cl,blocksize,ub,sizes,status=get_chblocks_mix!(m,I,supp0,ssupp,lt,fbasis,gbasis,cql,cliques,cliquesize,blocks,cl,blocksize,ub,sizes,method="chordal",chor_alg=TS)
    end
    if status==1
        opt,supp0,supp1,measure,moment=blockcpop_mix(n,m,dg,rlorder,supp,coe,cliques,cql,cliquesize,I,ncc,blocks,cl,blocksize,numeq=numeq,mix=true,QUIET=QUIET,solve=solve,solution=solution,extra_sos=extra_sos)
    else
        println("No higher CS-TSSOS hierarchy!")
    end
    if solution==true
        sol=approx_sol(moment,n,cliques,cql,cliquesize)
    else
        sol=nothing
    end
    data.supp0=supp0
    data.blocks=blocks
    data.cl=cl
    data.blocksize=blocksize
    data.ub=ub
    data.sizes=sizes
    return opt,sol,data
end

function blockupop_mix(n,d,supp,coe,cliques,cql,cliquesize,blocks,cl,blocksize;mix=true,QUIET=false,solve=true,solution=true,extra_sos=true)
    basis=Array{SparseMatrixCSC{UInt8,UInt32}}(undef,cql)
    col=Int[1]
    row=Int[]
    nz=UInt8[]
    if mix==false
        for i=1:cql
            basis[i]=sparse_basis(cliques[i],n,d)
            ssupp=sparse_basis(cliques[i],n,2*d)
            append!(col,ssupp.colptr[2:end].+(col[end]-1))
            append!(row,ssupp.rowval)
            append!(nz,ssupp.nzval)
        end
    else
        for i=1:cql
            basis[i]=sparse_basis(cliques[i],n,d)
            tcol=[1;basis[i].colptr]
            trow=basis[i].rowval
            tnz=basis[i].nzval
            for j=1:cl[i]
                for k=1:blocksize[i][j]
                    ind1=blocks[i][j][k]
                    t=ind1==1 ? 2 : k
                    for r=t:blocksize[i][j]
                        @inbounds ind2=blocks[i][j][r]
                        @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                        append!(col,col[end]+length(bi_row))
                        append!(row,bi_row)
                        append!(nz,bi_nz)
                    end
                end
            end
        end
        col0=copy(col)
        row0=copy(row)
        nz0=copy(nz)
        if extra_sos==true||solution==true
            for i=1:cql
                ssupp=sparse_basis(cliques[i],n,2)
                append!(col,ssupp.colptr[2:end].+(col[end]-1))
                append!(row,ssupp.rowval)
                append!(nz,ssupp.nzval)
            end
        end
    end
    supp1=SparseMatrixCSC(n,UInt32(length(col)),[1;col],row,nz)
    supp1=Array(supp1)
    supp1=unique(supp1,dims=2)
    supp1=sortslices(supp1,dims=2)
    supp1=sparse(supp1)
    if mix==true
        supp0=SparseMatrixCSC(n,UInt32(length(col0)),[1;col0],row0,nz0)
        supp0=Array(supp0)
        supp0=unique(supp0,dims=2)
        supp0=sparse(supp0)
    else
        supp0=supp1
    end
    moment=nothing
    objv=nothing
    if solve==true
        lsupp1=supp1.n
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:lsupp1]
        if mix==false
            pos0=Vector{Symmetric{VariableRef}}(undef, cql)
            for i=1:cql
                lb=basis[i].n+1
                pos0[i]=@variable(model, [1:lb, 1:lb], PSD)
                tcol=[1;basis[i].colptr]
                trow=basis[i].rowval
                tnz=basis[i].nzval
                for j=1:lb
                    for k=j:lb
                        bi_row,bi_nz=splus(trow[tcol[j]:(tcol[j+1]-1)],tnz[tcol[j]:(tcol[j+1]-1)],trow[tcol[k]:(tcol[k+1]-1)],tnz[tcol[k]:(tcol[k+1]-1)])
                        Locb=bfind_sparse(supp1,bi_row,bi_nz)
                        if j==k
                            @inbounds cons[Locb]+=pos0[i][j,k]
                        else
                            @inbounds cons[Locb]+=2*pos0[i][j,k]
                        end
                    end
                end
            end
        else
            pos1=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, cql)
            for i=1:cql
                if extra_sos==true||solution==true
                    pos0=Vector{Symmetric{VariableRef}}(undef, cql)
                    lb=cliquesize[i]+1
                    pos0[i]=@variable(model, [1:lb, 1:lb], PSD)
                    tcol=[1;basis[i].colptr[1:lb]]
                    trow=basis[i].rowval[1:lb-1]
                    tnz=basis[i].nzval[1:lb-1]
                    for j=1:lb
                        for k=j:lb
                            bi_row,bi_nz=splus(trow[tcol[j]:(tcol[j+1]-1)],tnz[tcol[j]:(tcol[j+1]-1)],trow[tcol[k]:(tcol[k+1]-1)],tnz[tcol[k]:(tcol[k+1]-1)])
                            Locb=bfind_sparse(supp1,bi_row,bi_nz)
                            if j==k
                                @inbounds cons[Locb]+=pos0[i][j,k]
                            else
                                @inbounds cons[Locb]+=2*pos0[i][j,k]
                            end
                        end
                    end
                end
                pos1[i]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i])
                tcol=[1;basis[i].colptr]
                trow=basis[i].rowval
                tnz=basis[i].nzval
                for k=1:cl[i]
                    if blocksize[i][k]==1
                       pos1[i][k]=@variable(model, lower_bound=0)
                       bi_row=trow[tcol[blocks[i][k][1]]:(tcol[blocks[i][k][1]+1]-1)]
                       bi_nz=2*tnz[tcol[blocks[i][k][1]]:(tcol[blocks[i][k][1]+1]-1)]
                       Locb=bfind_sparse(supp1,bi_row,bi_nz)
                       @inbounds cons[Locb]+=pos1[i][k]
                    else
                       pos1[i][k]=@variable(model, [1:blocksize[i][k], 1:blocksize[i][k]], PSD)
                       for j=1:blocksize[i][k]
                           ind1=blocks[i][k][j]
                           for r=j:blocksize[i][k]
                               @inbounds ind2=blocks[i][k][r]
                               @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                               Locb=bfind_sparse(supp1,bi_row,bi_nz)
                               if j==r
                                   @inbounds cons[Locb]+=pos1[i][k][j,r]
                               else
                                   @inbounds cons[Locb]+=2*pos1[i][k][j,r]
                               end
                           end
                       end
                    end
                end
            end
        end
        bc=zeros(lsupp1,1)
        for i=1:supp.n
            Locb=bfind_sparse(supp1,supp.rowval[supp.colptr[i]:(supp.colptr[i+1]-1)],supp.nzval[supp.colptr[i]:(supp.colptr[i+1]-1)])
            if Locb==0
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing
            else
               bc[Locb]=coe[i]
            end
        end
        @variable(model, lower)
        cons[1]+=lower
        @constraint(model, con[i=1:lsupp1], cons[i]==bc[i])
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
            moment=Vector{Union{Float64, Symmetric{Float64}, Array{Float64,2}}}(undef, cql)
            dual_var=-dual.(con)
            for i=1:cql
                lb=cliquesize[i]
                if lb==1
                    Locb=bfind_sparse(supp1,cliques[i],[1])
                    moment[i]=dual_var[Locb]
                else
                    tcol=UInt32[l for l=1:lb+1]
                    trow=cliques[i]
                    tnz=ones(UInt8, lb)
                    moment[i]=zeros(Float64,lb,lb)
                    for j=1:lb
                        for k=j:lb
                            @inbounds bi_row,bi_nz=splus(trow[tcol[j]:(tcol[j+1]-1)],tnz[tcol[j]:(tcol[j+1]-1)],trow[tcol[k]:(tcol[k+1]-1)],tnz[tcol[k]:(tcol[k+1]-1)])
                            Locb=bfind_sparse(supp1,bi_row,bi_nz)
                            moment[i][j,k]=dual_var[Locb]
                        end
                    end
                    moment[i]=Symmetric(moment[i],:U)
                end
            end
        end
    end
    return objv,supp0,moment
end

function blockcpop_mix(n,m,dg,rlorder,supp,coe,cliques,cql,cliquesize,I,ncc,blocks,cl,blocksize;numeq=0,mix=true,QUIET=false,solve=true,solution=false,cons_label=false,extra_sos=true,small=true)
    fbasis=Array{SparseMatrixCSC{UInt8,UInt32}}(undef,cql)
    gbasis=Array{SparseMatrixCSC{UInt8,UInt32}}(undef,m)
    col=Int[1]
    row=Int[]
    nz=UInt8[]
    if mix==false
        for i=1:cql
            fbasis[i]=sparse_basis(cliques[i],n,rlorder[i])
            ssupp=sparse_basis(cliques[i],n,2*rlorder[i])
            append!(col,ssupp.colptr[2:end].+(col[end]-1))
            append!(row,ssupp.rowval)
            append!(nz,ssupp.nzval)
            for j in I[i]
                # if j<=m-numeq
                    gbasis[j]=sparse_basis(cliques[i],n,rlorder[i]-ceil(Int, dg[j]/2))
                # else
                #     gbasis[j]=sparse_basis(cliques[i],n,2*rlorder[i]-dg[j])
                # end
                ssupp=sparse_basis(cliques[i],n,2*(rlorder[i]-ceil(Int, dg[j]/2)))
                tcol=[1;ssupp.colptr]
                trow=ssupp.rowval
                tnz=ssupp.nzval
                for t=1:ssupp.n+1
                    for s=1:supp[j+1].n
                        @inbounds bi_row,bi_nz=splus(trow[tcol[t]:(tcol[t+1]-1)],tnz[tcol[t]:(tcol[t+1]-1)],supp[j+1].rowval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)],supp[j+1].nzval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)])
                        append!(col,col[end]+length(bi_row))
                        append!(row,bi_row)
                        append!(nz,bi_nz)
                    end
                end
            end
        end
    else
        for i=1:cql
            fbasis[i]=sparse_basis(cliques[i],n,rlorder[i])
            tcol=[1;fbasis[i].colptr]
            trow=fbasis[i].rowval
            tnz=fbasis[i].nzval
            for j=1:cl[i][1]
                for k=1:blocksize[i][1][j]
                    ind1=blocks[i][1][j][k]
                    t=ind1==1 ? 2 : k
                    for r=t:blocksize[i][1][j]
                        @inbounds ind2=blocks[i][1][j][r]
                        @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                        append!(col,col[end]+length(bi_row))
                        append!(row,bi_row)
                        append!(nz,bi_nz)
                    end
                end
            end
        end
        if small==true||extra_sos==true
            col0=copy(col)
            row0=copy(row)
            nz0=copy(nz)
        end
        for i=1:cql
            if extra_sos==true||solution==true
                ssupp=sparse_basis(cliques[i],n,2)
                append!(col,ssupp.colptr[2:end].+(col[end]-1))
                append!(row,ssupp.rowval)
                append!(nz,ssupp.nzval)
            end
            for k=1:length(I[i])
                j=I[i][k]
                # if j<=m-numeq
                    gbasis[j]=sparse_basis(cliques[i],n,rlorder[i]-ceil(Int, dg[j]/2))
                    tcol=[1;gbasis[j].colptr]
                    trow=gbasis[j].rowval
                    tnz=gbasis[j].nzval
                    for l=1:cl[i][k+1]
                        for t=1:blocksize[i][k+1][l]
                            ind1=blocks[i][k+1][l][t]
                            for r=t:blocksize[i][k+1][l]
                                ind2=blocks[i][k+1][l][r]
                                for s=1:supp[j+1].n
                                    @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                                    @inbounds bi_row,bi_nz=splus(bi_row,bi_nz,supp[j+1].rowval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)],supp[j+1].nzval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)])
                                    if length(bi_row)!=0
                                       append!(col,col[end]+length(bi_row))
                                       append!(row,bi_row)
                                       append!(nz,bi_nz)
                                    end
                                end
                            end
                        end
                    end
                # else
                #     gbasis[j]=sparse_basis(cliques[i],n,2*rlorder[i]-dg[j])
                #     tcol=[1;gbasis[j].colptr]
                #     trow=gbasis[j].rowval
                #     tnz=gbasis[j].nzval
                #     for l in blocks[i][k+1][1]
                #         for s=1:supp[j+1].n
                #             @inbounds bi_row,bi_nz=splus(trow[tcol[l]:(tcol[l+1]-1)],tnz[tcol[l]:(tcol[l+1]-1)],supp[j+1].rowval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)],supp[j+1].nzval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)])
                #             if length(bi_row)!=0
                #                append!(col,col[end]+length(bi_row))
                #                append!(row,bi_row)
                #                append!(nz,bi_nz)
                #             end
                #         end
                #     end
                # end
            end
        end
    end
    for i ∈ ncc
        append!(col,supp[i+1].colptr[2:end].+(col[end]-1))
        append!(row,supp[i+1].rowval)
        append!(nz,supp[i+1].nzval)
    end
    supp1=SparseMatrixCSC(n,UInt32(length(col)),[1;col],row,nz)
    supp1=Array(supp1)
    supp1=unique(supp1,dims=2)
    supp1=sortslices(supp1,dims=2)
    supp1=sparse(supp1)
    if (small==true||extra_sos==true)&&mix==true
        supp0=SparseMatrixCSC(n,UInt32(length(col0)),[1;col0],row0,nz0)
        supp0=Array(supp0)
        supp0=unique(supp0,dims=2)
        supp0=sparse(supp0)
    else
        supp0=supp1
    end
    objv=nothing
    dual_var=nothing
    moment=nothing
    if solve==true
        lsupp1=supp1.n
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:lsupp1]
        if mix==false
            pos0=Vector{Symmetric{VariableRef}}(undef, cql)
            for i=1:cql
                lb=fbasis[i].n+1
                pos0[i]=@variable(model, [1:lb, 1:lb], PSD) # sigma_0
                tcol=[1;fbasis[i].colptr]
                trow=fbasis[i].rowval
                tnz=fbasis[i].nzval
                for j=1:lb
                    for k=j:lb
                        @inbounds bi_row,bi_nz=splus(trow[tcol[j]:(tcol[j+1]-1)],tnz[tcol[j]:(tcol[j+1]-1)],trow[tcol[k]:(tcol[k+1]-1)],tnz[tcol[k]:(tcol[k+1]-1)])
                        Locb=bfind_sparse(supp1,bi_row,bi_nz)
                        if j==k
                            @inbounds cons[Locb]+=pos0[i][j,k]
                        else
                            @inbounds cons[Locb]+=2*pos0[i][j,k]
                        end
                    end
                end
            end
        else
            pos1=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, cql)
            for i=1:cql
                if extra_sos==true||solution==true
                    pos0=Vector{Symmetric{VariableRef}}(undef, cql)
                    lb=cliquesize[i]+1
                    pos0[i]=@variable(model, [1:lb, 1:lb], PSD)
                    tcol=[1;fbasis[i].colptr[1:lb]]
                    trow=fbasis[i].rowval[1:lb-1]
                    tnz=fbasis[i].nzval[1:lb-1]
                    for t=1:lb
                        for r=t:lb
                            bi_row,bi_nz=splus(trow[tcol[t]:(tcol[t+1]-1)],tnz[tcol[t]:(tcol[t+1]-1)],trow[tcol[r]:(tcol[r+1]-1)],tnz[tcol[r]:(tcol[r+1]-1)])
                            Locb=bfind_sparse(supp1,bi_row,bi_nz)
                            if t==r
                                @inbounds cons[Locb]+=pos0[i][t,r]
                            else
                                @inbounds cons[Locb]+=2*pos0[i][t,r]
                            end
                        end
                    end
                end
                pos1[i]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i][1])
                tcol=[1;fbasis[i].colptr]
                trow=fbasis[i].rowval
                tnz=fbasis[i].nzval
                for l=1:cl[i][1]
                    if blocksize[i][1][l]==1
                       @inbounds pos1[i][l]=@variable(model, lower_bound=0)
                       @inbounds bi_row=trow[tcol[blocks[i][1][l][1]]:(tcol[blocks[i][1][l][1]+1]-1)]
                       @inbounds bi_nz=2*tnz[tcol[blocks[i][1][l][1]]:(tcol[blocks[i][1][l][1]+1]-1)]
                       Locb=bfind_sparse(supp1,bi_row,bi_nz)
                       @inbounds cons[Locb]+=pos1[i][l]
                    else
                       @inbounds bs=blocksize[i][1][l]
                       @inbounds pos1[i][l]=@variable(model, [1:bs, 1:bs], PSD)
                       for t=1:bs
                           @inbounds ind1=blocks[i][1][l][t]
                           for r=t:bs
                               @inbounds ind2=blocks[i][1][l][r]
                               @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                               Locb=bfind_sparse(supp1,bi_row,bi_nz)
                               if t==r
                                  @inbounds cons[Locb]+=pos1[i][l][t,r]
                               else
                                  @inbounds cons[Locb]+=2*pos1[i][l][t,r]
                               end
                           end
                       end
                    end
                end
            end
        end
        pos2=Vector{VariableRef}(undef, length(ncc))
        for k=1:length(ncc)
            i=ncc[k]
            if i<=m-numeq
                pos2[k]=@variable(model, lower_bound=0)
            else
                pos2[k]=@variable(model)
            end
            for j=1:supp[i+1].n
                Locb=bfind_sparse(supp1,supp[i+1].rowval[supp[i+1].colptr[j]:(supp[i+1].colptr[j+1]-1)],supp[i+1].nzval[supp[i+1].colptr[j]:(supp[i+1].colptr[j+1]-1)])
                cons[Locb]+=coe[i+1][j]*pos2[k]
            end
        end
        p=1
        if mix==false
            pos3=Vector{Union{VariableRef,Symmetric{VariableRef},Vector{VariableRef}}}(undef, m-length(ncc))
            for i=1:cql
                for j in I[i]
                    tcol=[1;gbasis[j].colptr]
                    trow=gbasis[j].rowval
                    tnz=gbasis[j].nzval
                    lb=gbasis[j].n+1
                    # if j<=m-numeq
                        if lb==1
                            if j<=m-numeq
                                pos3[p]=@variable(model, lower_bound=0)
                            else
                                pos3[p]=@variable(model)
                            end
                            for s=1:supp[j+1].n
                                Locb=bfind_sparse(supp1,supp[j+1].rowval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)],supp[j+1].nzval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)])
                                @inbounds cons[Locb]+=coe[j+1][s]*pos3[p]
                            end
                        else
                            if j<=m-numeq
                                pos3[p]=@variable(model, [1:lb, 1:lb], PSD)
                            else
                                pos3[p]=@variable(model, [1:lb, 1:lb], Symmetric)
                            end
                            for t=1:lb
                                for r=t:lb
                                    for s=1:supp[j+1].n
                                        @inbounds bi_row,bi_nz=splus(trow[tcol[t]:(tcol[t+1]-1)],tnz[tcol[t]:(tcol[t+1]-1)],trow[tcol[r]:(tcol[r+1]-1)],tnz[tcol[r]:(tcol[r+1]-1)])
                                        @inbounds bi_row,bi_nz=splus(bi_row,bi_nz,supp[j+1].rowval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)],supp[j+1].nzval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)])
                                        Locb=bfind_sparse(supp1,bi_row,bi_nz)
                                        if t==r
                                            @inbounds cons[Locb]+=coe[j+1][s]*pos3[p][t,r]
                                        else
                                            @inbounds cons[Locb],2*coe[j+1][s]*pos3[p][t,r]
                                        end
                                    end
                                end
                            end
                        end
                    # else
                    #     pos3[p]=@variable(model, [1:lb])
                    #     for t=1:lb
                    #         for s=1:supp[j+1].n
                    #             @inbounds bi_row,bi_nz=splus(trow[tcol[t]:(tcol[t+1]-1)],tnz[tcol[t]:(tcol[t+1]-1)],supp[j+1].rowval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)],supp[j+1].nzval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)])
                    #             Locb=bfind_sparse(supp1,bi_row,bi_nz)
                    #             @inbounds cons[Locb]+=coe[j+1][s]*pos3[p][t]
                    #         end
                    #     end
                    # end
                    p+=1
                end
            end
        else
            pos3=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m-length(ncc))
            for i=1:cql
                for k=1:length(I[i])
                    j=I[i][k]
                    tcol=[1;gbasis[j].colptr]
                    trow=gbasis[j].rowval
                    tnz=gbasis[j].nzval
                    # if j<=m-numeq
                        pos3[p]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef,cl[i][k+1])
                        for l=1:cl[i][k+1]
                            bs=blocksize[i][k+1][l]
                            if bs==1
                                if j<=m-numeq
                                    pos3[p][l]=@variable(model, lower_bound=0)
                                else
                                    pos3[p][l]=@variable(model)
                                end
                               for s=1:supp[j+1].n
                                   @inbounds ind1=blocks[i][k+1][l][1]
                                   @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],2*tnz[tcol[ind1]:(tcol[ind1+1]-1)],supp[j+1].rowval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)],supp[j+1].nzval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)])
                                   Locb=bfind_sparse(supp1,bi_row,bi_nz)
                                   @inbounds cons[Locb]+=coe[j+1][s]*pos3[p][l]
                               end
                            else
                                if j<=m-numeq
                                    pos3[p][l]=@variable(model, [1:bs, 1:bs], PSD)
                                else
                                    pos3[p][l]=@variable(model, [1:bs, 1:bs], Symmetric)
                                end
                                for t=1:bs
                                   ind1=blocks[i][k+1][l][t]
                                   for r=t:bs
                                       ind2=blocks[i][k+1][l][r]
                                       for s=1:supp[j+1].n
                                           @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                                           @inbounds bi_row,bi_nz=splus(bi_row,bi_nz,supp[j+1].rowval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)],supp[j+1].nzval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)])
                                           Locb=bfind_sparse(supp1,bi_row,bi_nz)
                                           if t==r
                                               @inbounds cons[Locb]+=coe[j+1][s]*pos3[p][l][t,r]
                                           else
                                               @inbounds cons[Locb]+=2*coe[j+1][s]*pos3[p][l][t,r]
                                           end
                                       end
                                   end
                               end
                            end
                        end
                    # else
                    #     pos3[p]=@variable(model, [1:blocksize[i][k+1][1]])
                    #     for r=1:blocksize[i][k+1][1]
                    #         t=blocks[i][k+1][1][r]
                    #         for s=1:supp[j+1].n
                    #             @inbounds bi_row,bi_nz=splus(trow[tcol[t]:(tcol[t+1]-1)],tnz[tcol[t]:(tcol[t+1]-1)],supp[j+1].rowval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)],supp[j+1].nzval[supp[j+1].colptr[s]:(supp[j+1].colptr[s+1]-1)])
                    #             Locb=bfind_sparse(supp1,bi_row,bi_nz)
                    #             @inbounds cons[Locb]+=coe[j+1][s]*pos3[p][r]
                    #         end
                    #     end
                    # end
                    p+=1
                end
            end
        end
        bc=zeros(lsupp1,1)
        for i=1:supp[1].n
            Locb=bfind_sparse(supp1,supp[1].rowval[supp[1].colptr[i]:(supp[1].colptr[i+1]-1)],supp[1].nzval[supp[1].colptr[i]:(supp[1].colptr[i+1]-1)])
            if Locb==0
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing,nothing,nothing
            else
               bc[Locb]=coe[1][i]
            end
        end
        @variable(model, lower)
        if cons_label==true||solution==true
            cons[1]+=lower
            @constraint(model, con[i=1:lsupp1], cons[i]==bc[i])
        else
            @constraint(model, cons[2:end].==bc[2:end])
            @constraint(model, cons[1]+lower==bc[1])
        end
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
            moment=Vector{Union{Float64, Symmetric{Float64}, Array{Float64,2}}}(undef, cql)
            dual_var=-dual.(con)
            for i=1:cql
                lb=cliquesize[i]
                if lb==1
                    Locb=bfind_sparse(supp1,cliques[i],[1])
                    moment[i]=dual_var[Locb]
                else
                    tcol=UInt32[l for l=1:lb+1]
                    trow=cliques[i]
                    tnz=ones(UInt8, lb)
                    moment[i]=zeros(Float64,lb,lb)
                    for j=1:lb
                        for k=j:lb
                            @inbounds bi_row,bi_nz=splus(trow[tcol[j]:(tcol[j+1]-1)],tnz[tcol[j]:(tcol[j+1]-1)],trow[tcol[k]:(tcol[k+1]-1)],tnz[tcol[k]:(tcol[k+1]-1)])
                            Locb=bfind_sparse(supp1,bi_row,bi_nz)
                            moment[i][j,k]=dual_var[Locb]
                        end
                    end
                    moment[i]=Symmetric(moment[i],:U)
                end
            end
        end
    end
    return objv,supp0,supp1,dual_var,moment
end

function get_blocks_mix(d,supp,cliques,cql,cliquesize;method="block",chor_alg="amd",merge=false)
    blocks=Vector{Vector{Vector{UInt16}}}(undef,cql)
    cl=Vector{UInt16}(undef,cql)
    ub=Vector{Vector{UInt16}}(undef,cql)
    sizes=Vector{Vector{UInt16}}(undef,cql)
    blocksize=Vector{Vector{Int}}(undef,cql)
    basis=Vector{Array{UInt8,2}}(undef,cql)
    for i=1:cql
        nvar=cliquesize[i]
        ssupp=zeros(UInt8,nvar,1)
        for j=1:supp.n
            if issubset(supp.rowval[supp.colptr[j]:(supp.colptr[j+1]-1)], cliques[i])
                bi=zeros(UInt8,nvar,1)
                for k=supp.colptr[j]:(supp.colptr[j+1]-1)
                    @inbounds locb=lbfind(cliques[i],cliquesize[i],supp.rowval[k])
                    @inbounds bi[locb]=supp.nzval[k]
                end
                ssupp=[ssupp bi]
            end
        end
        basis[i]=get_basis(nvar,d)
        if method=="block"
            blocks[i],cl[i],blocksize[i],ub[i],sizes[i]=get_blocks(nvar,ssupp,basis[i],QUIET=true)
        else
            blocks[i],cl[i],blocksize[i],ub[i],sizes[i]=get_cliques(nvar,ssupp,basis[i],QUIET=true,alg=chor_alg,merge=merge)
        end
    end
    return blocks,cl,blocksize,ub,sizes,basis
end

function get_hblocks_mix!(supp,basis,cliques,cql,cliquesize,blocks,cl,blocksize,ub,sizes;method="block",chor_alg="amd",merge=false)
    nub=Vector{Vector{UInt16}}(undef,cql)
    nsizes=Vector{Vector{UInt16}}(undef,cql)
    status=ones(UInt8,cql)
    for i=1:cql
        nvar=cliquesize[i]
        ssupp=zeros(UInt8,nvar,1)
        for j=1:supp.n
            if issubset(supp.rowval[supp.colptr[j]:(supp.colptr[j+1]-1)], cliques[i])
                bi=zeros(UInt8,nvar,1)
                for k=supp.colptr[j]:(supp.colptr[j+1]-1)
                    @inbounds locb=lbfind(cliques[i],cliquesize[i],supp.rowval[k])
                    @inbounds bi[locb]=supp.nzval[k]
                end
                ssupp=[ssupp bi]
            end
        end
        # ssupp=zeros(UInt8,nvar,Int(sum(blocksize[i].^2+blocksize[i])/2))
        # k=1
        # for s=1:cl[i]
        #     for j=1:blocksize[i][s]
        #         for r=j:blocksize[i][s]
        #             @inbounds bi=basis[i][:,blocks[i][s][j]]+basis[i][:,blocks[i][s][r]]
        #             @inbounds ssupp[:,k]=bi
        #             k+=1
        #         end
        #     end
        # end
        # ssupp=sortslices(ssupp,dims=2)
        # ssupp=unique(ssupp,dims=2)
        if method=="block"
            blocks[i],cl[i],blocksize[i],nub[i],nsizes[i],status[i]=get_hblocks(nvar,ssupp,basis[i],ub[i],sizes[i],QUIET=true)
        else
            blocks[i],cl[i],blocksize[i],nub[i],nsizes[i],status[i]=get_hcliques(nvar,ssupp,basis[i],ub[i],sizes[i],QUIET=true,alg=chor_alg,merge=merge)
        end
    end
    return blocks,cl,blocksize,nub,nsizes,maximum(status)
end

function assign_constraint(m,supp,cliques,cql,cliquesize;assign="first")
    I=[UInt16[] for i=1:cql]
    ncc=UInt16[]
    for i=2:m+1
        rind=unique(supp[i].rowval)
        if assign=="first"
            ind=findfirst(k->issubset(rind, cliques[k]), 1:cql)
            if ind!=nothing
                push!(I[ind], i-1)
            else
                push!(ncc, i-1)
            end
        else
            temp=UInt16[]
            for j=1:cql
                if issubset(rind, cliques[j])
                    push!(temp,j)
                end
            end
            if temp!=[]
                if assign=="min"
                    push!(I[temp[argmin(cliquesize[temp])]], i-1)
                else
                    push!(I[temp[argmax(cliquesize[temp])]], i-1)
                end
            else
                push!(ncc, i-1)
            end
        end
    end
    return I,ncc
end

function init_order(dg,I,cql;order="multi")
    rlorder=ones(Int,cql)
    if order=="multi"
        for i=1:cql
            if I[i]==[]
                rlorder[i]=1
            else
                rlorder[i]=ceil(Int, maximum(dg[I[i]])/2)
            end
        end
    else
        rlorder*=order
    end
    return rlorder
end

function get_cblocks_mix(dg,I,rlorder,m,supp,cliques,cql,cliquesize;method="block",chor_alg="amd",merge=false)
    blocks=Vector{Vector{Vector{Vector{UInt16}}}}(undef,cql)
    cl=Vector{Vector{UInt16}}(undef,cql)
    ub=Vector{Vector{UInt16}}(undef,cql)
    sizes=Vector{Vector{UInt16}}(undef,cql)
    blocksize=Vector{Vector{Vector{Int}}}(undef,cql)
    ssupp=Vector{Vector{Array{UInt8,2}}}(undef,cql)
    lt=Vector{Vector{UInt16}}(undef,cql)
    gbasis=Vector{Vector{Array{UInt8,2}}}(undef,cql)
    fbasis=Vector{Array{UInt8,2}}(undef,cql)
    lsupp=sum([supp[i].n for i=1:m+1])
    col=Int[1]
    row=Int[]
    nz=UInt8[]
    for i=1:m+1
        append!(col,supp[i].colptr[2:end].+(col[end]-1))
        append!(row,supp[i].rowval)
        append!(nz,supp[i].nzval)
    end
    for i=1:cql
        lc=length(I[i])
        blocks[i]=Vector{Vector{Vector{UInt16}}}(undef, lc+1)
        cl[i]=Vector{UInt16}(undef, lc+1)
        ub[i]=Vector{UInt16}(undef, lc+1)
        sizes[i]=Vector{UInt16}(undef, lc+1)
        blocksize[i]=Vector{Vector{Int}}(undef, lc+1)
        nvar=cliquesize[i]
        fsupp=zeros(UInt8,nvar,1)
        for j=1:lsupp
            if issubset(row[col[j]:(col[j+1]-1)], cliques[i])
                bi=zeros(UInt8,nvar,1)
                for k=col[j]:(col[j+1]-1)
                    @inbounds locb=lbfind(cliques[i],cliquesize[i],row[k])
                    @inbounds bi[locb]=nz[k]
                end
                fsupp=[fsupp bi]
            end
        end
        fsupp=unique(fsupp,dims=2)
        fbasis[i]=get_basis(cliquesize[i],rlorder[i])
        gbasis[i]=Vector{Array{UInt8,2}}(undef, lc)
        ssupp[i]=Vector{Array{UInt8,2}}(undef, lc+1)
        lt[i]=Vector{UInt16}(undef, lc+1)
        lt[i][1]=supp[1].n
        ssupp[i][1]=zeros(UInt8,nvar,supp[1].n)
        for s=1:lc
            t=I[i][s]
            # if t<=m-numeq
                gbasis[i][s]=get_basis(nvar,rlorder[i]-ceil(Int, dg[t]/2))
            # else
            #     gbasis[i][s]=get_basis(nvar,2*rlorder[i]-dg[t])
            # end
            ssupp[i][s+1]=zeros(UInt8,nvar,supp[t+1].n)
            lt[i][s+1]=supp[t+1].n
            for j=1:supp[t+1].n
                for k=supp[t+1].colptr[j]:(supp[t+1].colptr[j+1]-1)
                    @inbounds locb=lbfind(cliques[i],cliquesize[i],supp[t+1].rowval[k])
                    @inbounds ssupp[i][s+1][locb,j]=supp[t+1].nzval[k]
                end
            end
        end
        if method=="block"
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],ub[i],sizes[i]=get_cblocks(nvar,lc,fsupp,ssupp[i],lt[i],fbasis[i],gbasis[i],QUIET=true)
        else
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],ub[i],sizes[i]=get_ccliques(nvar,lc,fsupp,ssupp[i],lt[i],fbasis[i],gbasis[i],QUIET=true,alg=chor_alg,merge=merge)
        end
    end
    return blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis
end

function get_chblocks_mix!(m,I,supp,ssupp,lt,fbasis,gbasis,cql,cliques,cliquesize,blocks,cl,blocksize,ub,sizes;method="block",chor_alg="amd",merge=false)
    nub=Vector{Vector{UInt16}}(undef,cql)
    nsizes=Vector{Vector{UInt16}}(undef,cql)
    status=ones(UInt8,cql)
    for i=1:cql
        lc=length(I[i])
        nvar=cliquesize[i]
        fsupp=zeros(UInt8,nvar,1)
        for j=1:supp.n
            if issubset(supp.rowval[supp.colptr[j]:(supp.colptr[j+1]-1)], cliques[i])
                bi=zeros(UInt8,nvar,1)
                for k=supp.colptr[j]:(supp.colptr[j+1]-1)
                    @inbounds locb=lbfind(cliques[i],cliquesize[i],supp.rowval[k])
                    @inbounds bi[locb]=supp.nzval[k]
                end
                fsupp=[fsupp bi]
            end
        end
        # fsupp=zeros(UInt8,nvar,Int(sum(blocksize[i][1].^2+blocksize[i][1])/2))
        # k=1
        # for s=1:cl[i][1]
        #     for j=1:blocksize[i][1][s]
        #         for r=j:blocksize[i][1][s]
        #             @inbounds bi=fbasis[i][:,blocks[i][1][s][j]]+fbasis[i][:,blocks[i][1][s][r]]
        #             @inbounds fsupp[:,k]=bi
        #             k+=1
        #         end
        #     end
        # end
        # fsupp=sortslices(fsupp,dims=2)
        # fsupp=unique(fsupp,dims=2)
        if method=="block"
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],nub[i],nsizes[i],status[i]=get_chblocks!(nvar,lc,ssupp[i],lt[i],fbasis[i],gbasis[i],fsupp,blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],ub[i],sizes[i],QUIET=true)
        else
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],nub[i],nsizes[i],status[i]=get_chcliques!(nvar,lc,ssupp[i],lt[i],fbasis[i],gbasis[i],fsupp,blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],ub[i],sizes[i],QUIET=true,alg=chor_alg,merge=merge)
        end
    end
    return blocks,cl,blocksize,nub,nsizes,maximum(status)
end

function clique_decomp(n,supp;alg="amd")
    if alg=="greedy"
        G=CGraph()
        for i=1:n
            cadd_node!(G)
        end
    else
        A=zeros(UInt8,n,n)
    end
    for i = 1:supp.n
        if alg=="greedy"
            cadd_clique!(G,supp.rowval[supp.colptr[i]:(supp.colptr[i+1]-1)])
        else
            lcol=supp.colptr[i+1]-supp.colptr[i]
            A[supp.rowval[supp.colptr[i]:(supp.colptr[i+1]-1)],supp.rowval[supp.colptr[i]:(supp.colptr[i+1]-1)]]=ones(UInt8,lcol,lcol)
        end
    end
    if alg=="greedy"
        cliques,cql,cliquesize=chordal_extension(G, GreedyFillIn())
    else
        cliques,cql,cliquesize=cliquesFromSpMatD(A)
    end
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("------------------------------------------------------")
    return cliques,cql,cliquesize
end

function clique_cdecomp(n,m,dg,supp;order="multi",alg="amd")
    if alg=="greedy"
        G=CGraph()
        for i=1:n
            cadd_node!(G)
        end
    else
        A=zeros(UInt8,n,n)
    end
    for i=1:m+1
        if order=="multi"||i==1||order==ceil(Int, dg[i-1]/2)
            for j = 1:supp[i].n
                if alg=="greedy"
                    cadd_clique!(G,supp[i].rowval[supp[i].colptr[j]:(supp[i].colptr[j+1]-1)])
                else
                    lcol=supp[i].colptr[j+1]-supp[i].colptr[j]
                    A[supp[i].rowval[supp[i].colptr[j]:(supp[i].colptr[j+1]-1)],supp[i].rowval[supp[i].colptr[j]:(supp[i].colptr[j+1]-1)]]=ones(UInt8,lcol,lcol)
                end
            end
        else
            if alg=="greedy"
                cadd_clique!(G,unique(supp[i].rowval))
            else
                rind=unique(supp[i].rowval)
                lcol=length(rind)
                A[rind,rind]=ones(UInt8,lcol,lcol)
            end
        end
    end
    if alg=="greedy"
        cliques,cql,cliquesize=chordal_extension(G, GreedyFillIn())
    elseif alg=="amd"
        cliques,cql,cliquesize=cliquesFromSpMatD(A)
    else
        cliques,cql,cliquesize=max_cliques(A)
    end
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("------------------------------------------------------")
    return cliques,cql,cliquesize
end

function splus(arow,anz,brow,bnz)
    la=length(arow)
    lb=length(brow)
    row=UInt32[]
    nz=UInt8[]
    i=1
    j=1
    while i<=la&&j<=lb
        if arow[i]==brow[j]
            push!(row,arow[i])
            push!(nz,anz[i]+bnz[j])
            i+=1
            j+=1
        elseif arow[i]<brow[j]
            push!(row,arow[i])
            push!(nz,anz[i])
            i+=1
        else
            push!(row,brow[j])
            push!(nz,bnz[j])
            j+=1
        end
    end
    if i<=la&&j>lb
        append!(row,arow[i:end])
        append!(nz,anz[i:end])
    elseif i>la&&j<=lb
        append!(row,brow[j:end])
        append!(nz,bnz[j:end])
    end
    return row,nz
end

function sparse_basis(var,tvar,d)
    if d==0
        return SparseMatrixCSC(tvar,0,UInt32[1],UInt32[],UInt8[])
    else
        n=length(var)
        lb=binomial(n+d,d)-1
        col=UInt32[1;2]
        row=UInt32[var[1]]
        nz=UInt8[1]
        i=1
        while i<d+1
            if row[end]==var[n]&&nz[end]==i
               if i<d
                  push!(row,var[1])
                  push!(nz,i+1)
                  push!(col,col[end]+1)
               end
               i=i+1
           elseif nz[col[end-1]]==1
               j=1
               while row[col[end-1]]!=var[j]
                   j+=1
               end
               if col[end]-col[end-1]>1&&row[col[end-1]+1]==var[j+1]
                   append!(row,row[col[end-1]+1:col[end]-1])
                   append!(nz,nz[col[end-1]+1:col[end]-1])
                   nz[col[end]]+=1
                   push!(col,2*col[end]-col[end-1]-1)
               else
                   append!(row,row[col[end-1]:col[end]-1])
                   append!(nz,nz[col[end-1]:col[end]-1])
                   row[col[end]]=var[j+1]
                   push!(col,2*col[end]-col[end-1])
               end
           else
               if row[col[end-1]]==var[1]
                   if col[end]-col[end-1]>1&&row[col[end-1]+1]==var[2]
                       append!(row,row[col[end-1]:col[end]-1])
                       append!(nz,nz[col[end-1]:col[end]-1])
                       nz[col[end]]-=1
                       nz[col[end]+1]+=1
                       push!(col,2*col[end]-col[end-1])
                  else
                      append!(row,row[col[end-1]:col[end]-1])
                      insert!(row,col[end]+1,var[2])
                      append!(nz,nz[col[end-1]:col[end]-1])
                      nz[col[end]]-=1
                      insert!(nz,col[end]+1,1)
                      push!(col,2*col[end]-col[end-1]+1)
                  end
               else
                   j=1
                   while row[col[end-1]]!=var[j]
                       j+=1
                   end
                   if col[end]-col[end-1]>1&&row[col[end-1]+1]==var[j+1]
                       append!(row,row[col[end-1]:col[end]-1])
                       row[col[end]]=var[1]
                       append!(nz,nz[col[end-1]:col[end]-1])
                       nz[col[end]]-=1
                       nz[col[end]+1]+=1
                       push!(col,2*col[end]-col[end-1])
                   else
                       append!(row,row[col[end-1]:col[end]-1])
                       row[col[end]]=var[1]
                       insert!(row,col[end]+1,var[j+1])
                       append!(nz,nz[col[end-1]:col[end]-1])
                       nz[col[end]]-=1
                       insert!(nz,col[end]+1,1)
                       push!(col,2*col[end]-col[end-1]+1)
                   end
               end
            end
        end
        return SparseMatrixCSC(tvar,lb,col,row,nz)
    end
end

# function sort_sparse(s)
#     m=s.m
#     n=s.n
#     col=s.colptr
#     row=s.rowval
#     nz=s.nzval
#     i=2
#     while i<=n
#         corr_row=row[col[i]:(col[i+1]-1)]
#         corr_nz=nz[col[i]:(col[i+1]-1)]
#         j=i-1
#         while j>=1
#             pre_row=row[col[j]:(col[j+1]-1)]
#             pre_nz=nz[col[j]:(col[j+1]-1)]
#             comp=comp_sparse(corr_row,corr_nz,pre_row,pre_nz)
#             if comp==-1
#                 j-=1
#             else
#                 break
#             end
#         end
#         if j<i-1
#             for k=i:-1:j+2
#                 lprow=col[k]-col[k-1]
#                 row[(col[k+1]-lprow):(col[k+1]-1)]=row[(col[k-1]):(col[k]-1)]
#                 nz[(col[k+1]-lprow):(col[k+1]-1)]=nz[(col[k-1]):(col[k]-1)]
#                 col[k]=col[k+1]-lprow
#             end
#             row[(col[j+1]):(col[j+2]-1)]=corr_row
#             nz[(col[j+1]):(col[j+2]-1)]=corr_nz
#         end
#         i+=1
#     end
#     return SparseMatrixCSC(m,n,col,row,nz)
# end

function comp_sparse(corr_row,corr_nz,pre_row,pre_nz)
    i=1
    lc=length(corr_row)
    lp=length(pre_row)
    while i<=lc&&i<=lp
        if corr_row[i]>pre_row[i]
            return -1
        elseif corr_row[i]<pre_row[i]
            return 1
        elseif corr_nz[i]<pre_nz[i]
            return -1
        elseif corr_nz[i]>pre_nz[i]
            return 1
        else
            i+=1
        end
    end
    if lc>lp
        return 1
    elseif lc<lp
        return -1
    else
        return 0
    end
end

function bfind_sparse(s,row,nz)
    l=s.n
    if l==0
        return 0
    end
    if length(row)==0
        return 1
    end
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        mid_row=s.rowval[s.colptr[mid]:(s.colptr[mid+1]-1)]
        mid_nz=s.nzval[s.colptr[mid]:(s.colptr[mid+1]-1)]
        order=comp_sparse(mid_row,mid_nz,row,nz)
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

function approx_sol(moment,n,cliques,cql,cliquesize)
    qsol=Float64[]
    lcq=sum(cliquesize)
    A=zeros(lcq,n)
    q=1
    for k=1:cql
        cqs=cliquesize[k]
        if cqs==1
            append!(qsol, moment[k])
        else
            conn=var_block(moment[k],cqs)
            for i=1:length(conn)
                lconn=length(conn[i])
                if lconn==1&&abs(moment[k][conn[i][1],conn[i][1]])<=1e-5
                    push!(qsol, 0)
                else
                    if lconn==1
                        temp=[sqrt(abs(moment[k][conn[i][1],conn[i][1]]))]
                    else
                        F=eigen(Symmetric(moment[k][conn[i],conn[i]]), lconn:lconn)
                        temp=sqrt(F.values[1])*F.vectors[:,1]
                    end
                    for l=1:k-1
                        inter=intersect(cliques[l], cliques[k][conn[i]])
                        if inter!=[]
                            flag=0
                            for r in inter
                                if l==1
                                    ind1=lbfind(cliques[l],cliquesize[l],r)
                                else
                                    ind1=sum(cliquesize[s] for s=1:l-1)+lbfind(cliques[l],cliquesize[l],r)
                                end
                                ind2=lbfind(cliques[k][conn[i]],lconn,r)
                                if (qsol[ind1]>=1e-3&&temp[ind2]<=-1e-3)||(qsol[ind1]<=-1e-3&&temp[ind2]>=1e-3)
                                    temp=-temp
                                    flag=1
                                    break
                                elseif (qsol[ind1]>=1e-3&&temp[ind2]>=1e-3)||(qsol[ind1]<=-1e-3&&temp[ind2]<=-1e-3)
                                    flag=1
                                    break
                                end
                            end
                            if flag==1
                                break
                            end
                        end
                    end
                    append!(qsol, temp)
                end
            end
        end
        for j=1:cqs
            A[q,cliques[k][j]]=1
            q+=1
        end
    end
    return (A'*A)\(A'*qsol)
end

function var_block(A,lb)
    adj=[abs(A[i,j])<=1e-8 ? 0 : 1 for i=1:lb, j=1:lb]
    return connected_components(Graph(adj))
end

function seval(supp,coe,x)
    val=0
    col=supp.colptr
    row=supp.rowval
    nz=supp.nzval
    for i=1:supp.n
        temp=mapreduce(j->x[row[j]]^nz[j],*,col[i]:(col[i+1]-1),init=1)
        val+=coe[i]*temp
    end
    return val
end
