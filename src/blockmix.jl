-------------------------------------------------------------------
n=6
d=2
@polyvar x[1:n]
f=1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]

mon=monomials(f)
coe=coefficients(f)
lm=length(mon)
supp=zeros(UInt8,n,lm)
for i=1:lm
    for j=1:n
        supp[j,i]=MultivariatePolynomials.degree(mon[i],x[j])
    end
end
supp=sparse(supp)

cliques,cql,cliquesize=clique_decomp(n,supp)
mclique,lmc,blocks,cl,blocksize,ub,sizes,basis=get_blocks_mix(d,supp,cliques,cql,cliquesize,ts=4,method="block")
blocks,cl,blocksize,ub,sizes=get_hblocks_mix!(basis,mclique,lmc,cliquesize,blocks,cl,blocksize,ub,sizes,method="block")
objv=blockupop_mix(n,d,supp,coe,cliques,cql,cliquesize,mclique,lmc,blocks,cl,blocksize,ts=4,QUIET=true)

-------------------------------------------------------------------
n=6
m=2
@polyvar x[1:n]
f=1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[4]*x[5]*x[6]
g1=1-sum(x[1:3].^2)
g2=1-sum(x[3:6].^2)
pop=[f,g1,g2]
coe=Array{Vector{Float64}}(undef, m+1)
supp=Array{SparseMatrixCSC}(undef, m+1)
for k=1:m+1
    mon=monomials(pop[k])
    coe[k]=coefficients(pop[k])
    lt=length(mon)
    ssupp=zeros(UInt8,n,lt)
    for i=1:lt
        for j=1:n
            ssupp[j,i]=MultivariatePolynomials.degree(mon[i],x[j])
        end
    end
    supp[k]=sparse(ssupp)
end

rlorder=[2,1,1]

cliques,cql,cliquesize=clique_cdecomp(n,m,supp,rlorder)
mclique,I,ncc,lmc,blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis=get_cblocks_mix(rlorder,m,supp,cliques,cql,cliquesize;ts=4,method="block")
blocks,cl,blocksize,ub,sizes=get_chblocks_mix!(I,ssupp,lt,fbasis,gbasis,mclique,lmc,cliquesize,blocks,cl,blocksize,ub,sizes,method="clique")
objv=blockcpop_mix(n,rlorder,supp,coe,cliques,cql,cliquesize,mclique,I,ncc,lmc,blocks,cl,blocksize,numeq=0,ts=4,QUIET=true)

-------------------------------------------------------------------

function blockupop_mix(n,d,supp::SparseMatrixCSC,coe,cliques,cql,cliquesize,mclique,lmc,blocks,cl,blocksize;ts=20,QUIET=true)
    basis=Array{SparseMatrixCSC}(undef,cql)
    col=Int[1]
    row=Int[]
    nz=UInt8[]
    for i=1:cql
        basis[i]=sparse_basis(cliques[i],n,d)
        if cliquesize[i]<ts
            ssupp=sparse_basis(cliques[i],n,2*d)
            col=[col;ssupp.colptr[2:end].+(col[end]-1)]
            row=[row;ssupp.rowval]
            nz=[nz;ssupp.nzval]
        end
    end
    if lmc>0
        for i=1:lmc
            tcol=[1;basis[mclique[i]].colptr]
            trow=basis[mclique[i]].rowval
            tnz=basis[mclique[i]].nzval
            for j=1:cl[i]
                for k=1:blocksize[i][j]
                    ind1=blocks[i][j][k]
                    t=ind1==1 ? 2 : k
                    for r=t:blocksize[i][j]
                        @inbounds ind2=blocks[i][j][r]
                        @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                        col=[col;col[end]+length(bi_row)]
                        row=[row;bi_row]
                        nz=[nz;bi_nz]
                    end
                end
            end
        end
    end
    supp1=SparseMatrixCSC(n,length(col),[1;col],row,nz)
    supp1=copy(sort_sparse(supp1))
    lsupp1=supp1.n
    model=Model(with_optimizer(Mosek.Optimizer, QUIET=QUIET))
    cons=[AffExpr(0) for i=1:lsupp1]
    pos1=Vector{Symmetric{VariableRef}}(undef, cql-lmc)
#    gram1=Array{Any}(undef, cql-lmc)
    for i=1:cql
        if cliquesize[i]<ts
            lb=basis[i].n+1
            pos1[i]=@variable(model, [1:lb, 1:lb], PSD)
            tcol=[1;basis[i].colptr]
            trow=basis[i].rowval
            tnz=basis[i].nzval
            for j=1:lb
                for k=j:lb
                    bi_row,bi_nz=splus(trow[tcol[j]:(tcol[j+1]-1)],tnz[tcol[j]:(tcol[j+1]-1)],trow[tcol[k]:(tcol[k+1]-1)],tnz[tcol[k]:(tcol[k+1]-1)])
                    Locb=bfind_sparse(supp1,bi_row,bi_nz)
                    if j==k
                        @inbounds add_to_expression!(cons[Locb],pos1[i][j,k])
                    else
                        @inbounds add_to_expression!(cons[Locb],2,pos1[i][j,k])
                    end
                end
            end
        end
    end
    if lmc>0
        pos2=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, lmc)
        for k=1:lmc
            pos2[k]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k])
            tcol=[1;basis[mclique[k]].colptr]
            trow=basis[mclique[k]].rowval
            tnz=basis[mclique[k]].nzval
            for i=1:cl[k]
                if blocksize[k][i]==1
                   pos2[k][i]=@variable(model, lower_bound=0)
                   bi_row=trow[tcol[blocks[k][i][1]]:(tcol[blocks[k][i][1]+1]-1)]
                   bi_nz=2*tnz[tcol[blocks[k][i][1]]:(tcol[blocks[k][i][1]+1]-1)]
                   Locb=bfind_sparse(supp1,bi_row,bi_nz)
                   @inbounds add_to_expression!(cons[Locb],pos2[k][i])
                else
                   pos2[k][i]=@variable(model, [1:blocksize[k][i], 1:blocksize[k][i]], PSD)
                   for j=1:blocksize[k][i]
                       ind1=blocks[k][i][j]
                       for r=j:blocksize[k][i]
                           @inbounds ind2=blocks[k][i][r]
                           @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                           Locb=bfind_sparse(supp1,bi_row,bi_nz)
                           if j==r
                               @inbounds add_to_expression!(cons[Locb],pos2[k][i][j,r])
                           else
                               @inbounds add_to_expression!(cons[Locb],2,pos2[k][i][j,r])
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
           println(i)
           println("The monomial basis is not enough!")
           return nothing,nothing
        else
           bc[Locb]=coe[i]
        end
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
#    for i=1:cl
#        gram[i]=value.(pos[i])
#    end
    return objv
end

function blockcpop_mix(n,rlorder,supp,coe,cliques,cql,cliquesize,mclique,I,ncc,lmc,blocks,cl,blocksize;numeq=0,ts=20,QUIET=true)
    fbasis=Array{SparseMatrixCSC}(undef,cql)
    col=Int[1]
    row=Int[]
    nz=UInt8[]
    for i=1:cql
        fbasis[i]=sparse_basis(cliques[i],n,rlorder[1])
        if cliquesize[i]<ts
            ssupp=sparse_basis(cliques[i],n,2*rlorder[1])
            col=[col;ssupp.colptr[2:end].+(col[end]-1)]
            row=[row;ssupp.rowval]
            nz=[nz;ssupp.nzval]
        end
    end
    if lmc>0
        for i=1:lmc
            tcol=[1;fbasis[mclique[i]].colptr]
            trow=fbasis[mclique[i]].rowval
            tnz=fbasis[mclique[i]].nzval
            for j=1:cl[i][1]
                for k=1:blocksize[i][1][j]
                    ind1=blocks[i][1][j][k]
                    t=ind1==1 ? 2 : k
                    for r=t:blocksize[i][1][j]
                        @inbounds ind2=blocks[i][1][j][r]
                        @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                        col=[col;col[end]+length(bi_row)]
                        row=[row;bi_row]
                        nz=[nz;bi_nz]
                    end
                end
            end
        end
    end
    gbasis=Array{SparseMatrixCSC}(undef,m-length(ncc))
    for i ∈ ncc
        col=[col;supp[i+1].colptr[2:end].+(col[end]-1)]
        row=[row;supp[i+1].rowval]
        nz=[nz;supp[i+1].nzval]
    end
    for s=1:cql
        if cliquesize[s]>=ts
            t=findfirst(isequal(s), mclique)
            for i ∈ I[s]
                gbasis[i]=sparse_basis(cliques[s],n,rlorder[i+1])
                tcol=[1;gbasis[i].colptr]
                trow=gbasis[i].rowval
                tnz=gbasis[i].nzval
                for j=1:cl[t][i]
                    for k=1:blocksize[t][i][j]
                        ind1=blocks[t][i][j][k]
                        for r=k:blocksize[t][i][j]
                            ind2=blocks[t][i][j][r]
                            for p=1:supp[i+1].n
                                @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                                @inbounds bi_row,bi_nz=splus(bi_row,bi_nz,supp[i+1].rowval[supp[i+1].colptr[p]:(supp[i+1].colptr[p+1]-1)],supp[i+1].nzval[supp[i+1].colptr[p]:(supp[i+1].colptr[p+1]-1)])
                                if length(bi_row)!=0
                                   col=[col;col[end]+length(bi_row)]
                                   row=[row;bi_row]
                                   nz=[nz;bi_nz]
                                end
                            end
                        end
                    end
                end
            end
        else
            for i ∈ I[s]
                gbasis[i]=sparse_basis(cliques[s],n,rlorder[i+1])
                ssupp=sparse_basis(cliques[s],n,2*rlorder[i+1])
                tcol=[1;ssupp.colptr]
                trow=ssupp.rowval
                tnz=ssupp.nzval
                for j=1:ssupp.n+1
                    for p=1:supp[i+1].n
                        @inbounds bi_row,bi_nz=splus(trow[tcol[j]:(tcol[j+1]-1)],tnz[tcol[j]:(tcol[j+1]-1)],supp[i+1].rowval[supp[i+1].colptr[p]:(supp[i+1].colptr[p+1]-1)],supp[i+1].nzval[supp[i+1].colptr[p]:(supp[i+1].colptr[p+1]-1)])
                        col=[col;col[end]+length(bi_row)]
                        row=[row;bi_row]
                        nz=[nz;bi_nz]
                    end
                end
            end
        end
    end
    supp1=SparseMatrixCSC(n,length(col),[1;col],row,nz)
    supp1=copy(sort_sparse(supp1))
    lsupp1=supp1.n
    model=Model(with_optimizer(Mosek.Optimizer, QUIET=QUIET))
    cons=[AffExpr(0) for i=1:lsupp1]
    pos1=Vector{Symmetric{VariableRef}}(undef, cql-lmc)
#    gram1=Array{Any}(undef, cql-lmc)
    p=1
    for i=1:cql
        if cliquesize[i]<ts
            lb=fbasis[i].n+1
            pos1[p]=@variable(model, [1:lb, 1:lb], PSD)
            tcol=[1;fbasis[i].colptr]
            trow=fbasis[i].rowval
            tnz=fbasis[i].nzval
            for j=1:lb
                for k=j:lb
                    @inbounds bi_row,bi_nz=splus(trow[tcol[j]:(tcol[j+1]-1)],tnz[tcol[j]:(tcol[j+1]-1)],trow[tcol[k]:(tcol[k+1]-1)],tnz[tcol[k]:(tcol[k+1]-1)])
                    Locb=bfind_sparse(supp1,bi_row,bi_nz)
                    if j==k
                        @inbounds add_to_expression!(cons[Locb],pos1[p][j,k])
                    else
                        @inbounds add_to_expression!(cons[Locb],2,pos1[p][j,k])
                    end
                end
            end
            p+=1
        end
    end
    if lmc>0
        pos2=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, lmc)
        for k=1:lmc
            pos2[k]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k][1])
            tcol=[1;fbasis[mclique[k]].colptr]
            trow=fbasis[mclique[k]].rowval
            tnz=fbasis[mclique[k]].nzval
            for i=1:cl[k][1]
                if blocksize[k][1][i]==1
                   @inbounds pos2[k][i]=@variable(model, lower_bound=0)
                   @inbounds bi_row=trow[tcol[blocks[k][1][i][1]]:(tcol[blocks[k][1][i][1]+1]-1)]
                   @inbounds bi_nz=2*tnz[tcol[blocks[k][1][i][1]]:(tcol[blocks[k][1][i][1]+1]-1)]
                   Locb=bfind_sparse(supp1,bi_row,bi_nz)
                   @inbounds add_to_expression!(cons[Locb],pos2[k][i])
                else
                   @inbounds bs=blocksize[k][1][i]
                   @inbounds pos2[k][i]=@variable(model, [1:bs, 1:bs], PSD)
                   for j=1:bs
                       @inbounds ind1=blocks[k][1][i][j]
                       for r=j:bs
                           @inbounds ind2=blocks[k][1][i][r]
                           @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                           Locb=bfind_sparse(supp1,bi_row,bi_nz)
                           if j==r
                              @inbounds add_to_expression!(cons[Locb],pos2[k][i][j,r])
                           else
                              @inbounds add_to_expression!(cons[Locb],2,pos2[k][i][j,r])
                           end
                       end
                   end
                end
            end
        end
    end
    pos3=Vector{VariableRef}(undef, length(ncc))
    for i ∈ ncc
        if i<=m-numeq
            pos3[i]=@variable(model, lower_bound=0)
        else
            pos3[i]=@variable(model)
        end
        for j=1:supp[i+1].n
            Locb=bfind_sparse(supp1,supp[i+1].rowval[supp[i+1].colptr[j]:(supp[i+1].colptr[j+1]-1)],supp[i+1].nzval[supp[i+1].colptr[j]:(supp[i+1].colptr[j+1]-1)])
            add_to_expression!(cons[Locb],coe[i+1][j],pos3[i])
        end
    end
    bcon=sum([length(I[mclique[i]]) for i=1:lmc])
    pos4=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, bcon)
    nbcon=m-length(ncc)-bcon
    pos5=Vector{Symmetric{VariableRef}}(undef, bcon)
    p=1
    q=1
    for s=1:cql
        if cliquesize[s]>=ts
            t=findfirst(isequal(s), mclique)
            for i ∈ I[s]
                tcol=[1;gbasis[i].colptr]
                trow=gbasis[i].rowval
                tnz=gbasis[i].nzval
                pos4[p]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef,cl[t][i])
                for j=1:cl[t][i]
                    bs=blocksize[t][i][j]
                    if bs==1
                        if i<=m-numeq
                           pos4[p][j]=@variable(model, lower_bound=0)
                        else
                           pos4[p][j]=@variable(model)
                        end
                        for k=1:supp[i+1].n
                            @inbounds ind1=blocks[t][i][j][1]
                            @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],2*tnz[tcol[ind1]:(tcol[ind1+1]-1)],supp[i+1].rowval[supp[i+1].colptr[k]:(supp[i+1].colptr[k+1]-1)],supp[i+1].nzval[supp[i+1].colptr[k]:(supp[i+1].colptr[k+1]-1)])
                            Locb=bfind_sparse(supp1,bi_row,bi_nz)
                            @inbounds add_to_expression!(cons[Locb],coe[i+1][k],pos4[p][j])
                        end
                    else
                        if i<=m-numeq
                           pos4[p][j]=@variable(model, [1:bs, 1:bs], PSD)
                        else
                           pos4[p][j]=@variable(model, [1:bs, 1:bs], Symmetric)
                        end
                        for k=1:bs
                            ind1=blocks[t][i][j][k]
                            for r=k:bs
                                ind2=blocks[t][i][j][r]
                                for u=1:supp[i+1].n
                                    @inbounds bi_row,bi_nz=splus(trow[tcol[ind1]:(tcol[ind1+1]-1)],tnz[tcol[ind1]:(tcol[ind1+1]-1)],trow[tcol[ind2]:(tcol[ind2+1]-1)],tnz[tcol[ind2]:(tcol[ind2+1]-1)])
                                    @inbounds bi_row,bi_nz=splus(bi_row,bi_nz,supp[i+1].rowval[supp[i+1].colptr[u]:(supp[i+1].colptr[u+1]-1)],supp[i+1].nzval[supp[i+1].colptr[u]:(supp[i+1].colptr[u+1]-1)])
                                    Locb=bfind_sparse(supp1,bi_row,bi_nz)
                                    if k==r
                                        @inbounds add_to_expression!(cons[Locb],coe[i+1][u],pos4[p][j][k,r])
                                    else
                                        @inbounds add_to_expression!(cons[Locb],2*coe[i+1][u],pos4[p][j][k,r])
                                    end
                                end
                            end
                        end
                    end
                end
                p+=1
            end
        else
            for i ∈ I[s]
                tcol=[1;gbasis[i].colptr]
                trow=gbasis[i].rowval
                tnz=gbasis[i].nzval
                lb=gbasis[i].n+1
                if i<=m-numeq
                   pos5[q]=@variable(model, [1:lb, 1:lb], PSD)
                else
                   pos5[q]=@variable(model, [1:lb, 1:lb], Symmetric)
                end
                for j=1:lb
                    for k=j:lb
                        for r=1:supp[i+1].n
                            @inbounds bi_row,bi_nz=splus(trow[tcol[j]:(tcol[j+1]-1)],tnz[tcol[j]:(tcol[j+1]-1)],trow[tcol[k]:(tcol[k+1]-1)],tnz[tcol[k]:(tcol[k+1]-1)])
                            @inbounds bi_row,bi_nz=splus(bi_row,bi_nz,supp[i+1].rowval[supp[i+1].colptr[r]:(supp[i+1].colptr[r+1]-1)],supp[i+1].nzval[supp[i+1].colptr[r]:(supp[i+1].colptr[r+1]-1)])
                            Locb=bfind_sparse(supp1,bi_row,bi_nz)
                            if j==k
                                @inbounds add_to_expression!(cons[Locb],coe[i+1][r],pos5[q][j,k])
                            else
                                @inbounds add_to_expression!(cons[Locb],2*coe[i+1][r],pos5[q][j,k])
                            end
                        end
                    end
                end
                q+=1
            end
        end
    end
    bc=zeros(lsupp1,1)
    for i=1:supp[1].n
        Locb=bfind_sparse(supp1,supp[1].rowval[supp[1].colptr[i]:(supp[1].colptr[i+1]-1)],supp[1].nzval[supp[1].colptr[i]:(supp[1].colptr[i+1]-1)])
        if Locb==0
           println("The monomial basis is not enough!")
           return nothing,nothing
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
#    for i=1:cl
#        gram[i]=value.(pos[i])
#    end
    return objv
end

function get_blocks_mix(d,supp::SparseMatrixCSC,cliques,cql,cliquesize;ts=20,method="block")
    mclique=UInt8[]
    for i=1:cql
        if cliquesize[i]>=ts
            push!(mclique,i)
        end
    end
    lmc=length(mclique)
    blocks=Vector{Vector{Vector{UInt16}}}(undef,lmc)
    cl=Vector{UInt16}(undef,lmc)
    ub=Vector{Vector{UInt16}}(undef,lmc)
    sizes=Vector{Vector{UInt16}}(undef,lmc)
    blocksize=Vector{Vector{Int}}(undef,lmc)
    basis=Vector{Array{UInt8,2}}(undef,lmc)
    println("---------------------------------------------------")
    println("The block sizes for varible cliques of a large size:")
    for i=1:lmc
        ind=mclique[i]
        nvar=cliquesize[ind]
        ssupp=zeros(UInt8,nvar,1)
        for j=1:supp.n
            if issubset(supp.rowval[supp.colptr[j]:(supp.colptr[j+1]-1)], cliques[ind])
                bi=zeros(UInt8,nvar,1)
                for k=supp.colptr[j]:(supp.colptr[j+1]-1)
                    @inbounds locb=lbfind(cliques[ind],cliquesize[ind],supp.rowval[k])
                    @inbounds bi[locb]=supp.nzval[k]
                end
                ssupp=[ssupp bi]
            end
        end
        basis[i]=get_basis(nvar,d)
        if method=="block"
            blocks[i],cl[i],blocksize[i],ub[i],sizes[i]=get_blocks(nvar,ssupp,basis[i])
        else
            blocks[i],cl[i],blocksize[i],ub[i],sizes[i]=get_cliques(nvar,ssupp,basis[i])
        end
    end
    return mclique,lmc,blocks,cl,blocksize,ub,sizes,basis
end

function get_hblocks_mix!(basis,mclique,lmc,cliquesize,blocks,cl,blocksize,ub,sizes;method="block")
    nub=Vector{Vector{UInt16}}(undef,lmc)
    nsizes=Vector{Vector{UInt16}}(undef,lmc)
    println("---------------------------------------------------")
    println("The block sizes for varible cliques of a large size:")
    for i=1:lmc
        ind=mclique[i]
        nvar=cliquesize[ind]
        ssupp=zeros(UInt8,nvar,Int(sum(blocksize[i].^2+blocksize[i])/2))
        k=1
        for s=1:cl[i]
            for j=1:blocksize[i][s]
                for r=j:blocksize[i][s]
                    @inbounds bi=basis[i][:,blocks[i][s][j]]+basis[i][:,blocks[i][s][r]]
                    @inbounds ssupp[:,k]=bi
                    k+=1
                end
            end
        end
        ssupp=sortslices(ssupp,dims=2)
        ssupp=unique(ssupp,dims=2)
        if method=="block"
            blocks[i],cl[i],blocksize[i],nub[i],nsizes[i],status=get_hblocks!(nvar,ssupp,basis[i],ub[i],sizes[i])
        else
            blocks[i],cl[i],blocksize[i],nub[i],nsizes[i],status=get_hcliques!(nvar,ssupp,basis[i],ub[i],sizes[i])
        end
    end
    return blocks,cl,blocksize,nub,nsizes
end

function get_cblocks_mix(rlorder,m,supp,cliques,cql,cliquesize;ts=20,method="block")
    mclique=UInt8[]
    for i=1:cql
        if cliquesize[i]>=ts
            push!(mclique,i)
        end
    end
    lmc=length(mclique)
    I=[UInt16[] for i=1:cql]
    ncc=UInt16[]
    for i=2:m+1
        if rlorder[i]==0
            push!(ncc, i-1)
        else
            rind=unique(supp[i].rowval)
            ind=1
            while !issubset(rind, cliques[ind])
                ind+=1
            end
            push!(I[ind], i-1)
        end
    end
    blocks=Vector{Vector{Vector{Vector{UInt16}}}}(undef,lmc)
    cl=Vector{Vector{UInt16}}(undef,lmc)
    ub=Vector{Vector{UInt16}}(undef,lmc)
    sizes=Vector{Vector{UInt16}}(undef,lmc)
    blocksize=Vector{Vector{Vector{Int}}}(undef,lmc)
    ssupp=Vector{Vector{Array{UInt8,2}}}(undef,lmc)
    lt=Vector{Vector{UInt16}}(undef,lmc)
    gbasis=Vector{Vector{Array{UInt8,2}}}(undef,lmc)
    fbasis=Vector{Array{UInt8,2}}(undef,lmc)
    lsupp=sum([supp[i].n for i=1:m+1])
    col=Int[1]
    row=Int[]
    nz=UInt8[]
    for i=1:m+1
        col=[col;supp[i].colptr[2:end].+(col[end]-1)]
        row=[row;supp[i].rowval]
        nz=[nz;supp[i].nzval]
    end
    println("---------------------------------------------------")
    println("The block sizes for varible cliques of a large size:")
    for i=1:lmc
        lc=length(I[mclique[i]])
        blocks[i]=Vector{Vector{Vector{UInt16}}}(undef, lc+1)
        cl[i]=Vector{UInt16}(undef, lc+1)
        ub[i]=Vector{UInt16}(undef, lc+1)
        sizes[i]=Vector{UInt16}(undef, lc+1)
        blocksize[i]=Vector{Vector{Int}}(undef, lc+1)
        ind=mclique[i]
        nvar=cliquesize[ind]
        fsupp=zeros(UInt8,nvar,1)
        for j=1:lsupp
            if issubset(row[col[j]:(col[j+1]-1)], cliques[ind])
                bi=zeros(UInt8,nvar,1)
                for k=col[j]:(col[j+1]-1)
                    @inbounds locb=lbfind(cliques[ind],cliquesize[ind],row[k])
                    @inbounds bi[locb]=nz[k]
                end
                fsupp=[fsupp bi]
            end
        end
        fsupp=unique(fsupp,dims=2)
        fbasis[i]=get_basis(cliquesize[ind],rlorder[1])
        gbasis[i]=Vector{Array{UInt8,2}}(undef, lc)
        ssupp[i]=Vector{Array{UInt8,2}}(undef, lc+1)
        lt[i]=Vector{UInt16}(undef, lc+1)
        lt[i][1]=supp[1].n
        ssupp[i][1]=zeros(UInt8,nvar,supp[1].n)
        for s=1:lc
            t=I[mclique[i]][s]
            gbasis[i][s]=get_basis(nvar,rlorder[t+1])
            ssupp[i][s+1]=zeros(UInt8,nvar,supp[t+1].n)
            lt[i][s+1]=supp[t+1].n
            for j=1:supp[t+1].n
                for k=supp[t+1].colptr[j]:(supp[t+1].colptr[j+1]-1)
                    @inbounds locb=lbfind(cliques[ind],cliquesize[ind],supp[t+1].rowval[k])
                    @inbounds ssupp[i][s+1][locb,j]=supp[t+1].nzval[k]
                end
            end
        end
        if method=="block"
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],ub[i],sizes[i]=get_cblocks(nvar,lc,fsupp,ssupp[i],lt[i],fbasis[i],gbasis[i])
        else
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],ub[i],sizes[i]=get_ccliques(nvar,lc,fsupp,ssupp[i],lt[i],fbasis[i],gbasis[i])
        end
    end
    return mclique,I,ncc,lmc,blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis
end

function get_chblocks_mix!(I,ssupp,lt,fbasis,gbasis,mclique,lmc,cliquesize,blocks,cl,blocksize,ub,sizes;method="block")
    nub=Vector{Vector{UInt16}}(undef,lmc)
    nsizes=Vector{Vector{UInt16}}(undef,lmc)
    println("---------------------------------------------------")
    println("The block sizes for varible cliques of a large size:")
    for i=1:lmc
        lc=length(I[mclique[i]])
        ind=mclique[i]
        nvar=cliquesize[ind]
        fsupp=zeros(UInt8,nvar,Int(sum(blocksize[i][1].^2+blocksize[i][1])/2))
        k=1
        for s=1:cl[i][1]
            for j=1:blocksize[i][1][s]
                for r=j:blocksize[i][1][s]
                    @inbounds bi=fbasis[i][:,blocks[i][1][s][j]]+fbasis[i][:,blocks[i][1][s][r]]
                    @inbounds fsupp[:,k]=bi
                    k+=1
                end
            end
        end
        fsupp=sortslices(fsupp,dims=2)
        fsupp=unique(fsupp,dims=2)
        if method=="block"
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],nub[i],nsizes[i],status=get_chblocks!(nvar,lc,ssupp[i],lt[i],fbasis[i],gbasis[i],fsupp,ub[i],sizes[i])
        else
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],nub[i],nsizes[i],status=get_chcliques!(nvar,lc,ssupp[i],lt[i],fbasis[i],gbasis[i],fsupp,ub[i],sizes[i])
        end
    end
    return blocks,cl,blocksize,nub,nsizes
end

function clique_decomp(n,supp::SparseMatrixCSC)
    A=zeros(UInt8,n,n)
    for i = 1:supp.n
        lcol=supp.colptr[i+1]-supp.colptr[i]
        A[supp.rowval[supp.colptr[i]:(supp.colptr[i+1]-1)],supp.rowval[supp.colptr[i]:(supp.colptr[i+1]-1)]]=ones(UInt8,lcol,lcol)
    end
    cliques,cql,cliquesize=cliquesFromSpMatD(A)
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("The clique sizes of varibles:\n$uc\n$sizes")
    return cliques,cql,cliquesize
end

function clique_cdecomp(n,m,supp,rlorder)
    A=zeros(UInt8,n,n)
    for i=1:m+1
        if i==1||rlorder[i]==0
            for j = 1:supp[i].n
                lcol=supp[i].colptr[j+1]-supp[i].colptr[j]
                A[supp[i].rowval[supp[i].colptr[j]:(supp[i].colptr[j+1]-1)],supp[i].rowval[supp[i].colptr[j]:(supp[i].colptr[j+1]-1)]]=ones(UInt8,lcol,lcol)
            end
        else
            rind=unique(supp[i].rowval)
            lcol=length(rind)
            A[rind,rind]=ones(UInt8,lcol,lcol)
        end
    end
    cliques,cql,cliquesize=cliquesFromSpMatD(A)
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("The clique sizes of varibles:\n$uc\n$sizes")
    return cliques,cql,cliquesize
end

function splus(arow,anz,brow,bnz)
    la=length(arow)
    lb=length(brow)
    row=UInt16[]
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
    n=length(var)
    lb=binomial(n+d,d)-1
    col=UInt16[1;2]
    row=UInt16[var[1]]
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

function sort_sparse(s::SparseMatrixCSC)
    m=s.m
    n=s.n
    col=s.colptr
    row=s.rowval
    nz=s.nzval
    i=2
    while i<=n
        corr_row=row[col[i]:(col[i+1]-1)]
        corr_nz=nz[col[i]:(col[i+1]-1)]
        j=i-1
        comp=0
        while j>=1
            pre_row=row[col[j]:(col[j+1]-1)]
            pre_nz=nz[col[j]:(col[j+1]-1)]
            comp=comp_sparse(corr_row,corr_nz,pre_row,pre_nz)
            if comp==-1
                j-=1
            else
                break
            end
        end
        if j>=1&&comp==0
            deleteat!(row,col[i]:(col[i+1]-1))
            deleteat!(nz,col[i]:(col[i+1]-1))
            deleteat!(col,i)
            col[i:end].-=length(corr_row)
            n-=1
            i-=1
            while j>=2
                j-=1
                pre_row=row[col[j]:(col[j+1]-1)]
                pre_nz=nz[col[j]:(col[j+1]-1)]
                if comp_sparse(corr_row,corr_nz,pre_row,pre_nz)==0
                    deleteat!(row,col[j]:(col[j+1]-1))
                    deleteat!(nz,col[j]:(col[j+1]-1))
                    deleteat!(col,j)
                    col[j:end].-=length(corr_row)
                    n-=1
                    i-=1
                else
                    break
                end
            end
        elseif j<i-1
            for k=i:-1:j+2
                lprow=col[k]-col[k-1]
                row[(col[k+1]-lprow):(col[k+1]-1)]=row[(col[k-1]):(col[k]-1)]
                nz[(col[k+1]-lprow):(col[k+1]-1)]=nz[(col[k-1]):(col[k]-1)]
                col[k]=col[k+1]-lprow
            end
            row[(col[j+1]):(col[j+2]-1)]=corr_row
            nz[(col[j+1]):(col[j+2]-1)]=corr_nz
        else
        end
        i+=1
    end
    return SparseMatrixCSC(m,n,col,row,nz)
end

function comp_sparse(corr_row,corr_nz,pre_row,pre_nz)
    if sum(corr_nz)<sum(pre_nz)
        return -1
    elseif sum(corr_nz)>sum(pre_nz)
        return 1
    else
        i=1
        lc=length(corr_row)
        lp=length(pre_row)
        while i<=lc&&i<=lp
            if corr_row[end+1-i]>pre_row[end+1-i]
                return 1
            elseif corr_row[end+1-i]<pre_row[end+1-i]
                return -1
            elseif corr_nz[end+1-i]<pre_nz[end+1-i]
                return -1
            elseif corr_nz[end+1-i]>pre_nz[end+1-i]
                return 1
            else
                i+=1
            end
        end
    end
    return 0
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
