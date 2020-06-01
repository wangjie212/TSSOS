using DynamicPolynomials
using MultivariatePolynomials
using SparseArrays
using TSSOS

# unconstrained optimization using the CS-TSSOS hierarchy
# the Broyden banded function
n=20
@polyvar x[1:n]
jset=Array{Any}(undef, n)
f=x[1]
for i=1:n
    jset[i]=max(1,i-5):min(n,i+1)
    jset[i]=setdiff(jset[i],i)
    f0=(1+x[jset[i][1]])*x[jset[i][1]]
    ljset=length(jset[i])
    for j=2:ljset
        f0+=(1+x[jset[i][j]])*x[jset[i][j]]
    end
    global f+=(x[i]*(2+5*x[i]^2)+1-f0)^2
end
f-=x[1]
@time begin
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
d=3
blocks,cl,blocksize,ub,sizes,basis,status=get_blocks_mix(d,supp,nothing,cliques,cql,cliquesize,nothing,nothing,TS="MD")
objv,supp0,moment=blockupop_mix(n,d,supp,coe,cliques,cql,cliquesize,blocks,cl,blocksize)
end
@time begin
objv,supp0,moment=blockupop_mix(n,d,supp,coe,cliques,cql,cliquesize,nothing,nothing,nothing,mix=false)
end

# constrained optimization using the chordal-TSSOS hierarchy
# the modified Rosenbrock function
l=2
p=20
n=l*p
@polyvar x[1:n]
f=x[1]+1
for i=2:n
    global f+=100*(x[i]-x[i-1]^2)^2+(1-x[i])^2
end
f-=x[1]
pop=[f]
for i=1:l
    push!(pop,1-sum(x[(i-1)*p+1:i*p].^2))
end
m=length(pop)-1
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
dg=2*ones(Int,m)
order=2 # the relaxation order
@time begin
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,order,TS="MD",extra_sos=false)
end
@time begin
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,order,TS=false)
end

# the Broyden tridiagonal function
l=2
p=20
n=l*p
@polyvar x[1:n]
f=((3-2*x[1])*x[1]-2*x[2]+1)^2
for i=2:n-1
    global f+=((3-2*x[i])*x[i]-x[i-1]-2*x[i+1]+1)^2
end
f+=((3-2*x[n])*x[n]-x[n-1]+1)^2
pop=[f]
for i=1:l
    push!(pop,1-sum(x[(i-1)*p+1:i*p].^2))
end
m=length(pop)-1
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
dg=2*ones(Int,m)
@time begin
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,2,TS="MD",extra_sos=false)
end
@time begin
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,order,TS=false)
end

# the chained Wood function
l=2
p=20
n=l*p
@polyvar x[1:n]
f=x[1]+1
for i=1:2:n-3
    global f+=100*(x[i+1]-x[i]^2)^2+(1-x[i])^2+90*(x[i+3]-x[i+2]^2)^2+(1-x[i+2])^2+10*(x[i+1]+x[i+3]-2)^2+0.1*(x[i+1]-x[i+3])^2
end
f-=x[1]
pop=[f]
for i=1:l
    push!(pop,1-sum(x[(i-1)*p+1:i*p].^2))
end
m=length(pop)-1
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
dg=2*ones(Int,m)
@time begin
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,2,TS="MD",extra_sos=false)
end
@time begin
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,order,TS=false)
end

# Max-Cut problems
using LightGraphs
using DelimitedFiles
using MetaGraphs

io=open("E:\\Programs\\blockpop\\MaxCut\\g20.txt", "r")
data_cells, header_cells=readdlm(io, ' ', Int, '\n', header=true)
close(io)
lV=parse(Int, header_cells[1])
lE=parse(Int, header_cells[2])
G=MetaGraph(lV)
for i=1:lE
    add_edge!(G, data_cells[i,1], data_cells[i,2])
    set_prop!(G, data_cells[i,1], data_cells[i,2], :weight, data_cells[i,3])
end
n=lV
m=lV
numeq=lV
coe=Array{Vector{Float64}}(undef, m+1)
supp=Array{SparseMatrixCSC}(undef, m+1)
coe[1]=Float64[0]
col=UInt32[1, 1]
row=UInt32[]
nz=UInt8[]
for edge in edges(G)
    weight=get_prop(G, edge, :weight)
    coe[1][1]-=0.5*weight
    push!(coe[1], 0.5*weight)
    push!(col, col[end]+2)
    push!(row, src(edge), dst(edge))
    push!(nz, 1, 1)
end
supp[1]=SparseMatrixCSC(lV, lE+1, col, row, nz)
for i=1:lV
    coe[i+1]=[1;-1]
    col=[1;1;2]
    row=[i]
    nz=[2]
    supp[i+1]=SparseMatrixCSC(lV, 2, col, row, nz)
end
dg=2*ones(Int, lV)

# Shor relaxation
@time begin
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,1,numeq=numeq,TS=false)
end

# the 2-nd CS-TSSOS hierarchy
@time begin
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,2,numeq=numeq,TS="block")
end

# AC-OPF problems
using PowerModels
include("E:\\Programs\\blockpop\\OPF\\modelopf.jl")
cd("E:\\Programs\\blockpop\\OPF\\pglib")

opf_data = parse_file("pglib_opf_case5_pjm.m")

model=pop_opf_two(opf_data,normal=true)
n=model.n
m=model.m
numeq=model.numeq
nbus=model.nbus
ng=model.ng
nb=model.nb
supp=model.supp
coe=model.coe
dg=model.dg
# factor=maximum(abs.(coe[1]))
# coe[1]=coe[1]./factor

# Shor relaxation
@time begin
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,1,numeq=numeq,TS=false)
end

model=pop_opf(opf_data,vmc="quadratic",gen_model="two",normal=true)
n=model.n
m=model.m
numeq=model.numeq
nbus=model.nbus
ng=model.ng
nb=model.nb
supp=model.supp
coe=model.coe
dg=model.dg
# factor=maximum(abs.(coe[1]))
# coe[1]=coe[1]./factor

# the 2-nd CS-TSSOS hierarchy with multi order
@time begin
cliques,cql,cliquesize=clique_cdecomp(n,m,dg,supp,order="multi",alg="MF",minimize=true)
I,ncc=assign_constraint(m,supp,cliques,cql,cliquesize,assign="min")
rlorder=init_order(dg,I,cql,order="multi")
blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis,status=get_cblocks_mix!(dg,I,rlorder,m,supp,nothing,nothing,nothing,nothing,nothing,cliques,cql,cliquesize,nothing,nothing,nothing,nothing,nothing,TS="block")
opt,supp0,_,_,_=blockcpop_mix(n,m,dg,rlorder,supp,coe,cliques,cql,cliquesize,I,ncc,blocks,cl,blocksize,numeq=numeq,mix=true,QUIET=false,solve=true,solution=false,extra_sos=true)
end

# the 2-nd CSSOS hierarchy with multi order
@time begin
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,"multi",numeq=numeq,TS=false)
end

# the 2-nd CS-TSSOS hierarchy with uniform order
# @time begin
# cliques,cql,cliquesize=clique_opf_four(m,nbus,nb,supp,vmc="quadratic",alg="MF",minimize=true)
# I,ncc=assign_constraint(m,supp,cliques,cql,cliquesize,assign="min")
# rlorder=init_order(dg,I,cql,order=2)
# blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis,status=get_cblocks_mix!(dg,I,rlorder,m,supp,nothing,nothing,nothing,nothing,nothing,cliques,cql,cliquesize,nothing,nothing,nothing,nothing,nothing,TS="MD")
# opt,supp0,supp1,measure,moment=blockcpop_mix(n,m,dg,rlorder,supp,coe,cliques,cql,cliquesize,I,ncc,blocks,cl,blocksize,numeq=numeq,mix=true,QUIET=false,solve=true,solution=false,extra_sos=true)
# end

# the 2-nd CSSOS hierarchy with uniform order
# @time begin
# opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,2,numeq=numeq,TS=false)
# end
