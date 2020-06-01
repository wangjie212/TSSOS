# using DynamicPolynomials
# using MultivariatePolynomials
# using SparseArrays
# using TSSOS

# unconstrained optimization using the TSSOS hierarchy
@polyvar x[1:3]
f=1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
opt,sol,data=blockupop_first(f,x,newton=true,TS="block",solution=true)
opt,sol,data=blockupop_higher!(data,reducebasis=true,TS="MD",solution=true)

# constrained optimization using the TSSOS hierarchy
@polyvar x[1:3]
f=1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
g=1-x[1]^2-2*x[2]^2
h=x[2]^2+x[3]^2-1
pop=[f,g,h]
opt,sol,data=blockcpop_first(pop,x,2,numeq=1,TS="block",solution=true)
opt,sol,data=blockcpop_higher!(data,reducebasis=true,TS="MD",solution=true)

# unconstrained optimization using the CS-TSSOS hierarchy
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
blocks,cl,blocksize,ub,sizes,basis,status=get_blocks_mix(d,supp,nothing,cliques,cql,cliquesize,nothing,nothing,TS="block")
objv,supp0,moment=blockupop_mix(n,d,supp,coe,cliques,cql,cliquesize,blocks,cl,blocksize)
blocks,cl,blocksize,ub,sizes,basis,status=get_blocks_mix(d,supp0,basis,cliques,cql,cliquesize,ub,sizes,TS="MD")
objv,supp0,moment=blockupop_mix(n,d,supp,coe,cliques,cql,cliquesize,blocks,cl,blocksize)

# constrained optimization using the CS-TSSOS hierarchy
n=6
@polyvar x[1:n]
f=1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
pop=[f,1-sum(x[1:3].^2),1-sum(x[3:6].^2)]
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
dg=zeros(Int,m)
for i=1:m
    dg[i]=maxdegree(pop[i+1])
end
order=2
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,order,numeq=0,TS="block")
opt,sol,data=cs_tssos_higher!(data,TS="MD")
