using DynamicPolynomials
using MultivariatePolynomials
using SparseArrays
# using TSSOS

# unconstrained optimization using the TSSOS hierarchy
@polyvar x[1:3]
f=1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
opt,sol,data=tssos_first(f,x,newton=true,reducebasis=true,TS="MD",solution=true)
opt,sol,data=tssos_higher!(data,TS="MD",solution=true)

# constrained optimization using the TSSOS hierarchy
@polyvar x[1:3]
f=1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
g=1-x[1]^2-2*x[2]^2
h=x[2]^2+x[3]^2-1
pop=[f,g,h]
opt,sol,data=tssos_first(pop,x,2,numeq=1,reducebasis=true,TS="MD",solution=true)
opt,sol,data=tssos_higher!(data,TS="MD",solution=true)

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
blocks,cl,blocksize,ub,sizes,basis,status=get_blocks_mix(d,supp,cliques,cql,cliquesize,TS="MD")
objv,supp0,moment=blockupop_mix(n,d,supp,coe,cliques,cql,cliquesize,blocks,cl,blocksize)
blocks,cl,blocksize,ub,sizes,basis,status=get_blocks_mix(d,supp0,cliques,cql,cliquesize,basis=basis,ub=ub,sizes=sizes,TS="MD")
objv,supp0,moment=blockupop_mix(n,d,supp,coe,cliques,cql,cliquesize,blocks,cl,blocksize)

# constrained optimization using the CS-TSSOS hierarchy
n=6
@polyvar x[1:n]
f=1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
pop=[f,1-sum(x[1:3].^2),1-sum(x[3:6].^2)]
d=2
opt,sol,data=cs_tssos_first(pop,x,d,TS="MD",solution=true)
opt,sol,data=cs_tssos_higher!(data,TS="MD",solution=true)
