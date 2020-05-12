using DynamicPolynomials
using MultivariatePolynomials
using SparseArrays
using TSSOS

# unconstrained optimization using the TSSOS hierarchy
@polyvar x[1:2]
f=x[1]^4+x[2]^4-x[1]*x[2]
opt,sol,data=blockupop_first(f,x,newton=false,method="block")
opt,sol,data=blockupop_higher!(data,method="chordal")

# constrained optimization using the TSSOS hierarchy
@polyvar x[1:2]
f=x[1]^4+x[2]^4-x[1]*x[2]+1
g=1-x[1]^2-2*x[2]^2
h=x[1]^2+x[2]^2-1
pop=[f,g,h]
opt,sol,data=blockcpop_first(pop,x,2,numeq=1,method="block")
opt,sol,data=blockcpop_higher!(data,method="chordal")

# unconstrained optimization using the CS-TSSOS hierarchy
# n=6
# d=2
# @polyvar x[1:n]
# f=1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
# mon=monomials(f)
# coe=coefficients(f)
# lm=length(mon)
# supp=zeros(UInt8,n,lm)
# for i=1:lm
#     for j=1:n
#         supp[j,i]=MultivariatePolynomials.degree(mon[i],x[j])
#     end
# end
# supp=sparse(supp)
# cliques,cql,cliquesize=clique_decomp(n,supp)
# blocks,cl,blocksize,ub,sizes,basis=get_blocks_mix(d,supp,cliques,cql,cliquesize,method="chordal")
# blocks,cl,blocksize,ub,sizes,status=get_hblocks_mix!(supp0,basis,cliques,cql,cliquesize,blocks,cl,blocksize,ub,sizes,method="chordal",chor_alg="greedy")
# objv,supp0,moment=blockupop_mix(n,d,supp,coe,cliques,cql,cliquesize,blocks,cl,blocksize)

# constrained optimization using the CS-TSSOS hierarchy
n=6
@polyvar x[1:n]
f=1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
pop=[f,1-sum(x[1:3].^2),1-sum(x[1:4].^2)]
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
order=2
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,order,numeq=0,TS="block")
opt,sol,data=cs_tssos_higher!(data,TS="greedy")
