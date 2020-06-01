# TSSOS
TSSOS is a sparse polynomial optimization tool based on block moment-SOS hierarchies. To use TSSOS in Julia, run
```Julia
pkg> add https://github.com/wangjie212/TSSOS
 ```

## Dependencies
- Julia
- MOSEK
- JuMP

TSSOS has been tested on WINDOW 10, Julia 1.2, JuMP 0.21 and MOSEK 8.1.
## Usage
### Unconstrained polynomial optimization problems
The unconstrained polynomial optimization problem formulizes as
```
Inf{f(x): x∈R^n}
```
where f is a polynomial with variables x1,...,xn and of degree d.

Taking f=1+x1^4+x2^4+x3^4+x1\*x2\*x3+x2 as an example, to exetute the first TSSOS hierarchy, run
```Julia
using TSSOS
using DynamicPolynomials
@polyvar x[1:3]
f=1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
opt,sol,data=blockupop_first(f,x,TS="block")
```
By default, the monomial basis computed by the Newton polytope method is used. If one sets newton=false in the input,
```Julia
opt,sol,data=blockupop_first(f,x,newton=false,TS="block")
```
then the standard monomial basis will be used.

Two vectors will be outputed. The first vector is the size of blocks and the second vector is the number of blocks with size corresponding to the first vector.

To exetute higher TSSOS hierarchies, repeatedly run

```Julia
opt,sol,data=blockupop_higher!(data,TS="block")
```

Options:  
nb: specify the first nb variables to be binary variables (satisfying xi^2=xi)  
newton: true (using the monomial basis computed by the Newton polytope method), false  
TS (term sparsity): "block" (using the TSSOS hierarchy), "MD" or "MF" (using the chordal-TSSOS hierarchy), false (without term sparsity)  
solution: true (extract a solution), false (don't extract a solution)

### Constrained polynomial optimization problems
The constrained polynomial optimization problem formulizes as
```
Inf{f(x): x∈K}
```
where f is a polynomial and K is the basic semi-algebraic set
```
K={x∈R^n: g_j(x)>=0, j=1,...,m-numeq, g_j(x)=0, j=m-numeq+1,...,m},
```
for some polynomials g_j, j=1,...,m.

Taking f=1+x1^4+x2^4+x3^4+x1\*x2\*x3+x2 and K={x∈R^2: g_1=1-x1^2-2\*x2^2>=0, g_2=x2^2+x3^2-1=0} as an example, to exetute the first TSSOS hierarchy, run

```Julia
@polyvar x[1:3]
f=1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
g_1=1-x[1]^2-2*x[2]^2
g_2=x[2]^2+x[3]^2-1
pop=[f,g_1,g_2]
d=2 # the relaxation order of Lasserre hierarchy
opt,sol,data=blockcpop_first(pop,x,d,numeq=1,TS="block")
```

To exetute higher TSSOS hierarchies, repeatedly run

```Julia
julia> opt,sol,data=blockcpop_higher!(data,TS="block")
```

Options:  
nb: specify the first nb variables to be binary variables (satisfying xi^2=xi)  
TS: "block" (using the TSSOS hierarchy), "MD" or "MF" (using the chordal-TSSOS hierarchy), false (without term sparsity)  
solution: true (extract a solution), false (don't extract a solution)

One can also exploit correlative sparsity and term sparsity simultaneously, which is called the CS-TSSOS hierarchy.

```Julia
using SparseArrays
using MultivariatePolynomials
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
order=2 # the relaxation order of Lasserre hierarchy
opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,order,numeq=0,TS="block")
opt,sol,data=cs_tssos_higher!(data,TS="block")
```
Options:  
nb: specify the first nb variables to be binary variables (satisfying xi^2=xi)  
CS (correlative sparsity): "MD" or "MF" (an approximately minimum chordal extension), "NC" (no chordal extension)  
TS: "block" (the maximal chordal extension), "MD" or "MF" (an approximately minimum chordal extension), false (without term sparsity)  
order: d (the relaxation order of Lasserre hierarchy), "multi" (applying the lowest relaxation orders for each variable clique)  
extra_sos: true (adding a first-order moment matrix for each variable clique), false  
solution: true (extract a solution), false (don't extract a solution)


## References
[1] [TSSOS: A Moment-SOS hierarchy that exploits term sparsity](https://arxiv.org/abs/1912.08899)  
[2] [Chordal-TSSOS: a moment-SOS hierarchy that exploits term sparsity with chordal extension](https://arxiv.org/abs/2003.03210)  
[3] [CS-TSSOS: Correlative and term sparsity for large-scale polynomial optimization](https://arXiv:2005.02828)

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): jwang@laas.fr
