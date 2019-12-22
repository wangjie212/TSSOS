# TSSOS
TSSOS is a sparse polynomial optimization tool based on block Moment-SOS hierarchies. The Julia version of TSSOS provides a usage based on the Julia language. To use the Matlab version of TSSOS, one shoud use the *matlab* branch. Generally, the Julia version is more efficient.

To use the Julia version of TSSOS, run
```Julia
pkg> add https://github.com/wangjie212/TSSOS
 ```

## Dependencies
- Julia
- MOSEK 8.1

The Julia version of TSSOS has been tested on WINDOW 10, Julia 1.2.0, and MOSEK 8.1.
## Usage
### Unconstrained polynomial optimization problems
The unconstrained polynomial optimization problem formulizes as
```
Inf{f(x): x\in R^n}
```
where f is a polynomial with variables x1,...,xn and of degree d.

Taking f=x1^4+x2^4-x1\*x2 as an example, to exetute the first block hierarchy, run
```Julia
julia> using TSSOS
julia> using TypedPolynomials
julia> using MultivariatePolynomials
julia> n=2;d=4
@polyvar x[1:2]
f=x[1]^4+x[2]^4-x[1]*x[2]
mon=monomials(f)
coe=coefficients(f)
lm=length(mon)
supp=zeros(UInt8,n,lm)
for i=1:lm
    for j=1:n
        supp[j,i]=degree(mon[i],x[j])
    end
end
julia> opt,data,status=blockupop_first(n,d,supp,coe)
```
By default, a monomial basis computed by the Newton polytope method will be used. If we set the key newton=0 in the input,
```Julia
julia> opt,data,status=blockupop_first(n,d,supp,coe,newton=0)
```
then the standard monomial basis will be used.

Two vectors will be outputed. The first vector is the size of blocks and the second vector is the number of blocks of size corresponding to the first vector.

In most cases, the first block hierarchy already obtains the same optimum as the dense Moment-SOS relaxation.

To exetute higher block hierarchies, repeatedly run

```Julia
julia> opt,data,status=blockupop_higher!(n,data)
```

### Constrained polynomial optimization problems
The constrained polynomial optimization problem formulizes as
```
inf{f(x): x\in K}
```
where f is a polynomial and K is the basic semi-algebraic set
```
K={x\in R^n: g_j(x)>=0, j=1,...,m},
```
for some polynomials g_j, j=1,...,m.

Taking f=x1^4+x2^4-x1\*x2 and g_1=1-x1^2-2\*x2^2 as an example, to exetute the first block hierarchy, run

```Julia
julia> n=2;m=1
d=2 # the order of Lasserre's hierarchy
dg=[2] # the degree vector of {g_j}
@polyvar x[1:2]
f=x[1]^4+x[2]^4-x[1]*x[2]
g_1=1-x[1]^2-2*x[2]^2;
pop=[f,g_1]
coe=Array{Any}(undef, m+1)
mon=Array{Any}(undef, m+1)
ssupp=Array{Any}(undef, m+1)
lt=zeros(Int,1,m+1)
for k=1:m+1
    mon[k]=monomials(pop[k])
    coe[k]=coefficients(pop[k])
    lt[k]=length(mon[k])
    ssupp[k]=zeros(UInt8,n,lt[k])
    for i=1:lt[k]
        for j=1:n
            ssupp[k][j,i]=degree(mon[k][i],x[j])
        end
    end
end
supp=ssupp[1]
for i=2:m+1
    global supp=[supp ssupp[i]]
end
supp=unique(supp,dims=2)
julia> opt,data,status=blockcpop_first(n,m,d,dg,supp,ssupp,coe,lt)
```

In most cases, the first block hierarchy already obtains the same optimum as the dense Moment-SOS relaxation.

To exetute higher block hierarchies, repeatedly run

```Julia
julia> opt,data,status=blockcpop_higher!(n,m,data)
```

## Reference
For more details about TSSOS, please refer to [TSSOS: A Moment-SOS hierarchy that exploits term sparsity](https://arxiv.org/abs/1912.08899). If there are any problems, you can contact Jie Wang: wangjie212@mails.ucas.ac.cn.
