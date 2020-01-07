# TSSOS
TSSOS is a sparse polynomial optimization tool based on block Moment-SOS hierarchies. The Julia version of TSSOS provides a usage based on the Julia language. To use the Matlab version of TSSOS, one shoud use the *matlab* branch. We recommend the Julia version. Generally, the Julia version is more efficient.

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
julia> @polyvar x[1:2]
f=x[1]^4+x[2]^4-x[1]*x[2]
julia> opt,data=blockupop_first(f,x)
```
By default, a monomial basis computed by the Newton polytope method will be used. If we set the key newton=0 in the input,
```Julia
julia> opt,data=blockupop_first(f,x,newton=0)
```
then the standard monomial basis will be used.

Two vectors will be outputed. The first vector is the size of blocks and the second vector is the number of blocks of size corresponding to the first vector.

You can use the option method='clique' which may be more efficient.

In most cases, the first block hierarchy already obtains the same optimum as the dense Moment-SOS relaxation.

To exetute higher block hierarchies, repeatedly run

```Julia
julia> opt,data=blockupop_higher!(data)
```

### Constrained polynomial optimization problems
The constrained polynomial optimization problem formulizes as
```
Inf{f(x): x\in K}
```
where f is a polynomial and K is the basic semi-algebraic set
```
K={x\in R^n: g_j(x)>=0, j=1,...,m, h_k(x)=0, k=1,...,l},
```
for some polynomials g_j, j=1,...,m and h_k, k=1,...,l.

Taking f=x1^4+x2^4-x1\*x2 and g_1=2-x1^2-2\*x2^2, h_1=x1^2-x2^2-1 as an example, to exetute the first block hierarchy, run

```Julia
julia> @polyvar x[1:2]
f=x[1]^4+x[2]^4-x[1]*x[2]
g_1=2-x[1]^2-2*x[2]^2
h_1=x[1]^2-x[2]^2-1
pop=[f,g_1,h_1]
d=2 # the order of Lasserre's hierarchy
julia> opt,data=blockcpop_first(pop,x,d,numeq=1)
```

You can also use the option method='clique' which may be more efficient.

In most cases, the first block hierarchy already obtains the same optimum as the dense Moment-SOS relaxation.

To exetute higher block hierarchies, repeatedly run

```Julia
julia> opt,data=blockcpop_higher!(data,numeq=1)
```

## Reference
For more details about TSSOS, please refer to [TSSOS: A Moment-SOS hierarchy that exploits term sparsity](https://arxiv.org/abs/1912.08899). If there are any problems, you can contact Jie Wang: wangjie212@mails.ucas.ac.cn.
