# NCTSSOS
NCTSSOS is a sparse non-commutative polynomial optimization tool based on block moment-SOHS hierarchies. To use NCTSSOS in Julia, run
```Julia
pkg> add https://github.com/wangjie212/NCTSSOS
 ```

## Dependencies
- MultivariatePolynomials
- Julia
- MOSEK
- JuMP

NCTSSOS has been tested on WINDOW 10, Julia 1.2, JuMP 0.21 and MOSEK 8.1.
## Usage
### Unconstrained nc polynomial optimization problems
Taking f=1+x1^4+x2^4+x3^4+x1\*x2\*x3+x2 as an example, to exetute the first level of the NCTSSOS hierarchy, run
```Julia
using NCTSSOS
using DynamicPolynomials
@ncpolyvar x[1:3]
obj = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
opt,data = nctssos_first(obj,x,TS="MD",obj="eigen")
```

Two vectors will be outputed. The first vector is the sizes of blocks and the second vector is the numbers of blocks with sizes corresponding to the first vector.

To exetute higher levels of the NCTSSOS hierarchy, repeatedly run

```Julia
opt,data = nctssos_higher!(data,TS="MD")
```

Options:   
obj: "eigen" (implement the eigenvalue optimization), "trace" (implement the trace optimization)  
TS (term sparsity): "block" (using the maxmial chordal extension), "MD" or "MF" (using approximately minimum chordal extention), false (without term sparsity)  

### Constrained nc polynomial optimization problems
Taking f=2-x1^2+x1\*x2^2\*x1-x2^2 and g_1=4-x1^2-x2^2, g_2=x1\*x2+x2\*x1-2 as an example, to exetute the first level of the NCTSSOS hierarchy, run

```Julia
@ncpolyvar x[1:2]
obj=2-x[1]^2+x[1]*x[2]^2*x[1]-x[2]^2
ineq=[4-x[1]^2-x[2]^2, x[1]*x[2]+x[2]*x[1]-2]
eq=[]
pop=[obj; ineq; eq]
d=2 # the relaxation order
opt,data=nctssos_first(pop,x,d,TS="MD",obj="eigen")
```

To exetute higher levels of the NCTSSOS hierarchy, repeatedly run

```Julia
opt,data = nctssos_higher!(data,TS="MD")
```

Options:  
obj: "eigen" (implement the eigenvalue optimization), "trace" (implement the trace optimization)  
TS: "block" (using the maxmial chordal extension), "MD" or "MF" (using approximately minimum chordal extention), false (without term sparsity)  

To use combined correlative-term sparsity, run
```Julia
opt,data=cs_nctssos_first(pop,x,n,d,TS="block",obj="eigen")
```
and
```Julia
opt,data=cs_nctssos_higher!(data,TS="MD")
```

## References
[1] [Exploiting Term-Sparsity in Noncommutative Polynomial Optimization]

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): jwang@laas.fr
