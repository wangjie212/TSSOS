# TSSOS
TSSOS is a sparse polynomial optimization tool based on block Moment-SOS hierarchies. The Matlab version of TSSOS provides a usage based on MATLAB. Note that this Matlab version may be inefficient for large problems since it is written using high-level YALMIP code. For the efficency purpose, one should use the Julia version of TSSOS (the *master* branch).
## Dependencies
- MATLAB
- YALMIP
- an SDP solver (recommend MOSEK)

The Matlab version of TSSOS has been tested on WINDOW 10, MATLAB R2016a and MOSEK 9.
## Usage
### Unconstrained polynomial optimization problems
The unconstrained polynomial optimization problem formulizes as
```
Inf{f(x): x\in R^n}
```
where f is a polynomial with variables x1,...,xn and of degree d.

Taking f=x1^4+x2^4-x1\*x2 as an example, to exetute the first block hierarchy, run

```matlab
>>syms x1,x2;
>>f=x1^4+x2^4-x1*x2;
>>n=2;
>>d=4;
>>[opt,data,status]=blockpop_uncons_first(f,n,d,'mosek')
```

By default, a monomial basis computed by the Newton polytope method will be used. If we add 'newton',0 to the input,

```matlab
>>[opt,data,status]=blockpop_uncons_first(f,n,d,'mosek','newton',0)
```

then the standard monomial basis will be used.

Two lines will be outputed. The first line is the size of blocks and the second line is the number of blocks of size corresponding to the first line.

In most cases, the first block hierarchy already obtains the same optimum as the dense Moment-SOS relaxation.

To exetute higher block hierarchies, repeatedly run

```matlab
>>[opt,data,status]=blockpop_uncons_higher(n,data,'mosek')
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

```matlab
>>syms x1,x2;
>>f=x1^4+x2^4-x1*x2;
>>g_1=1-x1^2-2*x2^2;
>>n=2;
>>m=1;
>>d=2; % the order of Lasserre's hierarchy
>>dg=[2]；% the degree vector of {g_j}
>>[opt,data,status]=blockpop_cons_first(f,g,n,m,d,dg,'mosek')
```

In most cases, the first block hierarchy already obtains the same optimum as the dense Moment-SOS relaxation.

To exetute higher block hierarchies, repeatedly run

```matlab
>>[opt,data,status]=blockpop_cons_higher(n,m,data,'mosek')
```

## Reference
For more details about TSSOS, please refer to [TSSOS: A Moment-SOS hierarchy that exploits term sparsity](https://arxiv.org/abs/1912.08899). If there are any problems, you can contact Jie Wang: wangjie212@mails.ucas.ac.cn.
