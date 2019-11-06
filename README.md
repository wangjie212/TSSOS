# BlockPOP
***
BlockPOP is a sparse polynomial optimization tool based on blocking Moment-SOS hierarchies. The Matlab version of BlockPOP provides a usage based on MATLAB. Note that this Matlab version may be inefficient for large problems since it is written using high-level YALMIP code. For the efficency purpose, one should use the Julia version of BlockPOP.
***
## Dependencies
- MATLAB
- YALMIP
- a SDP solver like MOSEK

The Matlab version of BlockPOP has been tested on MATLAB R2016a.
***
## Usage
### Unconstrained polynomial optimization problems
$$
\inf_{\mathbf{x}}\{f(\mathbf{x}): \x\in\R^n\}
$$
where f is a polynomial with variables x1,...,xn and of degree d.

Taking f=x1^4+x2^4-x1\*x2 as an example, to exetute the first blocking hierarchy, run

```matlab
>>syms x1,x2;
>>f=x1^4+x2^4-x1*x2;
>>n=2;
>>d=4;
>>[opt,data,status]=blockpop_uncons_first(f,n,d,'mosek')
```

By default, a monomial basis computed by the Newton polytope method will be used. If we add 'newton',0 to the input, then the standard monomial basis will be used.

To exetute higher blocking hierarchies, repeatedly run

```matlab
>>[opt,data,status]=blockpop_uncons_higher(n,data,'mosek')
```

### Constrained polynomial optimization problems
$$
\inf_{\mathbf{x}}\{f(\mathbf{x}) : \mathbf{x}\in\mathbf{K}\}
$$
where f is a polynomial and $\mathbf{K}$ is the basic semi-algebraic set
$$
\mathbf{K} = \{\x\in\R^{n} : g_j(\x)\ge 0, j = 1,\ldots,m\},
$$
for some polynomials $g_j, j = 1,\ldots,m.$

Taking f=x1^4+x2^4-x1\*x2 and g=1-x1^2-2\*x2^2 as an example, to exetute the first blocking hierarchy, run

```matlab
>>syms x1,x2;
>>f=x1^4+x2^4-x1*x2;
>>g=1-x1^2-2*x2^2;
>>n=2;
>>m=1;
>>d=2; % the order of Lasserre's hierarchy
>>dg=[2]；
>>[opt,data,status]=blockpop_cons_first(f,g,n,m,d,dg,'mosek')
```


To exetute higher blocking hierarchies, repeatedly run

```matlab
>>[opt,data,status]=blockpop_cons_higher(n,m,data,'mosek')
```
