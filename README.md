# BlockPOP
BlockPOP is a sparse polynomial optimization tool based on blocking Moment-SOS hierarchies. The Julia version of BlockPOP provides a usage based on the Julia language. To use the Matlab version of BlockPOP, one shoud use the *matlab* branch. Generally, the Julia version is more efficient.
## Dependencies
- Julia
- MATLAB
- MOSEK

Since the Julia version of BlockPOP calls MATLAB to handle polynomials, one needs to add the path of MATLAB to the environment variable PATH. The Julia version of BlockPOP has been tested on WINDOW 10, Julia 1.2.0, MATLAB R2016a and MOSEK 8.1.
## Usage
### Unconstrained polynomial optimization problems
The unconstrained polynomial optimization problem formulizes as
```
Inf{f(x): x\in R^n}
```
where f is a polynomial with variables x1,...,xn and of degree d.

Taking f=x1^4+x2^4-x1\*x2 as an example, to exetute the first blocking hierarchy, run
```Julia
julia> using BlockPOP
julia> n=2;d=4
# call MATLAB
ms=MSession()
mat"x = sym('x',[1 $n]);
poly=x(1)^4+x(2)^4-x(1)*x(2);
[coe, terms] = coeffs(poly,x);
lt=length(terms);
supp=zeros($n,lt);
for i=1:lt
    for j=1:$n
        supp(j,i)=feval(symengine,'degree',terms(i),x(j));
    end
end
coe=double(coe)"
coe=jarray(get_mvariable(ms,:coe))
supp=jarray(get_mvariable(ms,:supp))
supp=convert(Array{UInt8},supp)
julia> opt,data,status=blockupop_first(n,d,supp,coe)
```
By default, a monomial basis computed by the Newton polytope method will be used. If we set the key newton=0 in the input,
```Julia
julia> opt,data,status=blockupop_first(n,d,supp,coe,newton=0)
```
then the standard monomial basis will be used.

Two vectors will be outputed. The first vector is the size of blocks and the second vector is the number of blocks of size corresponding to the first vector.

In most cases, the first blocking hierarchy already obtains the same optimum as the dense Moment-SOS relaxation.

To exetute higher blocking hierarchies, repeatedly run

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

Taking f=x1^4+x2^4-x1\*x2 and g_1=1-x1^2-2\*x2^2 as an example, to exetute the first blocking hierarchy, run

```Julia
julia> n=2;m=1
d=2 # the order of Lasserre's hierarchy
dg=[2] # the degree vector of {g_j}
# call MATLAB
ms=MSession()
mat"x = sym('x',[1 $n]);
f=x(1)^4+x(2)^4-x(1)*x(2);
g_1=1-x(1)^2-2*x(2)^2;
pop=[f,g_1];
coe=cell(1,$m+1);
terms=cell(1,$m+1);
ssupp=cell(1,$m+1);
supp=[];
lt=zeros(1,$m+1);
for k=1:$m+1
    [coe{k}, terms{k}] = coeffs(pop(k),x);
    lt(k)=length(terms{k});
    ssupp{k}=zeros($n,lt(k));
    for i=1:lt(k)
        for j=1:$n
            ssupp{k}(j,i)=feval(symengine,'degree',terms{k}(i),x(j));
        end
    end
    supp=[supp ssupp{k}];
end
for k=1:$m+1
    coe{k}=double(coe{k});
end"
coe=jarray(get_mvariable(ms,:coe))
supp=jarray(get_mvariable(ms,:supp))
supp=convert(Array{UInt8},supp)
supp=unique(supp,dims=2)
ssupp=jarray(get_mvariable(ms,:ssupp))
for k=1:m+1
    ssupp[k]=convert(Array{UInt8},ssupp[k])
end
lt=jarray(get_mvariable(ms,:lt))
lt=convert(Array{UInt32},lt)
julia> opt,data,status=blockcpop_first(n,m,d,dg,supp,ssupp,coe,lt)
```

In most cases, the first blocking hierarchy already obtains the same optimum as the dense Moment-SOS relaxation.

To exetute higher blocking hierarchies, repeatedly run

```Julia
julia> opt,data,status=blockcpop_higher!(n,m,data)
```

## Contact
If there is any problems, you can contact: wangjie212@mails.ucas.ac.cn.
