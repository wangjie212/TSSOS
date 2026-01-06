# Inclusion constants for free spectrahedra

## Four dichotomic qubit measurements

The following script allows to compute the bound for g=4 displayed in Section 7.2 from the paper [arXiv:2512.17706](https://arxiv.org/abs/2512.17706)

```Julia
using TSSOS
using DynamicPolynomials
using LinearAlgebra
using COSMO

g = 4
nj = 3
@polyvar a[1:g, 1:nj]
@polyvar x[1:g, 1:nj]
f = sum(a.*x)

epsg = [collect(c) for c in Iterators.product(fill([-1,1], g)...)]
vepsfun = r -> 1 - sum(j ->  (sum(epsg[r].*x[:,j]))^2, 1:nj)
ceps = map(vepsfun, 1:2^g)
va = i -> 1 - sum(a[i, :].^2)
ca = map(va, 1:g)
pop = [-f; ceps; ca]
r = 3

model = Model(optimizer_with_attributes(COSMO.Optimizer), "eps_abs"=>1e-6, "eps_rel"=>1e-6, "max_iter"=>1e7)
opt,sol,data = cs_tssos(pop, [a[:];x[:]], r, CS=false, TS="block", model=model)

model = Model(optimizer_with_attributes(COSMO.Optimizer), "eps_abs"=>1e-7, "eps_rel"=>1e-7, "max_iter"=>1e7)
opt,sol,data = cs_tssos(data, TS="block", model=model)
@show opt # 1.8029 
```
