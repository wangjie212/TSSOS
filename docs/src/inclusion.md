# Inclusion constants for free spectrahedra

## Four dichotomic qubit measurements

The following script allows to compute the bound for g=4 displayed in Section 7.2 from the paper [arXiv:2512.](https://arxiv.org/abs/2512.)

```Julia
using TSSOS
using DynamicPolynomials
using LinearAlgebra
using COSMO
settings = cosmo_para()
settings.eps_abs = 1e-6 # absolute residual tolerance
settings.eps_rel = 1e-6 # relative residual tolerance
settings.max_iter = 1e7 # maximum number of iterations

g=4
nj=3
@polyvar a[1:g,1:nj]
@polyvar x[1:g,1:nj]
f = sum(a.*x);

epsg = [collect(c) for c in Iterators.product(fill([-1,1],g)...)]
vepsfun = r -> 1 - sum(j ->  (sum(epsg[r].*x[:,j]))^2,1:nj);
ceps = map(vepsfun, 1:2^g);
va = i -> 1 - sum(a[i,:].^2);
ca = map(va,1:g)
pop = [-f;ceps;ca];
r = 3;

opt,sol,data = cs_tssos_first(pop, [a[:];x[:]], r, CS=false, TS="block", solver="COSMO", cosmo_setting=settings);

settings.eps_abs = 1e-7 # absolute residual tolerance
settings.eps_rel = 1e-7 # relative residual tolerance
settings.max_iter = 1e7 # maximum number of iterations

opt,sol,data = cs_tssos_higher!(data, TS="block", cosmo_setting=settings);
opt # 1.8029 
```
