# Polynomial Optimization

Polynomial optimization concerns minimizing a polynomial subject to a tuple of polynomial inequality constraints and equality constraints, which in general takes the form:

$$\mathrm{inf}_{\mathbf{x}\in\mathbb{R}^n}\ f(\mathbf{x})\ \text{ s.t. }\ g_1(\mathbf{x})\ge0,\ldots,g_m(\mathbf{x})\ge0,h_1(\mathbf{x})=0,\ldots,h_{\ell}(\mathbf{x})=0,$$

where $f,g_1,\ldots,g_m,h_1,\ldots,h_{\ell}\in\mathbb{R}[\mathbf{x}]$ are polynomials in variables $\mathbf{x}$.

To illustrate how to solve a polynomial optimization problem with TSSOS, let us consider $f=1+x_1^4+x_2^4+x_3^4+x_1x_2x_3+x_2$ and $g=1-x_1^2-2x_2^2, h=x_2^2+x_3^2-1$.

```Julia
@polyvar x[1:3]
f = 1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2]*x[3] + x[2]
g = 1 - x[1]^2 - 2*x[2]^2
h = x[2]^2 + x[3]^2 - 1
pop = [f, g, h]
d = 2 # set the relaxation order
opt,sol,data = tssos_first(pop, x, d, numeq=1, TS="block", solution=true) # compute the first TS step of the TSSOS hierarchy
opt,sol,data = tssos_higher!(data, TS="block", solution=true) # compute higher TS steps of the TSSOS hierarchy
```

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
nb | Specify the first **nb** variables to be $\pm1$ binary variables | 0
numeq | Specify the last **numeq** constraints to be equality constraints | 0
GroebnerBasis | Work in the quotient ring by computing a Gröbner basis | true
basis | Use customized monomial bases | []
reducebasis | Reduce the monomial bases | false
TS | Types of chordal extensions used in term sparsity iterations: "block"(maximal chordal extension), "signsymmetry" (sign symmetries), "MD" (approximately smallest chordal extension), false (invalidating term sparsity exploitation) | "block"
normality | Impose normality condtions | false
merge | Merge overlapping PSD blocks | false
md | Parameter for tunning the merging strength | 3
MomentOne | add a first-order moment PSD constraint | false
solver | Specify an SDP solver: "Mosek" or "COSMO" | "Mosek"
cosmo\_setting | Parameters for the COSMO solver: cosmo\_para(eps\_abs, eps\_rel, max\_iter, time\_limit) | cosmo\_para(1e-5, 1e-5, 1e4, 0)
mosek\_setting | Parameters for the Mosek solver: mosek\_para(tol\_pfeas, tol\_dfeas, tol\_relgap, time\_limit, num\_threads) | mosek\_para(1e-8, 1e-8, 1e-8, -1, 0)
QUIET | Silence the output| false
solve | Solve the SDP relaxation | true
dualize | Solve the dual SDP problem | false
Gram | Output Gram matrices | false
solution | Extract optimal solutions | false
rtol | tolerance for rank | 1e-2
gtol | tolerance for global optimality gap | 1e-2
ftol | tolerance for feasibility | 1e-3

## Correlative sparsity
The following is an example where one exploits correlative sparsity and term sparsity simultaneously.

```Julia
using DynamicPolynomials
n = 6
@polyvar x[1:n]
f = 1 + sum(x.^4) + x[1]*x[2]*x[3] + x[3]*x[4]*x[5] + x[3]*x[4]*x[6]+x[3]*x[5]*x[6] + x[4]*x[5]*x[6]
pop = [f, 1-sum(x[1:3].^2), 1-sum(x[3:6].^2)]
order = 2 # set the relaxation order
opt,sol,data = cs_tssos_first(pop, x, order, numeq=0, TS="MD", solution=true) # compute the first TS step of the CS-TSSOS hierarchy
opt,sol,data = cs_tssos_higher!(data, TS="MD", solution=true) # compute higher TS steps of the CS-TSSOS hierarchy
```

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
nb | Specify the first **nb** variables to be $\pm1$ binary variables | 0
numeq | Specify the last **numeq** constraints to be equality constraints | 0
CS | Types of chordal extensions in exploiting correlative sparsity: "MF" (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation) | "MF"
basis | Use customized monomial bases | []
hbasis | Use customized monomial bases associated with equality constraints | []
cliques | Use customized variable cliques | []
TS | Types of chordal extensions used in term sparsity iterations: "block"(maximal chordal extension), "signsymmetry" (sign symmetries), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations) | "block"
merge | Merge overlapping PSD blocks | false
md | Parameter for tunning the merging strength | 3
MomentOne | add a first-order moment PSD constraint for each variable clique | false
solver | Specify an SDP solver: "Mosek" or "COSMO" | "Mosek"
cosmo\_setting | Parameters for the COSMO solver: cosmo\_para(eps\_abs, eps\_rel, max\_iter, time\_limit) | cosmo\_para(1e-5, 1e-5, 1e4, 0)
mosek\_setting | Parameters for the Mosek solver: mosek\_para(tol\_pfeas, tol\_dfeas, tol\_relgap, time\_limit, num\_threads) | mosek\_para(1e-8, 1e-8, 1e-8, -1, 0)
QUIET | Silence the output| false
solve | Solve the SDP relaxation | true
dualize | Solve the dual SDP problem | false
Gram | Output Gram matrices | false
solution | Extract an optimal solution | false
rtol | tolerance for rank | 1e-2
gtol | tolerance for global optimality gap | 1e-2
ftol | tolerance for feasibility | 1e-3

## Compute a locally optimal solution
Moment-SOS relaxations provide lower bounds on the optimum of the polynomial optimization problem. As the complementary side, one could compute a locally optimal solution which provides an upper bound on the optimum of the polynomial optimization problem. The upper bound is useful in evaluating the quality (tightness) of those lower bounds provided by moment-SOS relaxations. In TSSOS, for a given polynomial optimization problem, a locally optimal solution could be obtained via the nonlinear programming solver [Ipopt](https://github.com/jump-dev/Ipopt.jl):

```Julia
obj,sol,status = local_solution(data.n, data.m, data.supp, data.coe, numeq=data.numeq, startpoint=rand(data.n))
```

## Methods
```@docs
tssos_first
tssos_higher!
cs_tssos_first
cs_tssos_higher!
local_solution
refine_sol
extract_solutions
extract_solutions_robust
```
