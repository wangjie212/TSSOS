# TSSOS
TSSOS aims to provide a user-friendly and efficient tool for solving optimization problems with polynomials, which is based on the structured moment-SOS hierarchy. To use TSSOS in Julia, run
```Julia
pkg> add https://github.com/wangjie212/TSSOS
 ```

 | **Documentation** |
 |:-----------------:|
 | [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://wangjie212.github.io/TSSOS/dev) |

## Dependencies
- [Julia](https://julialang.org/)
- [JuMP](https://github.com/jump-dev/JuMP.jl)
- [Mosek](https://www.mosek.com/) or [COSMO](https://github.com/oxfordcontrol/COSMO.jl)
- [CliqueTrees](https://github.com/AlgebraicJulia/CliqueTrees.jl)

TSSOS has been tested on Ubuntu and Windows.

## Usage
### Unconstrained polynomial optimization
An unconstrained polynomial optimization problem could be formulized as

$$\mathrm{inf}_{\mathbf{x}\in\mathbb{R}^n}\ f(\mathbf{x}),$$

where $f\in\mathbb{R}[\mathbf{x}]$ is a polynomial.

Taking $f=1+x_1^4+x_2^4+x_3^4+x_1x_2x_3+x_2$ as an example, to compute the first TS step of the TSSOS hierarchy, run

```Julia
using TSSOS
using DynamicPolynomials
@polyvar x[1:3]
f = 1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2]*x[3] + x[2]
opt,sol,data = tssos_first(f, x, TS="MD")
```
By default, the monomial basis computed by the Newton polytope method is used. If one sets newton=false:

```Julia
opt,sol,data = tssos_first(f, x, newton=false, TS="MD")
```
then the standard monomial basis will be used.

To compute higher TS steps of the TSSOS hierarchy, repeatedly run

```Julia
opt,sol,data = tssos_higher!(data, TS="MD")
```

Options  
**nb**: specify the first nb variables to be $\pm1$ binary variables   
**TS**: "block" by default (maximal chordal extension), "signsymmetry" (sign symmetries), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations)  
**solution**: true (extract optimal solutions), false  

Output  
**basis**: monomial basis  
**blocksize**: sizes of blocks  
**blocks**: block structrue  
**GramMat**: Gram matrices (set Gram=true)  
**flag**: 0 if global optimality is certified; 1 otherwise  

### Constrained polynomial optimization
A constrained polynomial optimization problem could be formulized as

$$\mathrm{inf}_{\mathbf{x}\in\mathbf{K}}\ f(\mathbf{x}),$$

where $f\in\mathbb{R}[\mathbf{x}]$ is a polynomial and $\mathbf{K}$ is the basic semialgebraic set

$$\mathbf{K}\coloneqq\lbrace \mathbf{x}\in\mathbb{R}^n \mid g_i(\mathbf{x})\ge0, i=1,\ldots,m,\ h_j(\mathbf{x})=0, j=1,\ldots,\ell\rbrace,$$

for some polynomials $g_i,h_j\in\mathbb{R}[\mathbf{x}]$.

Taking $f=1+x_1^4+x_2^4+x_3^4+x_1x_2x_3+x_2$ and $\mathbf{K}\coloneqq\lbrace \mathbf{x}\in\mathbb{R}^2 \mid g(\mathbf{x})=1-x_1^2-2x_2^2\ge0, h(\mathbf{x})=x_2^2+x_3^2-1=0\rbrace$ as an example, to compute the first TS step of the TSSOS hierarchy, run

```Julia
@polyvar x[1:3]
f = 1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2]*x[3] + x[2]
g = 1 - x[1]^2 - 2*x[2]^2
h = x[2]^2 + x[3]^2 - 1
pop = [f, g, h]
d = 2 # set the relaxation order
opt,sol,data = tssos_first(pop, x, d, numeq=1, TS="block", solution=true)
```

To compute higher TS steps of the TSSOS hierarchy, repeatedly run

```Julia
opt,sol,data = tssos_higher!(data, TS="MD", solution=true)
```

Options  
**nb**: specify the first nb variables to be $\pm1$ binary variables  
**TS**: "block" by default (maximal chordal extension), "signsymmetry" (sign symmetries), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations)  
**normality**: true (impose normality condtions), false   
**GroebnerBasis**: true (work in the quotient ring by computing a Gröbner basis), false  
**solution**: true (extract optimal solutions), false  

One could also exploit correlative sparsity and term sparsity simultaneously.

```Julia
using DynamicPolynomials
n = 6
@polyvar x[1:n]
f = 1 + sum(x.^4) + x[1]*x[2]*x[3] + x[3]*x[4]*x[5] + x[3]*x[4]*x[6] + x[3]*x[5]*x[6] + x[4]*x[5]*x[6]
pop = [f, 1 - sum(x[1:3].^2), 1 - sum(x[3:6].^2)]
order = 2 # set the relaxation order
opt,sol,data = cs_tssos_first(pop, x, order, numeq=0, TS="MD", solution=true) # compute the first TS step of the CS-TSSOS hierarchy
opt,sol,data = cs_tssos_higher!(data, TS="block", solution=true) # compute higher TS steps of the CS-TSSOS hierarchy
```

Options  
**nb**: specify the first nb variables to be $\pm1$ binary variables  
**CS**: "MF" by default (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation)   
**TS**: "block" by default (maximal chordal extension), "signsymmetry" (sign symmetries), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations)   
**MomentOne**: true (add a first-order moment PSD constraint for each variable clique), false  
**solution**: true (extract an approximately optimal solution), false  

One may set solver="Mosek" or solver="COSMO" to specify the SDP solver invoked by TSSOS. By default, the solver is Mosek.

The parameters of COSMO could be tuned by

```Julia
settings = cosmo_para()
settings.eps_abs = 1e-5 # absolute residual tolerance
settings.eps_rel = 1e-5 # relative residual tolerance
settings.max_iter = 1e4 # maximum number of iterations
settings.time_limit = 1e4 # limit of running time
```
and run for instance tssos_first(..., cosmo_setting=settings)

The parameters of Mosek could be tuned by

```Julia
settings = mosek_para()
settings.tol_pfeas = 1e-8 # primal feasibility tolerance
settings.tol_dfeas = 1e-8 # dual feasibility tolerance
settings.tol_relgap = 1e-8 # relative primal-dual gap tolerance
settings.time_limit = 1e4 # limit of running time
settings.num_threads = 0 # number of threads available for Mosek
```
and run for instance tssos_first(..., mosek_setting=settings)

Output  
**basis**: monomial basis  
**blocksize**: sizes of blocks  
**blocks**: block structrue  
**GramMat**: Gram matrices (set Gram=true)  
**moment**: moment matrices  
**flag**: 0 if global optimality is certified; 1 otherwise  

### Exploiting symmetries
TSSOS also supports exploiting symmetries for polynomial optimization problems (thanks to [SymbolicWedderburn.jl](https://github.com/kalmarek/SymbolicWedderburn.jl)).

```Julia
using DynamicPolynomials
using PermutationGroups
@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
order = 2 # set the relaxation order
G = PermGroup([perm"(1,2,3)", perm"(1,2)"]) # define the permutation group acting on variables
opt,basis,Gram = tssos_symmetry(pop, x, order, G)
```

Options  
**numeq**: number of equality constraints  

Output  
**basis**: symmetry adapted basis
**Gram**: Gram matrices

## The AC-OPF problem
Check out `example/runopf.jl` and `example/modelopf.jl`.

## Sum-of-squares optimization
TSSOS supports more general [sum-of-squares optimization](https://en.wikipedia.org/wiki/Sum-of-squares_optimization) (including polynomial optimization as a special case):

$$\mathrm{inf}_{\mathbf{y}\in\mathbb{R}^n}\ \mathbf{c}^{\intercal}\mathbf{y}$$

$$\mathrm{s.t.}\ a_{k0}+y_1a_{k1}+\cdots+y_na_{kn}\in\mathrm{SOS},\ k=1,\ldots,m.$$

where $\mathbf{c}\in\mathbb{R}^n$ and $a_{ki}\in\mathbb{R}[\mathbf{x}]$ are polynomials. SOS constraints could be handled with the routine **add_psatz!**:

```Julia
info = add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order, TS="block", SO=1, GroebnerBasis=false)
```
where **nonneg** is a nonnegative polynomial constrained to admit a Putinar's style SOS representation on the semialgebraic set defined by **ineq_cons** and **eq_cons**, and **SO** is the sparse order.

The following is a simple exmaple.

$$\mathrm{sup}\ \lambda$$

$$\mathrm{s.t.}\ x_1^2 + x_1x_2 + x_2^2 + x_2x_3 + x_3^2 - \lambda(x_1^2+x_2^2+x_3^2)=\sigma+\tau_1(x_1^2+x_2^2+y_1^2-1)+\tau_2(x_2^2+x_3^2+y_2^2-1),$$

$$\sigma\in\mathrm{SOS},\deg(\sigma)\le2d,\ \tau_1,\tau_2\in\mathbb{R}[\mathbf{x}],\deg(\tau_1),\deg(\tau_2)\le2d-2.$$

```Julia
using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using TSSOS

@polyvar x[1:3]
f = x[1]^2 + x[1]*x[2] + x[2]^2 + x[2]*x[3] + x[3]^2
d = 2 # set the relaxation order
@polyvar y[1:2]
h = [x[1]^2 + x[2]^2 + y[1]^2-1, x[2]^2 + x[3]^2 + y[2]^2 - 1]
model = Model(optimizer_with_attributes(Mosek.Optimizer))
@variable(model, lower)
nonneg = f - lower*sum(x.^2)
info = add_psatz!(model, nonneg, [x; y], [], h, d, TS="block", GroebnerBasis=true)
@objective(model, Max, lower)
optimize!(model)
```
Check out `example/sosprogram.jl` for a more complicated example.

## Compute a locally optimal solution
Moment-SOS relaxations provide lower bounds on the optimum of the polynomial optimization problem. As the complementary side, one could compute a locally optimal solution which provides an upper bound on the optimum of the polynomial optimization problem. The upper bound is useful in evaluating the quality (tightness) of those lower bounds provided by moment-SOS relaxations. In TSSOS, for a given polynomial optimization problem, a locally optimal solution could be obtained via the nonlinear programming solver [Ipopt](https://github.com/jump-dev/Ipopt.jl):

```Julia
obj,sol,status = local_solution(data.n, data.m, data.supp, data.coe, numeq=data.numeq, startpoint=rand(data.n))
```

## Complex polynomial optimization
TSSOS also supports solving complex polynomial optimization. See [Exploiting Sparsity in Complex Polynomial Optimization](https://arxiv.org/abs/2103.12444) for more details.

A general complex polynomial optimization problem could be formulized as

$$\mathrm{inf}_{\mathbf{z}\in\mathbf{K}}\ f(\mathbf{z},\bar{\mathbf{z}}),$$

with

$$\mathbf{K}\coloneqq\lbrace \mathbf{z}\in\mathbb{C}^n \mid g_i(\mathbf{z},\bar{\mathbf{z}})\ge0, i=1,\ldots,m,\ h_j(\mathbf{z},\bar{\mathbf{z}})=0, j=1,\ldots,\ell\rbrace,$$

where $\bar{\mathbf{z}}$ stands for the conjugate of $\mathbf{z}:=(z_1,\ldots,z_n)$, and $f, g_i, i=1,\ldots,m, h_j, j=1,\ldots,\ell$ are real-valued complex polynomials satisfying $\bar{f}=f$ and $\bar{g}_j=g_j$.

Consider the following example:

$$\mathrm{inf}\ 3-|z_1|^2-0.5\mathbf{i}z_1\bar{z}_2^2+0.5\mathbf{i}z_2^2\bar{z}_1$$

$$\mathrm{s.t.}\ z_2+\bar{z}_2\ge0,|z_1|^2-0.25z_1^2-0.25\bar{z}_1^2=1,|z_1|^2+|z_2|^2=3,\mathbf{i}z_2-\mathbf{i}\bar{z}_2=0.$$

```Julia
using DynamicPolynomials
n = 2 # set the number of complex variables
@complex_polyvar z[1:n]
f = 3 - x[1]*conj(x[1]) - 0.5im*x[1]*conj(x[2])^2 + 0.5im*x[2]^2*conj(x[1])
g1 = x[2] + conj(x[2])
g2 = x[1]*conj(x[1]) - 0.25*x[1]^2 - 0.25*conj(x[1])^2 - 1
g3 = x[1]*conj(x[1]) + x[2]*conj(x[2]) - 3
g4 = im*x[2] - im*conj(x[2])
pop = [f, g1, g2, g3, g4]
order = 2 # set the relaxation order
opt,sol,data = complex_tssos_first(pop, z, order, numeq=3, TS="block", solution=true) # no correlative sparsity
opt,sol,data = complex_cs_tssos_first(pop, z, order, numeq=3, TS="block", solution=true)
```
Options  
**nb**: specify the first nb complex variables to be of unit norm (satisfying $|z_i|=1$)  
**CS**: "MF" by default (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation)   
**TS**: "block" by default (maximal chordal extension), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations)  
**ConjugateBasis**: include conjugate variables in monomial bases (false by default)  
**normality**: specify the normal order   
**MomentOne**: true (add a first-order moment PSD constraint for each variable clique), false  

## Sums of rational functions optimization
The sum-of-rational-functions optimization problem could be formulized as

$$\mathrm{inf}_{\mathbf{x}\in\mathbf{K}}\ \sum\_{i=1}^N\frac{p_i(\mathbf{x})}{q_i(\mathbf{x})},$$

where $p_i,q_i\in\mathbb{R}[\mathbf{x}]$ are polynomials and $\mathbf{K}$ is the basic semialgebraic set

$$\mathbf{K}\coloneqq\lbrace \mathbf{x}\in\mathbb{R}^n \mid g_i(\mathbf{x})\ge0, i=1,\ldots,m,\ h_j(\mathbf{x})=0, j=1,\ldots,\ell\rbrace,$$

for some polynomials $g_i,h_j\in\mathbb{R}[\mathbf{x}]$.

Taking $\frac{p_1}{q_1}=\frac{x^2+y^2-yz}{1+2x^2+y^2+z^2}$, $\frac{p_2}{q_2}=\frac{y^2+x^2z}{1+x^2+2y^2+z^2}$, $\frac{p_3}{q_3}=\frac{z^2-x+y}{1+x^2+y^2+2z^2}$, and $\mathbf{K}\coloneqq\lbrace \mathbf{x}\in\mathbb{R}^2 \mid g=1-x^2-y^2-z^2\ge0\rbrace$ as an example, run

```Julia
@polyvar x y z
p = [x^2 + y^2 - y*z, y^2 + x^2*z, z^2 - x + y] # define the vector of denominators
q = [1 + 2x^2 + y^2 + z^2, 1 + x^2 + 2y^2 + z^2, 1 + x^2 + y^2 + 2z^2] # define the vector of numerator
g = [1 - x^2 - y^2 - z^2]
d = 2 # set the relaxation order
opt = SumOfRatios(p, q, g, [], [x;y;z], d, QUIET=true, SignSymmetry=true) # Without correlative sparsity
opt = SparseSumOfRatios(p, q, g, [], [x;y;z], d, QUIET=true, SignSymmetry=true) # With correlative sparsity
```

Options  
**SignSymmetry**: true (exploit sign symmetries), false
**GroebnerBasis**: true (work in the quotient ring by computing a Gröbner basis), false

## Polynomial matrix optimization
The polynomial matrix optimization problem aims to minimize the smallest eigenvalue of a polynomial matrix subject to a tuple of polynomial matrix inequalties (PMIs), which could be formulized as

$$\mathrm{inf}_{\mathbf{x}\in\mathbf{K}}\ \lambda\_{\mathrm{min}}(F(\mathbf{x})),$$

where $F\in\mathbb{S}[\mathbf{x}]^p$ is a $p\times p$ symmetric polynomial matrix and $\mathbf{K}$ is the basic semialgebraic set

$$\mathbf{K}\coloneqq\lbrace \mathbf{x}\in\mathbb{R}^n \mid G_j(\mathbf{x})\succeq0, j=1,\ldots,m\rbrace,$$

for some symmetric polynomial matrices $G_j\in\mathbb{S}[\mathbf{x}]^{q_j}, j=1,\ldots,m$. Note that when $p=1$, $\lambda_{\min}(F(\mathbf{x}))=F(\mathbf{x})$. More generally, one may consider

$$\mathrm{inf}_{\mathbf{y}\in\mathbb{R}^t}\ \mathbf{c}^{\intercal}\mathbf{y}$$

$$\mathrm{s.t.}\ F_{0}(\mathbf{x})+y_1F_{1}(\mathbf{x})+\cdots+y_tF_{t}(\mathbf{x})\succeq0 \textrm{ on } K,$$

where $F_i\in\mathbb{S}[\mathbf{x}]^{p}, i=0,1,\ldots,t$ are a tuple of symmetric polynomial matrices.

The following is a simple exmaple.

```Julia
using DynamicPolynomials
using TSSOS

@polyvar x[1:5]
F = [x[1]^4 x[1]^2 - x[2]*x[3] x[3]^2 - x[4]*x[5] x[1]*x[4] x[1]*x[5];
x[1]^2 - x[2]*x[3] x[2]^4 x[2]^2 - x[3]*x[4] x[2]*x[4] x[2]*x[5];
x[3]^2 - x[4]*x[5] x[2]^2 - x[3]*x[4] x[3]^4 x[4]^2 - x[1]*x[2] x[5]^2 - x[3]*x[5];
x[1]*x[4] x[2]*x[4] x[4]^2 - x[1]*x[2] x[4]^4 x[4]^2 - x[1]*x[3];
x[1]*x[5] x[2]*x[5] x[5]^2 - x[3]*x[5] x[4]^2 - x[1]*x[3] x[5]^4]
G = Vector{Matrix{Polynomial{true, Int}}}(undef, 2)
G[1] = [1 - x[1]^2 - x[2]^2 x[2]*x[3]; x[2]*x[3] 1 - x[3]^2]
G[2] = [1 - x[4]^2 x[4]*x[5]; x[4]*x[5] 1 - x[5]^2]
@time opt,sol,data = tssos_first(F, G, x, 3, TS="MD") # compute the first TS step of the TSSOS hierarchy
@time opt,sol,data = tssos_higher!(data, TS="MD") # compute higher TS steps of the TSSOS hierarchy
```

Options  
**CS**: "MF" by default (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation)   
**TS**: "block" by default (maximal chordal extension), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations)  

For more examples, please check out `example/pmi.jl`.

## Tips for modelling polynomial optimization problems
- When possible, explictly include a sphere/ball constraint (or multi-sphere/multi-ball constraints).
- When the feasible set is unbounded, try the homogenization technique introduced in [Homogenization for polynomial optimization with unbounded sets](https://link.springer.com/article/10.1007/s10107-022-01878-5).
- Scale the coefficients of the polynomial optimization problem to $[-1, 1]$.
- Scale the variables so that they take values in $[-1, 1]$ or $[0, 1]$.
- Try to include more (redundant) inequality constraints.

## Christoffel-Darboux Kernels 
Given a measure $\mu$ supported on $\Omega \in \mathbb{R}^n$, and a degree $d \in \mathbb{N}^{*}$, one can define, whenever the moment matrix $M_d^{\mu}$ is invertible,  the Christoffel-Darboux kernel of order $d$:

$$
K_d^\mu: \mathbb{R}^n\times\mathbb{R}^n \to \mathbb{R}, \quad (x,y)\mapsto v_d(x)^\top (M_d^{\mu})^{-1} v_d(y),
$$

where $v_d$ is the standard monomial basis of order $d$.

Then, for any $\mathbf{x} \in \mathbf{R}^n$, the Christoffel polynomial $\Lambda_d^{\mu}$ is an SOS polynomial of degree $2d$ defined as $\Lambda_d^{\mu}(\mathbf{x}):=K_d^\mu(x,x)$.

In practice, when solving a particular POP instance via moment-SOS relaxations, we have access to a sequence of pseudo-moments $y$, to which we can associate the Christoffel polynomial $\Lambda_d^{y}$. 
The sublevel sets of $\Lambda_d^{y}$ provide an approximation of the support of a measure that is concentrated on the global minimizers of the original POP.

For more information on how to construct Christoffel polynomials and how to use them to strengthen the bounds from (correlatively-sparse) Moment-SOS relaxations, visit [CDK_Bound_Strengthening](https://github.com/SoDvc2226/CDK_Bound_Strengthening) or `docs/src/technique.md`.

## Non-commutative polynomial optimization
Visit [NCTSSOS](https://github.com/wangjie212/NCTSSOS)

## Joint spetral radii
Visit [SparseJSR](https://github.com/wangjie212/SparseJSR)

## MATLAB implementation
Visit [SPOT](https://github.com/ComputationalRobotics/SPOT/tree/main)

## References
[1] [TSSOS: A Moment-SOS hierarchy that exploits term sparsity](https://arxiv.org/abs/1912.08899)  
[2] [Chordal-TSSOS: a moment-SOS hierarchy that exploits term sparsity with chordal extension](https://arxiv.org/abs/2003.03210)  
[3] [CS-TSSOS: Correlative and term sparsity for large-scale polynomial optimization](https://arXiv:2005.02828)  
[4] [TSSOS: a Julia library to exploit sparsity for large-scale polynomial optimization](https://arxiv.org/abs/2103.00915)  
[5] [Sparse polynomial optimization: theory and practice](https://arxiv.org/abs/2208.11158)  
[6] [Strengthening Lasserre's Hierarchy in Real and Complex Polynomial Optimization](https://arxiv.org/abs/2404.07125)  
[7] [Exploiting Sign Symmetries in Minimizing Sums of Rational Functions](https://arxiv.org/abs/2405.09419)  
[8] [Leveraging Christoffel-Darboux Kernels to Strenghten Moment-SOS Relaxations](https://arxiv.org/pdf/2501.14281)

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn  

[Victor Magron](https://homepages.laas.fr/vmagron/): vmagron@laas.fr

[Srećko Ðurašinović](https://github.com/SoDvc2226): srecko001@e.ntu.edu.sg
