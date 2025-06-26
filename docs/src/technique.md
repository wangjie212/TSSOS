# Techniques

## Homogenization


## Christoffel-Darboux Kernels

Here is the list of functions implemented in TSSOS that are based on computing Christoffel-Darboux kernels associated to a given POP:

```@docs
run_H1
run_H1CS
run_H2
run_H2CS
construct_CDK
construct_marginal_CDK
construct_CDK_cs
construct_marginal_CDK_cs
```

We provide below some illustrations. 
Let us consider the following POP (Example 3.2 from [[Sparse polynomial optimization: theory and practice]](https://arxiv.org/abs/2208.11158)):
```julia
n = 6
@polyvar x[1:n]
f = x[2]*x[5] + x[3]*x[6] - x[2]*x[3] - x[5]*x[6] + x[1]*(-x[1] + x[2] + x[3] - x[4] + x[5] + x[6])
g = [(6.36 - x[i]) * (x[i] - 4) for i in 1:6]
pop = [f, g...]
d = 1
```
We can solve the problem using the dense moment-SOS hierarchy:
```julia
opt, sol, data = cs_tssos_first(pop, x, d, TS=false, CS=false, solution=true)
```
Afterwards, one can try strenghtening the bound via **H1**:
```julia
N = 5
eps = 0.05
dc = 1
gap_tol = 0.1
resultH1 = run_H1(pop,x,d,dc,N,eps,gap_tol)
```
or via **H2**
```julia
dc = 1
local_sol = sol
tau = 1.1
resultH2 = run_H2(pop,x,d,dc,local_sol,tau)
```
Alternatively, the problem can be solved using the correlatively sparse moment-SOS hierarchy:
```julia
opt, sol, data = cs_tssos_first(pop, x, d, TS=false, solution=true)
```
Afterwards, the bound can be strengthened either via **H1CS**:
```julia
N = 5
eps = 0.05
dc = 1
gap_tol = 0.1
resultH1CS = run_H1CS(pop,x,d,dc,N,eps,gap_tol)
```
or via **H2CS**
```julia
dc = 1
local_sol = sol
tau = 1.1
resultH2CS = run_H2CS(pop,x,d,dc,local_sol,tau)
```

Moreover, here is how different Christoffel polynomials can be constructed using the output of the: 
- dense Moment-SOS relaxation of order 2
```julia
d = 2
opt, sol, data = cs_tssos_first(pop, x, d, TS=false, CS=false, solution=true, Mommat=true)

k = 4
dc = 1
CDK_order1 = construct_CDK(x, dc, data.moment[1])  # Constructs multivariate Christoffel polynomial of order dc=1 (quadratic CDK)
CDK_4_order1 = construct_marginal_CDK(x, k, dc, data.moment[1])  # Constructs marginal Christoffel polynomial, associated to x_4, of order dc=1 

dc = 2
CDK_order2 = construct_CDK(x, dc, data.moment[1])  # Construct multivariate Christoffel polynomial of order dc=2 (quartic CDK)
CDK_4_order2 = construct_marginal_CDK(x, k, dc, data.moment[1])  # Constructs marginal Christoffel polynomial, associated to x_4, of order dc=2 

```
- sparse Moment-SOS relaxation of order 2
```julia
d = 2
opt, sol, data = cs_tssos_first(pop, x, d, TS=false, solution=true, Mommat=true)

k = 4
dc = 1
CDK_sparse_order1 = construct_CDK_cs(x, dc, data.moment, data.cliques)  # Constructs multivariate Christoffel polynomial of order dc=1 for each identified clique.
CDK_sparse_4_order1 = construct_marginal_CDK_cs(x, k, dc, data.moment, data.cliques)  # Constructs marginal Christoffel polynomial, associated to x_4, of order dc=1 

dc = 2
CDK_sparse_order2 = construct_CDK_cs(x, dc, data.moment, data.cliques)  # Constructs multivariate Christoffel polynomial of order dc=2 for each identified clique.
CDK_sparse_4_order2 = construct_marginal_CDK_cs(x, k, dc, data.moment, data.cliques)  # Constructs marginal Christoffel polynomial, associated to x_4, of order dc=2

```

## Tips for modelling polynomial optimization problems
- When possible, explictly include a sphere/ball constraint (or multi-sphere/multi-ball constraints).
- When the feasible set is unbounded, try the homogenization technique introduced in [Homogenization for polynomial optimization with unbounded sets](https://link.springer.com/article/10.1007/s10107-022-01878-5).
- Scale the coefficients of the polynomial optimization problem to $[-1, 1]$.
- Scale the variables so that they take values in $[-1, 1]$ or $[0, 1]$.
- Try to include more (redundant) inequality constraints.
