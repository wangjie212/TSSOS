# TSSOS
TSSOS is a polynomial optimization tool based on the sparsity adapted moment-SOS hierarchies. To use TSSOS in Julia, run
```Julia
pkg> add https://github.com/wangjie212/TSSOS
 ```

 | **Documentation** |
 |:-----------------:|
 | [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://wangjie212.github.io/TSSOS/dev) |

## Dependencies
- [Julia](https://julialang.org/)
- [JuMP](https://github.com/jump-dev/JuMP.jl)
- [MOSEK](https://www.mosek.com/) or [COSMO](https://github.com/oxfordcontrol/COSMO.jl)

TSSOS has been tested on Ubuntu and Windows.
## Usage
### Unconstrained polynomial optimization problems
The unconstrained polynomial optimization problem formulizes as
```
Inf{f(x): x∈R^n}
```
where f is a polynomial with variables x1,...,xn and of degree d.

Taking f=1+x1^4+x2^4+x3^4+x1\*x2\*x3+x2 as an example, to execute the first level of the TSSOS hierarchy, run
```Julia
using TSSOS
using DynamicPolynomials
@polyvar x[1:3]
f = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
opt,sol,data = tssos_first(f, x, TS="MD")
```
By default, the monomial basis computed by the Newton polytope method is used. If one sets newton=false in the input,
```Julia
opt,sol,data = tssos_first(f, x, newton=false, TS="MD")
```
then the standard monomial basis will be used.

Two vectors will be output. The first vector includes the sizes of PSD blocks and the second vector includes the number of PSD blocks with sizes corresponding to the first vector.

To execute higher levels of the TSSOS hierarchy, repeatedly run

```Julia
opt,sol,data = tssos_higher!(data, TS="MD")
```

Options:  
nb: specify the first nb variables to be binary variables (satisfying xi^2=1)  
newton: true (use the monomial basis computed by the Newton polytope method), false  
TS (term sparsity): "block" (using the maximal chordal extension), "MD" (using approximately smallest chordal extensions), false (without term sparsity)  
solution: true (extract an (approximate optimal) solution), false (don't extract an (approximate optimal) solution)

### Constrained polynomial optimization problems
The constrained polynomial optimization problem formulizes as
```
Inf{f(x): x∈K}
```
where f is a polynomial and K is the basic semi-algebraic set
```
K={x∈R^n: g_j(x)>=0, j=1,...,m-numeq, g_j(x)=0, j=m-numeq+1,...,m},
```
for some polynomials g_j, j=1,...,m.

Taking f=1+x1^4+x2^4+x3^4+x1\*x2\*x3+x2 and K={x∈R^2: g_1=1-x1^2-2\*x2^2>=0, g_2=x2^2+x3^2-1=0} as an example, to execute the first level of the TSSOS hierarchy, run

```Julia
@polyvar x[1:3]
f = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
g_1 = 1-x[1]^2-2*x[2]^2
g_2 = x[2]^2+x[3]^2-1
pop = [f, g_1, g_2]
d = 2 # the relaxation order
opt,sol,data = tssos_first(pop, x, d, numeq=1, TS="MD")
```

To execute higher levels of the TSSOS hierarchy, repeatedly run

```Julia
opt,sol,data = tssos_higher!(data, TS="MD")
```

Options:  
nb: specify the first nb variables to be binary variables (satisfying xi^2=1)  
TS: "block" by default (using the maximal chordal extension), "MD" (using approximately smallest chordal extensions), false (without term sparsity)  
solution: true (extract an (approximate optimal) solution), false (don't extract an (approximate optimal) solution)

One can also exploit correlative sparsity and term sparsity simultaneously, which is called the CS-TSSOS hierarchy.

```Julia
using DynamicPolynomials
n = 6
@polyvar x[1:n]
f = 1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
pop = [f, 1-sum(x[1:3].^2), 1-sum(x[1:4].^2)]
order = 2 # the relaxation order
opt,sol,data = cs_tssos_first(pop, x, order, numeq=0, TS="MD")
opt,sol,data = cs_tssos_higher!(data, TS="MD")
```
Options:  
nb: specify the first nb variables to be binary variables (satisfying xi^2=1)  
CS (correlative sparsity): "MF" by default (generating an approximately smallest chordal extension), "NC" (without chordal extension), false (without correlative sparsity)   
TS: "block" (using the maximal chordal extension), "MD" (using approximately smallest chordal extensions), false (without term sparsity)  
order: d (the relaxation order), "min" (using the lowest relaxation order for each variable clique)  
MomentOne: true (adding a first-order moment matrix for each variable clique), false  
solution: true (extract an (approximate optimal) solution), false (don't extract an (approximate optimal) solution)

## The AC-OPF problem
See the file runopf.jl in examples.

## Complex polynomial optimization problems
TSSOS also supports solving complex polynomial optimization via sparsity adapted complex moment-SOHS hierarchy. See [Exploiting Sparsity in Complex Polynomial Optimization](https://arxiv.org/abs/2103.12444) for more details.

The complex polynomial optimization problem formulizes as
```
Inf{f(z,z'): z∈K}
```
with
```
K={z∈C^n: g_j(z,z')>=0, j=1,...,m-numeq, g_j(z,z')=0, j=m-numeq+1,...,m},
```
where ' stands for conjugate and f, g_j, j=1,...,m are real-valued polynomials satisfying f'=f and g_j'=g_j.

We use Vector{UInt16}[[1;2], [2;3]] to represent z_1z_2z_2'z_3'. Consider the example Inf{3-|z_1|^2-0.5iz_1z_2'^2+0.5iz_2^2z_1': z_2+z_2'>=0, |z_1|^2-0.25z_1^2-0.25z_1'^2=1, |z_1|^2+|z_2|^2=3, iz_2-iz_2'=0}.

```Julia
n = 2 # the number of complex variables
supp = Vector{Vector{Vector{UInt16}}}[[[[], []], [[1], [1]], [[1], [2;2]], [[2;2], [1]]],
[[[2], []], [[], [2]]], [[[], []], [[1], [1]], [[1;1], []], [[], [1;1]]],
[[[], []], [[1], [1]], [[2], [2]]], [[[2], []], [[], [2]]]]
coe = [[3;-1;-0.5im;0.5im], [1;1], [-1;1;-0.25;-0.25], [-3;1;1], [im;-im]]
order = 2 # the relaxation order
opt,sol,data = cs_tssos_first(supp, coe, n, order, numeq=3, TS="block")
```
Options are as above except that "solution" is not provided.

## Non-commutative polynomial optimization problems
Visit [NCTSSOS](https://github.com/wangjie212/NCTSSOS)

## Analysis of sparse dynamical systems
Visit [SparseDynamicSystem](https://github.com/wangjie212/SparseDynamicSystem)

## Joint spetral radii
Visit [SparseJSR](https://github.com/wangjie212/SparseJSR)

## References
[1] [TSSOS: A Moment-SOS hierarchy that exploits term sparsity](https://arxiv.org/abs/1912.08899)  
[2] [Chordal-TSSOS: a moment-SOS hierarchy that exploits term sparsity with chordal extension](https://arxiv.org/abs/2003.03210)  
[3] [CS-TSSOS: Correlative and term sparsity for large-scale polynomial optimization](https://arXiv:2005.02828)  
[4] [TSSOS: a Julia library to exploit sparsity for large-scale polynomial optimization](https://arxiv.org/abs/2103.00915)

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): jwang@laas.fr  
[Victor Magron](https://homepages.laas.fr/vmagron/): vmagron@laas.fr
