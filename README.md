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
- [Mosek](https://www.mosek.com/) or [COSMO](https://github.com/oxfordcontrol/COSMO.jl)
- [ChordalGraph](https://github.com/wangjie212/ChordalGraph)

TSSOS has been tested on Ubuntu and Windows.
## Usage
### Unconstrained polynomial optimization
The unconstrained polynomial optimization problem formulizes as
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
By default, the monomial basis computed by the Newton polytope method is used. If one sets newton=false in the input,
```Julia
opt,sol,data = tssos_first(f, x, newton=false, TS="MD")
```
then the standard monomial basis will be used.

Two vectors will be output. The first vector includes the sizes of PSD blocks and the second vector includes the number of PSD blocks with sizes corresponding to the first vector.

To compute higher TS steps of the TSSOS hierarchy, repeatedly run

```Julia
opt,sol,data = tssos_higher!(data, TS="MD")
```

Options:  
**nb**: specify the first nb variables to be binary variables (satisfying $x_i^2=1$)  
**newton**: true (using the monomial basis computed by the Newton polytope method), false  
**TS (term sparsity)**: "block" (using the maximal chordal extension), "MD" (using approximately smallest chordal extensions), false (without term sparsity)  
**solution**: true (extracting an (approximate optimal) solution), false  

Output:  
**basis**: monomial basis  
**cl**: numbers of blocks  
**blocksize**: sizes of blocks  
**blocks**: block structrue  
**GramMat**: Gram matrices (you need to set Gram=true)  
**flag**: 0 if global optimality is certified; 1 otherwise  

### Constrained polynomial optimization
The constrained polynomial optimization problem formulizes as
$$\mathrm{inf}_{\mathbf{x}\in\mathbf{K}}\ f(\mathbf{x}),$$
where $f\in\mathbb{R}[\mathbf{x}]$ is a polynomial and $\mathbf{K}$ is the basic semialgebraic set
$$\mathbf{K}=\lbrace \mathbf{x}\in\mathbb{R}^n \mid g_j(\mathbf{x})\ge0, j=1,\ldots,m-numeq, \,g_j(\mathbf{x})=0, j=m-numeq+1,\ldots,m\rbrace,$$
for some polynomials $g_j\in\mathbb{R}[\mathbf{x}], j=1,\ldots,m$.

Taking $f=1+x_1^4+x_2^4+x_3^4+x_1x_2x_3+x_2$ and $\mathbf{K}=\lbrace \mathbf{x}\in\mathbb{R}^2 \mid g_1=1-x_1^2-2x_2^2\ge0, g_2=x_2^2+x_3^2-1=0\rbrace$ as an example, to compute the first TS step of the TSSOS hierarchy, run

```Julia
@polyvar x[1:3]
f = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
g_1 = 1-x[1]^2-2*x[2]^2
g_2 = x[2]^2+x[3]^2-1
pop = [f, g_1, g_2]
d = 2 # set the relaxation order
opt,sol,data = tssos_first(pop, x, d, numeq=1, TS="MD")
```

To compute higher TS steps of the TSSOS hierarchy, repeatedly run

```Julia
opt,sol,data = tssos_higher!(data, TS="MD")
```

Options:  
**nb**: specify the first nb variables to be binary variables (satisfying $x_i^2=1$)  
**TS**: "block" by default (using the maximal chordal extension), "MD" (using approximately smallest chordal extensions), false (without term sparsity)  
**quotient**: true (working in the quotient ring by computing Gröbner basis), false  
**solution**: true (extracting an (approximate optimal) solution), false  

One can also exploit correlative sparsity and term sparsity simultaneously, which is called the CS-TSSOS hierarchy.

```Julia
using DynamicPolynomials
n = 6
@polyvar x[1:n]
f = 1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
pop = [f, 1-sum(x[1:3].^2), 1-sum(x[1:4].^2)]
order = 2 # set the relaxation order
opt,sol,data = cs_tssos_first(pop, x, order, numeq=0, TS="MD")
opt,sol,data = cs_tssos_higher!(data, TS="MD")
```
Options:  
**nb**: specify the first nb variables to be binary variables (satisfying $x_i^2=1$)  
**CS (correlative sparsity)**: "MF" by default (generating an approximately smallest chordal extension), "NC" (without chordal extension), false (without correlative sparsity)   
**TS**: "block" (using the maximal chordal extension), "MD" (using approximately smallest chordal extensions), false (without term sparsity)  
**order**: d (relaxation order), "min" (using the lowest relaxation order for each variable clique)  
**MomentOne**: true (adding a first-order moment matrix for each variable clique), false  
**solution**: true (extracting an (approximate optimal) solution), false  

You may set solver="Mosek" or solver="COSMO" to specify the SDP solver invoked by TSSOS. By default, the solver is Mosek.

You can tune the parameters of COSMO via

```
settings = cosmo_para()
settings.eps_abs = 1e-5 # absolute residual tolerance
settings.eps_rel = 1e-5 # relative residual tolerance
settings.max_iter = 1e4 # maximum number of iterations
```
and run for instance tssos_first(..., cosmo_setting=settings)

Output:  
**basis**: monomial basis  
**cl**: numbers of blocks  
**blocksize**: sizes of blocks  
**blocks**: block structrue  
**GramMat**: Gram matrices (you need to set Gram=true)  
**moment**: moment matrices (you need to set Mommat=true)  
**flag**: 0 if global optimality is certified; 1 otherwise  

## The AC-OPF problem
See the file runopf.jl as well as modelopf.jl in example.

## Sum-of-squares optimization
TSSOS supports more general [sum-of-squares optimization](https://en.wikipedia.org/wiki/Sum-of-squares_optimization) (including polynomial optimization as a special case):
$$
\begin{equation}
\begin{cases}
\mathrm{inf}_{\mathbb{y}\in\mathbb{R}^n}&\mathbb{c}^{\intercal}\mathbb{y}\\
\mathrm{s.t.}&a_{k0}+y_1a_{k1}+\cdots+y_na_{kn}\in\mathrm{SOS}, \quad k=1,\ldots,m.
\end{cases}\notag
\end{equation}
$$
where $\mathbb{c}\in\mathbb{R}^n$ and $a_{ki}\in\mathbb{R}[\mathbf{x}]$ are polynomials. The SOS constraints can be handled with the routine **add_psatz!**:

```Julia
model,info = add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order, TS="block", SO=1, Groebnerbasis=false)
```
where **nonneg** is a nonnegative polynomial constrained to be a Putinar's style SOS on the semialgebraic set defined by **ineq_cons** and **eq_cons**, and **SO** is the sparse order.

The following is a simple exmaple.
$$
\begin{equation}
\begin{cases}
\mathrm{sup}&\lambda\\
\mathrm{s.t.}&x_1^2 + x_1x_2 + x_2^2 + x_2x_3 + x_3^2 - \lambda(x_1^2+x_2^2+x_3^2)=\sigma+\tau_1(x_1^2+x_2^2+y_1^2-1)+\tau_2(x_2^2+x_3^2+y_2^2-1),\\
&\sigma\in\mathrm{SOS},\deg(\sigma)\le2d,\,\tau_1,\tau_2\in\mathbb{R}[\mathbf{x}],\deg(\tau_1),\deg(\tau_2)\le2d-2.
\end{cases}\notag
\end{equation}
$$

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
h = [x[1]^2+x[2]^2+y[1]^2-1, x[2]^2+x[3]^2+y[2]^2-1]
model = Model(optimizer_with_attributes(Mosek.Optimizer))
@variable(model, lower)
nonne = f - lower*sum(x.^2)
model,info = add_psatz!(model, nonne, [x; y], [], h, d, TS="block", Groebnerbasis=true)
@objective(model, Max, lower)
optimize!(model)
```
See the file example/sosprogram.jl for a more complicated example.

## Complex polynomial optimization
TSSOS also supports solving complex polynomial optimization via the sparsity adapted complex moment-HSOS hierarchies. See [Exploiting Sparsity in Complex Polynomial Optimization](https://arxiv.org/abs/2103.12444) for more details.

The complex polynomial optimization problem formulizes as
$$\mathrm{inf}_{\mathbf{z}\in\mathbf{K}}\ f(\mathbf{z},\bar{\mathbf{z}}),$$
with
$$\mathbf{K}=\lbrace \mathbf{z}\in\mathbb{C}^n \mid g_j(\mathbf{z},\bar{\mathbf{z}})\ge0, j=1,\ldots,m-numeq, \,g_j(\mathbf{z},\bar{\mathbf{z}})=0, j=m-numeq+1,\ldots,m\rbrace,$$
where $\bar{\mathbf{z}}$ stands for the conjugate of $\mathbf{z}:=(z_1,\ldots,z_n)$, and $f, g_j, j=1,\ldots,m$ are real-valued polynomials satisfying $\bar{f}=f$ and $\bar{g}_j=g_j$.

In TSSOS, we use $x_i$ to represent the complex variable $z_i$ and use $x_{n+i}$ to represent its conjugate $\bar{z}_i$. Consider the example
$$
\begin{equation}
\begin{cases}
\mathrm{inf}&3-|z_1|^2-0.5\mathbf{i}z_1\bar{z}_2^2+0.5\mathbf{i}z_2^2\bar{z}_1\\
\mathrm{s.t.}&z_2+\bar{z}_2\ge0,\\
&|z_1|^2-0.25z_1^2-0.25\bar{z}_1^2=1,\\
&|z_1|^2+|z_2|^2=3,\\
&\mathbf{i}z_2-\mathbf{i}\bar{z}_2=0.
\end{cases}\notag
\end{equation}
$$
It can be represented as
$$
\begin{equation}
\begin{cases}
\mathrm{inf}&3-x_1x_3-0.5\mathbf{i}x_1x_4^2+0.5\mathbf{i}x_2^2x_3\\
\mathrm{s.t.}&x_2+x_4\ge0,\\
&x_1x_3-0.25x_1^2-0.25x_3^2=1,\\
&x_1x_3+x_2x_4=3,\\
&\mathbf{i}x_2-\mathbf{i}x_4=0.
\end{cases}\notag
\end{equation}
$$

```Julia
using DynamicPolynomials
n = 2 # set the number of complex variables
@polyvar x[1:2n]
f = 3 - x[1]*x[3] - 0.5im*x[1]*x[4]^2 + 0.5im*x[2]^2*x[3]
g1 = x[2] + x[4]
g2 = x[1]*x[3] - 0.25*x[1]^2 - 0.25 x[3]^2 - 1
g3 = x[1]*x[3] + x[2]*x[4] - 3
g4 = im*x[2] - im*x[4]
pop = [f, g1, g2, g3, g4]
order = 2 # set the relaxation order
opt,sol,data = cs_tssos_first(pop, x, n, order, numeq=3, TS="block")
```
Options:  
**nb**: specify the first nb complex variables to be of unit norm (satisfying $|z_i|=1$)  
**CS (correlative sparsity)**: "MF" by default (generating an approximately smallest chordal extension), "NC" (without chordal extension), false (without correlative sparsity)   
**TS**: "block" (using the maximal chordal extension), "MD" (using approximately smallest chordal extensions), false (without term sparsity)  
**order**: d (relaxation order), "min" (using the lowest relaxation order for each variable clique)  
**MomentOne**: true (adding a first-order moment matrix for each variable clique), false  
**ipart**: true (with complex moment matrices), false (with real moment matrices)

## Non-commutative polynomial optimization
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
[5] [Sparse polynomial optimization: theory and practice](https://arxiv.org/abs/2208.11158)

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn  
[Victor Magron](https://homepages.laas.fr/vmagron/): vmagron@laas.fr
