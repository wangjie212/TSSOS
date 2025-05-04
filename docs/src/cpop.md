# Complex Polynomial Optimization

A general complex polynomial optimization problem could be formulized as

$$\mathrm{inf}_{\mathbf{z}\in\mathbf{K}}\ f(\mathbf{z},\bar{\mathbf{z}}),$$

with

$$\mathbf{K}\coloneqq\lbrace \mathbf{z}\in\mathbb{C}^n \mid g_i(\mathbf{z},\bar{\mathbf{z}})\ge0, i=1,\ldots,m,\ h_j(\mathbf{z},\bar{\mathbf{z}})=0, j=1,\ldots,\ell\rbrace,$$

where $\bar{\mathbf{z}}$ stands for the conjugate of $\mathbf{z}:=(z_1,\ldots,z_n)$, and $f, g_i, i=1,\ldots,m, h_j, j=1,\ldots,\ell$ are real-valued complex polynomials satisfying $\bar{f}=f$ and $\bar{g}_i=g_i$, $\bar{h}_j=h_j$.

In TSSOS, we use $x_i$ to represent the complex variable $z_i$ and use $x_{n+i}$ to represent its conjugate $\bar{z}_i$. Let us consider the following example:

$$\mathrm{inf}\ 3-|z_1|^2-0.5\mathbf{i}z_1\bar{z}_2^2+0.5\mathbf{i}z_2^2\bar{z}_1$$

$$\mathrm{s.t.}\ z_2+\bar{z}_2\ge0,|z_1|^2-0.25z_1^2-0.25\bar{z}_1^2=1,|z_1|^2+|z_2|^2=3,\mathbf{i}z_2-\mathbf{i}\bar{z}_2=0,$$

which is represented in TSSOS as

$$\mathrm{inf}\ 3-x_1x_3-0.5\mathbf{i}x_1x_4^2+0.5\mathbf{i}x_2^2x_3$$

$$\mathrm{s.t.}\ x_2+x_4\ge0,x_1x_3-0.25x_1^2-0.25x_3^2=1,x_1x_3+x_2x_4=3,\mathbf{i}x_2-\mathbf{i}x_4=0.$$

```Julia
using DynamicPolynomials
n = 2 # set the number of complex variables
@polyvar x[1:2n]
f = 3 - x[1]*x[3] - 0.5im*x[1]*x[4]^2 + 0.5im*x[2]^2*x[3]
g1 = x[2] + x[4]
h1 = x[1]*x[3] - 0.25*x[1]^2 - 0.25 x[3]^2 - 1
h2 = x[1]*x[3] + x[2]*x[4] - 3
h3 = im*x[2] - im*x[4]
pop = [f, g, h1, h2, h3]
order = 2 # set the relaxation order
opt,sol,data = cs_tssos_first(pop, x, n, order, numeq=3) # compute the first TS step of the CS-TSSOS hierarchy
opt,sol,data = cs_tssos_higher!(data) # compute higher TS steps of the CS-TSSOS hierarchy
```

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
nb | Specify the first **nb** complex variables to be of unit norm | 0
numeq | Specify the last **numeq** constraints to be equality constraints | 0
CS | Types of chordal extensions in exploiting correlative sparsity: "MF" (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation) | "MF"
cliques | Use customized variable cliques | []
TS | Types of chordal extensions used in term sparsity iterations: "block"(maximal chordal extension), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations) | "block"
ConjugateBasis | include conjugate variables in monomial bases | false
normality | Impose normality condtions of order **normality** | 1
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
solution | Extract an approximately optimal solution | false
tol | Tolerance for certifying global optimality | 1e-4

### References

1. [Exploiting Sparsity in Complex Polynomial Optimization](https://link.springer.com/article/10.1007/s10957-021-01975-z), Jie Wang and Victor Magron, 2021.