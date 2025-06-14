# Polynomial Matrix Optimization

The polynomial matrix optimization problem aims to minimize the smallest eigenvalue of a polynomial matrix subject to a tuple of polynomial matrix inequalties (PMIs), which could be formulized as

$$\mathrm{inf}_{\mathbf{x}\in\mathbf{K}}\ \lambda_{\mathrm{min}}(F(\mathbf{x})),$$

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

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
CS | Types of chordal extensions in exploiting correlative sparsity: "MF" (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation) | "MF"
cliques | Use customized variable cliques | []
TS | Types of chordal extensions used in term sparsity iterations: "block"(maximal chordal extension), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations) | "block"
QUIET | Silence the output| false
solve | Solve the SDP relaxation | true
Gram | Output Gram matrices | false
solution | Extract optimal solutions | false

### References

1. [Sparse Polynomial Matrix Optimization](https://arxiv.org/abs/2411.15479), Jared Miller, Jie Wang, and Feng Guo, 2024.