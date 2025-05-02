# Sum-Of-Rational-Functions Optimization

The sum-of-rational-functions optimization problem could be formulized as

$$\mathrm{inf}_{\mathbf{x}\in\mathbf{K}}\ \sum_{i=1}^N\frac{p_i(\mathbf{x})}{q_i(\mathbf{x})},$$

where $p_i,q_i\in\mathbb{R}[\mathbf{x}]$ are polynomials and $\mathbf{K}$ is the basic semialgebraic set

$$\mathbf{K}\coloneqq\lbrace \mathbf{x}\in\mathbb{R}^n \mid g_i(\mathbf{x})\ge0, i=1,\ldots,m,\ h_j(\mathbf{x})=0, j=1,\ldots,\ell\rbrace,$$

for some polynomials $g_i,h_j\in\mathbb{R}[\mathbf{x}]$.

The following is a simple example.

```Julia
@polyvar x y z
p = [x^2 + y^2 - y*z, y^2 + x^2*z, z^2 - x + y] # define the vector of denominators
q = [1 + 2x^2 + y^2 + z^2, 1 + x^2 + 2y^2 + z^2, 1 + x^2 + y^2 + 2z^2] # define the vector of numerator
g = [1 - x^2 - y^2 - z^2]
d = 2 # set the relaxation order
opt = SumOfRatios(p, q, g, [], [x;y;z], d, QUIET=true, SignSymmetry=true) # Without correlative sparsity
opt = SparseSumOfRatios(p, q, g, [], [x;y;z], d, QUIET=true, SignSymmetry=true) # With correlative sparsity
```

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
SignSymmetry | Exploit sign symmetries | true
GroebnerBasis | Work in the quotient ring by computing a Gröbner basis | false

## Methods
```@docs
SumOfRatios
SparseSumOfRatios
```

### References

1. [Exploiting Sign Symmetries in Minimizing Sums of Rational Functions](https://arxiv.org/abs/2405.09419), Feng Guo, Jie Wang, and Jianhao Zheng, 2024.