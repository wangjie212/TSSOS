# Sum-Of-Squares Optimization

A general [sum-of-squares optimization](https://en.wikipedia.org/wiki/Sum-of-squares_optimization) (including polynomial optimization as a special case) problem takes the form:

$$\mathrm{inf}_{\mathbf{y}\in\mathbb{R}^n}\ \mathbf{c}^{\intercal}\mathbf{y}$$

$$\mathrm{s.t.}\ a_{k0}+y_1a_{k1}+\cdots+y_na_{kn}\in\mathrm{SOS},\ k=1,\ldots,m.$$

where $\mathbf{c}\in\mathbb{R}^n$ and $a_{ki}\in\mathbb{R}[\mathbf{x}]$ are polynomials. In TSSOS, SOS constraints could be handled with the routine **add_psatz!**:

```Julia
info = add_psatz!(model, nonneg, vars, ineq_cons, eq_cons, order, TS="block", SO=1, Groebnerbasis=false)
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
h = [x[1]^2+x[2]^2+y[1]^2-1, x[2]^2+x[3]^2+y[2]^2-1]
model = Model(optimizer_with_attributes(Mosek.Optimizer))
@variable(model, lower)
nonneg = f - lower*sum(x.^2)
info = add_psatz!(model, nonneg, [x; y], [], h, d, TS="block", Groebnerbasis=true)
@objective(model, Max, lower)
optimize!(model)
```

### Keyword arguments
Argument | Description | Default value
--- | :--- | :---
CS | Types of chordal extensions in exploiting correlative sparsity: "MF" (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation) | "MF"
cliques | Use customized variable cliques | []
TS | Types of chordal extensions used in term sparsity iterations: "block"(maximal chordal extension), "signsymmetry" (sign symmetries), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations) | "block"
QUIET | Silence the output| false
SO | Specify the sparse order | 1
Groebnerbasis | Work in the quotient ring by computing a Gröbner basis | false

## Methods
```@docs
add_psatz!
add_poly!
add_SOS!
```