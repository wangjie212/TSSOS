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

## Image of a semialgebraic set by a polynomial map
Example 1, Section 6.1 from https://arxiv.org/pdf/1507.06143. 
The goal is to approximate the image of the two-dimensional unit ball $\mathbf{S} =$ { $\mathbf{x} \in \mathbb{R}^2 : x_1^2 + x_2^2 \leq 1$ } under the polynomial application $f(\mathbf{x})=(x_1+x_1 x_2, x_2-x_1^3)/2$. 
We choose $\mathbf{B}=\mathbf{S}$ since it can be checked that $f(\mathbf{S}) \subset \mathbf{B}$. 


```Julia
using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using SpecialFunctions
using TSSOS
using Plots

# Moments of the Lebesgue measure on the unit ball
function momball(a)
n,m = size(a);
y = zeros(m);
for k = 1:m
  if all(.!Bool.(rem.(a[:,k],2)))
    y[k]=prod(gamma.((a[:,k].+1)/2))/gamma(1+(n+sum(a[:,k]))/2);
 end
end
return y
end

# Half-degree of polynomials w(x) and v(x)
d = 5;

n = 2
@polyvar x[1:n] 
f = [x[1]+x[1]*x[2]; x[2]-x[1]^3 ] * 0.5;
gS = 1 - x[1]^2 - x[2]^2;
gB = gS;
model = Model(optimizer_with_attributes(Mosek.Optimizer))
#set_optimizer_attribute(model, MOI.Silent(), true)

v, vc, vb = add_poly!(model, x, 2*d);
w, wc, wb = add_poly!(model, x, 2*d);
vf=subs(v,x[1]=>f[1],x[2]=>f[2]);
dv = Int(ceil(maxdegree(vf)/2));

# Constraints
info1 = add_psatz!(model, vf, x, [gS], [], dv, QUIET=false, CS=false, TS=false, Groebnerbasis=false); # v o f >= 0 on S
info2 = add_psatz!(model, w-1-v, x, [gB], [], d, QUIET=false, CS=false, TS=false, Groebnerbasis=false); # w >= v + 1 on B
info3 = add_psatz!(model, w, x, [gB], [], d, QUIET=false, CS=false, TS=false, Groebnerbasis=false); # w >= 0 on B

supp = TSSOS.get_basis(n, 2*d)

# Lebesgue moments on B

moment = momball(supp);
@objective(model, Min, moment'*wc); # minimization of int w d_lambda


optimize!(model)
status = termination_status(model)
if status != MOI.OPTIMAL
    println("termination status: $status")
    status = primal_status(model)
    println("solution status: $status")
end
objv = objective_value(model);
wp = value.(wc)'*wb;

# Plot the superlevel set of w-1
x1 = range(-1, 1, length=1000);
x2 = range(-1, 1, length=1000);
hw(x1,x2) = if x1^2+x2^2<=1.0 wp(x1,x2) else 0.0 end;
zw = @. hw(x1', x2);
p=contour(x1,x2,zw, level=[1], color=[:white,:gray], levels=1, cbar=false,grid=false,fill=true)

# Sample the image set f(S)
N=10^5;
X = randn(2,N);
X=mapslices(c -> rand(1)[1]*c/sqrt(sum(c.^2)), X; dims=1);
f1=mapslices(c->f[1](c), X; dims=1)[1,:];
f2=mapslices(c->f[2](c), X; dims=1)[1,:];
scatter!(p,f1,f2,mc=:black,legend=false);

# Draw the unit circle
t=range(0,2*pi,length=100);
xt=cos.(t);
yt=sin.(t);
plot!(p,xt,yt,color=:black,legend=false, ylimits=(-1,1), xlimits=(-1,1), aspect_ratio=:equal);
```

<p align="center">
<img src="https://homepages.laas.fr/vmagron/files/tssos/image5.png" alt="Image approximation" style="width:40%; border:0;">
</p>

The black dots correspond to the image set of the points obtained by uniform sampling of $\mathbf{S}$ under $f$. 
The outer approximation obtained at the 5-th relaxation order is represented in light gray.
