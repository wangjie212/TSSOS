using DynamicPolynomials
using TSSOS
using JuMP
using MosekTools
using PermutationGroups
using Test


@testset begin

# Unconstrained Polynomial Optimization
@polyvar x[1:3]
f = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
opt,sol,data = tssos_first(f, x, newton=true, reducebasis=true, TS="block", solution=true, Gram=true, QUIET=true)
@test opt ≈ 0.4752748 atol = 1e-6
opt,sol,data = tssos_first(f, x, newton=false, reducebasis=false, TS="block", solution=true, Gram=true, QUIET=true)
@test opt ≈ 0.4752748 atol = 1e-6
opt,sol,data = tssos_higher!(data, TS="MD", solution=true, Gram=true, QUIET=true)
@test opt ≈ 0.4752748 atol = 1e-6

# Constrained Polynomial Optimization
@polyvar x[1:3]
f = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
g = 0.5-x[1]^2-2*x[2]^2
h = x[1]^2+x[2]^2+x[3]^2-1
pop = [f, g, h]
opt,sol,data = tssos_first(pop, x, 2, numeq=1, GroebnerBasis=false, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ 0.9555474 atol = 1e-6
opt,sol,data = tssos_first(pop, x, 2, numeq=1, GroebnerBasis=true, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ 0.9555474 atol = 1e-6
opt,sol,data = tssos_higher!(data, TS="block", solution=true, Gram=true, QUIET=true)
@test opt ≈ 0.9555474 atol = 1e-6
opt,sol,data = tssos_first(pop, x, 2, numeq=1, GroebnerBasis=true, TS="MD", Gram=true, solution=true, QUIET=true)
@test opt ≈ 0.9555474 atol = 1e-6
opt,sol,data = tssos_first(pop, x, 2, numeq=1, GroebnerBasis=false, TS=false, Gram=true, solution=true, QUIET=true)
@test opt ≈ 0.9555474 atol = 1e-6
opt,sol,data = tssos_first(pop, x, 2, numeq=1, GroebnerBasis=false, TS="signsymmetry", Gram=true, solution=true, QUIET=true)
@test opt ≈ 0.9555474 atol = 1e-6

@polyvar x[1:2]
pop = [-x[1]-x[2], 2x[1]^4-8x[1]^3+8x[1]^2+2-x[2], 4x[1]^4-32x[1]^3+88x[1]^2-96x[1]+36-x[2], (3-x[1])*x[1], (4-x[2])*x[2]]
opt,sol,data = tssos_first(pop, x, 3, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ -5.508013 atol = 1e-6

@polyvar x[1:3]
pop = [x[1]^2*x[2]^2*(x[1]^2+x[2]^2-3x[3]^2)+x[3]^6, 1-x[1]^2-x[2]^2-x[3]^2]
opt,sol,data = tssos_first(pop, x, 3, TS=false, Gram=true, solution=true, QUIET=true)
@test opt ≈ -0.0045964 atol = 1e-6

@polyvar x[1:4]
pop = [3-x[1]^2-x[3]^2+x[1]*x[2]^2+2x[2]*x[3]*x[4]-x[1]*x[4]^2, x[2], x[1]^2+3x[3]^2-2, x[4], x[1]^2+x[2]^2+x[3]^2+x[4]^2-3]
opt,sol,data = cs_tssos_first(pop, x, 2, numeq=3, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ -0.4142135 atol = 1e-6
opt,sol,data = cs_tssos_higher!(data, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ -0.4142135 atol = 1e-6

@polyvar x[1:6]
f = -25*(x[1]-2)^2-(x[2]-2)^2-(x[3]-1)^2-(x[4]-4)^2-(x[5]-1)^2-(x[6]-4)^2
pop = [f, (x[3]-3)^2+x[4]-4, (x[5]-3)^2+x[6]-4, 2-x[1]+3*x[2], 2+x[1]-x[2], 
6-x[1]-x[2], x[1]+x[2]-2, (5-x[3])*(x[3]-1), (6-x[4])*x[4], (5-x[5])*(x[5]-1),
(10-x[6])*x[6], x[1], x[2]]
d = 2
opt,sol,data = tssos_first(pop, x, d, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ -310 atol = 1e-3
opt,sol,data = cs_tssos_first(pop, x, d, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ -310 atol = 1e-3

@polyvar x[1:6]
f = 1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
pop = [f, 1-sum(x[1:3].^2), 1-sum(x[3:6].^2)]
d = 2
opt,sol,data = cs_tssos_first(pop, x, d, numeq=1, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ 0.71742533 atol = 1e-6
opt,sol,data = cs_tssos_higher!(data, TS="MD", solution=true, Gram=true, QUIET=true)
@test opt ≈ 0.71742533 atol = 1e-6
opt,sol,data = cs_tssos_first(pop, x, d, numeq=1, TS="signsymmetry", Gram=true, solution=true, QUIET=true)
@test opt ≈ 0.71742533 atol = 1e-6
opt,sol,data = cs_tssos_first(pop, x, d, numeq=1, TS=false, Gram=true, solution=true, QUIET=true)
@test opt ≈ 0.71742533 atol = 1e-6

# Solving systems of polynomial equations
@polyvar x[1:3]
pop = [1, 5x[1]^9-6x[1]^5*x[2]+x[1]*x[2]^4+2x[1]*x[3], -2x[1]^6*x[2]+2x[1]^2*x[2]^3+2x[2]*x[3], x[1]^2+x[2]^2-0.265625]
opt,sol,data = tssos_first(pop, x, 5, numeq=3, TS=false, GroebnerBasis=false, Gram=true, solution=true, QUIET=true)
@test opt ≈ 1 atol = 1e-6

@polyvar x[1:3]
pop = [1, x[1]^2-2x[1]*x[3]+5, x[1]*x[2]^2+x[2]*x[3]+1, 3x[2]^2-8x[1]*x[3]]
opt,sol,data = tssos_first(pop, x, 3, numeq=3, TS=false, Gram=true, solution=true, QUIET=true)
@test opt ≈ 1 atol = 1e-6

@polyvar x[1:6]
h1 = 2*x[6]^2 + 2x[5]^2 + 2x[4]^2 + 2x[3]^2 + 2x[2]^2 + x[1]^2 - x[1]
h2 = x[6]*x[5] + x[5]*x[4] + 2x[4]*x[3] + 2x[3]*x[2] + 2x[2]*x[1] - x[2]
h3 = 2x[6]*x[4] + 2x[5]*x[3] + 2x[4]*x[2] + x[2]^2 + 2x[3]*x[1] - x[3]
h4 = 2x[6]*x[3] + 2x[5]*x[2] + 2x[3]*x[2] + 2x[4]*x[1] - x[4]
h5 = x[3]^2 + 2x[6]*x[1] + 2x[5]*x[1] + 2x[4]*x[1] - x[5]
h6 = 2x[6] + 2x[5] + 2x[4] + 2x[3] + 2x[2] + x[1] - 1
opt,sol,data = tssos_first([1;h1;h2;h3;h4;h5;h6], x, 3, numeq=6, TS=false, Gram=true, solution=true, QUIET=true, rtol=1e-6)
@test opt ≈ 1 atol = 1e-6
opt,sol,data = tssos_first([1;x[1]-0.5;h1;h2;h3;h4;h5;h6], x, 2, numeq=6, TS=false, Gram=true, solution=true, QUIET=true, rtol=1e-4)
@test opt ≈ 1 atol = 1e-6

# SOS Programming
n = 3
@polyvar x[1:3]
f = [(x[1]^2+x[2]^2-1/4)*x[1], (x[2]^2+x[3]^2-1/4)*x[2], (x[2]^2+x[3]^2-1/4)*x[3]]
g = [1-x[1]^2, 1-x[2]^2, 1-x[3]^2]
d = 3
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
v, vc, vb = add_poly!(model, x, 2d-2)
w, wc, wb = add_poly!(model, x, 2d)
Lv = v - sum(f .* differentiate(v, x))
info1 = add_psatz!(model, Lv, x, g, [], d, TS="block", SO=1, constrs="con1")
info2 = add_psatz!(model, w, x, g, [], d, TS="block", SO=1)
info3 = add_psatz!(model, w-v-1, x, g, [], d, TS="block", SO=1)
moment = get_moment(wb, -ones(n), ones(n))
@objective(model, Min, sum(moment.*wc))
optimize!(model)
objv = objective_value(model)
@test objv ≈ 3.437648 atol = 1e-6

@polyvar x[1:3]
f = x[1]^2 + x[1]*x[2] + x[2]^2 + x[2]*x[3] + x[3]^2
d = 2 # set the relaxation order
@polyvar y[1:2]
h = [x[1]^2 + x[2]^2 + y[1]^2-1, x[2]^2 + x[3]^2 + y[2]^2 - 1]
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
@variable(model, lower)
nonneg = f - lower*sum(x.^2)
add_psatz!(model, nonneg, [x; y], [], h, d, TS="block", GroebnerBasis=true)
@objective(model, Max, lower)
optimize!(model)
opt = objective_value(model)
@test opt ≈ 0.29289321 atol = 1e-6

# Sum-Of-Ratio Optimization
@polyvar x y z
p = [x^2+y^2-y*z, y^2+x^2*z, z^2-x+y]
q = [1+2x^2+y^2+z^2, 1+x^2+2y^2+z^2, 1+x^2+y^2+2z^2]
g = [1-x^2-y^2-z^2]
d = 4
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
h1 = add_poly!(model, [x;y;z], 2d-2, signsymmetry=get_signsymmetry([p[2]; q[2]; g], [x;y;z]))[1]
h2 = add_poly!(model, [x;y;z], 2d-2, signsymmetry=get_signsymmetry([p[3]; q[3]; g], [x;y;z]))[1]
c = @variable(model)
add_psatz!(model, p[1]+(h1+h2-c)*q[1], [x;y;z], g, [], d, QUIET=true, TS="block", SO=1, GroebnerBasis=false)
add_psatz!(model, p[2]-h1*q[2], [x;y;z], g, [], d, QUIET=true, TS="block", SO=1, GroebnerBasis=false)
add_psatz!(model, p[3]-h2*q[3], [x;y;z], g, [], d, QUIET=true, TS="block", SO=1, GroebnerBasis=false)
@objective(model, Max, c)
optimize!(model)
objv = objective_value(model)
@test objv ≈ -0.346510 atol = 1e-6

# Polynomial Matrix Optimization
@polyvar x[1:2]
Q = [1/sqrt(2) -1/sqrt(3) 1/sqrt(6); 0 1/sqrt(3) 2/sqrt(6); 1/sqrt(2) 1/sqrt(3) -1/sqrt(6)]
F = Q*[-x[1]^2-x[2]^2 0 0; 0 -1/4*(x[1]+1)^2-1/4*(x[2]-1)^2 0; 0 0 -1/4*(x[1]-1)^2-1/4*(x[2]+1)^2]*Q'
G = [1-4x[1]*x[2] x[1]; x[1] 4-x[1]^2-x[2]^2]
opt,sol,data = tssos_first(F, [G], x, 2, TS="block", QUIET=true, Gram=true, Moment=true)
@test opt ≈ -4 atol = 1e-6
opt,sol,data = tssos_higher!(data, QUIET=true, Gram=true)
@test opt ≈ -4 atol = 1e-6

@polyvar x[1:3]
F = [1+x[1]^2 x[1]^2 x[2]^2; x[1]^2 1 x[3]^2; x[2]^2 x[3]^2 1]
G = [2-x[1]^2 1 x[1]+x[2]+x[3]; 1 2-x[2]^2 1; x[1]+x[2]+x[3] 1 2-x[3]^2]
opt,sol,data = tssos_first(F, [G], x, 2, TS=false, solution=true, Gram=true, QUIET=true)
@test opt ≈ -0.67070657 atol = 1e-6

@polyvar x[1:3]
F = [1+x[1]^2 x[1]^2 x[2]^2; x[1]^2 1 x[3]^2; x[2]^2 x[3]^2 1]
G = [2-x[1]^2 1 x[1]+x[2]+x[3]; 1 2-x[2]^2 1; x[1]+x[2]+x[3] 1 2-x[3]^2]
opt,sol,data = tssos_first(F, [G], x, 2, TS="block", Gram=true, QUIET=true)
@test opt ≈ -0.6823037 atol = 1e-6
opt,sol,data = tssos_higher!(data, QUIET=true)
@test opt ≈ -0.67070657 atol = 1e-6

# Exploiting Symmetry
@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry(pop, x, 2, G, QUIET=true)
@test opt ≈ -1.3987174 atol = 1e-6

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry(pop, x, 2, G, QUIET=true)
@test opt ≈ -1.3987174 atol = 1e-6

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry(pop, x, 2, G, numeq=1, QUIET=true)
@test opt ≈ -1.3987174 atol = 1e-6

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f; 1 .- x.^2]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry(pop, x, 2, G, numeq=3, SymmetricConstraint=false, QUIET=true)
@test opt ≈ 0 atol = 1e-6

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f; 1 .- x.^2]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry(pop, x, 2, G, SymmetricConstraint=false, QUIET=true)
@test opt ≈ -1.4174111 atol = 1e-6

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry_first(pop, x, 2, G, SymmetricConstraint=true, QUIET=true)
@test opt ≈ -1.3987174 atol = 1e-6
opt,data = tssos_symmetry_higher!(data, QUIET=true)
@test opt ≈ -1.3987174 atol = 1e-6

d = 2
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
lambda = @variable(model)
info = add_psatz_symmetry!(model, f - lambda, x, [1 - sum(x.^2)], [], d, G, SymmetricConstraint=true, TS="block", SO=2)
@objective(model, Max, lambda)
optimize!(model)
opt = objective_value(model)
@test opt ≈ -1.3987174 atol = 1e-6

# using Chebyshev basis
@polyvar x[1:3]
f = x[1]*x[2]*x[3]
g = [1-x[1]^2, 1-x[2]^2, 1-x[3]^2]
d = 2
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
@variable(model, lower)
add_psatz_cheby!(model, f-lower, x, g, [], d, TS="block", SO=1)
@objective(model, Max, lower)
optimize!(model)
opt = objective_value(model)
@test opt ≈ -1 atol = 1e-6

# sparse dynamic system
n = 3
@polyvar x[1:n]
f = [(x[1]^2+x[2]^2-1/4)*x[1], (x[2]^2+x[3]^2-1/4)*x[2], (x[2]^2+x[3]^2-1/4)*x[3]]
g = [1-x[1]^2, 1-x[2]^2, 1-x[3]^2]
d = 3
vsupp,vblocks,wsupp,wblocks,status = get_dynamic_sparsity(f, g, x, d, TS=["block","block"], SO=[2,1])
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
v, vc, vb = add_poly!(model, x, vsupp)
w, wc, wb = add_poly!(model, x, wsupp)
Lv = v - sum(f .* differentiate(v, x))
info1 = add_psatz!(model, Lv, x, g, [], d, blocks=[vblocks], constrs="con1")
info2 = add_psatz!(model, w, x, g, [], d, blocks=[wblocks], constrs="con2")
info3 = add_psatz!(model, w-v-1, x, g, [], d, blocks=[wblocks], constrs="con3")
moment = get_moment(wb, -ones(n), ones(n))
@objective(model, Min, sum(moment.*wc))
optimize!(model)
opt = objective_value(model)
@test opt ≈ 3.437648 atol = 1e-6

end
