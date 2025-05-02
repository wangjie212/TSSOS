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
opt,sol,data = tssos_first(f, x, newton=true, reducebasis=true, TS="block", solution=false, QUIET=true)
@test opt ≈ 0.4752748 atol = 1e-6
opt,sol,data = tssos_higher!(data, TS="MD", solution=true, QUIET=true)
@test opt ≈ 0.4752748 atol = 1e-6

# Constrained Polynomial Optimization
@polyvar x[1:3]
f = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
g = 1-x[1]^2-2*x[2]^2
h = x[2]^2+x[3]^2-1
pop = [f, g, h]
opt,sol,data = tssos_first(pop, x, 2, numeq=1, GroebnerBasis=false, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ 0.68619257 atol = 1e-6
opt,sol,data = tssos_first(pop, x, 2, numeq=1, GroebnerBasis=true, TS="block", solution=true, QUIET=true)
@test opt ≈ 0.68619257 atol = 1e-6
opt,sol,data = tssos_higher!(data, TS="MF", solution=true, QUIET=true)
@test opt ≈ 0.68619257 atol = 1e-6
opt,sol,data = tssos_first(pop, x, 2, numeq=1, GroebnerBasis=false, TS=false, solution=true, QUIET=true)
@test opt ≈ 0.68619257 atol = 1e-6

n = 6
@polyvar x[1:n]
f = 1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
pop = [f, 1-sum(x[1:3].^2), 1-sum(x[3:6].^2)]
d = 2
opt,sol,data = cs_tssos_first(pop, x, d, numeq=1, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ 0.71742533 atol = 1e-6
opt,sol,data = cs_tssos_higher!(data, TS="MD", solution=true, QUIET=true)
@test opt ≈ 0.71742533 atol = 1e-6

# SOS Programming
n = 3
@polyvar x[1:n]
f = [(x[1]^2+x[2]^2-1/4)*x[1], (x[2]^2+x[3]^2-1/4)*x[2], (x[2]^2+x[3]^2-1/4)*x[3]]
g = [1-x[1]^2, 1-x[2]^2, 1-x[3]^2]
d = 3
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
v, vc, vb = add_poly!(model, x, 2d-2)
w, wc, wb = add_poly!(model, x, 2d)
Lv = v - sum(f .* differentiate(v, x))
info1 = add_psatz!(model, Lv, x, g, [], d, TS=false, SO=1, constrs="con1")
info2 = add_psatz!(model, w, x, g, [], d, TS=false, SO=1)
info3 = add_psatz!(model, w-v-1, x, g, [], d, TS=false, SO=1)
moment = get_moment(wb, -ones(n), ones(n))
@objective(model, Min, sum(moment.*wc))
optimize!(model)
objv = objective_value(model)
@test objv ≈ 3.437648 atol = 1e-6

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
add_psatz!(model, p[1]+(h1+h2-c)*q[1], [x;y;z], g, [], d, QUIET=true, CS=false, TS=false, SO=1, GroebnerBasis=false)
add_psatz!(model, p[2]-h1*q[2], [x;y;z], g, [], d, QUIET=true, CS=false, TS="block", SO=1, GroebnerBasis=false)
add_psatz!(model, p[3]-h2*q[3], [x;y;z], g, [], d, QUIET=true, CS=false, TS="block", SO=1, GroebnerBasis=false)
@objective(model, Max, c)
optimize!(model)
objv = objective_value(model)
@test objv ≈ -0.346510 atol = 1e-6

# Polynomial Matrix Optimization
@polyvar x[1:2]
Q = [1/sqrt(2) -1/sqrt(3) 1/sqrt(6); 0 1/sqrt(3) 2/sqrt(6); 1/sqrt(2) 1/sqrt(3) -1/sqrt(6)]
F = Q*[-x[1]^2-x[2]^2 0 0; 0 -1/4*(x[1]+1)^2-1/4*(x[2]-1)^2 0; 0 0 -1/4*(x[1]-1)^2-1/4*(x[2]+1)^2]*Q'
G = [1-4x[1]*x[2] x[1]; x[1] 4-x[1]^2-x[2]^2]
opt,data = tssos_first(F, [G], x, 2, TS="block", QUIET=true)
@test opt ≈ -4 atol = 1e-6
opt,data = tssos_higher!(data, QUIET=true)
@test opt ≈ -4 atol = 1e-6

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

# Exploiting Symmetry
@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt = tssos_symmetry(pop, x, 2, G, QUIET=true)[1]
@test opt ≈ -1.3987174 atol = 1e-6

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt = tssos_symmetry(pop, x, 2, G, QUIET=true)[1]
@test opt ≈ -1.3987174 atol = 1e-6

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt = tssos_symmetry(pop, x, 2, G, numeq=1, QUIET=true)[1]
@test opt ≈ -1.3987174 atol = 1e-6

end
