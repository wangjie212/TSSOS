using DynamicPolynomials
using TSSOS
using JuMP
using MosekTools

@polyvar x[1:3]
f = x[1]*x[2]*x[3]
g = [1-x[1]^2, 1-x[2]^2, 1-x[3]^2]

opt,sol,data = tssos_first([f; g], x, 2, TS="block", QUIET=false)
# opt,sol,data = tssos_higher!(data, TS="block", QUIET=false)

d = 2
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), false)
@variable(model, lower)
info = add_psatz_cheby!(model, f-lower, x, g, [], d, TS="block", SO=1)
@objective(model, Max, lower)
optimize!(model)
objv = objective_value(model)
@show objv
