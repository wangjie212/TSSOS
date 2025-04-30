using JuMP
using MosekTools
using DynamicPolynomials
using TSSOS

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
objv = objective_value(model)
@show objv
## objv = 3.437648
