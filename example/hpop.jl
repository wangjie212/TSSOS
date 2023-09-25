using JuMP
using MosekTools
using DynamicPolynomials
using TSSOS

n = 3
@polyvar x[1:n]
f = x[1]^2 + x[1]*x[2] + x[2]^2 + x[2]*x[3] + x[3]^2
d = 2
@polyvar y[1:2]
h = [x[1]^2+x[2]^2+y[1]^2-1, x[2]^2+x[3]^2+y[2]^2-1]
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
@variable(model, lower)
nonne = f - lower*sum(x.^2)
model,info = add_psatz!(model, nonne, [x; y], [], h, d, QUIET=true, CS=true, TS="block", Groebnerbasis=true)
@objective(model, Max, lower)
optimize!(model)
status = termination_status(model)
if status != MOI.OPTIMAL
    println("termination status: $status")
    status = primal_status(model)
    println("solution status: $status")
end
obj = objective_value(model)
println(obj)
