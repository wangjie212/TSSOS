using DynamicPolynomials
using TSSOS
using JuMP
using MosekTools
using MultivariatePolynomials

s = 1
t = 6
@polyvar x[1:2]
@polyvar y
model = Model(optimizer_with_attributes(Mosek.Optimizer))
sos1 = add_SOS!(model, y, s)
sos2 = add_SOS!(model, y, s-1)
sos3 = add_SOS!(model, x[1], t)
sos4 = add_SOS!(model, x[1], t-1)
sos5 = add_SOS!(model, x[1], t)
sos6 = add_SOS!(model, x[1], t-1)
p = (sos1 + sos2*(1 - y^2))*(((x[1] - y)^2 - 2)^2 + ((x[2] - y)^2 - 2)^2 - 2)
q1 = sos3 + sos4*(1 - x[1]^2)
q2 = subs(q1, x[1]=>x[2])
q3 = sos5 + sos6*(1 - x[1]^2)
q4 = subs(q3, x[1]=>x[2])
msos1 = add_SOSMatrix!(model, x, 2, t)[1]
msos2 = add_SOSMatrix!(model, x, 2, t-1)[1]
msos3 = add_SOSMatrix!(model, x, 2, t-1)[1]
@variable(model, τ)
@constraint(model, coefficients(0.5*sum(MultivariatePolynomials.coefficient(p, y^i, [y])*(2/(i+1)) for i = 0:2:2s+4) + 1 + 0.5*(q1 + q2 - q3 - q4) - msos1[1,2] - msos2[1,2]*(1-x[1]^2) - msos3[1,2]*(1-x[2]^2)) .== 0)
@constraint(model, coefficients(τ - msos1[1,1] - msos2[1,1]*(1-x[1]^2) - msos3[1,1]*(1-x[2]^2)) .== 0)
@constraint(model, coefficients(((x[1] - x[2])^2 - 2)^2 - msos1[2,2] - msos2[2,2]*(1-x[1]^2) - msos3[2,2]*(1-x[2]^2)) .== 0)
@objective(model, Min, 0.5*τ + coefficients(sum(MultivariatePolynomials.coefficient(q3 - 0.1*q1, x[1]^i, [x[1]])*(2/(i+1)) for i = 0:2:2t))[1])
optimize!(model)
optimum = objective_value(model)
@show optimum
