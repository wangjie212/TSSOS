using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using TSSOS

function basis(x)
    basis = Monomial{true}[1]
    push!(basis, x...)
    for i = 1:length(x), j = i:length(x)
        push!(basis, x[i]*x[j])
    end
    return basis
end

n = 4
@polyvar x[1:n]
@polyvar y[1:3]
b = basis(x[1:n])
f = rand(length(b))'*b*y[1] + rand(length(b))'*b*y[2] + rand(length(b))'*b*y[3] + rand(length(b))'*b
g1 = rand(length(b))'*b*y[1] + rand(length(b))'*b*y[2] + rand(length(b))'*b*y[3] + rand(length(b))'*b
g2 = 1 - sum(x.^2)
g3 = rand(length(b))'*b - y[1]
g4 = y[1] + rand(length(b))'*b
g5 = rand(length(b))'*b - y[2]
g6 = y[2] + rand(length(b))'*b
g7 = rand(length(b))'*b - y[3]
g8 = y[3] + rand(length(b))'*b

opt,sol,data = tssos_first([f;g1;g2;g3;g4;g5;g6;g7;g8], [x;y], 3, TS=false, QUIET=true)

d = 3
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
s0 = add_SOS!(model, x, d)
s1 = add_SOS!(model, x, d-1)
s2 = add_SOS!(model, x, d-1)
s3 = add_SOS!(model, x, d-1)
s4 = add_SOS!(model, x, d-1)
s5 = add_SOS!(model, x, d-1)
s6 = add_SOS!(model, x, d-1)
s7 = add_SOS!(model, x, d-1)
s8 = add_SOS!(model, x, d-1)
@variable(model, lower)
@objective(model, Max, lower)
@constraint(model, coefficients(f - lower - s0 - s1*g1 - s2*g2 - s3*g3 - s4*g4 - s5*g5 - s6*g6 - s7*g7 - s8*g8) .== 0)
optimize!(model)
status = termination_status(model)
if status != MOI.OPTIMAL
    println("termination status: $status")
    status = primal_status(model)
    println("solution status: $status")
end
objv = objective_value(model)
@show objv
