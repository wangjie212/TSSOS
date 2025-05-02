using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using TSSOS

## Example 2
@polyvar x[1:5]
Φ = [1 - x[1]^4 - x[2]^4 + 0.1*x[3]^4, 10*x[3]^4 - x[1]^4 - x[2]^4]
ψ = [4*x[4]^2*(x[1]^2 + x[2]^2) - (sum(x[1:4].^2) - x[5]^2)^2, 6 - x[4], x[4] - 4, 1 - x[5], x[5] - 0.5]
@polyvar z # homogenization variable
Φ = homogenize.(Φ, z)
ψ = homogenize.(ψ, z)
d = 2 # relaxation order
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), false)
# generate an interpolant polynomial
# hc = @variable(model, [1:3])
# hb = x[1:3].^2
# h = 1 + hc'*hb
h,hc,hb = add_poly!(model, x[1:3], 2)
@constraint(model, sum(hc)==1)
hh = homogenize(h, z)
add_psatz!(model, hh, [x[1:3]; z], [Φ; z], [1-z^2-sum(x[1:3].^2)], d, QUIET=true, CS=false, TS=false)
add_psatz!(model, -hh, [x; z], [ψ; z], [1-z^2-sum(x.^2)], d, QUIET=true, CS=false, TS=false)
optimize!(model)
status = termination_status(model)
if status != MOI.OPTIMAL
    println("termination status: $status")
    status = primal_status(model)
    println("solution status: $status")
end

# retrieve the interpolant polynomial
p = value.(hc)'*hb
@show p