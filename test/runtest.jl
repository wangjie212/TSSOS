using DynamicPolynomials
using TSSOS

# unconstrained optimization using the TSSOS hierarchy
@polyvar x[1:3]
f = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
opt,sol,data = tssos_first(f, x, newton=true, reducebasis=true, TS="block", solution=false, QUIET=true)
# optimum = 0.475274778453039

opt,sol,data = tssos_higher!(data, TS="MD", solution=true, QUIET=true)
# optimum = 0.4752747815682939

# constrained optimization using the TSSOS hierarchy
@polyvar x[1:3]
f = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
g = 1-x[1]^2-2*x[2]^2
h = x[2]^2+x[3]^2-1
pop = [f, g, h]
opt,sol,data = tssos_first(pop, x, 2, numeq=1, quotient=false, TS="block", Gram=true, solution=true, QUIET=true)
opt,sol,data = tssos_first(pop, x, 2, numeq=1, quotient=true, TS="block", solution=true, QUIET=true)
# optimum = 0.6861925732538557

opt,sol,data = tssos_higher!(data, TS="MD", solution=true, QUIET=true)
# optimum = 0.6861925611033903

opt,sol,data = tssos_first(pop, x, 2, numeq=1, quotient=false, TS=false, solution=true, QUIET=true)
# sol = extract_solutions(pop, x, 2, opt, data.moment[1], numeq=1)
# sol = extract_solutions_robust(3, 2, data.moment[1])

# constrained optimization using the CS-TSSOS hierarchy
n = 6
@polyvar x[1:n]
f = 1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
pop = [f, 1-sum(x[1:3].^2), 1-sum(x[3:6].^2)]
d = 2
opt,sol,data = cs_tssos_first(pop, x, d, numeq=1, TS="block", Gram=true, solution=true, QUIET=true)
# optimum = 0.7174253329322346

opt,sol,data = cs_tssos_higher!(data, TS="MD", solution=true, QUIET=true)
# optimum = 0.7174253409780793

using JuMP
using MosekTools

@polyvar x[1:3]
f = x[1]^2 + x[1]*x[2] + x[2]^2 + x[2]*x[3] + x[3]^2
d = 2 # set the relaxation order
@polyvar y[1:2]
h = [x[1]^2 + x[2]^2 + y[1]^2-1, x[2]^2 + x[3]^2 + y[2]^2 - 1]
model = Model(optimizer_with_attributes(Mosek.Optimizer))
@variable(model, lower)
nonneg = f - lower*sum(x.^2)
model,info = add_psatz!(model, nonneg, [x; y], [], h, d, TS="block", Groebnerbasis=true)
@objective(model, Max, lower)
optimize!(model)
optimum = objective_value(model)
@show optimum
# optimum = 0.29289321228548393

println("Run successfully!")
