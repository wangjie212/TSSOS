using DynamicPolynomials
using TSSOS

# unconstrained optimization using the TSSOS hierarchy
@polyvar x[1:3]
f = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
opt,sol,data = tssos_first(f, x, newton=true, reducebasis=true, TS="block", solution=true, QUIET=true)
# optimum = 0.475274778453039

opt,sol,data = tssos_higher!(data, TS="MD", solution=true, QUIET=true)
# optimum = 0.4752747815682939

# constrained optimization using the TSSOS hierarchy
@polyvar x[1:3]
f = 1+x[1]^4+x[2]^4+x[3]^4+x[1]*x[2]*x[3]+x[2]
g = 1-x[1]^2-2*x[2]^2
h = x[2]^2+x[3]^2-1
pop = [f, g, h]
opt,sol,data = tssos_first(pop, x, 2, numeq=1, TS="block", solution=true, QUIET=true)
# optimum = 0.6861925732538557

opt,sol,data = tssos_higher!(data, TS="MD", solution=true, QUIET=true)
# optimum = 0.6861925611033903

# constrained optimization using the CS-TSSOS hierarchy
n = 6
@polyvar x[1:n]
f = 1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
pop = [f, 1-sum(x[1:3].^2), 1-sum(x[3:6].^2)]
d = 2
opt,sol,data = cs_tssos_first(pop, x, d, numeq=1, TS="block", solution=true, QUIET=true)
# optimum = 0.7174253329322346

opt,sol,data = cs_tssos_higher!(data, TS="MD", solution=true, QUIET=true)
# optimum = 0.7174253409780793

println("Run successfully!")
