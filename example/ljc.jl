using TSSOS
using DynamicPolynomials

N = 5
@polyvar x[1:N]
@polyvar y[1:N]
@polyvar z[1:N]
@polyvar τ[1:Int(N*(N-1)/2)]

f = 0
for i = 1:N-1, j = i+1:N
    f += τ[(i-1)*N - Int(i*(i+1)/2) + j]^6 - τ[(i-1)*N - Int(i*(i+1)/2) + j]^3
end
pop = [f]
for i = 1:N-1, j = i+1:N
    push!(pop, 7 - τ[(i-1)*N - Int(i*(i+1)/2) + j]^2 - x[i]^2 - x[j]^2 - y[i]^2 - y[j]^2 - z[i]^2 - z[j]^2)
end
push!(pop, 3N - sum(x.^2) - sum(y.^2) - sum(z.^2))
for i = 1:N-1, j = i+1:N
    push!(pop, τ[(i-1)*N - Int(i*(i+1)/2) + j]*((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2) - 1)
end
opt,sol,data = cs_tssos_first(pop, [x;y;z;τ], 3, numeq=Int(N*(N-1)/2), TS="MD")


using JuMP
using MosekTools

N = 2
@polyvar x[1:N]
@polyvar y[1:N]
@polyvar z[1:N]
@polyvar τ[1:Int(N*(N-1)/2)]
R = [[(x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2 for j = i+1:N] for i = 1:N-1]
ineqcons = Poly{Int}[]
eqcons = Poly{Int}[]
for i = 1:N-1, j = i+1:N
    push!(ineqcons, 7 - τ[(i-1)*N - Int(i*(i+1)/2) + j]^2 - x[i]^2 - x[j]^2 - y[i]^2 - y[j]^2 - z[i]^2 - z[j]^2)
    push!(eqcons, ((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2)^3 - τ[(i-1)*N - Int(i*(i+1)/2) + j]^2)
end
push!(ineqcons, 3N - sum(x.^2) - sum(y.^2) - sum(z.^2))

model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), false)
λ = @variable(model)
H = [[add_poly!(model, [x[i];x[j];y[i];y[j];z[i];z[j];τ[(i-1)*N-Int(i*(i+1)/2)+j]], 2)[1] for j = i+1:N] for i = 1:N-1]
add_psatz!(model, sum(sum(h) for h in H) - λ, [x;y;z;τ], ineqcons, eqcons, 3, CS=true, TS="MD", SO=1, QUIET=false)
for i = 1:N-1, j = i+1:N
    ind = (i-1)*N - Int(i*(i+1)/2) + j
    add_psatz!(model, 1 - τ[ind]^2 - τ[ind]^4*H[i][j-i], [x[i];x[j];y[i];y[j];z[i];z[j];τ[(i-1)*N-Int(i*(i+1)/2)+j]], [ineqcons[ind]], [eqcons[ind]], 3, TS="MD", SO=1, QUIET=false)
end
@objective(model, Max, λ)
optimize!(model)
objv = objective_value(model)
@show objv
