using TSSOS
using DynamicPolynomials

N = 20
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
    push!(pop, τ[(i-1)*N - Int(i*(i+1)/2) + j]*((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2) - 1)
end

opt,sol,data = cs_tssos_first(pop, [x;y;z;τ], 3, numeq=length(pop)-1)
