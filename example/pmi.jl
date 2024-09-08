using DynamicPolynomials
using TSSOS
using LinearAlgebra

## Inf mineig(F(x)) s.t. G1(x) >= 0, ..., Gm(x) >= 0
@polyvar x[1:2]
Q = [1/sqrt(2) -1/sqrt(3) 1/sqrt(6); 0 1/sqrt(3) 2/sqrt(6); 1/sqrt(2) 1/sqrt(3) -1/sqrt(6)]
F = Q*[-x[1]^2-x[2]^2 0 0; 0 -1/4*(x[1]+1)^2-1/4*(x[2]-1)^2 0; 0 0 -1/4*(x[1]-1)^2-1/4*(x[2]+1)^2]*Q'
G = [1-4x[1]*x[2] x[1]; x[1] 4-x[1]^2-x[2]^2]
opt,data = tssos_first(F, [G], x, 2, TS="block")
# opt,data = tssos_higher!(data)

@polyvar x[1:2]
m = 20
A = rand(m, m)
A = (A+A')/2
B = rand(m, m)
B = (B+B')/2
F = (1-x[1]^2-x[2]^2)*I(m) + (x[1]*x[2]-x[1]^2)*A + (2x[1]^2*x[2]^2-x[1]*x[2]-2x[2]^2)*B
G = [1-x[1]^2-x[2]^2]
@time opt,data = tssos_first(F, [G], x, 2, TS="block", QUIET=true)

## Inf b'*λ s.t. F0 + λ1*F1 + ... λt*Ft >=0 on {x ∈ R^n | G1(x) >= 0, ..., Gm(x) >= 0}
@polyvar x[1:3]
F = Vector{Matrix{Polynomial{true, Float64}}}(undef, 3)
F[1] = sum(x.^2)*[x[2]^4 0 0; 0 x[3]^4 0; 0 0 x[1]^4]
F[2] = sum(x.^2)*[0 x[1]^2*x[2]^2 0; x[1]^2*x[2]^2 0 0; 0 0 0]
F[3] = sum(x.^2)*[x[1]^4 0 0; 0 x[2]^4 x[2]^2*x[3]^2; 0 x[2]^2*x[3]^2 x[3]^4]
G = [1 - sum(x.^2)]
@time opt,data = LinearPMI_first([-10, 1], F, [G], x, 3, TS="block", QUIET=true)
@time opt,data = LinearPMI_higher!(data, TS="block", QUIET=true)

## polynomial matrix optimization with correlative sparsity
@polyvar x[1:3]
F = [x[1]*x[2] + x[2]*x[3] 2 + x[1] + x[3]; 2 + x[1] + x[3] 2*x[2]^2]
G = [[1 - x[1]^2 - x[2]^2], [1 - x[2]^2 - x[3]^2]]
@time opt,data = cs_tssos_first(F, G, x, 1, TS=false)
# @time opt,data = cs_tssos_higher!(data)


@polyvar x[1:3]
F = Matrix{Polynomial{true, Float64}}(undef, 1, 1)
F[1,1] = x[1]^4 + x[1]^2*x[3] + x[2]^4 + x[3]^4 + 1
G = Vector{Matrix{Polynomial{true, Float64}}}(undef, 1)
G[1] = Matrix{Polynomial{true, Float64}}(undef, 1, 1)
G[1][1,1] = - x[1]*x[2] + x[3]
@time opt,data = tssos_first(F, G, x, 2, TS="block")


@polyvar x[1:5]
F = [x[1]^4 x[1]^2-x[2]*x[3] x[3]^2-x[4]*x[5] 0 0;
x[1]^2-x[2]*x[3] x[2]^4 x[2]^2-x[3]*x[4] 0 0;
x[3]^2-x[4]*x[5] x[2]^2-x[3]*x[4] x[3]^4 x[4]^2-x[1]*x[2] x[5]^2-x[3]*x[4];
0 0 x[4]^2-x[1]*x[2] x[4]^4 x[4]^2-x[1]*x[3];
0 0 x[5]^2-x[3]*x[4] x[4]^2-x[1]*x[3] x[5]^4]
G = Vector{Matrix{Polynomial{true, Float64}}}(undef, 2)
G[1] = [1-x[1]^2-x[2]^2 x[2]*x[3]; x[2]*x[3] 1-x[3]^2]
G[2] = [1-x[4]^2 x[4]*x[5]; x[4]*x[5] 1-x[5]^2]
@time opt,data = tssos_first(F, G, x, 4, TS="MD", QUIET=true)
@time opt,data = tssos_higher!(data, TS="MD", QUIET=true)
println(maximum(maximum.([maximum.(data.blocksize[i]) for i = 1:data.cql])))


n = 5
r = 5
@polyvar x[1:n]
F = [sum(x[k]^2 for k = 1:n-2) sum(x[k]*x[k+1] for k = 1:n-1) 1;
sum(x[k]*x[k+1] for k = 1:n-1) sum(x[k]^2 for k = 2:n-1)  sum(x[k]*x[k+2] for k = 1:n-2);
1 sum(x[k]*x[k+2] for k = 1:n-2) sum(x[k]^2 for k = 3:n)]
G = Vector{Matrix{Polynomial{true, Float64}}}(undef, n-2)
for k = 1:n-2
    G[k] = [1-x[k]^2-x[k+1]^2 x[k+1]+0.5; x[k+1]+0.5 1-x[k+2]^2]
end
@time opt,data = cs_tssos_first(F, G, x, r, TS=false, QUIET=true)
# @time opt,data = cs_tssos_higher!(data, TS="block")
@time opt,data = cs_tssos_first(F, G, x, r, CS=false, TS=false, QUIET=true)


@polyvar x[1:5]
F = [x[1]^4+x[2]^4+1 x[1]*x[3]; x[1]*x[3] x[3]^4+x[4]^4+x[5]^4+0.5]
G = Vector{Matrix{Polynomial{true, Float64}}}(undef, 1)
G[1] = [1-x[1]^2 x[1]*x[2] x[1]*x[3] 0 0; x[1]*x[2] 1-x[2]^2 x[2]*x[3] 0 0; x[1]*x[3] x[2]*x[3] 1-x[3]^2 x[3]*x[4] x[3]*x[5]; 0 0 x[3]*x[4] 1-x[4]^2 x[4]*x[5]; 0 0 x[3]*x[5] x[4]*x[5] 1-x[5]^2]
@time opt,data = cs_tssos_first(F, G, x, 2, TS=false, QUIET=true)
# @time opt,data = tssos_higher!(data, TS="block", QUIET=true)
println(maximum(maximum.([maximum.(data.blocksize[i]) for i = 1:data.cql])))

@polyvar x[1:6]
F = [x[1]^4+x[2]^4+1 x[1]*x[3]; x[1]*x[3] x[3]^4+x[4]^4+x[5]^4+0.5]
G = Vector{Matrix{Polynomial{true, Float64}}}(undef, 2)
G[1] = [1-x[1]^2 x[1]*x[2] x[1]*x[3]; x[1]*x[2] 1-x[2]^2 x[2]*x[3]; x[1]*x[3] x[2]*x[3] x[6]^2]
G[2] = [1-x[3]^2-x[6]^2 x[3]*x[4] x[3]*x[5]; x[3]*x[4] 1-x[4]^2 x[4]*x[5]; x[3]*x[5] x[4]*x[5] 1-x[5]^2]
@time opt,data = cs_tssos_first(F, G, x, 2, TS=false, QUIET=true)
# @time opt,data = cs_tssos_higher!(data, TS="block", QUIET=true)
println(maximum(maximum.([maximum.(data.blocksize[i]) for i = 1:data.cql])))