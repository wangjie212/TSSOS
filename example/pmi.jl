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

## polynomial matrix optimization with correlative sparsity
@polyvar x[1:3]
F = [x[1]*x[2] + x[2]*x[3] 2 + x[1] + x[3]; 2 + x[1] + x[3] 2*x[2]^2]
G = [[1 - x[1]^2 - x[2]^2], [1 - x[2]^2 - x[3]^2]]
@time opt,data = cs_tssos_first(F, G, x, 1, TS=false)
# @time opt,data = cs_tssos_higher!(data)
