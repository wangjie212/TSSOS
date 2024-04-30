using DynamicPolynomials
using TSSOS
using LinearAlgebra

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
@time begin
opt,data = tssos_first(F, [G], x, 2, TS="block", QUIET=true)
println([Int.(data.blocksize[i]) for i = 1:2])
end
@time begin
opt,data = tssos_first(F, [G], x, 2, TS=false, QUIET=true)
println([Int.(data.blocksize[i]) for i = 1:2])
end