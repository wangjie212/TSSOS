using TSSOS
using DynamicPolynomials
using Test

@testset begin

supp = Vector{Vector{Vector{UInt16}}}[[[[], []], [[1], [1]], [[1], [2;2]], [[2;2], [1]]],
[[[2], []], [[], [2]]], [[[], []], [[1], [1]], [[1;1], []], [[], [1;1]]],
[[[], []], [[1], [1]], [[2], [2]]], [[[2], []], [[], [2]]]]
coe = Vector{ComplexF64}[[3;-1;-0.5im;0.5im], [1;1], [-1;1;-0.25;-0.25], [-3;1;1], [im;-im]]
opt,sol,data = cs_tssos_first(supp, coe, 2, 2, numeq=3, QUIET=true, TS=false, Gram=true)
@test opt ≈ 0.428174 atol = 1e-6

supp = Vector{Vector{Vector{UInt16}}}[[[[1;1], [2]], [[2], [1;1]], [[1], [2]], [[2], [1]], [[2], [3]], [[3], [2]]],
[[[], []], [[1], [1]], [[2],[2]]], [[[], []], [[2], [2]], [[3], [3]]],
[[[], []], [[], [2]], [[2], []]]]
coe = [[1;1;1;1;1;1], [1;-1;-1], [1;-1;-1], [1;1;1]]
opt,sol,data = cs_tssos_first(supp, coe, 3, 2, numeq=1, QUIET=true, TS="block", Gram=true)
@test opt ≈ -2.682188 atol = 1e-6

@polyvar z[1:4]
pop = [z[1]+z[2]+z[3]+z[4], (1+im)*z[1]^2*z[4]+(1-im)*z[2]*z[3]^2, 1-z[1]*z[3]-z[2]*z[4]]
opt,sol,data = cs_tssos_first(pop, z, 2, 2, numeq=1, QUIET=true, TS="block", ConjugateBasis=true, solution=true, Gram=true)
@test opt ≈ -2.7423247 atol = 1e-6
opt,sol,data = cs_tssos_first(pop, z, 2, 2, numeq=1, QUIET=true, TS="block", ConjugateBasis=true, normality=2, solution=true, Gram=true)
@test opt ≈ -2.7423247 atol = 1e-6
opt,sol,data = cs_tssos_first(pop, z, 2, 2, numeq=1, QUIET=true, TS="block", solution=true, Gram=true)
@test opt ≈ -2.7423247 atol = 1e-6
opt,sol,data = cs_tssos_higher!(data, QUIET=true, TS="block", solution=true, Gram=true)
@test opt ≈ -2.7423247 atol = 1e-6

@polyvar z[1:4]
pop = [z[1]*z[2] + z[3]*z[4], 1 - z[1]*z[3] - z[2]*z[4]]
opt,sol,data = cs_tssos_first(pop, z, 2, 2, numeq=1, QUIET=true, TS="block", solution=true, Gram=true)
@test opt ≈ -1 atol = 1e-6

@polyvar z1 z2 z3 z4 z5 z6
z = tuple(z1,z2,z3,z4,z5,z6)
f = - z3*z6
g1 = z1^2*z4^2 - z1^2 - z4^2 - 4*z3*z6 + 1
g2 = z2^2*z5^2 - z2^2 - z5^2 - 4*z3*z6 + 1
h1 = - z1*z4 - z2*z5 + 2/3
h2 = 3*z1*z2 + 3*z4*z5 - 2
h3 = 9*z1*z2*z4*z5 - 1
cpop = [f; g1; g2; h1; h2; h3]
opt,sol,data = cs_tssos_first(cpop, z, 3, 2, numeq=3, CS=false, TS="block", solution=true, QUIET=true, Gram=true)
@test opt ≈ -0.4444444 atol = 1e-6

n = 3
@polyvar z1 z2 z3 z4 z5 z6 z7 z8
z = tuple(z1,z2,z3,z4,z5,z6,z7,z8)
f = - z4*z8
g1 = 9 - 6*z1^3 + 12*z1^2*z2 + 12*z1^2*z3 - 6*z5^3 + 4*z1^3*z5^3 -
 8*z1^2*z2*z5^3 - 8*z1^2*z3*z5^3 + 12*z5^2*z6 -
   8*z1^3*z5^2*z6 + 16*z1^2*z2*z5^2*z6 + 16*z1^2*z3*z5^2*z6 +
 12*z5^2*z7 - 8*z1^3*z5^2*z7 + 16*z1^2*z2*z5^2*z7 +
   16*z1^2*z3*z5^2*z7 - 36*z4*z8
g2 = 9 + 12*z1*z2^2 - 6*z2^3 + 12*z2^2*z3 + 12*z5*z6^2 +
 16*z1*z2^2*z5*z6^2 - 8*z2^3*z5*z6^2 + 16*z2^2*z3*z5*z6^2 -
   6*z6^3 - 8*z1*z2^2*z6^3 + 4*z2^3*z6^3 - 8*z2^2*z3*z6^3 +
 12*z6^2*z7 + 16*z1*z2^2*z6^2*z7 - 8*z2^3*z6^2*z7 +
   16*z2^2*z3*z6^2*z7 - 36*z4*z8
g3 = 9 + 12*z1*z3^2 + 12*z2*z3^2 - 6*z3^3 + 12*z5*z7^2 +
 16*z1*z3^2*z5*z7^2 + 16*z2*z3^2*z5*z7^2 - 8*z3^3*z5*z7^2 +
   12*z6*z7^2 + 16*z1*z3^2*z6*z7^2 + 16*z2*z3^2*z6*z7^2 -
 8*z3^3*z6*z7^2 - 6*z7^3 - 8*z1*z3^2*z7^3 -
   8*z2*z3^2*z7^3 + 4*z3^3*z7^3 - 36*z4*z8
h1 = z1*z5 + z2*z6 + z3*z7 - n*(1/(n+1))^(2/n)
h2 = 2*z1*z2*z3 + 2*z5*z6*z7 + 1
h3 = 16*z1*z2*z3*z5*z6*z7 - 1
cpop = [f; g1; g2; g3; h1; h2; h3]
opt,sol,data = cs_tssos_first(cpop, z, n+1, 5, numeq=3, CS=false, TS="block", QUIET=true, solve=false)
opt,sol,data = cs_tssos_higher!(data, TS="block", QUIET=true, solve=false)
opt,sol,data = cs_tssos_higher!(data, TS="block", QUIET=true, Gram=true)
@test opt ≈ -0.5625 atol = 1e-6

N = 6
@polyvar z[1:2N+2]
f = z[N+1]^2 + z[2N+2]^2
cons = Vector{typeof(f)}(undef, N-2)
for k = 1:N-2
    cons[k] = z[N+1]^2 + z[2N+2]^2 - sum(z[i]*z[j+k]*z[j+N+1]*z[i+k+N+1] for i = 1:N-k, j = 1:N-k)
end
opt,sol,data = cs_tssos_first([f; cons], z, N+1, 3, nb=N+1, CS=false, TS="block", ConjugateBasis=true, QUIET=true, Gram=true)
@test opt ≈ 1 atol = 1e-6

end
