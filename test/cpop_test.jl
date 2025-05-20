using TSSOS
using DynamicPolynomials
using Test

@testset begin

@complex_polyvar z[1:2]
pop = [3-z[1]*conj(z[1])-0.5im*z[1]*conj(z[2])^2+0.5im*z[2]^2*conj(z[1]), z[2]+conj(z[2]), 
-1+z[1]*conj(z[1])-0.25*z[1]^2-0.25*conj(z[1])^2, -3+z[1]*conj(z[1])+z[2]*conj(z[2]), im*z[2]-im*conj(z[2])]
opt,sol,data = complex_tssos_first(pop, z, 2, numeq=3, QUIET=true, TS="block", solution=true, Gram=true)
@test opt ≈ -0.4165019 atol = 1e-6
opt,sol,data = complex_tssos_higher!(data, QUIET=true, TS="block", solution=true, Gram=true)
@test opt ≈ 0.4281746 atol = 1e-6

@complex_polyvar z[1:3]
pop = [z[1]^2*conj(z[2])+z[2]*conj(z[1])^2+z[1]*conj(z[2])+z[2]*conj(z[1])+z[2]*conj(z[3])+z[3]*conj(z[2]),
1-z[1]*conj(z[1])-z[2]*conj(z[2]), 1-z[2]*conj(z[2])-z[3]*conj(z[3]), 1+z[2]+conj(z[2])]
opt,sol,data = complex_cs_tssos_first(pop, z, 2, numeq=1, QUIET=true, TS="block", solution=true, Gram=true)
@test opt ≈ -2.682188 atol = 1e-6

@complex_polyvar z[1:3]
f = 0.5z[1]*conj(z[2])+0.5z[1]*conj(z[3])+0.5z[2]*conj(z[1])+0.25z[2]*conj(z[2])+0.25z[2]*conj(z[3])+0.5z[3]*conj(z[1])+
0.25z[3]*conj(z[2])+z[1]+z[2]+z[3]+conj(z[1])+conj(z[2])+conj(z[3])
pop = [f, 1-z[1]*conj(z[1]), 1-z[2]*conj(z[2]), 1-z[3]*conj(z[3])]
opt,sol,data = complex_cs_tssos_first(pop, z, 2, numeq=3, QUIET=true, TS=false, solution=true, Gram=true)
@test opt ≈ -3.75 atol = 1e-6

@complex_polyvar z[1:2]
pop = [z[1]+z[2]+conj(z[1])+conj(z[2]), (1+im)*z[1]^2*conj(z[2])+(1-im)*z[2]*conj(z[1])^2, 1 - z'*z]
opt,sol,data = complex_cs_tssos_first(pop, z, 2, numeq=1, QUIET=true, TS="block", ConjugateBasis=true, solution=true, Gram=true)
@test opt ≈ -2.7423247 atol = 1e-6
opt,sol,data = complex_cs_tssos_first(pop, z, 2, numeq=1, QUIET=true, TS="block", ConjugateBasis=true, normality=2, solution=true, Gram=true)
@test opt ≈ -2.7423247 atol = 1e-6
opt,sol,data = complex_cs_tssos_first(pop, z, 2, numeq=1, QUIET=true, TS="block", solution=true, Gram=true)
@test opt ≈ -2.7423247 atol = 1e-6
opt,sol,data = complex_cs_tssos_higher!(data, QUIET=true, TS="block", solution=true, Gram=true)
@test opt ≈ -2.7423247 atol = 1e-6

@complex_polyvar z[1:2]
pop = [z[1]*z[2] + conj(z[1])*conj(z[2]), 1 - z'*z]
opt,sol,data = complex_cs_tssos_first(pop, z, 2, numeq=1, QUIET=true, TS="block", solution=true, Gram=true)
@test opt ≈ -1 atol = 1e-6

@complex_polyvar x[1:3]
pop = [1, x[1]^2-2x[1]*x[3]+conj(x[1]^2)-2conj(x[1]*x[3])+10, im*(x[1]^2-2x[1]*x[3])-im*(conj(x[1]^2)-2conj(x[1]*x[3])), 
x[1]*x[2]^2+x[2]*x[3]+conj(x[1]*x[2]^2+x[2]*x[3])+2, im*(x[1]*x[2]^2+x[2]*x[3])-im*conj(x[1]*x[2]^2+x[2]*x[3])+2, 
3x[2]^2-8x[1]*x[3]+conj(3x[2]^2-8x[1]*x[3]), im*(3x[2]^2-8x[1]*x[3])-im*conj(3x[2]^2-8x[1]*x[3])]
opt,sol,data = complex_tssos_first(pop, x, 3, numeq=6, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ 1 atol = 1e-6

@complex_polyvar z[1:3]
pop = [1, z[1]^2+z[2]+z[3]+conj(z[1]^2+z[2]+z[3])+2, im*(z[1]^2+z[2]+z[3])-im*conj(z[1]^2+z[2]+z[3]),
z[2]^2+z[1]+z[3]+conj(z[2]^2+z[1]+z[3])+2, im*(z[2]^2+z[1]+z[3])-im*conj(z[2]^2+z[1]+z[3]),
z[3]^2+z[2]+z[1]+conj(z[3]^2+z[2]+z[1])+2, im*(z[3]^2+z[2]+z[1])-im*conj(z[3]^2+z[2]+z[1])]
opt,sol,data = complex_tssos_first(pop, z, 3, numeq=6, TS="block", Gram=true, solution=true, QUIET=true)
@test opt ≈ 1 atol = 1e-6

@complex_polyvar z1 z2 z3
z = tuple(z1,z2,z3)
f = - z3*conj(z3)
g1 = z1^2*conj(z1)^2 - z1^2 - conj(z1)^2 - 4*z3*conj(z3) + 1
g2 = z2^2*conj(z2)^2 - z2^2 - conj(z2)^2 - 4*z3*conj(z3) + 1
h1 = - z1*conj(z1) - z2*conj(z2) + 2/3
h2 = 3*z1*z2 + 3*conj(z1)*conj(z2) - 2
h3 = 9*z1*z2*conj(z1)*conj(z2) - 1
cpop = [f; g1; g2; h1; h2; h3]
opt,sol,data = complex_cs_tssos_first(cpop, z, 2, numeq=3, CS=false, TS="block", solution=true, QUIET=true, Gram=true)
@test opt ≈ -0.4444444 atol = 1e-6

n = 3
@complex_polyvar z1 z2 z3 z4
z = tuple(z1,z2,z3,z4)
f = - z4*conj(z4)
g1 = 9 - 6*z1^3 + 12*z1^2*z2 + 12*z1^2*z3 - 6*conj(z1)^3 + 4*z1^3*conj(z1)^3 -
 8*z1^2*z2*conj(z1)^3 - 8*z1^2*z3*conj(z1)^3 + 12*conj(z1)^2*conj(z2) -
   8*z1^3*conj(z1)^2*conj(z2) + 16*z1^2*z2*conj(z1)^2*conj(z2) + 16*z1^2*z3*conj(z1)^2*conj(z2) +
 12*conj(z1)^2*conj(z3) - 8*z1^3*conj(z1)^2*conj(z3) + 16*z1^2*z2*conj(z1)^2*conj(z3) +
   16*z1^2*z3*conj(z1)^2*conj(z3) - 36*z4*conj(z4)
g2 = 9 + 12*z1*z2^2 - 6*z2^3 + 12*z2^2*z3 + 12*conj(z1)*conj(z2)^2 +
 16*z1*z2^2*conj(z1)*conj(z2)^2 - 8*z2^3*conj(z1)*conj(z2)^2 + 16*z2^2*z3*conj(z1)*conj(z2)^2 -
   6*conj(z2)^3 - 8*z1*z2^2*conj(z2)^3 + 4*z2^3*conj(z2)^3 - 8*z2^2*z3*conj(z2)^3 +
 12*conj(z2)^2*conj(z3) + 16*z1*z2^2*conj(z2)^2*conj(z3) - 8*z2^3*conj(z2)^2*conj(z3) +
   16*z2^2*z3*conj(z2)^2*conj(z3) - 36*z4*conj(z4)
g3 = 9 + 12*z1*z3^2 + 12*z2*z3^2 - 6*z3^3 + 12*conj(z1)*conj(z3)^2 +
 16*z1*z3^2*conj(z1)*conj(z3)^2 + 16*z2*z3^2*conj(z1)*conj(z3)^2 - 8*z3^3*conj(z1)*conj(z3)^2 +
   12*conj(z2)*conj(z3)^2 + 16*z1*z3^2*conj(z2)*conj(z3)^2 + 16*z2*z3^2*conj(z2)*conj(z3)^2 -
 8*z3^3*conj(z2)*conj(z3)^2 - 6*conj(z3)^3 - 8*z1*z3^2*conj(z3)^3 -
   8*z2*z3^2*conj(z3)^3 + 4*z3^3*conj(z3)^3 - 36*z4*conj(z4)
h1 = z1*conj(z1) + z2*conj(z2) + z3*conj(z3) - n*(1/(n+1))^(2/n)
h2 = 2*z1*z2*z3 + 2*conj(z1)*conj(z2)*conj(z3) + 1
h3 = 16*z1*z2*z3*conj(z1)*conj(z2)*conj(z3) - 1
cpop = [f; g1; g2; g3; h1; h2; h3]
opt,sol,data = complex_cs_tssos_first(cpop, z, 5, numeq=3, CS=false, TS="block", QUIET=true, solve=false)
opt,sol,data = complex_cs_tssos_higher!(data, TS="block", QUIET=true, solve=false)
opt,sol,data = complex_cs_tssos_higher!(data, TS="block", QUIET=true, Gram=true)
@test opt ≈ -0.5625 atol = 1e-6

N = 6
@complex_polyvar z[1:N+1]
f = z[N+1]^2 + conj(z[N+1])^2
cons = Vector{typeof(f)}(undef, N-2)
for k = 1:N-2
    cons[k] = z[N+1]^2 + conj(z[N+1])^2 - sum(z[i]*z[j+k]*conj(z[j])*conj(z[i+k]) for i = 1:N-k, j = 1:N-k)
end
opt,sol,data = complex_cs_tssos_first([f; cons], z, 3, nb=N+1, CS=false, TS="block", ConjugateBasis=true, QUIET=true, Gram=true)
@test opt ≈ 1 atol = 1e-6

end
