using DynamicPolynomials

# Smale’s Mean Value Conjecture
n = 2
supp = Vector{Vector{Vector{UInt16}}}[[[[3], [3]]], [[[1;1], [1;1]], [[1;1], []], [[], [1;1]], [[3], [3]], [[], []]],
[[[2;2], [2;2]], [[2;2], []], [[], [2;2]], [[3], [3]], [[], []]],
[[[1;2], [1;2]], [[1;2], []], [[], [1;2]], [[], []]],
[[[1], [1]], [[2], [2]], [[3], [3]], [[], []]]]
coe = [[-1], [1;-1;-1;-4;1], [1;-1;-1;-4;1], [9;-3;-3;1], [-1;-1;-1;10/9]]

opt,sol,data = cs_tssos_first(supp, coe, 3, 3, numeq=2, CS=false, TS="block", ipart=false, QUIET=true)
opt,sol,data = cs_tssos_higher!(data, TS="block", ipart=false, QUIET=true)

n = 3
@polyvar z1 z2 z3 z4 z5 z6 z7 z8
z = tuple(z1,z2,z3,z4,z5,z6,z7,z8)
f = - z4*z8
g1 = 9 - 6*z1^3 + 12*z1^2*z2 + 12*z1^2*z3 - 6*z5^3 + 4*z1^3*z5^3 -
 8*z1^2*z2*z5^3 - 8*z1^2*z3*z5^3 + 12*z5^2*z6 -
   8*z1^3*z5^2*z6 + 16*z1^2*z2*z5^2*z6 + 16*z1^2*z3*z5^2*z6 +
 12*z5^2*z7 - 8*z1^3*z5^2*z7 + 16*z1^2*z2*z5^2*z7 +
   16*z1^2*z3*z5^2*z7 - 36z4*z8
g2 = 9 + 12*z1*z2^2 - 6*z2^3 + 12*z2^2*z3 + 12*z5*z6^2 +
 16*z1*z2^2*z5*z6^2 - 8*z2^3*z5*z6^2 + 16*z2^2*z3*z5*z6^2 -
   6*z6^3 - 8*z1*z2^2*z6^3 + 4*z2^3*z6^3 - 8*z2^2*z3*z6^3 +
 12*z6^2*z7 + 16*z1*z2^2*z6^2*z7 - 8*z2^3*z6^2*z7 +
   16*z2^2*z3*z6^2*z7 - 36z4*z8
g3 = 9 + 12*z1*z3^2 + 12*z2*z3^2 - 6*z3^3 + 12*z5*z7^2 +
 16*z1*z3^2*z5*z7^2 + 16*z2*z3^2*z5*z7^2 - 8*z3^3*z5*z7^2 +
   12*z6*z7^2 + 16*z1*z3^2*z6*z7^2 + 16*z2*z3^2*z6*z7^2 -
 8*z3^3*z6*z7^2 - 6*z7^3 - 8*z1*z3^2*z7^3 -
   8*z2*z3^2*z7^3 + 4*z3^3*z7^3 - 36z4*z8
h1 = z1*z5 + z2*z6 + z3*z7 + z4*z8 - n*(1/(n+1))^(2/n) - (n/(n+1))^2
h2 = 16*z1*z2*z3*z5*z6*z7 + 4*z1*z2*z3 + 4*z5*z6*z7 + 1
pop = [f; g1; g2; g3; h1; h2]

opt,sol,data = cs_tssos_first(pop, z, n+1, 7, numeq=2, CS=false, TS="block", ipart=false, QUIET=true)
opt,sol,data = cs_tssos_higher!(data, TS="block", ipart=false, QUIET=true)
println(sqrt(-opt))

# Mordell Inequality Conjecture
n = 2
@polyvar z1 z2 z3 z4
z = tuple(z1,z2,z3,z4)
f = 4*z1^3*z3^3 + 6*z1^2*z2*z3^3 - 6*z1*z2^2*z3^3 - 4*z2^3*z3^3 +
 6*z1^3*z3^2*z4 + 9*z1^2*z2*z3^2*z4 -
   9*z1*z2^2*z3^2*z4 - 6*z2^3*z3^2*z4 - 6*z1^3*z3*z4^2 -
 9*z1^2*z2*z3*z4^2 + 9*z1*z2^2*z3*z4^2 + 6*z2^3*z3*z4^2 -
   4*z1^3*z4^3 - 6*z1^2*z2*z4^3 + 6*z1*z2^2*z4^3 + 4*z2^3*z4^3
h1 = 2*z1*z3 + 2*z2*z4 + z1*z4 + z2*z3 -3
# f = z1^2*z2*z4^2*z5 - z1*z2^2*z4^2*z5 - z1^2*z3*z4^2*z5 +
#  z2^2*z3*z4^2*z5 + z1*z3^2*z4^2*z5 - z2*z3^2*z4^2*z5 -
#    z1^2*z2*z4*z5^2 + z1*z2^2*z4*z5^2 + z1^2*z3*z4*z5^2 -
#  z2^2*z3*z4*z5^2 - z1*z3^2*z4*z5^2 + z2*z3^2*z4*z5^2 -
#    z1^2*z2*z4^2*z6 + z1*z2^2*z4^2*z6 + z1^2*z3*z4^2*z6 -
#  z2^2*z3*z4^2*z6 - z1*z3^2*z4^2*z6 + z2*z3^2*z4^2*z6 +
#    z1^2*z2*z5^2*z6 - z1*z2^2*z5^2*z6 - z1^2*z3*z5^2*z6 +
#  z2^2*z3*z5^2*z6 + z1*z3^2*z5^2*z6 - z2*z3^2*z5^2*z6 +
#    z1^2*z2*z4*z6^2 - z1*z2^2*z4*z6^2 - z1^2*z3*z4*z6^2 +
#  z2^2*z3*z4*z6^2 + z1*z3^2*z4*z6^2 - z2*z3^2*z4*z6^2 -
#    z1^2*z2*z5*z6^2 + z1*z2^2*z5*z6^2 + z1^2*z3*z5*z6^2 -
#  z2^2*z3*z5*z6^2 - z1*z3^2*z5*z6^2 + z2*z3^2*z5*z6^2
# h1 = z1*z4 + z2*z5 + z3*z6 - 3
# h2 = z2*z4 + z3*z4 + z1*z5 + z3*z5 + z1*z6 + z2*z6 + 3
pop = [-f; h1]

opt,sol,data = cs_tssos_first(pop, z, n, 8, numeq=1, CS=false, TS="block", ipart=false, QUIET=true)
opt,sol,data = cs_tssos_higher!(data, TS="block", ipart=false, QUIET=true)
