using TSSOS
using DynamicPolynomials

supp = Vector{Vector{Vector{UInt16}}}[[[[], []], [[1], [1]], [[1], [2;2]], [[2;2], [1]]],
[[[2], []], [[], [2]]], [[[], []], [[1], [1]], [[1;1], []], [[], [1;1]]],
[[[], []], [[1], [1]], [[2], [2]]], [[[2], []], [[], [2]]]]
coe = [[3;-1;-0.5im;0.5im], [1;1], [-1;1;-0.25;-0.25], [-3;1;1], [im;-im]]
opt,sol,data = cs_tssos_first(supp, coe, 2, 2, numeq=3, QUIET=true, CS=false, TS=false)
# optimum = 0.428174

supp = Vector{Vector{Vector{UInt16}}}[[[[1;1], [2]], [[2], [1;1]], [[1], [2]], [[2], [1]], [[2], [3]], [[3], [2]]],
[[[], []], [[1], [1]], [[2],[2]]], [[[], []], [[2], [2]], [[3], [3]]],
[[[], []], [[], [2]], [[2], []]]]
coe = [[1;1;1;1;1;1], [1;-1;-1], [1;-1;-1], [1;1;1]]
opt,sol,data = cs_tssos_first(supp, coe, 3, 2, numeq=1, ipart=false, QUIET=true, CS=true, TS=false)
# optimum = -2.682188

@polyvar z[1:4]
pop = [z[1]+z[2]+z[3]+z[4], (1+im)*z[1]^2*z[4]+(1-im)*z[2]*z[3]^2, 1-z[1]*z[3]-z[2]*z[4]]
opt,sol,data = cs_tssos_first(pop, z, 2, 2, numeq=1, QUIET=true, CS=false, TS="block", solution=true)
# optimum = -2.74232
opt,sol,data = cs_tssos_higher!(data, QUIET=true, TS="block", solution=true)
# optimum = -2.74232

@polyvar z[1:4]
pop = [z[1]*z[2] + z[3]*z[4], 1 - z[1]*z[3] - z[2]*z[4]]
opt,sol,data = cs_tssos_first(pop, z, 2, 1, ipart=false, numeq=1, QUIET=true, CS=false, TS=false, solution=true)
# optimum = -1

println("Run successfully!")
