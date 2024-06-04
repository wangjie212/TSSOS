using TSSOS
using DynamicPolynomials

supp = Vector{Vector{Vector{UInt16}}}[[[[], []], [[1], [1]], [[1], [2;2]], [[2;2], [1]]],
[[[2], []], [[], [2]]], [[[], []], [[1], [1]], [[1;1], []], [[], [1;1]]],
[[[], []], [[1], [1]], [[2], [2]]], [[[2], []], [[], [2]]]]
coe = [[3;-1;-0.5im;0.5im], [1;1], [-1;1;-0.25;-0.25], [-3;1;1], [im;-im]]
opt,sol,data = cs_tssos_first(supp, coe, 2, 2, numeq=3, QUIET=true, CS=false, TS=false, Mommat=true)
# optimum = 0.15508925176378527

supp = Vector{Vector{Vector{UInt16}}}[[[[1;1], [2]], [[2], [1;1]], [[1], [2]], [[2], [1]], [[2], [3]], [[3], [2]]],
[[[], []], [[1], [1]], [[2],[2]]], [[[], []], [[2], [2]], [[3], [3]]],
[[[], []], [[], [2]], [[2], []]]]
coe = [[1;1;1;1;1;1], [1;-1;-1], [1;-1;-1], [1;1;1]]
opt,sol,data = cs_tssos_first(supp, coe, 3, 2, numeq=1, ipart=false, QUIET=true, CS=true, TS=false, Mommat=true)
# optimum = -2.6821884022760147

@polyvar z[1:4]
pop = [z[1]+z[2]+z[3]+z[4], (1+im)*z[1]^2*z[4]+(1-im)*z[2]*z[3]^2, 1-z[1]*z[3]-z[2]*z[4]]
opt,sol,data = cs_tssos_first(pop, z, 2, 4, numeq=1, QUIET=true, CS=false, TS="block", solution=true)
# optimum = -2.8284270957706616
opt,sol,data = cs_tssos_higher!(data, QUIET=true, TS="block", solution=true)
# optimum = -2.748735909698819

@polyvar z[1:4]
pop = [z[1]*z[2] + z[3]*z[4], 1 - z[1]*z[3] - z[2]*z[4]]
opt,sol,data = cs_tssos_first(pop, z, 2, 7, ipart=false, numeq=1, QUIET=true, CS=false, TS=false, solution=true)
# optimum = -1.072774066479055

println("Run successfully!")
