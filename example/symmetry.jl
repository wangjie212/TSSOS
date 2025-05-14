using PermutationGroups
using DynamicPolynomials
using TSSOS

@polyvar x[1:4]
f = sum(x) + sum(x.^2)
G = PermGroup([perm"(1,2,3,4)"]) # define the symmetry group
opt,data = tssos_symmetry([f], x, 1,  G)
# optimum = -1

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry([f], x, 2, G)
# optimum = -1.4174111

@polyvar x[1:3]
f = sum(x) + sum(x.^4) - 4*x[1]*x[2]*x[3]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry([f], x, 2, G)
# optimum = -2.1129138

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry(pop, x, 2, G)
# optimum = -1.3987174

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry(pop, x, 2, G, numeq=1)
# optimum = -1.3987174

@polyvar x[1:3]
f = sum(x.^2) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry_first(pop, x, 2, G, numeq=0, TS="block")
# optimum = 0
opt,data = tssos_symmetry_higher!(data, TS="block")
# optimum = 0
