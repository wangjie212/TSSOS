using TSSOS
using DynamicPolynomials

supp = Vector{Vector{Vector{UInt16}}}[[[[], []], [[1], [1]], [[1], [2;2]], [[2;2], [1]]],
[[[2], []], [[], [2]]], [[[], []], [[1], [1]], [[1;1], []], [[], [1;1]]],
[[[], []], [[1], [1]], [[2], [2]]], [[[2], []], [[], [2]]]]
coe = [[3;-1;-0.5im;0.5im], [1;1], [-1;1;-0.25;-0.25], [-3;1;1], [im;-im]]

@time begin
opt,sol,data = cs_tssos_first(supp, coe, 2, 3, numeq=3, QUIET=false, CS=false, TS=false, Gram=true, Mommat=true)
end

supp = [[[[1;1], [2]], [[2], [1;1]], [[1], [2]], [[2], [1]], [[2], [3]], [[3], [2]]],
[[[], []], [[1], [1]], [[2],[2]]], [[[], []], [[2], [2]], [[3], [3]]],
[[[], []], [[], [2]], [[2], []]]]
coe = [[1;1;1;1;1;1], [1;-1;-1], [1;-1;-1], [1;1;1]]

@time begin
opt,sol,data = cs_tssos_first(supp, coe, 3, 2, numeq=1, QUIET=false, CS=false, TS=false, Gram=true, Mommat=true)
end