include("D:\\Programs\\TSSOS\\src\\TSSOS.jl")
using .TSSOS
using DynamicPolynomials

b = 10
l = 1
n = (b-5)*l + 5
supp = Vector{Vector{Vector{Vector{UInt16}}}}(undef, l+1)
# coe = Vector{Vector{ComplexF64}}(undef, l+1)
coe = Vector{Vector{Float64}}(undef, l+1)
supp[1] = []
coe[1] = []
for i = 1:l,j = 1:2b
    if j < 3
        p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
        push!(supp[1], [p1, []], [[], p1])
    elseif j < 7
        p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
        p2 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
        push!(supp[1], [p1, p2], [p2, p1])
    elseif j < 13
        p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
        p2 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(2).*b))
        push!(supp[1], [p1, p2], [p2, p1])
    else
        p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(2).*b))
        p2 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(2).*b))
        push!(supp[1], [p1, p2], [p2, p1])
    end
    c1 = 2*rand(1)[1]-1
    # c2 = 2*rand(1)[1]-1
    # push!(coe[1], c1+c2*im, c1-c2*im)
    push!(coe[1], c1, c1)
end
supp[1],coe[1] = resort(supp[1],coe[1])
for i = 1:l
    supp[i+1] = [[[[], []]]; [[[j], [j]] for j=(b-5)*i-(b-6):(b-5)*i+5]]
    coe[i+1] = [1; -ones(b)]
end

@polyvar x[1:4]
pop = [3-x[1]^2-x[3]^2+x[3]*x[2]^2-2x[1]*x[2]*x[4]-x[3]*x[4]^2, x[2], x[1]^2+3x[3]^2-2, x[4], x[1]^2+x[2]^2+x[3]^2+x[4]^2-3]
# pop = [3-x[1]^2-x[3]^2+x[1]*x[2]^2+2x[2]*x[3]*x[4]-x[1]*x[4]^2, x[2], x[1]^2+3x[3]^2-2, x[4], x[1]^2+x[2]^2+x[3]^2+x[4]^2-3]
@time begin
opt,sol,data = tssos_first(pop, x, 2, numeq=3, TS=false, QUIET=true, quotient=false)
# opt,sol,data = tssos_higher!(data, TS="block", QUIET=true)
end
sol,ub,gap = refine_sol(opt, rand(4), data, QUIET=true)

supp = Vector{Vector{Vector{UInt16}}}[[[[], []], [[1], [1]], [[1], [2;2]], [[2;2], [1]]],
[[[2], []], [[], [2]]], [[[], []], [[1], [1]], [[1;1], []], [[], [1;1]]],
[[[2;2], []], [[], [2;2]], [[2], [2]]], [[[], []], [[1], [1]], [[2], [2]]]]
coe = [[3;-1;0.5;0.5], [1;1], [-1;1;-0.25;-0.25], [1;1;-2], [-3;1;1]]

@time begin
opt,sol,data = cs_tssos_first(supp, coe, 2, 4, numeq=3, QUIET=true, CS=false, TS=false, ipart=false)
end
@time begin
opt,sol,data = cs_tssos_first(supp, coe, n, 2, QUIET=true, ipart=false, TS="MD")
end
