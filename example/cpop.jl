using TSSOS
using DynamicPolynomials
using LinearAlgebra

b = 10
l = 1
n = (b-5)*l + 5
supp = Vector{Vector{Vector{Vector{UInt16}}}}(undef, l+1)
# coe = Vector{Vector{ComplexF64}}(undef, l+1)
coe = Vector{Vector{Float64}}(undef, l+1)
supp[1] = []
coe[1] = []
for i = 1:l,j = 1:2b
    # if j < 3
    #     p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
    #     push!(supp[1], [p1, []], [[], p1])
    # else
        # if j < 7
        p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
        p2 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
        push!(supp[1], [p1, p2], [p2, p1])
    # elseif j < 13
    #     p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
    #     p2 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(2).*b))
    #     push!(supp[1], [p1, p2], [p2, p1])
    # else
    #     p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(2).*b))
    #     p2 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(2).*b))
    #     push!(supp[1], [p1, p2], [p2, p1])
    # end
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

function bfind(A, l, a)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if ndims(A) == 2
            temp = A[:, mid]
        else
            temp = A[mid]
        end
        if temp == a
           return mid
        elseif temp < a
           low = mid+1
        else
           high = mid-1
        end
    end
    return 0
end

io = open("D:\\project\\cpop\\random\\supp_$n.txt", "w")
writedlm(io, convert(Vector{Vector{Vector{Int}}}, supp[1]), ';')
close(io)
io = open("D:\\project\\cpop\\random\\coe_$n.txt", "w")
writedlm(io, coe[1])
close(io)

@polyvar x[1:2]
pop = [-2x[1]*x[2]+x[1]^2+x[2]^2, x[1]*x[2]-x[1]-x[2]]
opt,sol,data = cs_tssos_first(pop, x, 1, 2, numeq=1, CS=false, TS=false, ipart=false, QUIET=true, Mommat=true)

@polyvar x[1:4]
# pop = [3-x[1]^2-x[3]^2+x[3]*x[2]^2-2x[1]*x[2]*x[4]-x[3]*x[4]^2, x[2], x[1]^2+3x[3]^2-2, x[4], x[1]^2+x[2]^2+x[3]^2+x[4]^2-3]
pop = [3-x[1]^2-x[3]^2+x[1]*x[2]^2+2x[2]*x[3]*x[4]-x[1]*x[4]^2, x[2], x[1]^2+3x[3]^2-2, x[4], x[1]^2+x[2]^2+x[3]^2+x[4]^2-3]
@time begin
opt,sol,data = cs_tssos_first(pop, x, 2, numeq=3, CS=false, TS=false, QUIET=true, Mommat=true)
end

supp = Vector{Vector{Vector{UInt16}}}[[[[], []], [[1], [1]], [[1], [2;2]], [[2;2], [1]]],
[[[2], []], [[], [2]]], [[[], []], [[1], [1]], [[1;1], []], [[], [1;1]]],
[[[2;2], []], [[], [2;2]], [[2], [2]]], [[[], []], [[1], [1]], [[2], [2]]]]
coe = [[3;-1;0.5;0.5], [1;1], [-1;1;-0.25;-0.25], [1;1;-2], [-3;1;1]]

@time begin
opt,sol,data = cs_tssos_first(supp, coe, b, 1, numeq=1, QUIET=true, CS=false, TS=false, ipart=false, Mommat=true)
end

M = Matrix{Float64}(data.Mmatrix[1][1])
F = eigen(M)
sol = F.vectors[2:end, end-1:end]
ceval(supp[1], coe[1], sol)

function ceval(supp, coe, sol)
    val = 0
    for i = 1:length(supp)
        val += coe[i]*prod(sol[supp[i][1],1] + sol[supp[i][1],2]*im)*prod(sol[supp[i][2],1] - sol[supp[i][2],2]*im)
    end
    return val
end

n = 300
@polyvar z[1:2n]
Q = randn(n, n)
Q = (Q+Q')/2
io = open("D:\\project\\cpop\\random\\Q_$n.txt", "w")
writedlm(io, Q)
close(io)
f = z[1:n]'*Q*z[n+1:end]
h = sum(z[i]*z[i+n] for i = 1:n) - 1
@time begin
opt,sol,data = cs_tssos_first([f; h], z, n, 1, numeq=1, QUIET=true, CS=false, TS=false, ipart=false)
end
