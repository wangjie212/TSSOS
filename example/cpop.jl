include("D:\\Programs\\TSSOS\\src\\TSSOS.jl")
using .TSSOS

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

@time begin
opt,sol,data = cs_tssos_first(supp, coe, n, 2, QUIET=true, TS="MD")
end
@time begin
opt,sol,data = cs_tssos_first(supp, coe, n, 2, QUIET=true, ipart=false, TS="MD")
end
