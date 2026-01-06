using TSSOS

function Broydenbanded(n::Int)
    supp = [UInt16[]]
    coe = Float64[n]
    for i = 1:n
        jset = max(1,i-5):min(n,i+1)
        jset = setdiff(jset,i)
        push!(supp, [i], [i;i], [i;i;i], [i;i;i;i], [i;i;i;i;i;i])
        push!(coe, 4, 4, 10, 20, 25)
        for j in jset
            push!(supp, [j], [j;j], sort([i;j]), sort([i;j;j]), sort([i;i;i;j]), sort([i;i;i;j;j]))
            push!(coe, -2, -2, -4, -4, -10, -10)
            for k in jset
                push!(supp, sort([j;k]), sort([j;j;k]), sort([j;k;k]), sort([j;j;k;k]))
                push!(coe, 1, 1, 1, 1)
            end
        end
    end
    return arrange(poly(supp, coe))
end

function genRosenbrock(n::Int)
    supp = [UInt16[]]
    coe = Float64[n]
    for i = 2:n
        push!(supp, [i-1;i-1;i-1;i-1], [i-1;i-1;i], [i], [i;i])
        push!(coe, 100, -200, -2, 101)
    end
    return arrange(poly(supp, coe))
end

function chainedWood(n::Int)
    supp = [UInt16[]]
    coe = Float64[21*n-41]
    for i = 1:2:n-3
        push!(supp, [i], [i;i], [i;i;i;i], [i;i;i+1], [i+1], [i+1;i+1], [i+1;i+3], [i+2], [i+2;i+2], [i+2;i+2;i+2;i+2], [i+2;i+2;i+3], [i+3], [i+3;i+3])
        push!(coe, -2, 1, 100, -200, -40, 110.1, 19.8, -2, 1, 90, -180, -40, 100.1)
    end
    return arrange(poly(supp, coe))
end

function Broydentridiagonal(n::Int)
    supp = [UInt16[]]
    coe = Float64[n]
    push!(supp, [1;1], [1;1;1;1], [2;2], [1;1;1], [1;2], [1], [1;1;2], [2])
    push!(coe, 5, 4, 4, -12, -12, 6, 8, -4)
    for i = 2:n-1
        push!(supp, [i;i], [i;i;i;i], [i-1;i-1], [i+1;i+1], [i;i;i], [i-1;i], [i;i+1], [i], [i-1;i;i], [i;i;i+1], [i-1;i+1], [i-1], [i+1])
        push!(coe, 5, 4, 1, 4, -12, -6, -12, 6, 4, 8, 4, -2, -4)
    end
    push!(supp, [n;n], [n;n;n;n], [n-1;n-1], [n;n;n], [n-1;n], [n], [n-1;n;n], [n-1])
    push!(coe, 5, 4, 1, -12, -6, 6, 4, -2)
    return arrange(poly(supp, coe))
end

# Chained singular polynomial
function Chainedsingular(n::Int)
    supp = Vector{UInt16}[]
    coe = Float64[]
    for i = 1:2:n-3
        push!(supp, [i;i], [i;i+1], [i+1;i+1], [i+2;i+2], [i+3;i+3], [i+2;i+3], [i+1;i+1;i+1;i+1], [i+1;i+1;i+2;i+2], [i+2;i+2;i+2;i+2], [i+1;i+1;i+1;i+2], [i+1;i+1;i+2;i+2], [i+1;i+2;i+2;i+2], [i;i;i;i], [i;i;i+3;i+3], [i+3;i+3;i+3;i+3], [i;i;i;i+3], [i;i;i+3;i+3], [i;i+3;i+3;i+3])
        push!(coe, 1, 20 ,100, 5, 5, -10, 1, 16, 16, -8, 8, -32, 10, 4000, 100000, -400, 2000, -40000)
    end
    return arrange(poly(supp, coe))
end

# randomly generated examples
l = 2
b = 15
n = (b-5)*l+5
supp = Vector{Vector{Vector{UInt16}}}(undef, l+1)
coe = Vector{Vector{Float64}}(undef, l+1)
supp[1] = Vector{Vector{UInt16}}[]
for i = 1:l, j = 1:b
    push!(supp[1], ceil.(UInt16, ((b-5)*i-(b-5)).+rand(4).*b))
end
coe[1] = 2*rand(b*l).-1
for i = 1:l
    supp[i+1] = [[[]]; [[j;j] for j=(b-5)*i-(b-6):(b-5)*i+5]]
    coe[i+1] = [1; -ones(b)]
end
