using DynamicPolynomials
using TSSOS
using Random

function cbasis(z)
    basis = Poly{Int}[1]
    for i = 1:length(z)
        push!(basis, z[i])
    end
    for i = 1:length(z), j = i:length(z)
        push!(basis, z[i]*z[j])
    end
    return basis
end

# minimizing a random complex quartic polynomial over the unit sphere
Random.seed!(1)
n = 10
@complex_polyvar z[1:n]
basis = cbasis(z)
P = randn(length(basis), length(basis))
Q = randn(length(basis), length(basis))
f = basis'*((P+P')/2+im*(Q-Q')/2)*basis
@time opt,sol,data = complex_tssos([f; z'*z - 1], z, 2, numeq=1, QUIET=true, solve=true, TS=false)
@time opt,sol,data = complex_tssos([f; z'*z - 1], z, 3, numeq=1, QUIET=true, solve=true, TS=false)

using LinearAlgebra
Random.seed!(1)
n = 1
@complex_polyvar x
basis = cbasis([x])
P = rand(length(basis), length(basis))
# Q = rand(length(basis), length(basis))
# A = (P+P')/2 + im*(Q-Q')/2 + 3*I(length(basis))
A = (P+P')/2 + I(length(basis))
A = [0 1 1; 1 1 0; 1 0 1]
f = basis'*A*basis + 0.5
@time opt,sol,data = complex_tssos([f], [x], 2, QUIET=true, normality=0, ConjugateBasis=false, TS=false)
@time opt,sol,data = complex_tssos([f], [x], 2, QUIET=true, ConjugateBasis=true, Gram=true, TS=false)

# minimizing a random complex quartic polynomial with unit-norm variables
Random.seed!(1)
n = 5
@complex_polyvar z[1:n]
basis = cbasis(z)
P = rand(length(basis), length(basis))
Q = rand(length(basis), length(basis))
pop = [basis'*((P+P')/2+im*(Q-Q')/2)*basis]
# pop = [basis'*((P+P')/2)*basis]
@time opt1,sol,data = complex_tssos(pop, z, 2, nb=n, QUIET=true, CS=false, ConjugateBasis=false)
@time opt2,sol,data = complex_tssos(pop, z, 2, nb=n, QUIET=true, CS=false, ConjugateBasis=true)

# minimizing large-scale randomly generated complex QCQPs
Random.seed!(1)
b = 10
l = 120
n = (b-5)*l + 5
pop = Vector{poly{ComplexF64}}(undef, l+1)
supp = Tuple{Vector{UInt16}, Vector{UInt16}}[]
coe = ComplexF64[]
for i = 1:l, j = 1:2b
    if j < 3
        p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
        push!(supp, tuple(p1, []), tuple([], p1))
    elseif j < 7
        p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
        p2 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
        push!(supp, tuple(p1, p2), tuple(p2, p1))
    elseif j < 13
        p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(1).*b))
        p2 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(2).*b))
        push!(supp, tuple(p1, p2), tuple(p2, p1))
    else
        p1 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(2).*b))
        p2 = sort(ceil.(UInt16, ((b-5)*i-(b-5)).+rand(2).*b))
        push!(supp, tuple(p1, p2), tuple(p2, p1))
    end
    c1 = 2*rand(1)[1]-1
    c2 = 2*rand(1)[1]-1
    push!(coe, c1+c2*im, c1-c2*im)
end
pop[1] = arrange(cpoly(supp, coe))
for i = 1:l
    pop[i+1] = cpoly(Tuple{Vector{UInt16}, Vector{UInt16}}[[tuple([], [])]; [tuple([j], [j]) for j=(b-5)*i-(b-6):(b-5)*i+5]], [1; -ones(b)])
end

@time opt,sol,data = complex_cs_tssos(pop, n, 2, numeq=l, TS="MD")

# AC-OPF problem
include("D:/Programs/TSSOS/example/modelopf.jl")
cd("D:/Programs/PolyOPF/pglib")
silence()

case = "pglib_opf_case30_as"
AC = 803.13
opfdata = parse_file(case * ".m")
model = pop_opf_com(opfdata, normal=true, AngleCons=true, LineLimit=true)
n = model.n
m = model.m
numeq = model.numeq
pop = model.pop
mc = maximum(abs.(pop[1].coe))
pop[1].coe = pop[1].coe./mc

t = @elapsed begin
opt,_,popd = complex_cs_tssos(pop, n, "min", numeq=numeq, CS="MF", TS="block", QUIET=true, MomentOne=false)
end
opt *= mc
maxc = maximum(popd.cliquesize) # maximal clique size
mb = 2*maximum(maximum.([maximum.(bs) for bs in popd.blocksize])) # maximal block size
gap = 100*(AC-opt)/AC # optimality gap
println("n = $n, m = $m")
println("mc = $maxc, opt = $opt, time = $t, mb = $mb, gap = $gap%")
