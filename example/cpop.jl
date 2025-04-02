using DynamicPolynomials
using TSSOS
using Random

function cbasis(z)
    basis = Monomial{true}[1]
    for i = 1:length(z)
        push!(basis, z[i])
    end
    # for i = 1:length(z), j = i:length(z)
    #     push!(basis, z[i]*z[j])
    # end
    return basis
end

# minimizing a random complex quartic polynomial over the unit sphere
Random.seed!(1)
n = 5
@polyvar z[1:2n]
basis1 = cbasis(z[1:n])
basis2 = cbasis(z[n+1:2n])
P = randn(length(basis1), length(basis1))
Q = randn(length(basis1), length(basis1))
f = basis2'*((P+P')/2+im*(Q-Q')/2)*basis1
h = sum(z[i]*z[i+n] for i = 1:n) - 1
@time begin
opt,sol,data = cs_tssos_first([f; h], z, n, 2, numeq=1, QUIET=false, solve=true, CS=false, TS=false)
end
# println(length(data.basis[1][1]))
@time begin
opt,sol,data = cs_tssos_first([f; h], z, n, 3, numeq=1, QUIET=false, solve=true, CS=false, TS=false)
end
# println(length(data.basis[1][1]))

# minimizing a random complex quartic polynomial with unit-norm variables
for i = 1:500
Random.seed!(i)
n = 3
@polyvar z[1:2n]
basis1 = cbasis(z[1:n])
basis2 = cbasis(z[n+1:2n])
P = rand(length(basis1), length(basis1))
Q = rand(length(basis1), length(basis1))
pop = [basis2'*((P+P')/2+im*(Q-Q')/2)*basis1]
# pop = [basis2'*((P+P')/2)*basis1]
@time begin
opt1,sol,data = cs_tssos_first(pop, z, n, 2, nb=n, QUIET=true, CS=false, TS=false)
end
# println(length(data.basis[1][1]))
@time begin
opt2,sol,data = cs_tssos_first(pop, z, n, 3, nb=n, QUIET=true, CS=false, TS=false)
end
if abs(opt1-opt2) > 1e-6
    println(i)
    break
end
end
# println(length(data.basis[1][1]))

# minimizing large-scale randomly generated complex QCQPs
Random.seed!(1)
b = 10
l = 120
n = (b-5)*l + 5
supp = Vector{Vector{Vector{Vector{UInt16}}}}(undef, l+1)
coe = Vector{Vector{ComplexF64}}(undef, l+1)
supp[1] = []
coe[1] = []
for i = 1:l, j = 1:2b
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
    c2 = 2*rand(1)[1]-1
    push!(coe[1], c1+c2*im, c1-c2*im)
end
supp[1],coe[1] = TSSOS.resort(supp[1],coe[1])
for i = 1:l
    supp[i+1] = [[[[], []]]; [[[j], [j]] for j=(b-5)*i-(b-6):(b-5)*i+5]]
    coe[i+1] = [1; -ones(b)]
end

@time begin
opt,sol,data = cs_tssos_first(supp, coe, n, 2, numeq=l, TS="MD")
end
# mb = 2*maximum(maximum.(data.sb)) # maximal block size
# println("n = $n, time = $time, mb = $mb")

# AC-OPF problem
include("D:/Programs/TSSOS/example/modelopf.jl")
cd("D:/Programs/PolyOPF/pglib")
silence()

case = "pglib_opf_case14_ieee"
AC = 2178.08
opfdata = parse_file(case * ".m")
model = pop_opf_com_vol(opfdata, normal=true, AngleCons=true, LineLimit=true)
n = model.n
m = model.m
numeq = model.numeq
supp = model.supp
coe = model.coe
mc = maximum(abs.(coe[1]))
coe[1] = coe[1]./mc

time = @elapsed begin
opt,sol,popd = cs_tssos_first(supp, coe, n, "min", numeq=numeq, tune=true, CS="MF", TS="block", MomentOne=false)
end
opt *= mc
maxc = maximum(popd.cliquesize) # maximal clique size
mb = 2*maximum(maximum.([maximum.(popd.blocksize[i]) for i = 1:popd.cql])) # maximal block size
gap = 100*(AC-opt)/AC # optimality gap
println("n = $n, m = $m")
println("mc = $maxc, opt = $opt, time = $time, mb = $mb, gap = $gap%")