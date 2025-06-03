using TSSOS
using DynamicPolynomials
using Random
using LinearAlgebra

# Minimizing a random complex quadratic polynomial with unit-norm variables
Random.seed!(1)
n = 10
@complex_polyvar z[1:n]
P = rand(n+1, n+1)
Q = rand(n+1, n+1)
pop = [[1; z]'*((P+P')/2+im*(Q-Q')/2)*[1; z]]
@time begin
opt,sol,data = complex_tssos_first(pop, z, 1, nb=n, QUIET=true, TS=false, normality=0)
end
@time begin
opt,sol,data = complex_tssos_first(pop, z, 2, nb=n, QUIET=true, TS=false, normality=0)
end
@time begin
opt,sol,data = complex_tssos_first(pop, z, 1, nb=n, QUIET=true, TS=false, normality=1)
end
# println(sum(eigvals(convert.(ComplexF64, data.moment[1][1])) .> 1e-4))

@polyvar x[1:2n]
rf = pop[1](z=>x[1:n]+im*x[n+1:2n])
rpop = [real.(coefficients(rf))'*monomials(rf)]
for i = 1:n
    push!(rpop, 1 - x[i]^2 - x[i+n]^2)
end
@time begin
opt,sol,data = tssos_first(rpop, x, 2, numeq=n, GroebnerBasis=true, QUIET=true, TS=false)
end

function basis(x)
    basis = Poly{Int}[1]
    push!(basis, x...)
    for i = 1:length(x), j = i:length(x)
        push!(basis, x[i]*x[j])
    end
    return basis
end

# Minimizing a random complex quartic polynomial on a unit sphere
Random.seed!(1)
n = 10
@complex_polyvar z[1:n]
cb = basis(z)
lcb = length(cb)
P = rand(lcb, lcb)
Q = rand(lcb, lcb)
pop = [cb'*((P+P')/2+im*(Q-Q')/2)*cb]
push!(pop, 1 - z'*z)
@time begin
opt,sol,data = complex_tssos_first(pop, z, 2, numeq=1, QUIET=true, TS=false, normality=0)
end
@time begin
opt,sol,data = complex_tssos_first(pop, z, 3, numeq=1, QUIET=true, TS=false, normality=0)
end
@time begin
opt,sol,data = complex_tssos_first(pop, z, 2, numeq=1, QUIET=true, TS=false, normality=1)
end

@polyvar x[1:2n]
rf = pop[1](z=>x[1:n]+im*x[n+1:2n])
rpop = [real.(coefficients(rf))'*monomials(rf), 1 - sum(x.^2)]
@time begin
opt,sol,data = tssos_first(rpop, x, 2, numeq=1, GroebnerBasis=false, QUIET=true, TS=false)
end
println(sum(eigvals(data.moment[1]) .> 1e-4))

# Minimizing a random complex quartic polynomial with CS on multi-spheres
Random.seed!(1)
l = 5
n = 4l + 2
@complex_polyvar z[1:n]
f = 0
for i = 1:l
    cb = basis(z[4i-3:4i+2])
    lcb = length(cb)
    P = rand(lcb, lcb)
    Q = rand(lcb, lcb)
    f += cb'*((P+P')/2+im*(Q-Q')/2)*cb
end
pop = [f]
for i = 1:l
    push!(pop, 1 - sum(z[4i-3:4i+2]'*z[4i-3:4i+2]))
end
@time begin
opt,sol,data = complex_tssos_first(pop, z, 2, numeq=l, QUIET=true, TS=false, normality=0)
end
@time begin
opt,sol,data = complex_tssos_first(pop, z, 3, numeq=l, QUIET=true, TS=false, normality=0)
end
@time begin
opt,sol,data = complex_tssos_first(pop, z, 2, numeq=l, QUIET=true, TS=false, normality=1)
end

@polyvar x[1:2n]
rf = pop[1](z=>x[1:n]+im*x[n+1:2n])
rpop = [real.(coefficients(rf))'*monomials(rf)]
for i = 1:l
    push!(rpop, 1 - sum(x[4i-3:4i+2].^2) - sum(x[n+4i-3:n+4i+2].^2))
end
@time begin
opt,sol,data = cs_tssos_first(rpop, x, 2, numeq=l, QUIET=true, TS=false, solution=false)
end
# println(maximum([sum(eigvals(data.moment[i][1]) .> 1e-4) for i = 1:l]))

# The AC-OPF problem
include("D:/Programs/TSSOS/example/modelopf.jl")
cd("D:/Programs/PolyOPF/pglib")
silence()

case = "pglib_opf_case30_as"
AC = 803.13
opfdata = parse_file(case * ".m")
model = pop_opf_com_vol(opfdata, normal=true, AngleCons=true, LineLimit=true)
n = model.n
m = model.m
numeq = model.numeq
supp = model.supp
coe = model.coe
mc = maximum(abs.(coe[1]))
coe[1] = coe[1]./mc

t = @elapsed begin
opt,_,popd = complex_cs_tssos_first(supp, coe, n, "min", numeq=numeq, QUIET=true, normality=0, CS="MF", TS="block")
end
opt *= mc
mb = maximum(maximum.([maximum.(popd.blocksize[i]) for i = 1:popd.cql])) # maximal block size
gap = (AC-opt)*100/AC # optimality gap
println("n = $n, m = $m")
println("opt = $opt, time = $t, mb = $mb, gap = $gap%")

t = @elapsed begin
opt,_,popd = complex_cs_tssos_first(supp, coe, n, "min", numeq=numeq, QUIET=true, normality=1, CS="MF", TS="block")
end
opt *= mc
mb = maximum(maximum.([maximum.(popd.blocksize[i]) for i = 1:popd.cql])) # maximal block size
gap = (AC-opt)*100/AC # optimality gap
println("n = $n, m = $m")
println("opt = $opt, time = $t, mb = $mb, gap = $gap%")
