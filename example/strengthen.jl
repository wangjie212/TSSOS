using TSSOS
using DynamicPolynomials
using Random
using LinearAlgebra

# Minimizing a random complex quadratic polynomial with unit-norm variables
Random.seed!(1)
n = 20
@polyvar z[1:2n]
P = rand(n+1, n+1)
Q = rand(n+1, n+1)
pop = [[1; z[n+1:2n]]'*((P+P')/2+im*(Q-Q')/2)*[1; z[1:n]]]
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 1, nb=n, QUIET=true, CS=false, TS=false, normality=0)
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 2, nb=n, QUIET=true, CS=false, TS=false, normality=0)
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 1, nb=n, QUIET=true, CS=false, TS=false, normality=1, mosek_setting=mosek_para(1e-8, 1e-8, 1e-8, -1, 0))
end

eigvals(convert.(ComplexF64, data.moment[1][1]))

@polyvar x[1:2n]
rf = pop[1](z[1:n]=>x[1:n]+im*x[n+1:2n], z[n+1:2n]=>x[1:n]-im*x[n+1:2n])
rpop = [real.(coefficients(rf))'*monomials(rf)]
for i = 1:n
    push!(rpop, 1 - x[i]^2 - x[i+n]^2)
end
@time begin
opt,sol,data = tssos_first(rpop, x, 2, numeq=n, Groebnerbasis=true, QUIET=true, TS=false)
end

function basis(x)
    basis = Monomial{true}[1]
    push!(basis, x...)
    for i = 1:length(x), j = i:length(x)
        push!(basis, x[i]*x[j])
    end
    return basis
end

# Minimizing a random complex quartic polynomial on a unit sphere
Random.seed!(1)
n = 5
@polyvar z[1:2n]
cb1 = basis(z[1:n])
cb2 = basis(z[n+1:2n])
lcb = length(cb1)
P = rand(lcb, lcb)
Q = rand(lcb, lcb)
pop = [cb2'*((P+P')/2+im*(Q-Q')/2)*cb1]
push!(pop, 1 - sum(z[1:n]'*z[n+1:2n]))
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 2, numeq=1, QUIET=true, CS=false, TS=false)
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 3, numeq=1, QUIET=true, CS=false, TS=false)
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 2, numeq=1, QUIET=true, CS=false, TS=false, normality=2)
end

eigvals(convert.(ComplexF64, data.moment[1][1]))

@polyvar x[1:2n]
rf = pop[1](z[1:n]=>x[1:n]+im*x[n+1:2n], z[n+1:2n]=>x[1:n]-im*x[n+1:2n])
rpop = [real.(coefficients(rf))'*monomials(rf), 1 - sum(x.^2)]
@time begin
opt,sol,data = tssos_first(rpop, x, 2, numeq=1, Groebnerbasis=false, QUIET=true, TS=false)
end

# Minimizing a random complex quartic polynomial with CS on multi-spheres
Random.seed!(2)
l = 5
n = 4l + 2
@polyvar z[1:2n]
f = 0
for i = 1:l
    cb1 = basis(z[4i-3:4i+2])
    cb2 = basis(z[n+4i-3:n+4i+2])
    lcb = length(cb1)
    P = rand(lcb, lcb)
    Q = rand(lcb, lcb)
    f += cb2'*((P+P')/2+im*(Q-Q')/2)*cb1
end
pop = [f]
for i = 1:l
    push!(pop, 1 - sum(z[4i-3:4i+2]'*z[n+4i-3:n+4i+2]))
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 2, numeq=l, QUIET=true, TS=false)
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 3, numeq=l, QUIET=true, TS=false)
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 2, numeq=l, QUIET=true, TS=false, normality=2)
end

@polyvar x[1:2n]
rf = pop[1](z[1:n]=>x[1:n]+im*x[n+1:2n], z[n+1:2n]=>x[1:n]-im*x[n+1:2n])
rpop = [real.(coefficients(rf))'*monomials(rf)]
for i = 1:l
    push!(rpop, 1 - sum(x[4i-3:4i+2].^2) - sum(x[n+4i-3:n+4i+2].^2))
end
@time begin
opt,sol,data = cs_tssos_first(rpop, x, 2, numeq=l, QUIET=true, TS=false, solution=false)
end

# The AC-OPF problem
include("D:/Programs/TSSOS/example/modelopf.jl")
cd("D:/Programs/PolyOPF/pglib")
silence()

case = "pglib_opf_case179_goc"
AC = 754270
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
opt,sol,popd = cs_tssos_first(supp, coe, n, "min", numeq=numeq, QUIET=true, normality=0, CS="MF", TS="block")
end
opt *= mc
mb = maximum(maximum.([maximum.(popd.blocksize[i]) for i = 1:popd.cql])) # maximal block size
gap = (AC-opt)*100/AC # optimality gap
println("n = $n, m = $m")
println("opt = $opt, time = $t, mb = $mb, gap = $gap%")

t = @elapsed begin
opt,sol,popd = cs_tssos_first(supp, coe, n, "min", numeq=numeq, QUIET=true, normality=1, CS="MF", TS="block")
end
opt *= mc
mb = maximum(maximum.([maximum.(popd.blocksize[i]) for i = 1:popd.cql])) # maximal block size
gap = (AC-opt)*100/AC # optimality gap
println("n = $n, m = $m")
println("opt = $opt, time = $t, mb = $mb, gap = $gap%")
