using TSSOS
using DynamicPolynomials
using Random

# Minimizing a random real quadratic polynomial with binary variables
Random.seed!(1)
n = 5
@polyvar x[1:n]
Q = rand(n+1, n+1)
Q = (Q+Q')/2
f = [1; x]'*Q*[1; x]
@time begin
opt,sol,data = cs_tssos_first([f], x, 1, nb=n, QUIET=true, CS=false, TS=false)
end
@time begin
opt,sol,data = cs_tssos_first([f], x, 2, nb=n, QUIET=true, CS=false, TS=false)
end
@time begin
opt,sol,data = cs_tssos_first([f], x, 1, nb=n, QUIET=true, CS=false, TS=false, normality=1)
end
@time begin
opt,sol,data = tssos_first([f], x, 1, nb=n, QUIET=true, TS=false, normality=1)
end

# Minimizing a random complex quadratic polynomial with unit-norm variables
Random.seed!(1)
n = 10
@polyvar z[1:2n]
P = rand(n+1, n+1)
Q = rand(n+1, n+1)
pop = [[1; z[n+1:2n]]'*((P+P')/2+im*(Q-Q')/2)*[1; z[1:n]]]
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 1, nb=n, QUIET=true, CS=false, TS=false)
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 2, nb=n, QUIET=true, CS=false, TS=false)
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 1, nb=n, QUIET=true, CS=false, TS=false, normality=1)
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
opt,sol,data = cs_tssos_first(pop, z, n, 2, numeq=1, QUIET=true, CS=false, TS=false, normality=1)
end

# Minimizing a random sparse complex quartic polynomial on a unit sphere
function sbasis(x)
    basis = Monomial{true}[1]
    for i = 1:length(x), j = i:length(x)
        push!(basis, x[i]*x[j])
    end
    return basis
end

Random.seed!(1)
n = 6
@polyvar z[1:2n]
cb1 = sbasis(z[1:n])
cb2 = sbasis(z[n+1:2n])
P = rand(length(cb1), length(cb1))
Q = rand(length(cb1), length(cb1))
pop = [cb2'*((P+P')/2+im*(Q-Q')/2)*cb1]
push!(pop, 1 - sum(z[1:n]'*z[n+1:2n]))
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 2, numeq=1, QUIET=true, CS=false, TS="block")
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 3, numeq=1, QUIET=true, CS=false, TS="block", solve=false)
opt,sol,data = cs_tssos_higher!(data, QUIET=true, TS="block")
end
@time begin
opt,sol,data = cs_tssos_first(pop, z, n, 2, numeq=1, QUIET=true, CS=false, TS="block", normality=1, NormalSparse=true)
end

# Minimizing a random complex quartic polynomial with CS on multi-spheres
Random.seed!(1)
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
opt,sol,data = cs_tssos_first(pop, z, n, 2, numeq=l, QUIET=true, TS=false, normality=1)
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

time = @elapsed begin
opt,sol,popd = cs_tssos_first(supp, coe, n, "min", numeq=numeq, QUIET=true, normality=0, tune=true, CS="MF", TS="block")
end
opt *= mc
mb = maximum(maximum.(popd.sb)) # maximal block size
gap = (AC-opt)*100/AC # optimality gap
println("n = $n, m = $m")
println("opt = $opt, time = $time, mb = $mb, gap = $gap%")

time = @elapsed begin
opt,sol,popd = cs_tssos_first(supp, coe, n, "min", numeq=numeq, QUIET=true, normality=1, tune=true, CS="MF", TS="block")
end
opt *= mc
mb = maximum(maximum.(popd.sb)) # maximal block size
gap = (AC-opt)*100/AC # optimality gap
println("n = $n, m = $m")
println("opt = $opt, time = $time, mb = $mb, gap = $gap%")
