using TSSOS
using DynamicPolynomials
using MultivariatePolynomials
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

# minimizing a random complex quadratic polynomial over the unit sphere
Random.seed!(1)
n = 10
@complex_polyvar z[1:n]
Q = rand(n+1, n+1)
Q = (Q+Q')/2
f = [1; z]'*Q*[1; z]
h = z'*z - 1
@time begin
opt,sol,data = complex_tssos_first([f; h], z, 1, numeq=1, QUIET=true, TS=false)
end

@polyvar x[1:2n]
rf = f(z=>x[1:n]+im*x[n+1:2n])
rpop = [real.(coefficients(rf))'*monomials(rf), 1-sum(x.^2)]
@time begin
opt,sol,data = tssos_first(rpop, x, 1, numeq=1, GroebnerBasis=false, QUIET=false, TS="block")
end

opt,sol = local_solution(2n, 1, data.supp, data.coe, numeq=1, startpoint=rand(2n), QUIET=true)
println(opt)

# minimizing a random complex quartic polynomial over the unit sphere
Random.seed!(1)
n = 20
@complex_polyvar z[1:n]
basis = cbasis(z)
P = rand(length(basis), length(basis))
f = basis'*((P+P')/2)*basis
h = z'*z - 1

@time begin
opt,sol,data = complex_tssos_first([f; h], z, 2, numeq=1, QUIET=true, TS=false)
end

@polyvar x[1:2n]
rf = f(z=>x[1:n]+im*x[n+1:2n])
rpop = [real.(coefficients(rf))'*monomials(rf), 1-sum(x.^2)]
@time begin
opt,sol,data = tssos_first(rpop, x, 2, numeq=1, GroebnerBasis=false, QUIET=true, TS="block")
end

opt,sol = local_solution(2n, 1, data.supp, data.coe, numeq=1, startpoint=rand(2n), QUIET=true)
println(opt)
