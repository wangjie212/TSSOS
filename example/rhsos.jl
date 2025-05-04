using TSSOS
using DynamicPolynomials
using MultivariatePolynomials
using Random

function cbasis(z)
    basis = Monomial{true}[1]
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
n = 50
@polyvar z[1:2n]
Q = rand(n+1, n+1)
Q = (Q+Q')/2
f = [1; z[n+1:end]]'*Q*[1; z[1:n]]
h = sum(z[i]*z[i+n] for i = 1:n) - 1
@time begin
opt,sol,data = cs_tssos_first([f; h], z, n, 1, numeq=1, QUIET=true, CS=false, TS=false)
end

@polyvar x[1:2n]
rf = f(z[1:n]=>x[1:n]+im*x[n+1:2n], z[n+1:2n]=>x[1:n]-im*x[n+1:2n])
rpop = [real.(coefficients(rf))'*monomials(rf), 1-sum(x.^2)]
@time begin
opt,sol,data = tssos_first(rpop, x, 1, numeq=1, GroebnerBasis=false, QUIET=true, TS="block")
end

opt,sol = local_solution(2n, 1, data.supp, data.coe, numeq=1, startpoint=rand(2n), QUIET=true)
println(opt)

# minimizing a random complex quartic polynomial over the unit sphere
Random.seed!(1)
n = 20
@polyvar z[1:2n]
basis1 = cbasis(z[1:n])
basis2 = cbasis(z[n+1:2n])
P = rand(length(basis1), length(basis1))
f = basis2'*((P+P')/2)*basis1
h = sum(z[i]*z[i+n] for i = 1:n) - 1

@time begin
opt,sol,data = cs_tssos_first([f; h], z, n, 2, numeq=1, QUIET=true, CS=false, TS=false)
end

@polyvar x[1:2n]
rf = f(z[1:n]=>x[1:n]+im*x[n+1:2n], z[n+1:2n]=>x[1:n]-im*x[n+1:2n])
rpop = [real.(coefficients(rf))'*monomials(rf), 1-sum(x.^2)]
@time begin
opt,sol,data = tssos_first(rpop, x, 2, numeq=1, GroebnerBasis=false, QUIET=true, TS="block")
end
    
opt,sol = local_solution(2n, 1, data.supp, data.coe, numeq=1, startpoint=rand(2n), QUIET=true)
println(opt)
