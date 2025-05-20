using TSSOS
using DynamicPolynomials

## Polyphase Code Waveform Design
## formulation with inequality constraints
N = 6
@complex_polyvar z[1:N+1]
f = z[N+1]^2 + conj(z[N+1])^2
cons = Vector{typeof(f)}(undef, N-2)
for k = 1:N-2
    cons[k] = z[N+1]^2 + conj(z[N+1])^2 - sum(z[i]*z[j+k]*conj(z[j]*z[i+k]) for i = 1:N-k, j = 1:N-k)
end

order = 4
@time begin
opt,sol,data = complex_tssos_first([f; cons], z, order, nb=N+1, TS="block", ConjugateBasis=true, normality=0, QUIET=false)
end
println(opt^0.5)

# writetofile="D:/project/ManiDSDP/polyphasecode4.sdpa"

## formulation with equality constraints
N = 8
@complex_polyvar z[1:2N-1]
f = z[N+1]^2 + conj(z[N+1])^2
cons = Vector{typeof(f)}(undef, N-2)
for k = 1:N-2
    cons[k] = z[N+1]^2 + conj(z[N+1])^2 - z[N+1+k]^2 - conj(z[N+1+k]^2) - 2 - sum(z[i]*z[j+k]*conj(z[j]*z[i+k]) for i = 1:N-k, j = 1:N-k)
    # z[N+1+k]*conj(z[N+1+k])
end

order = 3
@time begin
opt,sol,data = complex_tssos_first([f; cons], z, order, numeq=N-2, nb=2N-1, CS=false, TS="block", QUIET=false)
opt,sol,data = complex_tssos_higher!(data, TS="block", solve=true, QUIET=false)
end
println(opt^0.5)

# find a solution
mm = data.moment[1]
b = [acos(mm[N+3+2k][1,2]) for k=1:N-2]
# b[1] = b[8] = -b[1]
# b[2] = b[7] = -b[2]
# b[3] = b[6] = -b[3]
# b[4] = b[5] = -b[4]
A = diagm(2*ones(N-2))
A[1,2] = A[N-2,N-3] = -1
for k = 2:N-3
    A[k,k-1] = A[k,k+1] = -1
end
θ = A\b
θ = [0;θ;0]
sol = cos.(θ)+sin.(θ)*im
abs.([pop[j+1](z=>sol) for j=1:N-2])

# compute a local solution
N = 15
@polyvar x[1:N]
@polyvar y[1:N]
@polyvar t
pop = Vector{Poly}(undef, 2N-1)
pop[1] = t^2
for k = 1:N-2
    pop[k+1] = t^2 - sum(x[j]*x[j+k]+y[j]*y[j+k] for j=1:N-k)^2 - sum(x[j]*y[j+k]-y[j]*x[j+k] for j=1:N-k)^2
end
for k = 1:N
    pop[k+N-1] = 1 - x[k]^2 - y[k]^2
end

@time begin
opt,sol,data = tssos_first(pop, [x;y;t], 2, numeq=N, GroebnerBasis=false, TS=false, QUIET=true, solve=false)
# opt,sol,data = tssos_higher!(data, TS="block", QUIET=true)
end
# println(opt^0.5)

opt,sol = local_solution(data.n, data.m, data.supp, data.coe, numeq=N, startpoint=rand(2N+1), QUIET=false)
println(opt^0.5)

# Another model
N = 6
@complex_polyvar z[1:N]
f = sum(sum(z[i]*conj(z[i+j]) for i = 1:N-j)*sum(conj(z[i])*z[i+j] for i = 1:N-j) for j = 1:N-2)
order = 4
@time begin
opt,sol,data = complex_tssos_first([f], z, order, nb=N, TS="block", QUIET=false)
end

N = 12
@polyvar x[1:N]
@polyvar y[1:N]
f = sum(sum(x[j]*x[j+k]+y[j]*y[j+k] for j=1:N-k)^2 + sum(x[j]*y[j+k]-y[j]*x[j+k] for j=1:N-k)^2 for k=1:N-2)
cons = Vector{typeof(f)}(undef, N)
for k = 1:N
    cons[k] = 1 - x[k]^2 - y[k]^2
end
@time begin
opt,sol,data = tssos_first([f; cons], [x;y], 2, numeq=N, TS=false, GroebnerBasis=false, QUIET=true, solve=false)
# opt,sol,data = tssos_higher!(data, TS="block", QUIET=true)
end

opt,sol = local_solution(data.n, data.m, data.supp, data.coe, numeq=N, startpoint=rand(2N), QUIET=false)

# extract a pair of conjugate solutions
using LinearAlgebra
using LDLFactorizations
using DynamicPolynomials

n = 5
@complex_polyvar z[1:n]
# Q = rand(n, n)
# c = rand(n)
basis = cbasis(z)
Q = rand(length(basis), length(basis))
# ind = [Int.(floor.(rand(2)*length(basis)).+1) for i = 1:5]
# for item in ind
#     Q[item[1],item[2]] = rand(1)[1]
# end
Q = (Q+Q')/2
f = basis'*Q*basis
# f = z'*Q*z + c'*(z+conj(z))
# f = (1+im)*z[1]*z[2]*z[3] + (1+im)*z[2]*z[3]*z[4] + (1-im)*conj(z[1]*z[2]*z[3]) + (1-im)*conj(z[2]*z[3]*z[4])
h = 1 - z'*z
opt,sol,data = complex_tssos_first([f; h], z, 2, numeq=1, QUIET=true, TS=false)
M = Matrix{Float64}(data.moment[1][1][1:n+1,1:n+1])
F = ldl(M)
ind = findall(diag(F.D).>1e-4)
if length(ind) > 1
    temp = F.L[2:end, ind[2]]
    temp[ind[2]-1] = 1
    v = sqrt(M[ind[2],ind[2]]-F.L[ind[2],1]^2)
    sol = Matrix([F.L[2:end, 1] v*temp])
    obj = f(z=>sol[:,1]+im*sol[:,2])
    println([rank(M, 1e-4), abs(opt-obj)])
end
