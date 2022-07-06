using TSSOS

# N Polyphase Code Waveform Design, N+1 complex variables
# N bit Code Elements x[n], n = 1,...,N, the (N+1)th complex variables is the Object variable
N = 4   # N+1 complex variables
#  N           PSL           SDP Time(s)                   SDP Time(s)                      SDP Time(s)
#                           without ipart=false           with ipart=false                with ipart=false
#                           but All Equalities            and (N-2)Inequalities           and All Equalities
#  7   0.5219283102856175    6.4136603                     6.3560014                       1.4861946
#  8   0.6482710182019227    33.2900535                    5.2663592                       3.4292143
#  9   0.11192087452688078   66.5328702                    33.4952673                      9.724971
#  10  0.5804674373017424    437.9319477                   85.0143015                      26.1291904
#  11  0.49426371490943544   582.3490046                   294.8161438                     66.0632282
#  12  0.4928513896548784    3899.5969825                  1111.954345                     169.9979179

supp = Vector{Vector{Vector{Vector{UInt16}}}}(undef, N-1)
coe  = Vector{Vector{Float64}}(undef, N-1)

# Minimize the Object |x[N+1]|^2
supp[1] = [[[N+1],[N+1]]]
coe[1] = [1]

# sum(x[i]*x[i+k]') * sum(x[j]'*x[j+k]) <= x[N+1]*x[N+1]', i,j = 1,...,N-k, k = 1,...,N-2
for k = 1:N-2
    supp[k+1] = Vector{Vector{Vector{UInt16}}}(undef, (N-k)*(N-k)+1)
    coe[k+1] = Vector{Float64}(undef, (N-k)*(N-k)+1)
    for i = 1:N-k, j =1:N-k
        supp[k+1][j+(i-1)*(N-k)] = sort.([[i;(j+k)],[(i+k);j]])
        coe[k+1][j+(i-1)*(N-k)] = -1.0
    end
    supp[k+1][(N-k)*(N-k)+1] = [[N+1], [N+1]]
    coe[k+1][(N-k)*(N-k)+1] = 1.0
    # supp[k+1][(N-k)*(N-k)+2] = [[N+1+k], [N+1+k]]
    # coe[k+1][(N-k)*(N-k)+2] = -1.0
end

order = 3
@time begin
opt,sol,data = cs_tssos_first(supp, coe, N+1, order, CS=false, TS="block", ipart=false, solve=false, nb=N, QUIET=true)
opt,sol,data = cs_tssos_higher!(data, TS="block", ipart=false, QUIET=true, Mommat=true)
end
println(opt^0.5)

# compute a local solution
using DynamicPolynomials

N = 4
@polyvar x[1:N]
@polyvar y[1:N]
@polyvar t
pop = Vector{Polynomial{true,Float64}}(undef, 2N-1)
pop[1] = t^2
for k = 1:N-2
    pop[k+1] = t^2 - sum(x[j]*x[j+k]+y[j]*y[j+k] for j=1:N-k)^2 - sum(x[j]*y[j+k]-y[j]*x[j+k] for j=1:N-k)^2
end
for k = 1:N
    pop[k+N-1] = 1 - x[k]^2 - y[k]^2
end

@time begin
opt,sol,data = tssos_first(pop, [x;y;t], 2, numeq=N, TS="block", QUIET=true, solve=false, quotient=true)
# opt,sol,data = tssos_higher!(data, TS="block", QUIET=true)
end
# println(opt^0.5)
sol,ub,gap = refine_sol(1, rand(2N+1), data, QUIET=true)
println([ub^0.5])

@polyvar x[1:2N]
pop = Vector{Polynomial{true,Float64}}(undef, N-2)
for j = 1:N-2
    pop[j] = sum(x[i]*x[i+j+N] for i = 1:N-j)
end

z1 = 1
z2 = cos(18pi/31)+sin(18pi/31)*im
z3 = z2
z4 = 1
abs.([pop[j](x=>[z1;z2;z3;z4;conj(z1);conj(z2);conj(z3);conj(z4)]) for j=1:N-2])
