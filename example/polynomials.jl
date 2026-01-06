# Rosenbrock-Lerner polynomial
N = 60
b = 15
@polyvar x[1:N]
f_R = 10*sum((x[3:N-1].+x[2:N-2].-x[1:N-3].^2).^2)+sum((ntuple(i->1,N-3).-x[1:N-3].-x[4:N]).^2)
f_Q = sum([x[i]*i/j*x[j] for i=1:b,j=1:b])
f = f_R + f_Q

# Broyden tridiagonal polynomial
l = 5
p = 20
n = l*p
@polyvar x[1:n]
f = ((3-2*x[1])*x[1]-2*x[2]+1)^2
for i = 2:n-1
    f += ((3-2*x[i])*x[i]-x[i-1]-2*x[i+1]+1)^2
end
f += ((3-2*x[n])*x[n]-x[n-1]+1)^2
pop = [f]
for i = 1:l
    push!(pop, 1-sum(x[(i-1)*p+1:i*p].^2))
end

# Chained singular polynomial
l = 2
p = 15
n = l*p
@polyvar x[1:n]
f = (x[1]+10*x[2])^2+5*(x[3]-x[4])^2+(x[2]-2*x[3])^4+10*(x[1]-10*x[4])^4
for i = 3:2:n-3
    f += (x[i]+10*x[i+1])^2+5*(x[i+2]-x[i+3])^2+(x[i+1]-2*x[i+2])^4+10*(x[i]-10*x[i+3])^4
end
pop = [f]
for i = 1:l
    push!(pop, 1-sum(x[(i-1)*p+1:i*p].^2))
end

# Broyden banded polynomial
n = 400
@polyvar x[1:n]
f = 0
for i = 1:n
    jset = max(1,i-5):min(n,i+1)
    jset = setdiff(jset,i)
    f0 = sum([(1+x[jset[j]])*x[jset[j]] for j=1:length(jset)])
    f += (x[i]*(2+5*x[i]^2)+1-f0)^2
end

# chained Wood polynomial
l = 50
p = 20
n = l*p
@polyvar x[1:n]
f = 1
for i = 1:2:n-3
    f += 100*(x[i+1]-x[i]^2)^2+(1-x[i])^2+90*(x[i+3]-x[i+2]^2)^2+(1-x[i+2])^2+10*(x[i+1]+x[i+3]-2)^2+0.1*(x[i+1]-x[i+3])^2
end
pop = [f]
for i = 1:l
    push!(pop, 1-sum(x[(i-1)*p+1:i*p].^2))
end

# generalized Rosenbrock polynomial
l = 50
p = 20
n = l*p
@polyvar x[1:n]
f = 1
for i = 2:n
    f += 100*(x[i]-x[i-1]^2)^2+(1-x[i])^2
end
pop = [f]
for i = 1:l
    push!(pop, 1-sum(x[(i-1)*p+1:i*p].^2))
end
