using DynamicPolynomials
using TSSOS

# unconstrained optimization using the chordal-TSSOS hierarchy
# randomly generated examples of type I
@polyvar x1 x2 x3 x4 x5 x6 x7 x8
x = tuple(x1,x2,x3,x4,x5,x6,x7,x8)
f=include("E:\\Programs\\blocksos\\poly1\\F1.txt")
@time begin
opt,sol,data=blockupop_first(f,x,newton=true,TS="MD")
end
@time begin
opt,sol,data=blockupop_higher!(data,TS="MD")
end

# randomly generated examples of type II
@polyvar x1 x2 x3 x4 x5 x6 x7 x8
x = tuple(x1,x2,x3,x4,x5,x6,x7,x8)
f=include("E:\\Programs\\blocksos\\poly2\\G1.txt")
@time begin
opt,sol,data=blockupop_first(f,x,newton=false,TS="MD")
end
@time begin
opt,sol,data=blockupop_higher!(data,TS="MD")
end

# the Broyden banded function
n=6
@polyvar x[1:n]
jset=Array{Any}(undef, n)
f=x[1]
for i=1:n
    jset[i]=max(1,i-5):min(n,i+1)
    jset[i]=setdiff(jset[i],i)
    f0=(1+x[jset[i][1]])*x[jset[i][1]]
    ljset=length(jset[i])
    for j=2:ljset
        f0+=(1+x[jset[i][j]])*x[jset[i][j]]
    end
    global f+=(x[i]*(2+5*x[i]^2)+1-f0)^2
end
f-=x[1]
@time begin
opt,sol,data=blockupop_first(f,x,newton=false,TS="MD")
end

# the modified generalized Rosenbrock function
n=10
@polyvar x[1:n]
f=x[1]+1
for i=2:n
    global f+=100*(x[i]-x[i-1]^2)^2+(1-x[i])^2
end
for i=1:n
    for j=i+1:n
        global f+=x[i]^2*x[j]^2
    end
end
f-=x[1]
@time begin
opt,sol,data=blockupop_first(f,x,newton=false,TS="MD")
end

# the modified chained singular function
n=10
@polyvar x[1:n]
f=(x[1]+10*x[2])^2+5*(x[3]-x[4])^2+(x[2]-2*x[3])^4+10*(x[1]-10*x[4])^4
for i=3:2:n-3
    global f+=(x[i]+10*x[i+1])^2+5*(x[i+2]-x[i+3])^2+(x[i+1]-2*x[i+2])^4+10*(x[i]-10*x[i+3])^4
end
for i=1:n
    for j=i+1:n
        global f+=x[i]^2*x[j]^2
    end
end
@time begin
opt,sol,data=blockupop_first(f,x,newton=false,TS="MD")
end

# constrained optimization using the chordal-TSSOS hierarchy
# minimization of randomly generated examples over the unit ball
@polyvar x1 x2 x3 x4 x5 x6
x = tuple(x1,x2,x3,x4,x5,x6)
f=include("E:\\Programs\\blocksos\\cpoly\\H1.txt")
pop=[f,1-sum(x.^2)]
d=2 # the relaxation order
@time begin
opt,sol,data=blockcpop_first(pop,x,d,TS="MD")
end
@time begin
opt,sol,data=blockcpop_higher!(data,TS="MD")
end

# minimization of randomly generated examples over unit hypercube
@polyvar x1 x2 x3 x4 x5 x6
x = tuple(x1,x2,x3,x4,x5,x6)
f=include("E:\\Programs\\blocksos\\cpoly\\H1.txt")
pop=[f,1-x1^2,1-x2^2,1-x3^2,1-x4^2,1-x5^2,1-x6^2]
d=2 # the relaxation order
@time begin
opt,sol,data=blockcpop_first(pop,x,d,TS="MD")
end
@time begin
opt,sol,data=blockcpop_higher!(data,TS="MD")
end

# minimization of the Broyden tridiagonal function over the unit ball
n=10
@polyvar x[1:n]
f=((3-2*x[1])*x[1]-2*x[2]+1)^2
for i=2:n-1
    global f+=((3-2*x[i])*x[i]-x[i-1]-2*x[i+1]+1)^2
end
f+=((3-2*x[n])*x[n]-x[n-1]+1)^2
pop=[f,1-sum(x.^2)]
d=2 # the relaxation order
@time begin
opt,sol,data=blockcpop_first(pop,x,d,TS="MD")
end
@time begin
opt,sol,data=blockcpop_higher!(data,TS="MD")
end

# minimization of the generalized Rosenbrock function over the unit ball
n=10
@polyvar x[1:n]
f=x[1]+1
for i=2:n
    global f+=100*(x[i]-x[i-1]^2)^2+(1-x[i])^2
end
f-=x[1]
pop=[f,1-sum(x.^2)]
d=2 # the relaxation order
@time begin
opt,sol,data=blockcpop_first(pop,x,d,TS="MD")
end
@time begin
opt,sol,data=blockcpop_higher!(data,TS="MD")
end
