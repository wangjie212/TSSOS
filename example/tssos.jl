using DynamicPolynomials
using TSSOS

# unconstrained optimization using the TSSOS hierarchy
# randomly generated examples of type I
@polyvar x1 x2 x3 x4 x5 x6 x7 x8
x = tuple(x1,x2,x3,x4,x5,x6,x7,x8)
f=include("E:\\Programs\\blocksos\\poly1\\F1.txt")
@time begin
opt,sol,data=blockupop_first(f,x,newton=true,TS="block")
end
@time begin
opt,sol,data=blockupop_higher!(data,TS="block")
end

# randomly generated examples of type II
@polyvar x1 x2 x3 x4 x5 x6 x7 x8
x = tuple(x1,x2,x3,x4,x5,x6,x7,x8)
f=include("E:\\Programs\\blocksos\\poly2\\G1.txt")
@time begin
opt,sol,data=blockupop_first(f,x,newton=false,TS="block")
end
@time begin
opt,sol,data=blockupop_higher!(data,TS="block")
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
opt,sol,data=blockupop_first(f,x,newton=false,TS="block")
end

# the network system of type I
n=10
@polyvar x[1:n]
a=1.0.+rand(n)
b=(0.5.+rand(n,n))./n
f=x[1]
for i=1:n
   global f+=a[i]*(x[i]^2+x[i]^4)
end
for i=1:n
   for j=1:n
       global f-=b[i,j]*x[i]^2*x[j]^2
   end
end
f-=x[1]
@time begin
opt,sol,data=blockupop_first(f,x,newton=false,TS="MD")
end

# the network system of type II
n=50
@polyvar x[1:n]
a=0.5.+rand(n)
b=(0.5.+rand(n,n))./n
f=x[1]
for i=1:n
    global f+=a[i]*(1/2*x[i]^2-1/4*x[i]^4)
end
for i=1:n
    for j=1:n
        global f+=1/8*b[i,j]*(x[i]-x[j])^4
    end
end
f-=x[1]
pop=[f]
g=2
for i=1:n
    push!(pop,x[i]^2*(2-x[i]^2))
end
@time begin
opt,sol,data=blockcpop_first(pop,x,2,TS="MD")
end

# constrained optimization using the TSSOS hierarchy
# minimization over the unit ball
@polyvar x1 x2 x3 x4 x5 x6
x = tuple(x1,x2,x3,x4,x5,x6)
f=include("E:\\Programs\\blocksos\\cpoly\\H1.txt")
pop=[f,1-sum(x.^2)]
d=2 # the relaxation order
@time begin
opt,sol,data=blockcpop_first(pop,x,d,TS="block")
end
@time begin
opt,sol,data=blockcpop_higher!(data,TS="block")
end

# minimization over unit hypercube
@polyvar x1 x2 x3 x4 x5 x6
x = tuple(x1,x2,x3,x4,x5,x6)
f=include("E:\\Programs\\blocksos\\cpoly\\H1.txt")
pop=[f,1-x1^2,1-x2^2,1-x3^2,1-x4^2,1-x5^2,1-x6^2]
d=2 # the relaxation order
@time begin
opt,sol,data=blockcpop_first(pop,x,d,TS="block")
end
@time begin
opt,sol,data=blockcpop_higher!(data,TS="block")
end
