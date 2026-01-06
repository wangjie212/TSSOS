using DynamicPolynomials
using TSSOS
using Random

# Example 5.1
@polyvar x[1:4]
f = x[3]^2*(x[1]^2 + x[1]^4*x[2]^2 + x[3]^4 - 3x[1]^2*x[2]^2) + x[2]^8 + x[1]^2*x[2]^2*x[4]^2
# Homogenization without CS
solve_hpop(f, x, [], [], 5, QUIET=true, CS=false, TS="block")
# Homogenization with CS type 1
solve_hpop(f, x, [], [], 5, QUIET=true, type=1, TS="block", SO=2, nnhomovar=true)
# Homogenization with CS type 2
solve_hpop(f, x, [], [], 5, QUIET=true, type=2, TS="block")
# Homogenization with CS type 3
solve_hpop(f, x, [], [], 5, QUIET=true, type=3, TS="block")
# No homogenization with CS
opt,sol,data = cs_tssos([f], x, 5, TS="block", solution=true, QUIET=true)

#######################################################################
# Example 5.0
@polyvar x[1:6]
f = x[1]^6 + x[2]^6 + 1 + 3*x[1]^2*x[2]^2 - x[1]^4*(x[2]^2 + 1) - x[2]^4*(x[1]^2 + 1) - (x[1]^2 + x[2]^2) + x[3]^2*(x[1]^2 + x[5]^2 - 2)^2 +
(x[4] + x[5] + 1)^2*x[6]^2 + (x[6] - 1)^6
# Homogenization without CS
solve_hpop(f, x, [], [], 4, QUIET=true, CS=false, TS="block", SO=2, nnhomovar=true)
# Homogenization with CS type 1
solve_hpop(f, x, [], [], 5, QUIET=true, type=1, ε=1e-4, TS="block", SO=2, nnhomovar=true)
# Homogenization with CS type 2
solve_hpop(f, x, [], [], 5, QUIET=true, type=2, TS="block", SO=2, nnhomovar=true)
# Homogenization with CS type 3
solve_hpop(f, x, [], [], 5, QUIET=true, type=3, TS="block", SO=2, nnhomovar=true)
# No homogenization with CS
@time begin
opt,sol,data = cs_tssos([f], x, 7, TS="block", solution=true, QUIET=true)
end
# Find a local solution
obj,sol,status = local_solution(data.n, data.m, data.supp, data.coe, startpoint=rand(data.n))

#######################################################################
# Example 5.2
@polyvar x[1:10]
f1 = sum(x[1:4].^4) + (1-x[1])*(1-x[2])*(1-x[3])*(1-x[4]) + (x[1]-1)*(x[1]-x[2])*(x[1]-x[3])*(x[1]-x[4]) + (x[2]-1)*(x[2]-x[1])*(x[2]-x[3])*(x[2]-x[4]) + (x[3]-1)*(x[3]-x[1])*(x[3]-x[2])*(x[3]-x[4]) + (x[4]-1)*(x[4]-x[1])*(x[4]-x[2])*(x[4]-x[3])
f2 = sum(x[4:7].^4) + (1-x[4])*(1-x[5])*(1-x[6])*(1-x[7]) + (x[4]-1)*(x[4]-x[5])*(x[4]-x[6])*(x[4]-x[7]) + (x[5]-1)*(x[5]-x[4])*(x[5]-x[6])*(x[5]-x[7]) + (x[6]-1)*(x[6]-x[4])*(x[6]-x[5])*(x[6]-x[7]) + (x[7]-1)*(x[7]-x[4])*(x[7]-x[5])*(x[7]-x[6])
f3 = sum(x[7:10].^4) + (1-x[7])*(1-x[8])*(1-x[9])*(1-x[10]) + (x[7]-1)*(x[7]-x[8])*(x[7]-x[9])*(x[7]-x[10]) + (x[8]-1)*(x[8]-x[7])*(x[8]-x[9])*(x[8]-x[10]) + (x[9]-1)*(x[9]-x[7])*(x[9]-x[8])*(x[9]-x[10]) + (x[10]-1)*(x[10]-x[7])*(x[10]-x[8])*(x[10]-x[9])
f = f1 + f2 + f3
# Homogenization without CS
solve_hpop(f, x, [], [], 2, QUIET=true, CS=false, TS="block")
# Homogenization with CS type 1
solve_hpop(f, x, [], [], 2, QUIET=true, type=1, TS="block", SO=2, nnhomovar=true)
# Homogenization with CS type 2
solve_hpop(f, x, [], [], 2, QUIET=true, type=2, TS="block")
# Homogenization with CS type 3
solve_hpop(f, x, [], [], 2, QUIET=true, type=3, TS="block")
# No homogenization with CS
@time begin
opt,sol,data = cs_tssos([f], x, 4, TS="block", solution=true, QUIET=true)
end

#######################################################################
# Example 5.3
@polyvar x[1:20]
f1 = x[1]^2*(x[1]-1)^2 + x[2]^2*(x[2]-1)^2 + x[3]^2*(x[3]-1)^2 + 2*x[1]*x[2]*x[3]*(sum(x[1:3])-2) + 0.25*((x[1]-1)^2+(x[2]-1)^2+(x[3]-1)^2+(x[4]-1)^2) + (x[4]*x[5]-1)^2
f2 = x[1]^2*(x[1]-1)^2 + x[2]^2*(x[2]-1)^2 + x[6]^2*(x[6]-1)^2 + 2*x[1]*x[2]*x[6]*(sum(x[1:2])+x[6]-2) + 0.25*((x[1]-1)^2+(x[2]-1)^2+(x[6]-1)^2+(x[7]-1)^2) + (x[7]*x[8]-1)^2
f3 = x[1]^2*(x[1]-1)^2 + x[2]^2*(x[2]-1)^2 + x[9]^2*(x[9]-1)^2 + 2*x[1]*x[2]*x[9]*(sum(x[1:2])+x[9]-2) + 0.25*((x[1]-1)^2+(x[2]-1)^2+(x[9]-1)^2+(x[10]-1)^2) + (x[10]*x[11]-1)^2
f4 = x[1]^2*(x[1]-1)^2 + x[2]^2*(x[2]-1)^2 + x[12]^2*(x[12]-1)^2 + 2*x[1]*x[2]*x[12]*(sum(x[1:2])+x[12]-2) + 0.25*((x[1]-1)^2+(x[2]-1)^2+(x[12]-1)^2+(x[13]-1)^2) + (x[13]*x[14]-1)^2
f5 = x[1]^2*(x[1]-1)^2 + x[2]^2*(x[2]-1)^2 + x[15]^2*(x[15]-1)^2 + 2*x[1]*x[2]*x[15]*(sum(x[1:2])+x[15]-2) + 0.25*((x[1]-1)^2+(x[2]-1)^2+(x[15]-1)^2+(x[16]-1)^2) + (x[16]*x[17]-1)^2
f6 = x[1]^2*(x[1]-1)^2 + x[2]^2*(x[2]-1)^2 + x[18]^2*(x[18]-1)^2 + 2*x[1]*x[2]*x[18]*(sum(x[1:2])+x[18]-2) + 0.25*((x[1]-1)^2+(x[2]-1)^2+(x[18]-1)^2+(x[19]-1)^2) + (x[19]*x[20]-1)^2
f = f1 + f2 + f3 + f4 + f5 + f6
# Homogenization without CS
solve_hpop(f, x, [], [], 2, QUIET=true, CS=false, TS="block")
# Homogenization with CS type 1
solve_hpop(f, x, [], [], 2, QUIET=true, type=1, TS="block", SO=2, nnhomovar=true)
# Homogenization with CS type 2
solve_hpop(f, x, [], [], 3, QUIET=true, type=2, TS="block", SO=1, nnhomovar=true)
# Homogenization with CS type 3
solve_hpop(f, x, [], [], 3, QUIET=true, type=3, TS="block", SO=1, nnhomovar=true)
# No homogenization with CS
@time begin
opt,sol,data = cs_tssos([f], x, 2, TS="block", solution=false, QUIET=true)
end

#######################################################################
# Example 5.4
@polyvar x[1:5]
f = x[1]^2 + 3x[2]^2 - 2x[2]*x[3]^2 + x[3]^4 - x[2]*(x[4]^2 + x[5]^2)
g = [x[2]-1; x[1]^2-2*x[1]*x[2]-1; x[1]^2+2*x[1]*x[2]-1; x[2]-x[4]^2-x[5]^2]
# Homogenization without CS
solve_hpop(f, x, g, [], 4, QUIET=true, CS=false, TS="block", SO=2)
# Homogenization with CS type 1
solve_hpop(f, x, g, [], 4, QUIET=true, type=1, ε=1e-4, TS="block", SO=2)
# Homogenization with CS type 2
solve_hpop(f, x, g, [], 4, QUIET=true, type=2, TS="block", SO=2)
# Homogenization with CS type 3
solve_hpop(f, x, g, [], 4, QUIET=true, type=3, TS="block", SO=2)
# No homogenization with CS
@time begin
opt,sol,data = cs_tssos([f; g], x, 2, TS="block", solution=false, QUIET=true)
end
# Find a local solution
# obj,sol,status = local_solution(data.n, data.m, data.supp, data.coe, startpoint=rand(data.n))

#######################################################################
# Example 5.5
@polyvar x[1:7]
f = x[1]^4*x[2]^2 + x[2]^4*x[3]^2 + x[1]^2*x[3]^4 - 3*(x[1]*x[2]*x[3])^2 + x[2]^2 + x[7]^2*(sum(x[1:3].^2)) + 
x[4]^2*x[5]^2*(10 - x[6]^2) + x[7]^2*(x[4]^2 + 2x[5]^2 + 3x[6]^2)
g = [x[1]-x[2]*x[3], -x[2]+x[3]^2, 1-sum(x[4:6].^2)]
# Homogenization without CS
solve_hpop(f, x, g, [], 5, QUIET=true, CS=false, TS="block", SO=2, nnhomovar=true)
# Homogenization with CS type 1
solve_hpop(f, x, g, [], 5, QUIET=true, type=1, ε=1e-4, TS="block", SO=2, nnhomovar=true)
# Homogenization with CS type 2
solve_hpop(f, x, g, [], 5, QUIET=true, type=2, TS="block", SO=2, nnhomovar=true)
# Homogenization with CS type 3
solve_hpop(f, x, g, [], 5, QUIET=true, type=3, TS="block", SO=2, nnhomovar=true)
# No homogenization with CS
@time begin
opt,sol,data = cs_tssos([f; g], x, 3, TS="block", solution=true, QUIET=true)
end

#######################################################################
# Example 5.6
p = 3
@polyvar x[1:8*p+2]
f = 0
for i = 1:p
    f += (sum(x[8*i-7:8*i+2].^2)+1)^2 - 4*(x[8*i-7]^2*x[8*i-6]^2+x[8*i-6]^2*x[8*i-5]^2+x[8*i-5]^2*x[8*i-4]^2+x[8*i-4]^2*x[8*i-3]^2+x[8*i-3]^2*x[8*i-7]^2) - 4*(x[8*i-2]^2*x[8*i-1]^2+x[8*i-1]^2*x[8*i]^2+x[8*i]^2*x[8*i+1]^2+x[8*i+1]^2*x[8*i+2]^2+x[8*i+2]^2*x[8*i-2]^2) + 0.2*sum(x[8*i-7:8*i+2].^4)
end
g = [sum(x[8*i-7:8*i+2].^2) - 1 for i=1:p]
# Homogenization without CS
solve_hpop(f, x, g, [], 4, QUIET=true, CS=false, TS="block")
# Homogenization with CS type 1
solve_hpop(f, x, g, [], 4, QUIET=true, type=1, ε=1e-4, TS="block")
# Homogenization with CS type 2
solve_hpop(f, x, g, [], 4, QUIET=true, type=2, TS="block")
# Homogenization with CS type 3
solve_hpop(f, x, g, [], 4, QUIET=true, type=3, TS="block")
# No homogenization with CS
@time begin
opt,sol,data = cs_tssos([f; g], x, 4, TS="block", solution=true, QUIET=true)
end

#######################################################################
# Example 5.7
Random.seed!(0)
n = 40
u = 3
p = trunc(Int, n/u)
@polyvar x[1:n]
A = rand(u, u)
b = rand(u)
f = (x[1:u].^2)'*A'*A*(x[1:u].^2) + b'*(x[1:u].^2)
g = [sum(x[1:u].^4) - 1]
for i = 2:p-1
    A = rand(u+1, u+1)
    b = rand(u+1)
    f += (x[u*(i-1):u*i].^2)'*A'*A*(x[u*(i-1):u*i].^2) + b'*(x[u*(i-1):u*i].^2)
    g = [g; sum(x[u*(i-1):u*i].^4) - 1]
end
A = rand(n-u*(p-1)+1, n-u*(p-1)+1)
b = rand(n-u*(p-1)+1)
f += (x[u*(p-1):n].^2)'*A'*A*(x[u*(p-1):n].^2) + b'*(x[u*(p-1):n].^2)
g = [g; sum(x[u*(p-1):n].^4) - 1]
# Homogenization without CS
solve_hpop(f, x, g, [], 4, QUIET=true, CS=false, TS="block")
# Homogenization with CS type 1
solve_hpop(f, x, g, [], 4, QUIET=true, type=1, ε=1e-4, TS="block")
# Homogenization with CS type 2
solve_hpop(f, x, g, [], 4, QUIET=true, type=2, TS="block")
# Homogenization with CS type 3
solve_hpop(f, x, g, [], 4, QUIET=true, type=3, TS="block")
# No homogenization with CS
@time begin
opt,sol,data = cs_tssos([f; g], x, 2, TS="block", solution=false, QUIET=true)
end
# Find a local solution
obj,sol,status = local_solution(data.n, data.m, data.supp, data.coe, startpoint=rand(data.n))
