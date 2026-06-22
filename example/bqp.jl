using DynamicPolynomials
using TSSOS


n = 20
@polyvar x[1:n]
P = randn(n, n)
for i = 1:n
    P[i,i] = 0
end
f = x'*((P+P')/2)*x
opt,sol,data = tssos([f], x, 1, nb=n, QUIET=true, solve=true, TS=false)
opt,sol,data = tssos([f], x, 2, nb=n, QUIET=true, solve=true, TS=false)
opt,sol,data = tssos([f, f-opt], x, 1, nb=n, QUIET=true, solve=true, TS=false)