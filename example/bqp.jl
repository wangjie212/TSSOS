using DynamicPolynomials
using TSSOS


n = 10
@polyvar x[1:n]
P = randn(n, n)
for i = 1:n
    P[i,i] = 0
end
f = x'*((P+P')/2)*x
opt1,sol,data = tssos([f], x, 1, nb=n, QUIET=true, solve=true, TS=false)
opt2,sol,data = tssos([f], x, 2, nb=n, QUIET=true, solve=true, TS=false)
opt1,sol,data = tssos([f, f-opt1], x, 1, nb=n, QUIET=true, solve=true, TS=false)