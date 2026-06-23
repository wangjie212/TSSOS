using DynamicPolynomials
using TSSOS


n = 4
@polyvar x[1:n]
f = (sum(x) + 1)^2
opt1,sol,data = tssos([f], x, 2, nb=n, QUIET=true, solve=true, TS=false)
opt2,sol,data = tssos([f], x, 3, nb=n, QUIET=true, solve=true, TS=false)
opt1,sol,data = tssos([f, f-opt1], x, 2, nb=n, QUIET=true, solve=true, TS=false)