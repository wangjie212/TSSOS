using TSSOS
using DynamicPolynomials

# Smale’s Mean Value Conjecture
n = 2
@complex_polyvar z1 z2 z3
z = tuple(z1,z2,z3)
f = - z3*conj(z3)
g1 = z1^2*conj(z1)^2 - z1^2 - conj(z1)^2 - 4*z3*conj(z3) + 1
g2 = z2^2*conj(z2)^2 - z2^2 - conj(z2)^2 - 4*z3*conj(z3) + 1
h1 = - z1*conj(z1) - z2*conj(z2) + 2/3
h2 = 3*z1*z2 + 3*conj(z1)*conj(z2) - 2
h3 = 9*z1*z2*conj(z1)*conj(z2) - 1
cpop = [f; g1; g2; h1; h2; h3]

@time begin
opt,sol,data = complex_tssos_first(cpop, z, 2, numeq=3, TS="block", QUIET=true)
end
println(sqrt(-opt))

n = 3
@complex_polyvar z1 z2 z3 z4
z = tuple(z1,z2,z3,z4)
f = - z4*conj(z4)
g1 = 9 - 6*z1^3 + 12*z1^2*z2 + 12*z1^2*z3 - 6*conj(z1)^3 + 4*z1^3*conj(z1)^3 -
 8*z1^2*z2*conj(z1)^3 - 8*z1^2*z3*conj(z1)^3 + 12*conj(z1)^2*conj(z2) -
   8*z1^3*conj(z1)^2*conj(z2) + 16*z1^2*z2*conj(z1)^2*conj(z2) + 16*z1^2*z3*conj(z1)^2*conj(z2) +
 12*conj(z1)^2*conj(z3) - 8*z1^3*conj(z1)^2*conj(z3) + 16*z1^2*z2*conj(z1)^2*conj(z3) +
   16*z1^2*z3*conj(z1)^2*conj(z3) - 36*z4*conj(z4)
g2 = 9 + 12*z1*z2^2 - 6*z2^3 + 12*z2^2*z3 + 12*conj(z1)*conj(z2)^2 +
 16*z1*z2^2*conj(z1)*conj(z2)^2 - 8*z2^3*conj(z1)*conj(z2)^2 + 16*z2^2*z3*conj(z1)*conj(z2)^2 -
   6*conj(z2)^3 - 8*z1*z2^2*conj(z2)^3 + 4*z2^3*conj(z2)^3 - 8*z2^2*z3*conj(z2)^3 +
 12*conj(z2)^2*conj(z3) + 16*z1*z2^2*conj(z2)^2*conj(z3) - 8*z2^3*conj(z2)^2*conj(z3) +
   16*z2^2*z3*conj(z2)^2*conj(z3) - 36*z4*conj(z4)
g3 = 9 + 12*z1*z3^2 + 12*z2*z3^2 - 6*z3^3 + 12*conj(z1)*conj(z3)^2 +
 16*z1*z3^2*conj(z1)*conj(z3)^2 + 16*z2*z3^2*conj(z1)*conj(z3)^2 - 8*z3^3*conj(z1)*conj(z3)^2 +
   12*conj(z2)*conj(z3)^2 + 16*z1*z3^2*conj(z2)*conj(z3)^2 + 16*z2*z3^2*conj(z2)*conj(z3)^2 -
 8*z3^3*conj(z2)*conj(z3)^2 - 6*conj(z3)^3 - 8*z1*z3^2*conj(z3)^3 -
   8*z2*z3^2*conj(z3)^3 + 4*z3^3*conj(z3)^3 - 36*z4*conj(z4)
h1 = z1*conj(z1) + z2*conj(z2) + z3*conj(z3) - n*(1/(n+1))^(2/n)
h2 = 2*z1*z2*z3 + 2*conj(z1)*conj(z2)*conj(z3) + 1
h3 = 16*z1*z2*z3*conj(z1)*conj(z2)*conj(z3) - 1
cpop = [f; g1; g2; g3; h1; h2; h3]

@time begin
opt,sol,data = complex_tssos_first(cpop, z, 5, numeq=3, TS="block", QUIET=true, solve=false)
opt,sol,data = complex_tssos_higher!(data, TS="block", QUIET=true, solve=false)
opt,sol,data = complex_tssos_higher!(data, TS="block", QUIET=true)  
end
println(sqrt(-opt))

n = 4
@complex_polyvar z1 z2 z3 z4 z5
z = tuple(z1,z2,z3,z4,z5)
f = - z5*conj(z5)
g1 = 36 - 18*z1^4 + 30*z1^3*z2 + 30*z1^3*z3 - 60*z1^2*z2*z3 + 30*z1^3*z4 - 60*z1^2*z2*z4 - 60*z1^2*z3*z4 - 18*conj(z1)^4 +
  9*z1^4*conj(z1)^4 - 15*z1^3*z2*conj(z1)^4 - 15*z1^3*z3*conj(z1)^4 + 30*z1^2*z2*z3*conj(z1)^4 - 15*z1^3*z4*conj(z1)^4 + 30*z1^2*z2*z4*conj(z1)^4 +
  30*z1^2*z3*z4*conj(z1)^4 + 30*conj(z1)^3*conj(z2) - 15*z1^4*conj(z1)^3*conj(z2) + 25*z1^3*z2*conj(z1)^3*conj(z2) + 25*z1^3*z3*conj(z1)^3*conj(z2) -
  50*z1^2*z2*z3*conj(z1)^3*conj(z2) + 25*z1^3*z4*conj(z1)^3*conj(z2) - 50*z1^2*z2*z4*conj(z1)^3*conj(z2) - 50*z1^2*z3*z4*conj(z1)^3*conj(z2) + 30*conj(z1)^3*conj(z3) -
  15*z1^4*conj(z1)^3*conj(z3) + 25*z1^3*z2*conj(z1)^3*conj(z3) + 25*z1^3*z3*conj(z1)^3*conj(z3) - 50*z1^2*z2*z3*conj(z1)^3*conj(z3) + 25*z1^3*z4*conj(z1)^3*conj(z3) -
  50*z1^2*z2*z4*conj(z1)^3*conj(z3) - 50*z1^2*z3*z4*conj(z1)^3*conj(z3) - 60*conj(z1)^2*conj(z2)*conj(z3) + 30*z1^4*conj(z1)^2*conj(z2)*conj(z3) - 50*z1^3*z2*conj(z1)^2*conj(z2)*conj(z3) -
  50*z1^3*z3*conj(z1)^2*conj(z2)*conj(z3) + 100*z1^2*z2*z3*conj(z1)^2*conj(z2)*conj(z3) - 50*z1^3*z4*conj(z1)^2*conj(z2)*conj(z3) + 100*z1^2*z2*z4*conj(z1)^2*conj(z2)*conj(z3) +
  100*z1^2*z3*z4*conj(z1)^2*conj(z2)*conj(z3) + 30*conj(z1)^3*conj(z4) - 15*z1^4*conj(z1)^3*conj(z4) + 25*z1^3*z2*conj(z1)^3*conj(z4) + 25*z1^3*z3*conj(z1)^3*conj(z4) -
  50*z1^2*z2*z3*conj(z1)^3*conj(z4) + 25*z1^3*z4*conj(z1)^3*conj(z4) - 50*z1^2*z2*z4*conj(z1)^3*conj(z4) - 50*z1^2*z3*z4*conj(z1)^3*conj(z4) - 60*conj(z1)^2*conj(z2)*conj(z4) +
  30*z1^4*conj(z1)^2*conj(z2)*conj(z4) - 50*z1^3*z2*conj(z1)^2*conj(z2)*conj(z4) - 50*z1^3*z3*conj(z1)^2*conj(z2)*conj(z4) + 100*z1^2*z2*z3*conj(z1)^2*conj(z2)*conj(z4) -
  50*z1^3*z4*conj(z1)^2*conj(z2)*conj(z4) + 100*z1^2*z2*z4*conj(z1)^2*conj(z2)*conj(z4) + 100*z1^2*z3*z4*conj(z1)^2*conj(z2)*conj(z4) - 60*conj(z1)^2*conj(z3)*conj(z4) +
  30*z1^4*conj(z1)^2*conj(z3)*conj(z4) - 50*z1^3*z2*conj(z1)^2*conj(z3)*conj(z4) - 50*z1^3*z3*conj(z1)^2*conj(z3)*conj(z4) + 100*z1^2*z2*z3*conj(z1)^2*conj(z3)*conj(z4) -
  50*z1^3*z4*conj(z1)^2*conj(z3)*conj(z4) + 100*z1^2*z2*z4*conj(z1)^2*conj(z3)*conj(z4) + 100*z1^2*z3*z4*conj(z1)^2*conj(z3)*conj(z4) - 144*z5*conj(z5)
g2 = 36 + 30*z1*z2^3 - 18*z2^4 - 60*z1*z2^2*z3 + 30*z2^3*z3 - 60*z1*z2^2*z4 + 30*z2^3*z4 - 60*z2^2*z3*z4 +
  30*conj(z1)*conj(z2)^3 + 25*z1*z2^3*conj(z1)*conj(z2)^3 - 15*z2^4*conj(z1)*conj(z2)^3 - 50*z1*z2^2*z3*conj(z1)*conj(z2)^3 + 25*z2^3*z3*conj(z1)*conj(z2)^3 -
  50*z1*z2^2*z4*conj(z1)*conj(z2)^3 + 25*z2^3*z4*conj(z1)*conj(z2)^3 - 50*z2^2*z3*z4*conj(z1)*conj(z2)^3 - 18*conj(z2)^4 - 15*z1*z2^3*conj(z2)^4 + 9*z2^4*conj(z2)^4 +
  30*z1*z2^2*z3*conj(z2)^4 - 15*z2^3*z3*conj(z2)^4 + 30*z1*z2^2*z4*conj(z2)^4 - 15*z2^3*z4*conj(z2)^4 + 30*z2^2*z3*z4*conj(z2)^4 -
  60*conj(z1)*conj(z2)^2*conj(z3) - 50*z1*z2^3*conj(z1)*conj(z2)^2*conj(z3) + 30*z2^4*conj(z1)*conj(z2)^2*conj(z3) + 100*z1*z2^2*z3*conj(z1)*conj(z2)^2*conj(z3) -
  50*z2^3*z3*conj(z1)*conj(z2)^2*conj(z3) + 100*z1*z2^2*z4*conj(z1)*conj(z2)^2*conj(z3) - 50*z2^3*z4*conj(z1)*conj(z2)^2*conj(z3) + 100*z2^2*z3*z4*conj(z1)*conj(z2)^2*conj(z3) +
  30*conj(z2)^3*conj(z3) + 25*z1*z2^3*conj(z2)^3*conj(z3) - 15*z2^4*conj(z2)^3*conj(z3) - 50*z1*z2^2*z3*conj(z2)^3*conj(z3) + 25*z2^3*z3*conj(z2)^3*conj(z3) -
  50*z1*z2^2*z4*conj(z2)^3*conj(z3) + 25*z2^3*z4*conj(z2)^3*conj(z3) - 50*z2^2*z3*z4*conj(z2)^3*conj(z3) - 60*conj(z1)*conj(z2)^2*conj(z4) - 50*z1*z2^3*conj(z1)*conj(z2)^2*conj(z4) +
  30*z2^4*conj(z1)*conj(z2)^2*conj(z4) + 100*z1*z2^2*z3*conj(z1)*conj(z2)^2*conj(z4) - 50*z2^3*z3*conj(z1)*conj(z2)^2*conj(z4) + 100*z1*z2^2*z4*conj(z1)*conj(z2)^2*conj(z4) -
  50*z2^3*z4*conj(z1)*conj(z2)^2*conj(z4) + 100*z2^2*z3*z4*conj(z1)*conj(z2)^2*conj(z4) + 30*conj(z2)^3*conj(z4) + 25*z1*z2^3*conj(z2)^3*conj(z4) - 15*z2^4*conj(z2)^3*conj(z4) -
  50*z1*z2^2*z3*conj(z2)^3*conj(z4) + 25*z2^3*z3*conj(z2)^3*conj(z4) - 50*z1*z2^2*z4*conj(z2)^3*conj(z4) + 25*z2^3*z4*conj(z2)^3*conj(z4) -
  50*z2^2*z3*z4*conj(z2)^3*conj(z4) - 60*conj(z2)^2*conj(z3)*conj(z4) - 50*z1*z2^3*conj(z2)^2*conj(z3)*conj(z4) + 30*z2^4*conj(z2)^2*conj(z3)*conj(z4) +
  100*z1*z2^2*z3*conj(z2)^2*conj(z3)*conj(z4) - 50*z2^3*z3*conj(z2)^2*conj(z3)*conj(z4) + 100*z1*z2^2*z4*conj(z2)^2*conj(z3)*conj(z4) - 50*z2^3*z4*conj(z2)^2*conj(z3)*conj(z4) +
  100*z2^2*z3*z4*conj(z2)^2*conj(z3)*conj(z4) - 144*z5*conj(z5)
g3 = 36 - 60*z1*z2*z3^2 + 30*z1*z3^3 + 30*z2*z3^3 - 18*z3^4 - 60*z1*z3^2*z4 - 60*z2*z3^2*z4 + 30*z3^3*z4 -
  60*conj(z1)*conj(z2)*conj(z3)^2 + 100*z1*z2*z3^2*conj(z1)*conj(z2)*conj(z3)^2 - 50*z1*z3^3*conj(z1)*conj(z2)*conj(z3)^2 - 50*z2*z3^3*conj(z1)*conj(z2)*conj(z3)^2 +
  30*z3^4*conj(z1)*conj(z2)*conj(z3)^2 + 100*z1*z3^2*z4*conj(z1)*conj(z2)*conj(z3)^2 + 100*z2*z3^2*z4*conj(z1)*conj(z2)*conj(z3)^2 - 50*z3^3*z4*conj(z1)*conj(z2)*conj(z3)^2 +
  30*conj(z1)*conj(z3)^3 - 50*z1*z2*z3^2*conj(z1)*conj(z3)^3 + 25*z1*z3^3*conj(z1)*conj(z3)^3 + 25*z2*z3^3*conj(z1)*conj(z3)^3 - 15*z3^4*conj(z1)*conj(z3)^3 -
  50*z1*z3^2*z4*conj(z1)*conj(z3)^3 - 50*z2*z3^2*z4*conj(z1)*conj(z3)^3 + 25*z3^3*z4*conj(z1)*conj(z3)^3 + 30*conj(z2)*conj(z3)^3 - 50*z1*z2*z3^2*conj(z2)*conj(z3)^3 +
  25*z1*z3^3*conj(z2)*conj(z3)^3 + 25*z2*z3^3*conj(z2)*conj(z3)^3 - 15*z3^4*conj(z2)*conj(z3)^3 - 50*z1*z3^2*z4*conj(z2)*conj(z3)^3 - 50*z2*z3^2*z4*conj(z2)*conj(z3)^3 +
  25*z3^3*z4*conj(z2)*conj(z3)^3 - 18*conj(z3)^4 + 30*z1*z2*z3^2*conj(z3)^4 - 15*z1*z3^3*conj(z3)^4 - 15*z2*z3^3*conj(z3)^4 + 9*z3^4*conj(z3)^4 +
  30*z1*z3^2*z4*conj(z3)^4 + 30*z2*z3^2*z4*conj(z3)^4 - 15*z3^3*z4*conj(z3)^4 - 60*conj(z1)*conj(z3)^2*conj(z4) + 100*z1*z2*z3^2*conj(z1)*conj(z3)^2*conj(z4) -
  50*z1*z3^3*conj(z1)*conj(z3)^2*conj(z4) - 50*z2*z3^3*conj(z1)*conj(z3)^2*conj(z4) + 30*z3^4*conj(z1)*conj(z3)^2*conj(z4) + 100*z1*z3^2*z4*conj(z1)*conj(z3)^2*conj(z4) +
  100*z2*z3^2*z4*conj(z1)*conj(z3)^2*conj(z4) - 50*z3^3*z4*conj(z1)*conj(z3)^2*conj(z4) - 60*conj(z2)*conj(z3)^2*conj(z4) + 100*z1*z2*z3^2*conj(z2)*conj(z3)^2*conj(z4) -
  50*z1*z3^3*conj(z2)*conj(z3)^2*conj(z4) - 50*z2*z3^3*conj(z2)*conj(z3)^2*conj(z4) + 30*z3^4*conj(z2)*conj(z3)^2*conj(z4) + 100*z1*z3^2*z4*conj(z2)*conj(z3)^2*conj(z4) +
  100*z2*z3^2*z4*conj(z2)*conj(z3)^2*conj(z4) - 50*z3^3*z4*conj(z2)*conj(z3)^2*conj(z4) + 30*conj(z3)^3*conj(z4) - 50*z1*z2*z3^2*conj(z3)^3*conj(z4) + 25*z1*z3^3*conj(z3)^3*conj(z4) +
  25*z2*z3^3*conj(z3)^3*conj(z4) - 15*z3^4*conj(z3)^3*conj(z4) - 50*z1*z3^2*z4*conj(z3)^3*conj(z4) - 50*z2*z3^2*z4*conj(z3)^3*conj(z4) + 25*z3^3*z4*conj(z3)^3*conj(z4) - 144*z5*conj(z5)
g4 = 36 - 60*z1*z2*z4^2 - 60*z1*z3*z4^2 - 60*z2*z3*z4^2 + 30*z1*z4^3 + 30*z2*z4^3 + 30*z3*z4^3 - 18*z4^4 -
  60*conj(z1)*conj(z2)*conj(z4)^2 + 100*z1*z2*z4^2*conj(z1)*conj(z2)*conj(z4)^2 + 100*z1*z3*z4^2*conj(z1)*conj(z2)*conj(z4)^2 + 100*z2*z3*z4^2*conj(z1)*conj(z2)*conj(z4)^2 -
  50*z1*z4^3*conj(z1)*conj(z2)*conj(z4)^2 - 50*z2*z4^3*conj(z1)*conj(z2)*conj(z4)^2 - 50*z3*z4^3*conj(z1)*conj(z2)*conj(z4)^2 + 30*z4^4*conj(z1)*conj(z2)*conj(z4)^2 - 60*conj(z1)*conj(z3)*conj(z4)^2 +
  100*z1*z2*z4^2*conj(z1)*conj(z3)*conj(z4)^2 + 100*z1*z3*z4^2*conj(z1)*conj(z3)*conj(z4)^2 + 100*z2*z3*z4^2*conj(z1)*conj(z3)*conj(z4)^2 - 50*z1*z4^3*conj(z1)*conj(z3)*conj(z4)^2 -
  50*z2*z4^3*conj(z1)*conj(z3)*conj(z4)^2 - 50*z3*z4^3*conj(z1)*conj(z3)*conj(z4)^2 + 30*z4^4*conj(z1)*conj(z3)*conj(z4)^2 - 60*conj(z2)*conj(z3)*conj(z4)^2 +
  100*z1*z2*z4^2*conj(z2)*conj(z3)*conj(z4)^2 + 100*z1*z3*z4^2*conj(z2)*conj(z3)*conj(z4)^2 + 100*z2*z3*z4^2*conj(z2)*conj(z3)*conj(z4)^2 - 50*z1*z4^3*conj(z2)*conj(z3)*conj(z4)^2 -
  50*z2*z4^3*conj(z2)*conj(z3)*conj(z4)^2 - 50*z3*z4^3*conj(z2)*conj(z3)*conj(z4)^2 + 30*z4^4*conj(z2)*conj(z3)*conj(z4)^2 + 30*conj(z1)*conj(z4)^3 - 50*z1*z2*z4^2*conj(z1)*conj(z4)^3 -
  50*z1*z3*z4^2*conj(z1)*conj(z4)^3 - 50*z2*z3*z4^2*conj(z1)*conj(z4)^3 + 25*z1*z4^3*conj(z1)*conj(z4)^3 + 25*z2*z4^3*conj(z1)*conj(z4)^3 + 25*z3*z4^3*conj(z1)*conj(z4)^3 -
  15*z4^4*conj(z1)*conj(z4)^3 + 30*conj(z2)*conj(z4)^3 - 50*z1*z2*z4^2*conj(z2)*conj(z4)^3 - 50*z1*z3*z4^2*conj(z2)*conj(z4)^3 - 50*z2*z3*z4^2*conj(z2)*conj(z4)^3 +
  25*z1*z4^3*conj(z2)*conj(z4)^3 + 25*z2*z4^3*conj(z2)*conj(z4)^3 + 25*z3*z4^3*conj(z2)*conj(z4)^3 - 15*z4^4*conj(z2)*conj(z4)^3 + 30*conj(z3)*conj(z4)^3 -
  50*z1*z2*z4^2*conj(z3)*conj(z4)^3 - 50*z1*z3*z4^2*conj(z3)*conj(z4)^3 - 50*z2*z3*z4^2*conj(z3)*conj(z4)^3 + 25*z1*z4^3*conj(z3)*conj(z4)^3 +
  25*z2*z4^3*conj(z3)*conj(z4)^3 + 25*z3*z4^3*conj(z3)*conj(z4)^3 - 15*z4^4*conj(z3)*conj(z4)^3 - 18*conj(z4)^4 + 30*z1*z2*z4^2*conj(z4)^4 + 30*z1*z3*z4^2*conj(z4)^4 +
  30*z2*z3*z4^2*conj(z4)^4 - 15*z1*z4^3*conj(z4)^4 - 15*z2*z4^3*conj(z4)^4 - 15*z3*z4^3*conj(z4)^4 + 9*z4^4*conj(z4)^4 - 144*z5*conj(z5)
h1 = z1*conj(z1) + z2*conj(z2) + z3*conj(z3) + z4*conj(z4) - n*(1/(n+1))^(2/n)
h2 = 5*z1*z2*z3*z4 + 5*conj(z1)*conj(z2)*conj(z3)*conj(z4) - 2
h3 = 25*z1*z2*z3*z4*conj(z1)*conj(z2)*conj(z3)*conj(z4) - 1
cpop = [f; g1; g2; g3; g4; h1; h2; h3]

@time begin
opt,sol,data = complex_tssos_first(cpop, z, 4, numeq=3, TS="block", QUIET=true, normality=1, solve=false)
opt,sol,data = complex_tssos_higher!(data, TS="block", QUIET=true, solve=false)
opt,sol,data = complex_tssos_higher!(data, TS="block", QUIET=true, solve=false)
opt,sol,data = complex_tssos_higher!(data, TS="block", QUIET=true)  
end
println(sqrt(-opt))

n = 5
@complex_polyvar z1 z2 z3 z4 z5 z6
z = tuple(z1,z2,z3,z4,z5,z6)
f = - z6*conj(z6)
g1 = 900 + 360*z1^5 - 540*z1^4*z2 - 540*z1^4*z3 + 900*z1^3*z2*z3 - 540*z1^4*z4 + 900*z1^3*z2*z4 + 900*z1^3*z3*z4 -
  1800*z1^2*z2*z3*z4 - 540*z1^4*z5 + 900*z1^3*z2*z5 + 900*z1^3*z3*z5 - 1800*z1^2*z2*z3*z5 + 900*z1^3*z4*z5 -
  1800*z1^2*z2*z4*z5 - 1800*z1^2*z3*z4*z5 + 900*conj(z4)*conj(z5)*conj(z1)^3 + 360*z1^5*conj(z4)*conj(z5)*conj(z1)^3 - 540*z1^4*conj(z4)*conj(z5)*z2*conj(z1)^3 -
  540*z1^4*conj(z4)*conj(z5)*z3*conj(z1)^3 + 900*z1^3*conj(z4)*conj(z5)*z2*z3*conj(z1)^3 - 540*z1^4*conj(z4)*conj(z5)*z4*conj(z1)^3 +
  900*z1^3*conj(z4)*conj(z5)*z2*z4*conj(z1)^3 + 900*z1^3*conj(z4)*conj(z5)*z3*z4*conj(z1)^3 - 1800*z1^2*conj(z4)*conj(z5)*z2*z3*z4*conj(z1)^3 -
  540*z1^4*conj(z4)*conj(z5)*z5*conj(z1)^3 + 900*z1^3*conj(z4)*conj(z5)*z2*z5*conj(z1)^3 + 900*z1^3*conj(z4)*conj(z5)*z3*z5*conj(z1)^3 -
  1800*z1^2*conj(z4)*conj(z5)*z2*z3*z5*conj(z1)^3 + 900*z1^3*conj(z4)*conj(z5)*z4*z5*conj(z1)^3 - 1800*z1^2*conj(z4)*conj(z5)*z2*z4*z5*conj(z1)^3 -
  1800*z1^2*conj(z4)*conj(z5)*z3*z4*z5*conj(z1)^3 - 540*conj(z4)*conj(z1)^4 - 216*z1^5*conj(z4)*conj(z1)^4 - 540*conj(z5)*conj(z1)^4 - 216*z1^5*conj(z5)*conj(z1)^4 +
  324*z1^4*conj(z4)*z2*conj(z1)^4 + 324*z1^4*conj(z5)*z2*conj(z1)^4 + 324*z1^4*conj(z4)*z3*conj(z1)^4 + 324*z1^4*conj(z5)*z3*conj(z1)^4 -
  540*z1^3*conj(z4)*z2*z3*conj(z1)^4 - 540*z1^3*conj(z5)*z2*z3*conj(z1)^4 + 324*z1^4*conj(z4)*z4*conj(z1)^4 + 324*z1^4*conj(z5)*z4*conj(z1)^4 -
  540*z1^3*conj(z4)*z2*z4*conj(z1)^4 - 540*z1^3*conj(z5)*z2*z4*conj(z1)^4 - 540*z1^3*conj(z4)*z3*z4*conj(z1)^4 - 540*z1^3*conj(z5)*z3*z4*conj(z1)^4 +
  1080*z1^2*conj(z4)*z2*z3*z4*conj(z1)^4 + 1080*z1^2*conj(z5)*z2*z3*z4*conj(z1)^4 + 324*z1^4*conj(z4)*z5*conj(z1)^4 + 324*z1^4*conj(z5)*z5*conj(z1)^4 -
  540*z1^3*conj(z4)*z2*z5*conj(z1)^4 - 540*z1^3*conj(z5)*z2*z5*conj(z1)^4 - 540*z1^3*conj(z4)*z3*z5*conj(z1)^4 - 540*z1^3*conj(z5)*z3*z5*conj(z1)^4 +
  1080*z1^2*conj(z4)*z2*z3*z5*conj(z1)^4 + 1080*z1^2*conj(z5)*z2*z3*z5*conj(z1)^4 - 540*z1^3*conj(z4)*z4*z5*conj(z1)^4 - 540*z1^3*conj(z5)*z4*z5*conj(z1)^4 +
  1080*z1^2*conj(z4)*z2*z4*z5*conj(z1)^4 + 1080*z1^2*conj(z5)*z2*z4*z5*conj(z1)^4 + 1080*z1^2*conj(z4)*z3*z4*z5*conj(z1)^4 +
  1080*z1^2*conj(z5)*z3*z4*z5*conj(z1)^4 + 360*conj(z1)^5 + 144*z1^5*conj(z1)^5 - 216*z1^4*z2*conj(z1)^5 - 216*z1^4*z3*conj(z1)^5 +
  360*z1^3*z2*z3*conj(z1)^5 - 216*z1^4*z4*conj(z1)^5 + 360*z1^3*z2*z4*conj(z1)^5 + 360*z1^3*z3*z4*conj(z1)^5 - 720*z1^2*z2*z3*z4*conj(z1)^5 -
  216*z1^4*z5*conj(z1)^5 + 360*z1^3*z2*z5*conj(z1)^5 + 360*z1^3*z3*z5*conj(z1)^5 - 720*z1^2*z2*z3*z5*conj(z1)^5 + 360*z1^3*z4*z5*conj(z1)^5 -
  720*z1^2*z2*z4*z5*conj(z1)^5 - 720*z1^2*z3*z4*z5*conj(z1)^5 - 1800*conj(z4)*conj(z5)*conj(z1)^2*conj(z2) - 720*z1^5*conj(z4)*conj(z5)*conj(z1)^2*conj(z2) +
  1080*z1^4*conj(z4)*conj(z5)*z2*conj(z1)^2*conj(z2) + 1080*z1^4*conj(z4)*conj(z5)*z3*conj(z1)^2*conj(z2) - 1800*z1^3*conj(z4)*conj(z5)*z2*z3*conj(z1)^2*conj(z2) +
  1080*z1^4*conj(z4)*conj(z5)*z4*conj(z1)^2*conj(z2) - 1800*z1^3*conj(z4)*conj(z5)*z2*z4*conj(z1)^2*conj(z2) - 1800*z1^3*conj(z4)*conj(z5)*z3*z4*conj(z1)^2*conj(z2) +
  3600*z1^2*conj(z4)*conj(z5)*z2*z3*z4*conj(z1)^2*conj(z2) + 1080*z1^4*conj(z4)*conj(z5)*z5*conj(z1)^2*conj(z2) - 1800*z1^3*conj(z4)*conj(z5)*z2*z5*conj(z1)^2*conj(z2) -
  1800*z1^3*conj(z4)*conj(z5)*z3*z5*conj(z1)^2*conj(z2) + 3600*z1^2*conj(z4)*conj(z5)*z2*z3*z5*conj(z1)^2*conj(z2) - 1800*z1^3*conj(z4)*conj(z5)*z4*z5*conj(z1)^2*conj(z2) +
  3600*z1^2*conj(z4)*conj(z5)*z2*z4*z5*conj(z1)^2*conj(z2) + 3600*z1^2*conj(z4)*conj(z5)*z3*z4*z5*conj(z1)^2*conj(z2) + 900*conj(z4)*conj(z1)^3*conj(z2) +
  360*z1^5*conj(z4)*conj(z1)^3*conj(z2) + 900*conj(z5)*conj(z1)^3*conj(z2) + 360*z1^5*conj(z5)*conj(z1)^3*conj(z2) - 540*z1^4*conj(z4)*z2*conj(z1)^3*conj(z2) -
  540*z1^4*conj(z5)*z2*conj(z1)^3*conj(z2) - 540*z1^4*conj(z4)*z3*conj(z1)^3*conj(z2) - 540*z1^4*conj(z5)*z3*conj(z1)^3*conj(z2) + 900*z1^3*conj(z4)*z2*z3*conj(z1)^3*conj(z2) +
  900*z1^3*conj(z5)*z2*z3*conj(z1)^3*conj(z2) - 540*z1^4*conj(z4)*z4*conj(z1)^3*conj(z2) - 540*z1^4*conj(z5)*z4*conj(z1)^3*conj(z2) + 900*z1^3*conj(z4)*z2*z4*conj(z1)^3*conj(z2) +
  900*z1^3*conj(z5)*z2*z4*conj(z1)^3*conj(z2) + 900*z1^3*conj(z4)*z3*z4*conj(z1)^3*conj(z2) + 900*z1^3*conj(z5)*z3*z4*conj(z1)^3*conj(z2) -
  1800*z1^2*conj(z4)*z2*z3*z4*conj(z1)^3*conj(z2) - 1800*z1^2*conj(z5)*z2*z3*z4*conj(z1)^3*conj(z2) - 540*z1^4*conj(z4)*z5*conj(z1)^3*conj(z2) -
  540*z1^4*conj(z5)*z5*conj(z1)^3*conj(z2) + 900*z1^3*conj(z4)*z2*z5*conj(z1)^3*conj(z2) + 900*z1^3*conj(z5)*z2*z5*conj(z1)^3*conj(z2) +
  900*z1^3*conj(z4)*z3*z5*conj(z1)^3*conj(z2) + 900*z1^3*conj(z5)*z3*z5*conj(z1)^3*conj(z2) - 1800*z1^2*conj(z4)*z2*z3*z5*conj(z1)^3*conj(z2) -
  1800*z1^2*conj(z5)*z2*z3*z5*conj(z1)^3*conj(z2) + 900*z1^3*conj(z4)*z4*z5*conj(z1)^3*conj(z2) + 900*z1^3*conj(z5)*z4*z5*conj(z1)^3*conj(z2) -
  1800*z1^2*conj(z4)*z2*z4*z5*conj(z1)^3*conj(z2) - 1800*z1^2*conj(z5)*z2*z4*z5*conj(z1)^3*conj(z2) - 1800*z1^2*conj(z4)*z3*z4*z5*conj(z1)^3*conj(z2) -
  1800*z1^2*conj(z5)*z3*z4*z5*conj(z1)^3*conj(z2) - 540*conj(z1)^4*conj(z2) - 216*z1^5*conj(z1)^4*conj(z2) + 324*z1^4*z2*conj(z1)^4*conj(z2) + 324*z1^4*z3*conj(z1)^4*conj(z2) -
  540*z1^3*z2*z3*conj(z1)^4*conj(z2) + 324*z1^4*z4*conj(z1)^4*conj(z2) - 540*z1^3*z2*z4*conj(z1)^4*conj(z2) - 540*z1^3*z3*z4*conj(z1)^4*conj(z2) +
  1080*z1^2*z2*z3*z4*conj(z1)^4*conj(z2) + 324*z1^4*z5*conj(z1)^4*conj(z2) - 540*z1^3*z2*z5*conj(z1)^4*conj(z2) - 540*z1^3*z3*z5*conj(z1)^4*conj(z2) +
  1080*z1^2*z2*z3*z5*conj(z1)^4*conj(z2) - 540*z1^3*z4*z5*conj(z1)^4*conj(z2) + 1080*z1^2*z2*z4*z5*conj(z1)^4*conj(z2) + 1080*z1^2*z3*z4*z5*conj(z1)^4*conj(z2) -
  1800*conj(z4)*conj(z5)*conj(z1)^2*conj(z3) - 720*z1^5*conj(z4)*conj(z5)*conj(z1)^2*conj(z3) + 1080*z1^4*conj(z4)*conj(z5)*z2*conj(z1)^2*conj(z3) + 1080*z1^4*conj(z4)*conj(z5)*z3*conj(z1)^2*conj(z3) -
  1800*z1^3*conj(z4)*conj(z5)*z2*z3*conj(z1)^2*conj(z3) + 1080*z1^4*conj(z4)*conj(z5)*z4*conj(z1)^2*conj(z3) - 1800*z1^3*conj(z4)*conj(z5)*z2*z4*conj(z1)^2*conj(z3) -
  1800*z1^3*conj(z4)*conj(z5)*z3*z4*conj(z1)^2*conj(z3) + 3600*z1^2*conj(z4)*conj(z5)*z2*z3*z4*conj(z1)^2*conj(z3) + 1080*z1^4*conj(z4)*conj(z5)*z5*conj(z1)^2*conj(z3) -
  1800*z1^3*conj(z4)*conj(z5)*z2*z5*conj(z1)^2*conj(z3) - 1800*z1^3*conj(z4)*conj(z5)*z3*z5*conj(z1)^2*conj(z3) + 3600*z1^2*conj(z4)*conj(z5)*z2*z3*z5*conj(z1)^2*conj(z3) -
  1800*z1^3*conj(z4)*conj(z5)*z4*z5*conj(z1)^2*conj(z3) + 3600*z1^2*conj(z4)*conj(z5)*z2*z4*z5*conj(z1)^2*conj(z3) + 3600*z1^2*conj(z4)*conj(z5)*z3*z4*z5*conj(z1)^2*conj(z3) +
  900*conj(z4)*conj(z1)^3*conj(z3) + 360*z1^5*conj(z4)*conj(z1)^3*conj(z3) + 900*conj(z5)*conj(z1)^3*conj(z3) + 360*z1^5*conj(z5)*conj(z1)^3*conj(z3) - 540*z1^4*conj(z4)*z2*conj(z1)^3*conj(z3) -
  540*z1^4*conj(z5)*z2*conj(z1)^3*conj(z3) - 540*z1^4*conj(z4)*z3*conj(z1)^3*conj(z3) - 540*z1^4*conj(z5)*z3*conj(z1)^3*conj(z3) + 900*z1^3*conj(z4)*z2*z3*conj(z1)^3*conj(z3) +
  900*z1^3*conj(z5)*z2*z3*conj(z1)^3*conj(z3) - 540*z1^4*conj(z4)*z4*conj(z1)^3*conj(z3) - 540*z1^4*conj(z5)*z4*conj(z1)^3*conj(z3) + 900*z1^3*conj(z4)*z2*z4*conj(z1)^3*conj(z3) +
  900*z1^3*conj(z5)*z2*z4*conj(z1)^3*conj(z3) + 900*z1^3*conj(z4)*z3*z4*conj(z1)^3*conj(z3) + 900*z1^3*conj(z5)*z3*z4*conj(z1)^3*conj(z3) -
  1800*z1^2*conj(z4)*z2*z3*z4*conj(z1)^3*conj(z3) - 1800*z1^2*conj(z5)*z2*z3*z4*conj(z1)^3*conj(z3) - 540*z1^4*conj(z4)*z5*conj(z1)^3*conj(z3) -
  540*z1^4*conj(z5)*z5*conj(z1)^3*conj(z3) + 900*z1^3*conj(z4)*z2*z5*conj(z1)^3*conj(z3) + 900*z1^3*conj(z5)*z2*z5*conj(z1)^3*conj(z3) +
  900*z1^3*conj(z4)*z3*z5*conj(z1)^3*conj(z3) + 900*z1^3*conj(z5)*z3*z5*conj(z1)^3*conj(z3) - 1800*z1^2*conj(z4)*z2*z3*z5*conj(z1)^3*conj(z3) -
  1800*z1^2*conj(z5)*z2*z3*z5*conj(z1)^3*conj(z3) + 900*z1^3*conj(z4)*z4*z5*conj(z1)^3*conj(z3) + 900*z1^3*conj(z5)*z4*z5*conj(z1)^3*conj(z3) -
  1800*z1^2*conj(z4)*z2*z4*z5*conj(z1)^3*conj(z3) - 1800*z1^2*conj(z5)*z2*z4*z5*conj(z1)^3*conj(z3) - 1800*z1^2*conj(z4)*z3*z4*z5*conj(z1)^3*conj(z3) -
  1800*z1^2*conj(z5)*z3*z4*z5*conj(z1)^3*conj(z3) - 540*conj(z1)^4*conj(z3) - 216*z1^5*conj(z1)^4*conj(z3) + 324*z1^4*z2*conj(z1)^4*conj(z3) + 324*z1^4*z3*conj(z1)^4*conj(z3) -
  540*z1^3*z2*z3*conj(z1)^4*conj(z3) + 324*z1^4*z4*conj(z1)^4*conj(z3) - 540*z1^3*z2*z4*conj(z1)^4*conj(z3) - 540*z1^3*z3*z4*conj(z1)^4*conj(z3) +
  1080*z1^2*z2*z3*z4*conj(z1)^4*conj(z3) + 324*z1^4*z5*conj(z1)^4*conj(z3) - 540*z1^3*z2*z5*conj(z1)^4*conj(z3) - 540*z1^3*z3*z5*conj(z1)^4*conj(z3) +
  1080*z1^2*z2*z3*z5*conj(z1)^4*conj(z3) - 540*z1^3*z4*z5*conj(z1)^4*conj(z3) + 1080*z1^2*z2*z4*z5*conj(z1)^4*conj(z3) + 1080*z1^2*z3*z4*z5*conj(z1)^4*conj(z3) -
  1800*conj(z4)*conj(z1)^2*conj(z2)*conj(z3) - 720*z1^5*conj(z4)*conj(z1)^2*conj(z2)*conj(z3) - 1800*conj(z5)*conj(z1)^2*conj(z2)*conj(z3) - 720*z1^5*conj(z5)*conj(z1)^2*conj(z2)*conj(z3) +
  1080*z1^4*conj(z4)*z2*conj(z1)^2*conj(z2)*conj(z3) + 1080*z1^4*conj(z5)*z2*conj(z1)^2*conj(z2)*conj(z3) + 1080*z1^4*conj(z4)*z3*conj(z1)^2*conj(z2)*conj(z3) +
  1080*z1^4*conj(z5)*z3*conj(z1)^2*conj(z2)*conj(z3) - 1800*z1^3*conj(z4)*z2*z3*conj(z1)^2*conj(z2)*conj(z3) - 1800*z1^3*conj(z5)*z2*z3*conj(z1)^2*conj(z2)*conj(z3) +
  1080*z1^4*conj(z4)*z4*conj(z1)^2*conj(z2)*conj(z3) + 1080*z1^4*conj(z5)*z4*conj(z1)^2*conj(z2)*conj(z3) - 1800*z1^3*conj(z4)*z2*z4*conj(z1)^2*conj(z2)*conj(z3) -
  1800*z1^3*conj(z5)*z2*z4*conj(z1)^2*conj(z2)*conj(z3) - 1800*z1^3*conj(z4)*z3*z4*conj(z1)^2*conj(z2)*conj(z3) - 1800*z1^3*conj(z5)*z3*z4*conj(z1)^2*conj(z2)*conj(z3) +
  3600*z1^2*conj(z4)*z2*z3*z4*conj(z1)^2*conj(z2)*conj(z3) + 3600*z1^2*conj(z5)*z2*z3*z4*conj(z1)^2*conj(z2)*conj(z3) + 1080*z1^4*conj(z4)*z5*conj(z1)^2*conj(z2)*conj(z3) +
  1080*z1^4*conj(z5)*z5*conj(z1)^2*conj(z2)*conj(z3) - 1800*z1^3*conj(z4)*z2*z5*conj(z1)^2*conj(z2)*conj(z3) - 1800*z1^3*conj(z5)*z2*z5*conj(z1)^2*conj(z2)*conj(z3) -
  1800*z1^3*conj(z4)*z3*z5*conj(z1)^2*conj(z2)*conj(z3) - 1800*z1^3*conj(z5)*z3*z5*conj(z1)^2*conj(z2)*conj(z3) + 3600*z1^2*conj(z4)*z2*z3*z5*conj(z1)^2*conj(z2)*conj(z3) +
  3600*z1^2*conj(z5)*z2*z3*z5*conj(z1)^2*conj(z2)*conj(z3) - 1800*z1^3*conj(z4)*z4*z5*conj(z1)^2*conj(z2)*conj(z3) - 1800*z1^3*conj(z5)*z4*z5*conj(z1)^2*conj(z2)*conj(z3) +
  3600*z1^2*conj(z4)*z2*z4*z5*conj(z1)^2*conj(z2)*conj(z3) + 3600*z1^2*conj(z5)*z2*z4*z5*conj(z1)^2*conj(z2)*conj(z3) + 3600*z1^2*conj(z4)*z3*z4*z5*conj(z1)^2*conj(z2)*conj(z3) +
  3600*z1^2*conj(z5)*z3*z4*z5*conj(z1)^2*conj(z2)*conj(z3) + 900*conj(z1)^3*conj(z2)*conj(z3) + 360*z1^5*conj(z1)^3*conj(z2)*conj(z3) - 540*z1^4*z2*conj(z1)^3*conj(z2)*conj(z3) -
  540*z1^4*z3*conj(z1)^3*conj(z2)*conj(z3) + 900*z1^3*z2*z3*conj(z1)^3*conj(z2)*conj(z3) - 540*z1^4*z4*conj(z1)^3*conj(z2)*conj(z3) + 900*z1^3*z2*z4*conj(z1)^3*conj(z2)*conj(z3) +
  900*z1^3*z3*z4*conj(z1)^3*conj(z2)*conj(z3) - 1800*z1^2*z2*z3*z4*conj(z1)^3*conj(z2)*conj(z3) - 540*z1^4*z5*conj(z1)^3*conj(z2)*conj(z3) +
  900*z1^3*z2*z5*conj(z1)^3*conj(z2)*conj(z3) + 900*z1^3*z3*z5*conj(z1)^3*conj(z2)*conj(z3) - 1800*z1^2*z2*z3*z5*conj(z1)^3*conj(z2)*conj(z3) +
  900*z1^3*z4*z5*conj(z1)^3*conj(z2)*conj(z3) - 1800*z1^2*z2*z4*z5*conj(z1)^3*conj(z2)*conj(z3) - 1800*z1^2*z3*z4*z5*conj(z1)^3*conj(z2)*conj(z3) - 3600*z6*conj(z6)
g2 = 900 - 540*z1*z2^4 + 360*z2^5 + 900*z1*z2^3*z3 - 540*z2^4*z3 + 900*z1*z2^3*z4 - 540*z2^4*z4 - 1800*z1*z2^2*z3*z4 +
  900*z2^3*z3*z4 + 900*z1*z2^3*z5 - 540*z2^4*z5 - 1800*z1*z2^2*z3*z5 + 900*z2^3*z3*z5 - 1800*z1*z2^2*z4*z5 +
  900*z2^3*z4*z5 - 1800*z2^2*z3*z4*z5 - 1800*conj(z4)*conj(z5)*conj(z1)*conj(z2)^2 + 1080*z1*conj(z4)*conj(z5)*z2^4*conj(z1)*conj(z2)^2 -
  720*conj(z4)*conj(z5)*z2^5*conj(z1)*conj(z2)^2 - 1800*z1*conj(z4)*conj(z5)*z2^3*z3*conj(z1)*conj(z2)^2 + 1080*conj(z4)*conj(z5)*z2^4*z3*conj(z1)*conj(z2)^2 -
  1800*z1*conj(z4)*conj(z5)*z2^3*z4*conj(z1)*conj(z2)^2 + 1080*conj(z4)*conj(z5)*z2^4*z4*conj(z1)*conj(z2)^2 + 3600*z1*conj(z4)*conj(z5)*z2^2*z3*z4*conj(z1)*conj(z2)^2 -
  1800*conj(z4)*conj(z5)*z2^3*z3*z4*conj(z1)*conj(z2)^2 - 1800*z1*conj(z4)*conj(z5)*z2^3*z5*conj(z1)*conj(z2)^2 + 1080*conj(z4)*conj(z5)*z2^4*z5*conj(z1)*conj(z2)^2 +
  3600*z1*conj(z4)*conj(z5)*z2^2*z3*z5*conj(z1)*conj(z2)^2 - 1800*conj(z4)*conj(z5)*z2^3*z3*z5*conj(z1)*conj(z2)^2 + 3600*z1*conj(z4)*conj(z5)*z2^2*z4*z5*conj(z1)*conj(z2)^2 -
  1800*conj(z4)*conj(z5)*z2^3*z4*z5*conj(z1)*conj(z2)^2 + 3600*conj(z4)*conj(z5)*z2^2*z3*z4*z5*conj(z1)*conj(z2)^2 + 900*conj(z4)*conj(z5)*conj(z2)^3 -
  540*z1*conj(z4)*conj(z5)*z2^4*conj(z2)^3 + 360*conj(z4)*conj(z5)*z2^5*conj(z2)^3 + 900*z1*conj(z4)*conj(z5)*z2^3*z3*conj(z2)^3 - 540*conj(z4)*conj(z5)*z2^4*z3*conj(z2)^3 +
  900*z1*conj(z4)*conj(z5)*z2^3*z4*conj(z2)^3 - 540*conj(z4)*conj(z5)*z2^4*z4*conj(z2)^3 - 1800*z1*conj(z4)*conj(z5)*z2^2*z3*z4*conj(z2)^3 +
  900*conj(z4)*conj(z5)*z2^3*z3*z4*conj(z2)^3 + 900*z1*conj(z4)*conj(z5)*z2^3*z5*conj(z2)^3 - 540*conj(z4)*conj(z5)*z2^4*z5*conj(z2)^3 -
  1800*z1*conj(z4)*conj(z5)*z2^2*z3*z5*conj(z2)^3 + 900*conj(z4)*conj(z5)*z2^3*z3*z5*conj(z2)^3 - 1800*z1*conj(z4)*conj(z5)*z2^2*z4*z5*conj(z2)^3 +
  900*conj(z4)*conj(z5)*z2^3*z4*z5*conj(z2)^3 - 1800*conj(z4)*conj(z5)*z2^2*z3*z4*z5*conj(z2)^3 + 900*conj(z4)*conj(z1)*conj(z2)^3 + 900*conj(z5)*conj(z1)*conj(z2)^3 -
  540*z1*conj(z4)*z2^4*conj(z1)*conj(z2)^3 - 540*z1*conj(z5)*z2^4*conj(z1)*conj(z2)^3 + 360*conj(z4)*z2^5*conj(z1)*conj(z2)^3 + 360*conj(z5)*z2^5*conj(z1)*conj(z2)^3 +
  900*z1*conj(z4)*z2^3*z3*conj(z1)*conj(z2)^3 + 900*z1*conj(z5)*z2^3*z3*conj(z1)*conj(z2)^3 - 540*conj(z4)*z2^4*z3*conj(z1)*conj(z2)^3 - 540*conj(z5)*z2^4*z3*conj(z1)*conj(z2)^3 +
  900*z1*conj(z4)*z2^3*z4*conj(z1)*conj(z2)^3 + 900*z1*conj(z5)*z2^3*z4*conj(z1)*conj(z2)^3 - 540*conj(z4)*z2^4*z4*conj(z1)*conj(z2)^3 - 540*conj(z5)*z2^4*z4*conj(z1)*conj(z2)^3 -
  1800*z1*conj(z4)*z2^2*z3*z4*conj(z1)*conj(z2)^3 - 1800*z1*conj(z5)*z2^2*z3*z4*conj(z1)*conj(z2)^3 + 900*conj(z4)*z2^3*z3*z4*conj(z1)*conj(z2)^3 +
  900*conj(z5)*z2^3*z3*z4*conj(z1)*conj(z2)^3 + 900*z1*conj(z4)*z2^3*z5*conj(z1)*conj(z2)^3 + 900*z1*conj(z5)*z2^3*z5*conj(z1)*conj(z2)^3 -
  540*conj(z4)*z2^4*z5*conj(z1)*conj(z2)^3 - 540*conj(z5)*z2^4*z5*conj(z1)*conj(z2)^3 - 1800*z1*conj(z4)*z2^2*z3*z5*conj(z1)*conj(z2)^3 -
  1800*z1*conj(z5)*z2^2*z3*z5*conj(z1)*conj(z2)^3 + 900*conj(z4)*z2^3*z3*z5*conj(z1)*conj(z2)^3 + 900*conj(z5)*z2^3*z3*z5*conj(z1)*conj(z2)^3 -
  1800*z1*conj(z4)*z2^2*z4*z5*conj(z1)*conj(z2)^3 - 1800*z1*conj(z5)*z2^2*z4*z5*conj(z1)*conj(z2)^3 + 900*conj(z4)*z2^3*z4*z5*conj(z1)*conj(z2)^3 +
  900*conj(z5)*z2^3*z4*z5*conj(z1)*conj(z2)^3 - 1800*conj(z4)*z2^2*z3*z4*z5*conj(z1)*conj(z2)^3 - 1800*conj(z5)*z2^2*z3*z4*z5*conj(z1)*conj(z2)^3 - 540*conj(z4)*conj(z2)^4 -
  540*conj(z5)*conj(z2)^4 + 324*z1*conj(z4)*z2^4*conj(z2)^4 + 324*z1*conj(z5)*z2^4*conj(z2)^4 - 216*conj(z4)*z2^5*conj(z2)^4 - 216*conj(z5)*z2^5*conj(z2)^4 -
  540*z1*conj(z4)*z2^3*z3*conj(z2)^4 - 540*z1*conj(z5)*z2^3*z3*conj(z2)^4 + 324*conj(z4)*z2^4*z3*conj(z2)^4 + 324*conj(z5)*z2^4*z3*conj(z2)^4 -
  540*z1*conj(z4)*z2^3*z4*conj(z2)^4 - 540*z1*conj(z5)*z2^3*z4*conj(z2)^4 + 324*conj(z4)*z2^4*z4*conj(z2)^4 + 324*conj(z5)*z2^4*z4*conj(z2)^4 +
  1080*z1*conj(z4)*z2^2*z3*z4*conj(z2)^4 + 1080*z1*conj(z5)*z2^2*z3*z4*conj(z2)^4 - 540*conj(z4)*z2^3*z3*z4*conj(z2)^4 - 540*conj(z5)*z2^3*z3*z4*conj(z2)^4 -
  540*z1*conj(z4)*z2^3*z5*conj(z2)^4 - 540*z1*conj(z5)*z2^3*z5*conj(z2)^4 + 324*conj(z4)*z2^4*z5*conj(z2)^4 + 324*conj(z5)*z2^4*z5*conj(z2)^4 +
  1080*z1*conj(z4)*z2^2*z3*z5*conj(z2)^4 + 1080*z1*conj(z5)*z2^2*z3*z5*conj(z2)^4 - 540*conj(z4)*z2^3*z3*z5*conj(z2)^4 - 540*conj(z5)*z2^3*z3*z5*conj(z2)^4 +
  1080*z1*conj(z4)*z2^2*z4*z5*conj(z2)^4 + 1080*z1*conj(z5)*z2^2*z4*z5*conj(z2)^4 - 540*conj(z4)*z2^3*z4*z5*conj(z2)^4 - 540*conj(z5)*z2^3*z4*z5*conj(z2)^4 +
  1080*conj(z4)*z2^2*z3*z4*z5*conj(z2)^4 + 1080*conj(z5)*z2^2*z3*z4*z5*conj(z2)^4 - 540*conj(z1)*conj(z2)^4 + 324*z1*z2^4*conj(z1)*conj(z2)^4 -
  216*z2^5*conj(z1)*conj(z2)^4 - 540*z1*z2^3*z3*conj(z1)*conj(z2)^4 + 324*z2^4*z3*conj(z1)*conj(z2)^4 - 540*z1*z2^3*z4*conj(z1)*conj(z2)^4 +
  324*z2^4*z4*conj(z1)*conj(z2)^4 + 1080*z1*z2^2*z3*z4*conj(z1)*conj(z2)^4 - 540*z2^3*z3*z4*conj(z1)*conj(z2)^4 - 540*z1*z2^3*z5*conj(z1)*conj(z2)^4 +
  324*z2^4*z5*conj(z1)*conj(z2)^4 + 1080*z1*z2^2*z3*z5*conj(z1)*conj(z2)^4 - 540*z2^3*z3*z5*conj(z1)*conj(z2)^4 + 1080*z1*z2^2*z4*z5*conj(z1)*conj(z2)^4 -
  540*z2^3*z4*z5*conj(z1)*conj(z2)^4 + 1080*z2^2*z3*z4*z5*conj(z1)*conj(z2)^4 + 360*conj(z2)^5 - 216*z1*z2^4*conj(z2)^5 + 144*z2^5*conj(z2)^5 +
  360*z1*z2^3*z3*conj(z2)^5 - 216*z2^4*z3*conj(z2)^5 + 360*z1*z2^3*z4*conj(z2)^5 - 216*z2^4*z4*conj(z2)^5 - 720*z1*z2^2*z3*z4*conj(z2)^5 +
  360*z2^3*z3*z4*conj(z2)^5 + 360*z1*z2^3*z5*conj(z2)^5 - 216*z2^4*z5*conj(z2)^5 - 720*z1*z2^2*z3*z5*conj(z2)^5 + 360*z2^3*z3*z5*conj(z2)^5 -
  720*z1*z2^2*z4*z5*conj(z2)^5 + 360*z2^3*z4*z5*conj(z2)^5 - 720*z2^2*z3*z4*z5*conj(z2)^5 - 1800*conj(z4)*conj(z5)*conj(z2)^2*conj(z3) +
  1080*z1*conj(z4)*conj(z5)*z2^4*conj(z2)^2*conj(z3) - 720*conj(z4)*conj(z5)*z2^5*conj(z2)^2*conj(z3) - 1800*z1*conj(z4)*conj(z5)*z2^3*z3*conj(z2)^2*conj(z3) +
  1080*conj(z4)*conj(z5)*z2^4*z3*conj(z2)^2*conj(z3) - 1800*z1*conj(z4)*conj(z5)*z2^3*z4*conj(z2)^2*conj(z3) + 1080*conj(z4)*conj(z5)*z2^4*z4*conj(z2)^2*conj(z3) +
  3600*z1*conj(z4)*conj(z5)*z2^2*z3*z4*conj(z2)^2*conj(z3) - 1800*conj(z4)*conj(z5)*z2^3*z3*z4*conj(z2)^2*conj(z3) - 1800*z1*conj(z4)*conj(z5)*z2^3*z5*conj(z2)^2*conj(z3) +
  1080*conj(z4)*conj(z5)*z2^4*z5*conj(z2)^2*conj(z3) + 3600*z1*conj(z4)*conj(z5)*z2^2*z3*z5*conj(z2)^2*conj(z3) - 1800*conj(z4)*conj(z5)*z2^3*z3*z5*conj(z2)^2*conj(z3) +
  3600*z1*conj(z4)*conj(z5)*z2^2*z4*z5*conj(z2)^2*conj(z3) - 1800*conj(z4)*conj(z5)*z2^3*z4*z5*conj(z2)^2*conj(z3) + 3600*conj(z4)*conj(z5)*z2^2*z3*z4*z5*conj(z2)^2*conj(z3) -
  1800*conj(z4)*conj(z1)*conj(z2)^2*conj(z3) - 1800*conj(z5)*conj(z1)*conj(z2)^2*conj(z3) + 1080*z1*conj(z4)*z2^4*conj(z1)*conj(z2)^2*conj(z3) + 1080*z1*conj(z5)*z2^4*conj(z1)*conj(z2)^2*conj(z3) -
  720*conj(z4)*z2^5*conj(z1)*conj(z2)^2*conj(z3) - 720*conj(z5)*z2^5*conj(z1)*conj(z2)^2*conj(z3) - 1800*z1*conj(z4)*z2^3*z3*conj(z1)*conj(z2)^2*conj(z3) -
  1800*z1*conj(z5)*z2^3*z3*conj(z1)*conj(z2)^2*conj(z3) + 1080*conj(z4)*z2^4*z3*conj(z1)*conj(z2)^2*conj(z3) + 1080*conj(z5)*z2^4*z3*conj(z1)*conj(z2)^2*conj(z3) -
  1800*z1*conj(z4)*z2^3*z4*conj(z1)*conj(z2)^2*conj(z3) - 1800*z1*conj(z5)*z2^3*z4*conj(z1)*conj(z2)^2*conj(z3) + 1080*conj(z4)*z2^4*z4*conj(z1)*conj(z2)^2*conj(z3) +
  1080*conj(z5)*z2^4*z4*conj(z1)*conj(z2)^2*conj(z3) + 3600*z1*conj(z4)*z2^2*z3*z4*conj(z1)*conj(z2)^2*conj(z3) + 3600*z1*conj(z5)*z2^2*z3*z4*conj(z1)*conj(z2)^2*conj(z3) -
  1800*conj(z4)*z2^3*z3*z4*conj(z1)*conj(z2)^2*conj(z3) - 1800*conj(z5)*z2^3*z3*z4*conj(z1)*conj(z2)^2*conj(z3) - 1800*z1*conj(z4)*z2^3*z5*conj(z1)*conj(z2)^2*conj(z3) -
  1800*z1*conj(z5)*z2^3*z5*conj(z1)*conj(z2)^2*conj(z3) + 1080*conj(z4)*z2^4*z5*conj(z1)*conj(z2)^2*conj(z3) + 1080*conj(z5)*z2^4*z5*conj(z1)*conj(z2)^2*conj(z3) +
  3600*z1*conj(z4)*z2^2*z3*z5*conj(z1)*conj(z2)^2*conj(z3) + 3600*z1*conj(z5)*z2^2*z3*z5*conj(z1)*conj(z2)^2*conj(z3) - 1800*conj(z4)*z2^3*z3*z5*conj(z1)*conj(z2)^2*conj(z3) -
  1800*conj(z5)*z2^3*z3*z5*conj(z1)*conj(z2)^2*conj(z3) + 3600*z1*conj(z4)*z2^2*z4*z5*conj(z1)*conj(z2)^2*conj(z3) + 3600*z1*conj(z5)*z2^2*z4*z5*conj(z1)*conj(z2)^2*conj(z3) -
  1800*conj(z4)*z2^3*z4*z5*conj(z1)*conj(z2)^2*conj(z3) - 1800*conj(z5)*z2^3*z4*z5*conj(z1)*conj(z2)^2*conj(z3) + 3600*conj(z4)*z2^2*z3*z4*z5*conj(z1)*conj(z2)^2*conj(z3) +
  3600*conj(z5)*z2^2*z3*z4*z5*conj(z1)*conj(z2)^2*conj(z3) + 900*conj(z4)*conj(z2)^3*conj(z3) + 900*conj(z5)*conj(z2)^3*conj(z3) - 540*z1*conj(z4)*z2^4*conj(z2)^3*conj(z3) -
  540*z1*conj(z5)*z2^4*conj(z2)^3*conj(z3) + 360*conj(z4)*z2^5*conj(z2)^3*conj(z3) + 360*conj(z5)*z2^5*conj(z2)^3*conj(z3) + 900*z1*conj(z4)*z2^3*z3*conj(z2)^3*conj(z3) +
  900*z1*conj(z5)*z2^3*z3*conj(z2)^3*conj(z3) - 540*conj(z4)*z2^4*z3*conj(z2)^3*conj(z3) - 540*conj(z5)*z2^4*z3*conj(z2)^3*conj(z3) + 900*z1*conj(z4)*z2^3*z4*conj(z2)^3*conj(z3) +
  900*z1*conj(z5)*z2^3*z4*conj(z2)^3*conj(z3) - 540*conj(z4)*z2^4*z4*conj(z2)^3*conj(z3) - 540*conj(z5)*z2^4*z4*conj(z2)^3*conj(z3) -
  1800*z1*conj(z4)*z2^2*z3*z4*conj(z2)^3*conj(z3) - 1800*z1*conj(z5)*z2^2*z3*z4*conj(z2)^3*conj(z3) + 900*conj(z4)*z2^3*z3*z4*conj(z2)^3*conj(z3) +
  900*conj(z5)*z2^3*z3*z4*conj(z2)^3*conj(z3) + 900*z1*conj(z4)*z2^3*z5*conj(z2)^3*conj(z3) + 900*z1*conj(z5)*z2^3*z5*conj(z2)^3*conj(z3) -
  540*conj(z4)*z2^4*z5*conj(z2)^3*conj(z3) - 540*conj(z5)*z2^4*z5*conj(z2)^3*conj(z3) - 1800*z1*conj(z4)*z2^2*z3*z5*conj(z2)^3*conj(z3) -
  1800*z1*conj(z5)*z2^2*z3*z5*conj(z2)^3*conj(z3) + 900*conj(z4)*z2^3*z3*z5*conj(z2)^3*conj(z3) + 900*conj(z5)*z2^3*z3*z5*conj(z2)^3*conj(z3) -
  1800*z1*conj(z4)*z2^2*z4*z5*conj(z2)^3*conj(z3) - 1800*z1*conj(z5)*z2^2*z4*z5*conj(z2)^3*conj(z3) + 900*conj(z4)*z2^3*z4*z5*conj(z2)^3*conj(z3) +
  900*conj(z5)*z2^3*z4*z5*conj(z2)^3*conj(z3) - 1800*conj(z4)*z2^2*z3*z4*z5*conj(z2)^3*conj(z3) - 1800*conj(z5)*z2^2*z3*z4*z5*conj(z2)^3*conj(z3) + 900*conj(z1)*conj(z2)^3*conj(z3) -
  540*z1*z2^4*conj(z1)*conj(z2)^3*conj(z3) + 360*z2^5*conj(z1)*conj(z2)^3*conj(z3) + 900*z1*z2^3*z3*conj(z1)*conj(z2)^3*conj(z3) - 540*z2^4*z3*conj(z1)*conj(z2)^3*conj(z3) +
  900*z1*z2^3*z4*conj(z1)*conj(z2)^3*conj(z3) - 540*z2^4*z4*conj(z1)*conj(z2)^3*conj(z3) - 1800*z1*z2^2*z3*z4*conj(z1)*conj(z2)^3*conj(z3) +
  900*z2^3*z3*z4*conj(z1)*conj(z2)^3*conj(z3) + 900*z1*z2^3*z5*conj(z1)*conj(z2)^3*conj(z3) - 540*z2^4*z5*conj(z1)*conj(z2)^3*conj(z3) -
  1800*z1*z2^2*z3*z5*conj(z1)*conj(z2)^3*conj(z3) + 900*z2^3*z3*z5*conj(z1)*conj(z2)^3*conj(z3) - 1800*z1*z2^2*z4*z5*conj(z1)*conj(z2)^3*conj(z3) +
  900*z2^3*z4*z5*conj(z1)*conj(z2)^3*conj(z3) - 1800*z2^2*z3*z4*z5*conj(z1)*conj(z2)^3*conj(z3) - 540*conj(z2)^4*conj(z3) + 324*z1*z2^4*conj(z2)^4*conj(z3) -
  216*z2^5*conj(z2)^4*conj(z3) - 540*z1*z2^3*z3*conj(z2)^4*conj(z3) + 324*z2^4*z3*conj(z2)^4*conj(z3) - 540*z1*z2^3*z4*conj(z2)^4*conj(z3) +
  324*z2^4*z4*conj(z2)^4*conj(z3) + 1080*z1*z2^2*z3*z4*conj(z2)^4*conj(z3) - 540*z2^3*z3*z4*conj(z2)^4*conj(z3) - 540*z1*z2^3*z5*conj(z2)^4*conj(z3) +
  324*z2^4*z5*conj(z2)^4*conj(z3) + 1080*z1*z2^2*z3*z5*conj(z2)^4*conj(z3) - 540*z2^3*z3*z5*conj(z2)^4*conj(z3) + 1080*z1*z2^2*z4*z5*conj(z2)^4*conj(z3) -
  540*z2^3*z4*z5*conj(z2)^4*conj(z3) + 1080*z2^2*z3*z4*z5*conj(z2)^4*conj(z3) - 3600*z6*conj(z6)
g3 = 900 + 900*z1*z2*z3^3 - 540*z1*z3^4 - 540*z2*z3^4 + 360*z3^5 - 1800*z1*z2*z3^2*z4 + 900*z1*z3^3*z4 +
  900*z2*z3^3*z4 - 540*z3^4*z4 - 1800*z1*z2*z3^2*z5 + 900*z1*z3^3*z5 + 900*z2*z3^3*z5 - 540*z3^4*z5 -
  1800*z1*z3^2*z4*z5 - 1800*z2*z3^2*z4*z5 + 900*z3^3*z4*z5 - 1800*conj(z4)*conj(z5)*conj(z1)*conj(z3)^2 -
  1800*z1*conj(z4)*conj(z5)*z2*z3^3*conj(z1)*conj(z3)^2 + 1080*z1*conj(z4)*conj(z5)*z3^4*conj(z1)*conj(z3)^2 + 1080*conj(z4)*conj(z5)*z2*z3^4*conj(z1)*conj(z3)^2 -
  720*conj(z4)*conj(z5)*z3^5*conj(z1)*conj(z3)^2 + 3600*z1*conj(z4)*conj(z5)*z2*z3^2*z4*conj(z1)*conj(z3)^2 - 1800*z1*conj(z4)*conj(z5)*z3^3*z4*conj(z1)*conj(z3)^2 -
  1800*conj(z4)*conj(z5)*z2*z3^3*z4*conj(z1)*conj(z3)^2 + 1080*conj(z4)*conj(z5)*z3^4*z4*conj(z1)*conj(z3)^2 + 3600*z1*conj(z4)*conj(z5)*z2*z3^2*z5*conj(z1)*conj(z3)^2 -
  1800*z1*conj(z4)*conj(z5)*z3^3*z5*conj(z1)*conj(z3)^2 - 1800*conj(z4)*conj(z5)*z2*z3^3*z5*conj(z1)*conj(z3)^2 + 1080*conj(z4)*conj(z5)*z3^4*z5*conj(z1)*conj(z3)^2 +
  3600*z1*conj(z4)*conj(z5)*z3^2*z4*z5*conj(z1)*conj(z3)^2 + 3600*conj(z4)*conj(z5)*z2*z3^2*z4*z5*conj(z1)*conj(z3)^2 - 1800*conj(z4)*conj(z5)*z3^3*z4*z5*conj(z1)*conj(z3)^2 -
  1800*conj(z4)*conj(z5)*conj(z2)*conj(z3)^2 - 1800*z1*conj(z4)*conj(z5)*z2*z3^3*conj(z2)*conj(z3)^2 + 1080*z1*conj(z4)*conj(z5)*z3^4*conj(z2)*conj(z3)^2 +
  1080*conj(z4)*conj(z5)*z2*z3^4*conj(z2)*conj(z3)^2 - 720*conj(z4)*conj(z5)*z3^5*conj(z2)*conj(z3)^2 + 3600*z1*conj(z4)*conj(z5)*z2*z3^2*z4*conj(z2)*conj(z3)^2 -
  1800*z1*conj(z4)*conj(z5)*z3^3*z4*conj(z2)*conj(z3)^2 - 1800*conj(z4)*conj(z5)*z2*z3^3*z4*conj(z2)*conj(z3)^2 + 1080*conj(z4)*conj(z5)*z3^4*z4*conj(z2)*conj(z3)^2 +
  3600*z1*conj(z4)*conj(z5)*z2*z3^2*z5*conj(z2)*conj(z3)^2 - 1800*z1*conj(z4)*conj(z5)*z3^3*z5*conj(z2)*conj(z3)^2 - 1800*conj(z4)*conj(z5)*z2*z3^3*z5*conj(z2)*conj(z3)^2 +
  1080*conj(z4)*conj(z5)*z3^4*z5*conj(z2)*conj(z3)^2 + 3600*z1*conj(z4)*conj(z5)*z3^2*z4*z5*conj(z2)*conj(z3)^2 + 3600*conj(z4)*conj(z5)*z2*z3^2*z4*z5*conj(z2)*conj(z3)^2 -
  1800*conj(z4)*conj(z5)*z3^3*z4*z5*conj(z2)*conj(z3)^2 - 1800*conj(z4)*conj(z1)*conj(z2)*conj(z3)^2 - 1800*conj(z5)*conj(z1)*conj(z2)*conj(z3)^2 - 1800*z1*conj(z4)*z2*z3^3*conj(z1)*conj(z2)*conj(z3)^2 -
  1800*z1*conj(z5)*z2*z3^3*conj(z1)*conj(z2)*conj(z3)^2 + 1080*z1*conj(z4)*z3^4*conj(z1)*conj(z2)*conj(z3)^2 + 1080*z1*conj(z5)*z3^4*conj(z1)*conj(z2)*conj(z3)^2 +
  1080*conj(z4)*z2*z3^4*conj(z1)*conj(z2)*conj(z3)^2 + 1080*conj(z5)*z2*z3^4*conj(z1)*conj(z2)*conj(z3)^2 - 720*conj(z4)*z3^5*conj(z1)*conj(z2)*conj(z3)^2 - 720*conj(z5)*z3^5*conj(z1)*conj(z2)*conj(z3)^2 +
  3600*z1*conj(z4)*z2*z3^2*z4*conj(z1)*conj(z2)*conj(z3)^2 + 3600*z1*conj(z5)*z2*z3^2*z4*conj(z1)*conj(z2)*conj(z3)^2 - 1800*z1*conj(z4)*z3^3*z4*conj(z1)*conj(z2)*conj(z3)^2 -
  1800*z1*conj(z5)*z3^3*z4*conj(z1)*conj(z2)*conj(z3)^2 - 1800*conj(z4)*z2*z3^3*z4*conj(z1)*conj(z2)*conj(z3)^2 - 1800*conj(z5)*z2*z3^3*z4*conj(z1)*conj(z2)*conj(z3)^2 +
  1080*conj(z4)*z3^4*z4*conj(z1)*conj(z2)*conj(z3)^2 + 1080*conj(z5)*z3^4*z4*conj(z1)*conj(z2)*conj(z3)^2 + 3600*z1*conj(z4)*z2*z3^2*z5*conj(z1)*conj(z2)*conj(z3)^2 +
  3600*z1*conj(z5)*z2*z3^2*z5*conj(z1)*conj(z2)*conj(z3)^2 - 1800*z1*conj(z4)*z3^3*z5*conj(z1)*conj(z2)*conj(z3)^2 - 1800*z1*conj(z5)*z3^3*z5*conj(z1)*conj(z2)*conj(z3)^2 -
  1800*conj(z4)*z2*z3^3*z5*conj(z1)*conj(z2)*conj(z3)^2 - 1800*conj(z5)*z2*z3^3*z5*conj(z1)*conj(z2)*conj(z3)^2 + 1080*conj(z4)*z3^4*z5*conj(z1)*conj(z2)*conj(z3)^2 +
  1080*conj(z5)*z3^4*z5*conj(z1)*conj(z2)*conj(z3)^2 + 3600*z1*conj(z4)*z3^2*z4*z5*conj(z1)*conj(z2)*conj(z3)^2 + 3600*z1*conj(z5)*z3^2*z4*z5*conj(z1)*conj(z2)*conj(z3)^2 +
  3600*conj(z4)*z2*z3^2*z4*z5*conj(z1)*conj(z2)*conj(z3)^2 + 3600*conj(z5)*z2*z3^2*z4*z5*conj(z1)*conj(z2)*conj(z3)^2 - 1800*conj(z4)*z3^3*z4*z5*conj(z1)*conj(z2)*conj(z3)^2 -
  1800*conj(z5)*z3^3*z4*z5*conj(z1)*conj(z2)*conj(z3)^2 + 900*conj(z4)*conj(z5)*conj(z3)^3 + 900*z1*conj(z4)*conj(z5)*z2*z3^3*conj(z3)^3 - 540*z1*conj(z4)*conj(z5)*z3^4*conj(z3)^3 -
  540*conj(z4)*conj(z5)*z2*z3^4*conj(z3)^3 + 360*conj(z4)*conj(z5)*z3^5*conj(z3)^3 - 1800*z1*conj(z4)*conj(z5)*z2*z3^2*z4*conj(z3)^3 +
  900*z1*conj(z4)*conj(z5)*z3^3*z4*conj(z3)^3 + 900*conj(z4)*conj(z5)*z2*z3^3*z4*conj(z3)^3 - 540*conj(z4)*conj(z5)*z3^4*z4*conj(z3)^3 -
  1800*z1*conj(z4)*conj(z5)*z2*z3^2*z5*conj(z3)^3 + 900*z1*conj(z4)*conj(z5)*z3^3*z5*conj(z3)^3 + 900*conj(z4)*conj(z5)*z2*z3^3*z5*conj(z3)^3 -
  540*conj(z4)*conj(z5)*z3^4*z5*conj(z3)^3 - 1800*z1*conj(z4)*conj(z5)*z3^2*z4*z5*conj(z3)^3 - 1800*conj(z4)*conj(z5)*z2*z3^2*z4*z5*conj(z3)^3 +
  900*conj(z4)*conj(z5)*z3^3*z4*z5*conj(z3)^3 + 900*conj(z4)*conj(z1)*conj(z3)^3 + 900*conj(z5)*conj(z1)*conj(z3)^3 + 900*z1*conj(z4)*z2*z3^3*conj(z1)*conj(z3)^3 +
  900*z1*conj(z5)*z2*z3^3*conj(z1)*conj(z3)^3 - 540*z1*conj(z4)*z3^4*conj(z1)*conj(z3)^3 - 540*z1*conj(z5)*z3^4*conj(z1)*conj(z3)^3 - 540*conj(z4)*z2*z3^4*conj(z1)*conj(z3)^3 -
  540*conj(z5)*z2*z3^4*conj(z1)*conj(z3)^3 + 360*conj(z4)*z3^5*conj(z1)*conj(z3)^3 + 360*conj(z5)*z3^5*conj(z1)*conj(z3)^3 - 1800*z1*conj(z4)*z2*z3^2*z4*conj(z1)*conj(z3)^3 -
  1800*z1*conj(z5)*z2*z3^2*z4*conj(z1)*conj(z3)^3 + 900*z1*conj(z4)*z3^3*z4*conj(z1)*conj(z3)^3 + 900*z1*conj(z5)*z3^3*z4*conj(z1)*conj(z3)^3 +
  900*conj(z4)*z2*z3^3*z4*conj(z1)*conj(z3)^3 + 900*conj(z5)*z2*z3^3*z4*conj(z1)*conj(z3)^3 - 540*conj(z4)*z3^4*z4*conj(z1)*conj(z3)^3 - 540*conj(z5)*z3^4*z4*conj(z1)*conj(z3)^3 -
  1800*z1*conj(z4)*z2*z3^2*z5*conj(z1)*conj(z3)^3 - 1800*z1*conj(z5)*z2*z3^2*z5*conj(z1)*conj(z3)^3 + 900*z1*conj(z4)*z3^3*z5*conj(z1)*conj(z3)^3 +
  900*z1*conj(z5)*z3^3*z5*conj(z1)*conj(z3)^3 + 900*conj(z4)*z2*z3^3*z5*conj(z1)*conj(z3)^3 + 900*conj(z5)*z2*z3^3*z5*conj(z1)*conj(z3)^3 -
  540*conj(z4)*z3^4*z5*conj(z1)*conj(z3)^3 - 540*conj(z5)*z3^4*z5*conj(z1)*conj(z3)^3 - 1800*z1*conj(z4)*z3^2*z4*z5*conj(z1)*conj(z3)^3 -
  1800*z1*conj(z5)*z3^2*z4*z5*conj(z1)*conj(z3)^3 - 1800*conj(z4)*z2*z3^2*z4*z5*conj(z1)*conj(z3)^3 - 1800*conj(z5)*z2*z3^2*z4*z5*conj(z1)*conj(z3)^3 +
  900*conj(z4)*z3^3*z4*z5*conj(z1)*conj(z3)^3 + 900*conj(z5)*z3^3*z4*z5*conj(z1)*conj(z3)^3 + 900*conj(z4)*conj(z2)*conj(z3)^3 + 900*conj(z5)*conj(z2)*conj(z3)^3 +
  900*z1*conj(z4)*z2*z3^3*conj(z2)*conj(z3)^3 + 900*z1*conj(z5)*z2*z3^3*conj(z2)*conj(z3)^3 - 540*z1*conj(z4)*z3^4*conj(z2)*conj(z3)^3 - 540*z1*conj(z5)*z3^4*conj(z2)*conj(z3)^3 -
  540*conj(z4)*z2*z3^4*conj(z2)*conj(z3)^3 - 540*conj(z5)*z2*z3^4*conj(z2)*conj(z3)^3 + 360*conj(z4)*z3^5*conj(z2)*conj(z3)^3 + 360*conj(z5)*z3^5*conj(z2)*conj(z3)^3 -
  1800*z1*conj(z4)*z2*z3^2*z4*conj(z2)*conj(z3)^3 - 1800*z1*conj(z5)*z2*z3^2*z4*conj(z2)*conj(z3)^3 + 900*z1*conj(z4)*z3^3*z4*conj(z2)*conj(z3)^3 +
  900*z1*conj(z5)*z3^3*z4*conj(z2)*conj(z3)^3 + 900*conj(z4)*z2*z3^3*z4*conj(z2)*conj(z3)^3 + 900*conj(z5)*z2*z3^3*z4*conj(z2)*conj(z3)^3 -
  540*conj(z4)*z3^4*z4*conj(z2)*conj(z3)^3 - 540*conj(z5)*z3^4*z4*conj(z2)*conj(z3)^3 - 1800*z1*conj(z4)*z2*z3^2*z5*conj(z2)*conj(z3)^3 -
  1800*z1*conj(z5)*z2*z3^2*z5*conj(z2)*conj(z3)^3 + 900*z1*conj(z4)*z3^3*z5*conj(z2)*conj(z3)^3 + 900*z1*conj(z5)*z3^3*z5*conj(z2)*conj(z3)^3 +
  900*conj(z4)*z2*z3^3*z5*conj(z2)*conj(z3)^3 + 900*conj(z5)*z2*z3^3*z5*conj(z2)*conj(z3)^3 - 540*conj(z4)*z3^4*z5*conj(z2)*conj(z3)^3 - 540*conj(z5)*z3^4*z5*conj(z2)*conj(z3)^3 -
  1800*z1*conj(z4)*z3^2*z4*z5*conj(z2)*conj(z3)^3 - 1800*z1*conj(z5)*z3^2*z4*z5*conj(z2)*conj(z3)^3 - 1800*conj(z4)*z2*z3^2*z4*z5*conj(z2)*conj(z3)^3 -
  1800*conj(z5)*z2*z3^2*z4*z5*conj(z2)*conj(z3)^3 + 900*conj(z4)*z3^3*z4*z5*conj(z2)*conj(z3)^3 + 900*conj(z5)*z3^3*z4*z5*conj(z2)*conj(z3)^3 + 900*conj(z1)*conj(z2)*conj(z3)^3 +
  900*z1*z2*z3^3*conj(z1)*conj(z2)*conj(z3)^3 - 540*z1*z3^4*conj(z1)*conj(z2)*conj(z3)^3 - 540*z2*z3^4*conj(z1)*conj(z2)*conj(z3)^3 + 360*z3^5*conj(z1)*conj(z2)*conj(z3)^3 -
  1800*z1*z2*z3^2*z4*conj(z1)*conj(z2)*conj(z3)^3 + 900*z1*z3^3*z4*conj(z1)*conj(z2)*conj(z3)^3 + 900*z2*z3^3*z4*conj(z1)*conj(z2)*conj(z3)^3 -
  540*z3^4*z4*conj(z1)*conj(z2)*conj(z3)^3 - 1800*z1*z2*z3^2*z5*conj(z1)*conj(z2)*conj(z3)^3 + 900*z1*z3^3*z5*conj(z1)*conj(z2)*conj(z3)^3 +
  900*z2*z3^3*z5*conj(z1)*conj(z2)*conj(z3)^3 - 540*z3^4*z5*conj(z1)*conj(z2)*conj(z3)^3 - 1800*z1*z3^2*z4*z5*conj(z1)*conj(z2)*conj(z3)^3 -
  1800*z2*z3^2*z4*z5*conj(z1)*conj(z2)*conj(z3)^3 + 900*z3^3*z4*z5*conj(z1)*conj(z2)*conj(z3)^3 - 540*conj(z4)*conj(z3)^4 - 540*conj(z5)*conj(z3)^4 -
  540*z1*conj(z4)*z2*z3^3*conj(z3)^4 - 540*z1*conj(z5)*z2*z3^3*conj(z3)^4 + 324*z1*conj(z4)*z3^4*conj(z3)^4 + 324*z1*conj(z5)*z3^4*conj(z3)^4 +
  324*conj(z4)*z2*z3^4*conj(z3)^4 + 324*conj(z5)*z2*z3^4*conj(z3)^4 - 216*conj(z4)*z3^5*conj(z3)^4 - 216*conj(z5)*z3^5*conj(z3)^4 +
  1080*z1*conj(z4)*z2*z3^2*z4*conj(z3)^4 + 1080*z1*conj(z5)*z2*z3^2*z4*conj(z3)^4 - 540*z1*conj(z4)*z3^3*z4*conj(z3)^4 - 540*z1*conj(z5)*z3^3*z4*conj(z3)^4 -
  540*conj(z4)*z2*z3^3*z4*conj(z3)^4 - 540*conj(z5)*z2*z3^3*z4*conj(z3)^4 + 324*conj(z4)*z3^4*z4*conj(z3)^4 + 324*conj(z5)*z3^4*z4*conj(z3)^4 +
  1080*z1*conj(z4)*z2*z3^2*z5*conj(z3)^4 + 1080*z1*conj(z5)*z2*z3^2*z5*conj(z3)^4 - 540*z1*conj(z4)*z3^3*z5*conj(z3)^4 - 540*z1*conj(z5)*z3^3*z5*conj(z3)^4 -
  540*conj(z4)*z2*z3^3*z5*conj(z3)^4 - 540*conj(z5)*z2*z3^3*z5*conj(z3)^4 + 324*conj(z4)*z3^4*z5*conj(z3)^4 + 324*conj(z5)*z3^4*z5*conj(z3)^4 +
  1080*z1*conj(z4)*z3^2*z4*z5*conj(z3)^4 + 1080*z1*conj(z5)*z3^2*z4*z5*conj(z3)^4 + 1080*conj(z4)*z2*z3^2*z4*z5*conj(z3)^4 +
  1080*conj(z5)*z2*z3^2*z4*z5*conj(z3)^4 - 540*conj(z4)*z3^3*z4*z5*conj(z3)^4 - 540*conj(z5)*z3^3*z4*z5*conj(z3)^4 - 540*conj(z1)*conj(z3)^4 -
  540*z1*z2*z3^3*conj(z1)*conj(z3)^4 + 324*z1*z3^4*conj(z1)*conj(z3)^4 + 324*z2*z3^4*conj(z1)*conj(z3)^4 - 216*z3^5*conj(z1)*conj(z3)^4 +
  1080*z1*z2*z3^2*z4*conj(z1)*conj(z3)^4 - 540*z1*z3^3*z4*conj(z1)*conj(z3)^4 - 540*z2*z3^3*z4*conj(z1)*conj(z3)^4 + 324*z3^4*z4*conj(z1)*conj(z3)^4 +
  1080*z1*z2*z3^2*z5*conj(z1)*conj(z3)^4 - 540*z1*z3^3*z5*conj(z1)*conj(z3)^4 - 540*z2*z3^3*z5*conj(z1)*conj(z3)^4 + 324*z3^4*z5*conj(z1)*conj(z3)^4 +
  1080*z1*z3^2*z4*z5*conj(z1)*conj(z3)^4 + 1080*z2*z3^2*z4*z5*conj(z1)*conj(z3)^4 - 540*z3^3*z4*z5*conj(z1)*conj(z3)^4 - 540*conj(z2)*conj(z3)^4 -
  540*z1*z2*z3^3*conj(z2)*conj(z3)^4 + 324*z1*z3^4*conj(z2)*conj(z3)^4 + 324*z2*z3^4*conj(z2)*conj(z3)^4 - 216*z3^5*conj(z2)*conj(z3)^4 +
  1080*z1*z2*z3^2*z4*conj(z2)*conj(z3)^4 - 540*z1*z3^3*z4*conj(z2)*conj(z3)^4 - 540*z2*z3^3*z4*conj(z2)*conj(z3)^4 + 324*z3^4*z4*conj(z2)*conj(z3)^4 +
  1080*z1*z2*z3^2*z5*conj(z2)*conj(z3)^4 - 540*z1*z3^3*z5*conj(z2)*conj(z3)^4 - 540*z2*z3^3*z5*conj(z2)*conj(z3)^4 + 324*z3^4*z5*conj(z2)*conj(z3)^4 +
  1080*z1*z3^2*z4*z5*conj(z2)*conj(z3)^4 + 1080*z2*z3^2*z4*z5*conj(z2)*conj(z3)^4 - 540*z3^3*z4*z5*conj(z2)*conj(z3)^4 + 360*conj(z3)^5 +
  360*z1*z2*z3^3*conj(z3)^5 - 216*z1*z3^4*conj(z3)^5 - 216*z2*z3^4*conj(z3)^5 + 144*z3^5*conj(z3)^5 - 720*z1*z2*z3^2*z4*conj(z3)^5 +
  360*z1*z3^3*z4*conj(z3)^5 + 360*z2*z3^3*z4*conj(z3)^5 - 216*z3^4*z4*conj(z3)^5 - 720*z1*z2*z3^2*z5*conj(z3)^5 + 360*z1*z3^3*z5*conj(z3)^5 +
  360*z2*z3^3*z5*conj(z3)^5 - 216*z3^4*z5*conj(z3)^5 - 720*z1*z3^2*z4*z5*conj(z3)^5 - 720*z2*z3^2*z4*z5*conj(z3)^5 + 360*z3^3*z4*z5*conj(z3)^5 - 3600*z6*conj(z6)
g4 = 900 + 360*conj(z4)^5 - 540*conj(z4)^4*conj(z5) - 1800*z1*z2*z3*z4^2 - 720*z1*conj(z4)^5*z2*z3*z4^2 + 1080*z1*conj(z4)^4*conj(z5)*z2*z3*z4^2 +
  900*z1*z2*z4^3 + 360*z1*conj(z4)^5*z2*z4^3 - 540*z1*conj(z4)^4*conj(z5)*z2*z4^3 + 900*z1*z3*z4^3 + 360*z1*conj(z4)^5*z3*z4^3 -
  540*z1*conj(z4)^4*conj(z5)*z3*z4^3 + 900*z2*z3*z4^3 + 360*conj(z4)^5*z2*z3*z4^3 - 540*conj(z4)^4*conj(z5)*z2*z3*z4^3 - 540*z1*z4^4 -
  216*z1*conj(z4)^5*z4^4 + 324*z1*conj(z4)^4*conj(z5)*z4^4 - 540*z2*z4^4 - 216*conj(z4)^5*z2*z4^4 + 324*conj(z4)^4*conj(z5)*z2*z4^4 -
  540*z3*z4^4 - 216*conj(z4)^5*z3*z4^4 + 324*conj(z4)^4*conj(z5)*z3*z4^4 + 360*z4^5 + 144*conj(z4)^5*z4^5 - 216*conj(z4)^4*conj(z5)*z4^5 -
  1800*z1*z2*z4^2*z5 - 720*z1*conj(z4)^5*z2*z4^2*z5 + 1080*z1*conj(z4)^4*conj(z5)*z2*z4^2*z5 - 1800*z1*z3*z4^2*z5 -
  720*z1*conj(z4)^5*z3*z4^2*z5 + 1080*z1*conj(z4)^4*conj(z5)*z3*z4^2*z5 - 1800*z2*z3*z4^2*z5 - 720*conj(z4)^5*z2*z3*z4^2*z5 +
  1080*conj(z4)^4*conj(z5)*z2*z3*z4^2*z5 + 900*z1*z4^3*z5 + 360*z1*conj(z4)^5*z4^3*z5 - 540*z1*conj(z4)^4*conj(z5)*z4^3*z5 +
  900*z2*z4^3*z5 + 360*conj(z4)^5*z2*z4^3*z5 - 540*conj(z4)^4*conj(z5)*z2*z4^3*z5 + 900*z3*z4^3*z5 + 360*conj(z4)^5*z3*z4^3*z5 -
  540*conj(z4)^4*conj(z5)*z3*z4^3*z5 - 540*z4^4*z5 - 216*conj(z4)^5*z4^4*z5 + 324*conj(z4)^4*conj(z5)*z4^4*z5 - 540*conj(z4)^4*conj(z1) +
  900*conj(z4)^3*conj(z5)*conj(z1) + 1080*z1*conj(z4)^4*z2*z3*z4^2*conj(z1) - 1800*z1*conj(z4)^3*conj(z5)*z2*z3*z4^2*conj(z1) - 540*z1*conj(z4)^4*z2*z4^3*conj(z1) +
  900*z1*conj(z4)^3*conj(z5)*z2*z4^3*conj(z1) - 540*z1*conj(z4)^4*z3*z4^3*conj(z1) + 900*z1*conj(z4)^3*conj(z5)*z3*z4^3*conj(z1) - 540*conj(z4)^4*z2*z3*z4^3*conj(z1) +
  900*conj(z4)^3*conj(z5)*z2*z3*z4^3*conj(z1) + 324*z1*conj(z4)^4*z4^4*conj(z1) - 540*z1*conj(z4)^3*conj(z5)*z4^4*conj(z1) + 324*conj(z4)^4*z2*z4^4*conj(z1) -
  540*conj(z4)^3*conj(z5)*z2*z4^4*conj(z1) + 324*conj(z4)^4*z3*z4^4*conj(z1) - 540*conj(z4)^3*conj(z5)*z3*z4^4*conj(z1) - 216*conj(z4)^4*z4^5*conj(z1) +
  360*conj(z4)^3*conj(z5)*z4^5*conj(z1) + 1080*z1*conj(z4)^4*z2*z4^2*z5*conj(z1) - 1800*z1*conj(z4)^3*conj(z5)*z2*z4^2*z5*conj(z1) +
  1080*z1*conj(z4)^4*z3*z4^2*z5*conj(z1) - 1800*z1*conj(z4)^3*conj(z5)*z3*z4^2*z5*conj(z1) + 1080*conj(z4)^4*z2*z3*z4^2*z5*conj(z1) -
  1800*conj(z4)^3*conj(z5)*z2*z3*z4^2*z5*conj(z1) - 540*z1*conj(z4)^4*z4^3*z5*conj(z1) + 900*z1*conj(z4)^3*conj(z5)*z4^3*z5*conj(z1) -
  540*conj(z4)^4*z2*z4^3*z5*conj(z1) + 900*conj(z4)^3*conj(z5)*z2*z4^3*z5*conj(z1) - 540*conj(z4)^4*z3*z4^3*z5*conj(z1) + 900*conj(z4)^3*conj(z5)*z3*z4^3*z5*conj(z1) +
  324*conj(z4)^4*z4^4*z5*conj(z1) - 540*conj(z4)^3*conj(z5)*z4^4*z5*conj(z1) - 540*conj(z4)^4*conj(z2) + 900*conj(z4)^3*conj(z5)*conj(z2) +
  1080*z1*conj(z4)^4*z2*z3*z4^2*conj(z2) - 1800*z1*conj(z4)^3*conj(z5)*z2*z3*z4^2*conj(z2) - 540*z1*conj(z4)^4*z2*z4^3*conj(z2) +
  900*z1*conj(z4)^3*conj(z5)*z2*z4^3*conj(z2) - 540*z1*conj(z4)^4*z3*z4^3*conj(z2) + 900*z1*conj(z4)^3*conj(z5)*z3*z4^3*conj(z2) - 540*conj(z4)^4*z2*z3*z4^3*conj(z2) +
  900*conj(z4)^3*conj(z5)*z2*z3*z4^3*conj(z2) + 324*z1*conj(z4)^4*z4^4*conj(z2) - 540*z1*conj(z4)^3*conj(z5)*z4^4*conj(z2) + 324*conj(z4)^4*z2*z4^4*conj(z2) -
  540*conj(z4)^3*conj(z5)*z2*z4^4*conj(z2) + 324*conj(z4)^4*z3*z4^4*conj(z2) - 540*conj(z4)^3*conj(z5)*z3*z4^4*conj(z2) - 216*conj(z4)^4*z4^5*conj(z2) +
  360*conj(z4)^3*conj(z5)*z4^5*conj(z2) + 1080*z1*conj(z4)^4*z2*z4^2*z5*conj(z2) - 1800*z1*conj(z4)^3*conj(z5)*z2*z4^2*z5*conj(z2) +
  1080*z1*conj(z4)^4*z3*z4^2*z5*conj(z2) - 1800*z1*conj(z4)^3*conj(z5)*z3*z4^2*z5*conj(z2) + 1080*conj(z4)^4*z2*z3*z4^2*z5*conj(z2) -
  1800*conj(z4)^3*conj(z5)*z2*z3*z4^2*z5*conj(z2) - 540*z1*conj(z4)^4*z4^3*z5*conj(z2) + 900*z1*conj(z4)^3*conj(z5)*z4^3*z5*conj(z2) -
  540*conj(z4)^4*z2*z4^3*z5*conj(z2) + 900*conj(z4)^3*conj(z5)*z2*z4^3*z5*conj(z2) - 540*conj(z4)^4*z3*z4^3*z5*conj(z2) + 900*conj(z4)^3*conj(z5)*z3*z4^3*z5*conj(z2) +
  324*conj(z4)^4*z4^4*z5*conj(z2) - 540*conj(z4)^3*conj(z5)*z4^4*z5*conj(z2) + 900*conj(z4)^3*conj(z1)*conj(z2) - 1800*conj(z4)^2*conj(z5)*conj(z1)*conj(z2) -
  1800*z1*conj(z4)^3*z2*z3*z4^2*conj(z1)*conj(z2) + 3600*z1*conj(z4)^2*conj(z5)*z2*z3*z4^2*conj(z1)*conj(z2) + 900*z1*conj(z4)^3*z2*z4^3*conj(z1)*conj(z2) -
  1800*z1*conj(z4)^2*conj(z5)*z2*z4^3*conj(z1)*conj(z2) + 900*z1*conj(z4)^3*z3*z4^3*conj(z1)*conj(z2) - 1800*z1*conj(z4)^2*conj(z5)*z3*z4^3*conj(z1)*conj(z2) +
  900*conj(z4)^3*z2*z3*z4^3*conj(z1)*conj(z2) - 1800*conj(z4)^2*conj(z5)*z2*z3*z4^3*conj(z1)*conj(z2) - 540*z1*conj(z4)^3*z4^4*conj(z1)*conj(z2) +
  1080*z1*conj(z4)^2*conj(z5)*z4^4*conj(z1)*conj(z2) - 540*conj(z4)^3*z2*z4^4*conj(z1)*conj(z2) + 1080*conj(z4)^2*conj(z5)*z2*z4^4*conj(z1)*conj(z2) -
  540*conj(z4)^3*z3*z4^4*conj(z1)*conj(z2) + 1080*conj(z4)^2*conj(z5)*z3*z4^4*conj(z1)*conj(z2) + 360*conj(z4)^3*z4^5*conj(z1)*conj(z2) - 720*conj(z4)^2*conj(z5)*z4^5*conj(z1)*conj(z2) -
  1800*z1*conj(z4)^3*z2*z4^2*z5*conj(z1)*conj(z2) + 3600*z1*conj(z4)^2*conj(z5)*z2*z4^2*z5*conj(z1)*conj(z2) - 1800*z1*conj(z4)^3*z3*z4^2*z5*conj(z1)*conj(z2) +
  3600*z1*conj(z4)^2*conj(z5)*z3*z4^2*z5*conj(z1)*conj(z2) - 1800*conj(z4)^3*z2*z3*z4^2*z5*conj(z1)*conj(z2) + 3600*conj(z4)^2*conj(z5)*z2*z3*z4^2*z5*conj(z1)*conj(z2) +
  900*z1*conj(z4)^3*z4^3*z5*conj(z1)*conj(z2) - 1800*z1*conj(z4)^2*conj(z5)*z4^3*z5*conj(z1)*conj(z2) + 900*conj(z4)^3*z2*z4^3*z5*conj(z1)*conj(z2) -
  1800*conj(z4)^2*conj(z5)*z2*z4^3*z5*conj(z1)*conj(z2) + 900*conj(z4)^3*z3*z4^3*z5*conj(z1)*conj(z2) - 1800*conj(z4)^2*conj(z5)*z3*z4^3*z5*conj(z1)*conj(z2) -
  540*conj(z4)^3*z4^4*z5*conj(z1)*conj(z2) + 1080*conj(z4)^2*conj(z5)*z4^4*z5*conj(z1)*conj(z2) - 540*conj(z4)^4*conj(z3) + 900*conj(z4)^3*conj(z5)*conj(z3) +
  1080*z1*conj(z4)^4*z2*z3*z4^2*conj(z3) - 1800*z1*conj(z4)^3*conj(z5)*z2*z3*z4^2*conj(z3) - 540*z1*conj(z4)^4*z2*z4^3*conj(z3) +
  900*z1*conj(z4)^3*conj(z5)*z2*z4^3*conj(z3) - 540*z1*conj(z4)^4*z3*z4^3*conj(z3) + 900*z1*conj(z4)^3*conj(z5)*z3*z4^3*conj(z3) - 540*conj(z4)^4*z2*z3*z4^3*conj(z3) +
  900*conj(z4)^3*conj(z5)*z2*z3*z4^3*conj(z3) + 324*z1*conj(z4)^4*z4^4*conj(z3) - 540*z1*conj(z4)^3*conj(z5)*z4^4*conj(z3) + 324*conj(z4)^4*z2*z4^4*conj(z3) -
  540*conj(z4)^3*conj(z5)*z2*z4^4*conj(z3) + 324*conj(z4)^4*z3*z4^4*conj(z3) - 540*conj(z4)^3*conj(z5)*z3*z4^4*conj(z3) - 216*conj(z4)^4*z4^5*conj(z3) +
  360*conj(z4)^3*conj(z5)*z4^5*conj(z3) + 1080*z1*conj(z4)^4*z2*z4^2*z5*conj(z3) - 1800*z1*conj(z4)^3*conj(z5)*z2*z4^2*z5*conj(z3) +
  1080*z1*conj(z4)^4*z3*z4^2*z5*conj(z3) - 1800*z1*conj(z4)^3*conj(z5)*z3*z4^2*z5*conj(z3) + 1080*conj(z4)^4*z2*z3*z4^2*z5*conj(z3) -
  1800*conj(z4)^3*conj(z5)*z2*z3*z4^2*z5*conj(z3) - 540*z1*conj(z4)^4*z4^3*z5*conj(z3) + 900*z1*conj(z4)^3*conj(z5)*z4^3*z5*conj(z3) -
  540*conj(z4)^4*z2*z4^3*z5*conj(z3) + 900*conj(z4)^3*conj(z5)*z2*z4^3*z5*conj(z3) - 540*conj(z4)^4*z3*z4^3*z5*conj(z3) + 900*conj(z4)^3*conj(z5)*z3*z4^3*z5*conj(z3) +
  324*conj(z4)^4*z4^4*z5*conj(z3) - 540*conj(z4)^3*conj(z5)*z4^4*z5*conj(z3) + 900*conj(z4)^3*conj(z1)*conj(z3) - 1800*conj(z4)^2*conj(z5)*conj(z1)*conj(z3) -
  1800*z1*conj(z4)^3*z2*z3*z4^2*conj(z1)*conj(z3) + 3600*z1*conj(z4)^2*conj(z5)*z2*z3*z4^2*conj(z1)*conj(z3) + 900*z1*conj(z4)^3*z2*z4^3*conj(z1)*conj(z3) -
  1800*z1*conj(z4)^2*conj(z5)*z2*z4^3*conj(z1)*conj(z3) + 900*z1*conj(z4)^3*z3*z4^3*conj(z1)*conj(z3) - 1800*z1*conj(z4)^2*conj(z5)*z3*z4^3*conj(z1)*conj(z3) +
  900*conj(z4)^3*z2*z3*z4^3*conj(z1)*conj(z3) - 1800*conj(z4)^2*conj(z5)*z2*z3*z4^3*conj(z1)*conj(z3) - 540*z1*conj(z4)^3*z4^4*conj(z1)*conj(z3) +
  1080*z1*conj(z4)^2*conj(z5)*z4^4*conj(z1)*conj(z3) - 540*conj(z4)^3*z2*z4^4*conj(z1)*conj(z3) + 1080*conj(z4)^2*conj(z5)*z2*z4^4*conj(z1)*conj(z3) -
  540*conj(z4)^3*z3*z4^4*conj(z1)*conj(z3) + 1080*conj(z4)^2*conj(z5)*z3*z4^4*conj(z1)*conj(z3) + 360*conj(z4)^3*z4^5*conj(z1)*conj(z3) - 720*conj(z4)^2*conj(z5)*z4^5*conj(z1)*conj(z3) -
  1800*z1*conj(z4)^3*z2*z4^2*z5*conj(z1)*conj(z3) + 3600*z1*conj(z4)^2*conj(z5)*z2*z4^2*z5*conj(z1)*conj(z3) - 1800*z1*conj(z4)^3*z3*z4^2*z5*conj(z1)*conj(z3) +
  3600*z1*conj(z4)^2*conj(z5)*z3*z4^2*z5*conj(z1)*conj(z3) - 1800*conj(z4)^3*z2*z3*z4^2*z5*conj(z1)*conj(z3) + 3600*conj(z4)^2*conj(z5)*z2*z3*z4^2*z5*conj(z1)*conj(z3) +
  900*z1*conj(z4)^3*z4^3*z5*conj(z1)*conj(z3) - 1800*z1*conj(z4)^2*conj(z5)*z4^3*z5*conj(z1)*conj(z3) + 900*conj(z4)^3*z2*z4^3*z5*conj(z1)*conj(z3) -
  1800*conj(z4)^2*conj(z5)*z2*z4^3*z5*conj(z1)*conj(z3) + 900*conj(z4)^3*z3*z4^3*z5*conj(z1)*conj(z3) - 1800*conj(z4)^2*conj(z5)*z3*z4^3*z5*conj(z1)*conj(z3) -
  540*conj(z4)^3*z4^4*z5*conj(z1)*conj(z3) + 1080*conj(z4)^2*conj(z5)*z4^4*z5*conj(z1)*conj(z3) + 900*conj(z4)^3*conj(z2)*conj(z3) - 1800*conj(z4)^2*conj(z5)*conj(z2)*conj(z3) -
  1800*z1*conj(z4)^3*z2*z3*z4^2*conj(z2)*conj(z3) + 3600*z1*conj(z4)^2*conj(z5)*z2*z3*z4^2*conj(z2)*conj(z3) + 900*z1*conj(z4)^3*z2*z4^3*conj(z2)*conj(z3) -
  1800*z1*conj(z4)^2*conj(z5)*z2*z4^3*conj(z2)*conj(z3) + 900*z1*conj(z4)^3*z3*z4^3*conj(z2)*conj(z3) - 1800*z1*conj(z4)^2*conj(z5)*z3*z4^3*conj(z2)*conj(z3) +
  900*conj(z4)^3*z2*z3*z4^3*conj(z2)*conj(z3) - 1800*conj(z4)^2*conj(z5)*z2*z3*z4^3*conj(z2)*conj(z3) - 540*z1*conj(z4)^3*z4^4*conj(z2)*conj(z3) +
  1080*z1*conj(z4)^2*conj(z5)*z4^4*conj(z2)*conj(z3) - 540*conj(z4)^3*z2*z4^4*conj(z2)*conj(z3) + 1080*conj(z4)^2*conj(z5)*z2*z4^4*conj(z2)*conj(z3) -
  540*conj(z4)^3*z3*z4^4*conj(z2)*conj(z3) + 1080*conj(z4)^2*conj(z5)*z3*z4^4*conj(z2)*conj(z3) + 360*conj(z4)^3*z4^5*conj(z2)*conj(z3) - 720*conj(z4)^2*conj(z5)*z4^5*conj(z2)*conj(z3) -
  1800*z1*conj(z4)^3*z2*z4^2*z5*conj(z2)*conj(z3) + 3600*z1*conj(z4)^2*conj(z5)*z2*z4^2*z5*conj(z2)*conj(z3) - 1800*z1*conj(z4)^3*z3*z4^2*z5*conj(z2)*conj(z3) +
  3600*z1*conj(z4)^2*conj(z5)*z3*z4^2*z5*conj(z2)*conj(z3) - 1800*conj(z4)^3*z2*z3*z4^2*z5*conj(z2)*conj(z3) + 3600*conj(z4)^2*conj(z5)*z2*z3*z4^2*z5*conj(z2)*conj(z3) +
  900*z1*conj(z4)^3*z4^3*z5*conj(z2)*conj(z3) - 1800*z1*conj(z4)^2*conj(z5)*z4^3*z5*conj(z2)*conj(z3) + 900*conj(z4)^3*z2*z4^3*z5*conj(z2)*conj(z3) -
  1800*conj(z4)^2*conj(z5)*z2*z4^3*z5*conj(z2)*conj(z3) + 900*conj(z4)^3*z3*z4^3*z5*conj(z2)*conj(z3) - 1800*conj(z4)^2*conj(z5)*z3*z4^3*z5*conj(z2)*conj(z3) -
  540*conj(z4)^3*z4^4*z5*conj(z2)*conj(z3) + 1080*conj(z4)^2*conj(z5)*z4^4*z5*conj(z2)*conj(z3) - 1800*conj(z4)^2*conj(z1)*conj(z2)*conj(z3) +
  3600*z1*conj(z4)^2*z2*z3*z4^2*conj(z1)*conj(z2)*conj(z3) - 1800*z1*conj(z4)^2*z2*z4^3*conj(z1)*conj(z2)*conj(z3) - 1800*z1*conj(z4)^2*z3*z4^3*conj(z1)*conj(z2)*conj(z3) -
  1800*conj(z4)^2*z2*z3*z4^3*conj(z1)*conj(z2)*conj(z3) + 1080*z1*conj(z4)^2*z4^4*conj(z1)*conj(z2)*conj(z3) + 1080*conj(z4)^2*z2*z4^4*conj(z1)*conj(z2)*conj(z3) +
  1080*conj(z4)^2*z3*z4^4*conj(z1)*conj(z2)*conj(z3) - 720*conj(z4)^2*z4^5*conj(z1)*conj(z2)*conj(z3) + 3600*z1*conj(z4)^2*z2*z4^2*z5*conj(z1)*conj(z2)*conj(z3) +
  3600*z1*conj(z4)^2*z3*z4^2*z5*conj(z1)*conj(z2)*conj(z3) + 3600*conj(z4)^2*z2*z3*z4^2*z5*conj(z1)*conj(z2)*conj(z3) - 1800*z1*conj(z4)^2*z4^3*z5*conj(z1)*conj(z2)*conj(z3) -
  1800*conj(z4)^2*z2*z4^3*z5*conj(z1)*conj(z2)*conj(z3) - 1800*conj(z4)^2*z3*z4^3*z5*conj(z1)*conj(z2)*conj(z3) + 1080*conj(z4)^2*z4^4*z5*conj(z1)*conj(z2)*conj(z3) - 3600*z6*conj(z6)
g5 = 900 - 540*conj(z4)*conj(z5)^4 + 360*conj(z5)^5 - 1800*z1*z2*z3*z5^2 + 1080*z1*conj(z4)*conj(z5)^4*z2*z3*z5^2 - 720*z1*conj(z5)^5*z2*z3*z5^2 -
  1800*z1*z2*z4*z5^2 + 1080*z1*conj(z4)*conj(z5)^4*z2*z4*z5^2 - 720*z1*conj(z5)^5*z2*z4*z5^2 - 1800*z1*z3*z4*z5^2 +
  1080*z1*conj(z4)*conj(z5)^4*z3*z4*z5^2 - 720*z1*conj(z5)^5*z3*z4*z5^2 - 1800*z2*z3*z4*z5^2 + 1080*conj(z4)*conj(z5)^4*z2*z3*z4*z5^2 -
  720*conj(z5)^5*z2*z3*z4*z5^2 + 900*z1*z2*z5^3 - 540*z1*conj(z4)*conj(z5)^4*z2*z5^3 + 360*z1*conj(z5)^5*z2*z5^3 + 900*z1*z3*z5^3 -
  540*z1*conj(z4)*conj(z5)^4*z3*z5^3 + 360*z1*conj(z5)^5*z3*z5^3 + 900*z2*z3*z5^3 - 540*conj(z4)*conj(z5)^4*z2*z3*z5^3 +
  360*conj(z5)^5*z2*z3*z5^3 + 900*z1*z4*z5^3 - 540*z1*conj(z4)*conj(z5)^4*z4*z5^3 + 360*z1*conj(z5)^5*z4*z5^3 + 900*z2*z4*z5^3 -
  540*conj(z4)*conj(z5)^4*z2*z4*z5^3 + 360*conj(z5)^5*z2*z4*z5^3 + 900*z3*z4*z5^3 - 540*conj(z4)*conj(z5)^4*z3*z4*z5^3 +
  360*conj(z5)^5*z3*z4*z5^3 - 540*z1*z5^4 + 324*z1*conj(z4)*conj(z5)^4*z5^4 - 216*z1*conj(z5)^5*z5^4 - 540*z2*z5^4 +
  324*conj(z4)*conj(z5)^4*z2*z5^4 - 216*conj(z5)^5*z2*z5^4 - 540*z3*z5^4 + 324*conj(z4)*conj(z5)^4*z3*z5^4 - 216*conj(z5)^5*z3*z5^4 -
  540*z4*z5^4 + 324*conj(z4)*conj(z5)^4*z4*z5^4 - 216*conj(z5)^5*z4*z5^4 + 360*z5^5 - 216*conj(z4)*conj(z5)^4*z5^5 + 144*conj(z5)^5*z5^5 +
  900*conj(z4)*conj(z5)^3*conj(z1) - 540*conj(z5)^4*conj(z1) - 1800*z1*conj(z4)*conj(z5)^3*z2*z3*z5^2*conj(z1) + 1080*z1*conj(z5)^4*z2*z3*z5^2*conj(z1) -
  1800*z1*conj(z4)*conj(z5)^3*z2*z4*z5^2*conj(z1) + 1080*z1*conj(z5)^4*z2*z4*z5^2*conj(z1) - 1800*z1*conj(z4)*conj(z5)^3*z3*z4*z5^2*conj(z1) +
  1080*z1*conj(z5)^4*z3*z4*z5^2*conj(z1) - 1800*conj(z4)*conj(z5)^3*z2*z3*z4*z5^2*conj(z1) + 1080*conj(z5)^4*z2*z3*z4*z5^2*conj(z1) +
  900*z1*conj(z4)*conj(z5)^3*z2*z5^3*conj(z1) - 540*z1*conj(z5)^4*z2*z5^3*conj(z1) + 900*z1*conj(z4)*conj(z5)^3*z3*z5^3*conj(z1) - 540*z1*conj(z5)^4*z3*z5^3*conj(z1) +
  900*conj(z4)*conj(z5)^3*z2*z3*z5^3*conj(z1) - 540*conj(z5)^4*z2*z3*z5^3*conj(z1) + 900*z1*conj(z4)*conj(z5)^3*z4*z5^3*conj(z1) - 540*z1*conj(z5)^4*z4*z5^3*conj(z1) +
  900*conj(z4)*conj(z5)^3*z2*z4*z5^3*conj(z1) - 540*conj(z5)^4*z2*z4*z5^3*conj(z1) + 900*conj(z4)*conj(z5)^3*z3*z4*z5^3*conj(z1) - 540*conj(z5)^4*z3*z4*z5^3*conj(z1) -
  540*z1*conj(z4)*conj(z5)^3*z5^4*conj(z1) + 324*z1*conj(z5)^4*z5^4*conj(z1) - 540*conj(z4)*conj(z5)^3*z2*z5^4*conj(z1) + 324*conj(z5)^4*z2*z5^4*conj(z1) -
  540*conj(z4)*conj(z5)^3*z3*z5^4*conj(z1) + 324*conj(z5)^4*z3*z5^4*conj(z1) - 540*conj(z4)*conj(z5)^3*z4*z5^4*conj(z1) + 324*conj(z5)^4*z4*z5^4*conj(z1) +
  360*conj(z4)*conj(z5)^3*z5^5*conj(z1) - 216*conj(z5)^4*z5^5*conj(z1) + 900*conj(z4)*conj(z5)^3*conj(z2) - 540*conj(z5)^4*conj(z2) - 1800*z1*conj(z4)*conj(z5)^3*z2*z3*z5^2*conj(z2) +
  1080*z1*conj(z5)^4*z2*z3*z5^2*conj(z2) - 1800*z1*conj(z4)*conj(z5)^3*z2*z4*z5^2*conj(z2) + 1080*z1*conj(z5)^4*z2*z4*z5^2*conj(z2) -
  1800*z1*conj(z4)*conj(z5)^3*z3*z4*z5^2*conj(z2) + 1080*z1*conj(z5)^4*z3*z4*z5^2*conj(z2) - 1800*conj(z4)*conj(z5)^3*z2*z3*z4*z5^2*conj(z2) +
  1080*conj(z5)^4*z2*z3*z4*z5^2*conj(z2) + 900*z1*conj(z4)*conj(z5)^3*z2*z5^3*conj(z2) - 540*z1*conj(z5)^4*z2*z5^3*conj(z2) +
  900*z1*conj(z4)*conj(z5)^3*z3*z5^3*conj(z2) - 540*z1*conj(z5)^4*z3*z5^3*conj(z2) + 900*conj(z4)*conj(z5)^3*z2*z3*z5^3*conj(z2) - 540*conj(z5)^4*z2*z3*z5^3*conj(z2) +
  900*z1*conj(z4)*conj(z5)^3*z4*z5^3*conj(z2) - 540*z1*conj(z5)^4*z4*z5^3*conj(z2) + 900*conj(z4)*conj(z5)^3*z2*z4*z5^3*conj(z2) - 540*conj(z5)^4*z2*z4*z5^3*conj(z2) +
  900*conj(z4)*conj(z5)^3*z3*z4*z5^3*conj(z2) - 540*conj(z5)^4*z3*z4*z5^3*conj(z2) - 540*z1*conj(z4)*conj(z5)^3*z5^4*conj(z2) + 324*z1*conj(z5)^4*z5^4*conj(z2) -
  540*conj(z4)*conj(z5)^3*z2*z5^4*conj(z2) + 324*conj(z5)^4*z2*z5^4*conj(z2) - 540*conj(z4)*conj(z5)^3*z3*z5^4*conj(z2) + 324*conj(z5)^4*z3*z5^4*conj(z2) -
  540*conj(z4)*conj(z5)^3*z4*z5^4*conj(z2) + 324*conj(z5)^4*z4*z5^4*conj(z2) + 360*conj(z4)*conj(z5)^3*z5^5*conj(z2) - 216*conj(z5)^4*z5^5*conj(z2) -
  1800*conj(z4)*conj(z5)^2*conj(z1)*conj(z2) + 900*conj(z5)^3*conj(z1)*conj(z2) + 3600*z1*conj(z4)*conj(z5)^2*z2*z3*z5^2*conj(z1)*conj(z2) - 1800*z1*conj(z5)^3*z2*z3*z5^2*conj(z1)*conj(z2) +
  3600*z1*conj(z4)*conj(z5)^2*z2*z4*z5^2*conj(z1)*conj(z2) - 1800*z1*conj(z5)^3*z2*z4*z5^2*conj(z1)*conj(z2) + 3600*z1*conj(z4)*conj(z5)^2*z3*z4*z5^2*conj(z1)*conj(z2) -
  1800*z1*conj(z5)^3*z3*z4*z5^2*conj(z1)*conj(z2) + 3600*conj(z4)*conj(z5)^2*z2*z3*z4*z5^2*conj(z1)*conj(z2) - 1800*conj(z5)^3*z2*z3*z4*z5^2*conj(z1)*conj(z2) -
  1800*z1*conj(z4)*conj(z5)^2*z2*z5^3*conj(z1)*conj(z2) + 900*z1*conj(z5)^3*z2*z5^3*conj(z1)*conj(z2) - 1800*z1*conj(z4)*conj(z5)^2*z3*z5^3*conj(z1)*conj(z2) +
  900*z1*conj(z5)^3*z3*z5^3*conj(z1)*conj(z2) - 1800*conj(z4)*conj(z5)^2*z2*z3*z5^3*conj(z1)*conj(z2) + 900*conj(z5)^3*z2*z3*z5^3*conj(z1)*conj(z2) -
  1800*z1*conj(z4)*conj(z5)^2*z4*z5^3*conj(z1)*conj(z2) + 900*z1*conj(z5)^3*z4*z5^3*conj(z1)*conj(z2) - 1800*conj(z4)*conj(z5)^2*z2*z4*z5^3*conj(z1)*conj(z2) +
  900*conj(z5)^3*z2*z4*z5^3*conj(z1)*conj(z2) - 1800*conj(z4)*conj(z5)^2*z3*z4*z5^3*conj(z1)*conj(z2) + 900*conj(z5)^3*z3*z4*z5^3*conj(z1)*conj(z2) +
  1080*z1*conj(z4)*conj(z5)^2*z5^4*conj(z1)*conj(z2) - 540*z1*conj(z5)^3*z5^4*conj(z1)*conj(z2) + 1080*conj(z4)*conj(z5)^2*z2*z5^4*conj(z1)*conj(z2) -
  540*conj(z5)^3*z2*z5^4*conj(z1)*conj(z2) + 1080*conj(z4)*conj(z5)^2*z3*z5^4*conj(z1)*conj(z2) - 540*conj(z5)^3*z3*z5^4*conj(z1)*conj(z2) +
  1080*conj(z4)*conj(z5)^2*z4*z5^4*conj(z1)*conj(z2) - 540*conj(z5)^3*z4*z5^4*conj(z1)*conj(z2) - 720*conj(z4)*conj(z5)^2*z5^5*conj(z1)*conj(z2) + 360*conj(z5)^3*z5^5*conj(z1)*conj(z2) +
  900*conj(z4)*conj(z5)^3*conj(z3) - 540*conj(z5)^4*conj(z3) - 1800*z1*conj(z4)*conj(z5)^3*z2*z3*z5^2*conj(z3) + 1080*z1*conj(z5)^4*z2*z3*z5^2*conj(z3) -
  1800*z1*conj(z4)*conj(z5)^3*z2*z4*z5^2*conj(z3) + 1080*z1*conj(z5)^4*z2*z4*z5^2*conj(z3) - 1800*z1*conj(z4)*conj(z5)^3*z3*z4*z5^2*conj(z3) +
  1080*z1*conj(z5)^4*z3*z4*z5^2*conj(z3) - 1800*conj(z4)*conj(z5)^3*z2*z3*z4*z5^2*conj(z3) + 1080*conj(z5)^4*z2*z3*z4*z5^2*conj(z3) +
  900*z1*conj(z4)*conj(z5)^3*z2*z5^3*conj(z3) - 540*z1*conj(z5)^4*z2*z5^3*conj(z3) + 900*z1*conj(z4)*conj(z5)^3*z3*z5^3*conj(z3) - 540*z1*conj(z5)^4*z3*z5^3*conj(z3) +
  900*conj(z4)*conj(z5)^3*z2*z3*z5^3*conj(z3) - 540*conj(z5)^4*z2*z3*z5^3*conj(z3) + 900*z1*conj(z4)*conj(z5)^3*z4*z5^3*conj(z3) - 540*z1*conj(z5)^4*z4*z5^3*conj(z3) +
  900*conj(z4)*conj(z5)^3*z2*z4*z5^3*conj(z3) - 540*conj(z5)^4*z2*z4*z5^3*conj(z3) + 900*conj(z4)*conj(z5)^3*z3*z4*z5^3*conj(z3) - 540*conj(z5)^4*z3*z4*z5^3*conj(z3) -
  540*z1*conj(z4)*conj(z5)^3*z5^4*conj(z3) + 324*z1*conj(z5)^4*z5^4*conj(z3) - 540*conj(z4)*conj(z5)^3*z2*z5^4*conj(z3) + 324*conj(z5)^4*z2*z5^4*conj(z3) -
  540*conj(z4)*conj(z5)^3*z3*z5^4*conj(z3) + 324*conj(z5)^4*z3*z5^4*conj(z3) - 540*conj(z4)*conj(z5)^3*z4*z5^4*conj(z3) + 324*conj(z5)^4*z4*z5^4*conj(z3) +
  360*conj(z4)*conj(z5)^3*z5^5*conj(z3) - 216*conj(z5)^4*z5^5*conj(z3) - 1800*conj(z4)*conj(z5)^2*conj(z1)*conj(z3) + 900*conj(z5)^3*conj(z1)*conj(z3) +
  3600*z1*conj(z4)*conj(z5)^2*z2*z3*z5^2*conj(z1)*conj(z3) - 1800*z1*conj(z5)^3*z2*z3*z5^2*conj(z1)*conj(z3) + 3600*z1*conj(z4)*conj(z5)^2*z2*z4*z5^2*conj(z1)*conj(z3) -
  1800*z1*conj(z5)^3*z2*z4*z5^2*conj(z1)*conj(z3) + 3600*z1*conj(z4)*conj(z5)^2*z3*z4*z5^2*conj(z1)*conj(z3) - 1800*z1*conj(z5)^3*z3*z4*z5^2*conj(z1)*conj(z3) +
  3600*conj(z4)*conj(z5)^2*z2*z3*z4*z5^2*conj(z1)*conj(z3) - 1800*conj(z5)^3*z2*z3*z4*z5^2*conj(z1)*conj(z3) - 1800*z1*conj(z4)*conj(z5)^2*z2*z5^3*conj(z1)*conj(z3) +
  900*z1*conj(z5)^3*z2*z5^3*conj(z1)*conj(z3) - 1800*z1*conj(z4)*conj(z5)^2*z3*z5^3*conj(z1)*conj(z3) + 900*z1*conj(z5)^3*z3*z5^3*conj(z1)*conj(z3) -
  1800*conj(z4)*conj(z5)^2*z2*z3*z5^3*conj(z1)*conj(z3) + 900*conj(z5)^3*z2*z3*z5^3*conj(z1)*conj(z3) - 1800*z1*conj(z4)*conj(z5)^2*z4*z5^3*conj(z1)*conj(z3) +
  900*z1*conj(z5)^3*z4*z5^3*conj(z1)*conj(z3) - 1800*conj(z4)*conj(z5)^2*z2*z4*z5^3*conj(z1)*conj(z3) + 900*conj(z5)^3*z2*z4*z5^3*conj(z1)*conj(z3) -
  1800*conj(z4)*conj(z5)^2*z3*z4*z5^3*conj(z1)*conj(z3) + 900*conj(z5)^3*z3*z4*z5^3*conj(z1)*conj(z3) + 1080*z1*conj(z4)*conj(z5)^2*z5^4*conj(z1)*conj(z3) -
  540*z1*conj(z5)^3*z5^4*conj(z1)*conj(z3) + 1080*conj(z4)*conj(z5)^2*z2*z5^4*conj(z1)*conj(z3) - 540*conj(z5)^3*z2*z5^4*conj(z1)*conj(z3) +
  1080*conj(z4)*conj(z5)^2*z3*z5^4*conj(z1)*conj(z3) - 540*conj(z5)^3*z3*z5^4*conj(z1)*conj(z3) + 1080*conj(z4)*conj(z5)^2*z4*z5^4*conj(z1)*conj(z3) -
  540*conj(z5)^3*z4*z5^4*conj(z1)*conj(z3) - 720*conj(z4)*conj(z5)^2*z5^5*conj(z1)*conj(z3) + 360*conj(z5)^3*z5^5*conj(z1)*conj(z3) - 1800*conj(z4)*conj(z5)^2*conj(z2)*conj(z3) +
  900*conj(z5)^3*conj(z2)*conj(z3) + 3600*z1*conj(z4)*conj(z5)^2*z2*z3*z5^2*conj(z2)*conj(z3) - 1800*z1*conj(z5)^3*z2*z3*z5^2*conj(z2)*conj(z3) +
  3600*z1*conj(z4)*conj(z5)^2*z2*z4*z5^2*conj(z2)*conj(z3) - 1800*z1*conj(z5)^3*z2*z4*z5^2*conj(z2)*conj(z3) + 3600*z1*conj(z4)*conj(z5)^2*z3*z4*z5^2*conj(z2)*conj(z3) -
  1800*z1*conj(z5)^3*z3*z4*z5^2*conj(z2)*conj(z3) + 3600*conj(z4)*conj(z5)^2*z2*z3*z4*z5^2*conj(z2)*conj(z3) - 1800*conj(z5)^3*z2*z3*z4*z5^2*conj(z2)*conj(z3) -
  1800*z1*conj(z4)*conj(z5)^2*z2*z5^3*conj(z2)*conj(z3) + 900*z1*conj(z5)^3*z2*z5^3*conj(z2)*conj(z3) - 1800*z1*conj(z4)*conj(z5)^2*z3*z5^3*conj(z2)*conj(z3) +
  900*z1*conj(z5)^3*z3*z5^3*conj(z2)*conj(z3) - 1800*conj(z4)*conj(z5)^2*z2*z3*z5^3*conj(z2)*conj(z3) + 900*conj(z5)^3*z2*z3*z5^3*conj(z2)*conj(z3) -
  1800*z1*conj(z4)*conj(z5)^2*z4*z5^3*conj(z2)*conj(z3) + 900*z1*conj(z5)^3*z4*z5^3*conj(z2)*conj(z3) - 1800*conj(z4)*conj(z5)^2*z2*z4*z5^3*conj(z2)*conj(z3) +
  900*conj(z5)^3*z2*z4*z5^3*conj(z2)*conj(z3) - 1800*conj(z4)*conj(z5)^2*z3*z4*z5^3*conj(z2)*conj(z3) + 900*conj(z5)^3*z3*z4*z5^3*conj(z2)*conj(z3) +
  1080*z1*conj(z4)*conj(z5)^2*z5^4*conj(z2)*conj(z3) - 540*z1*conj(z5)^3*z5^4*conj(z2)*conj(z3) + 1080*conj(z4)*conj(z5)^2*z2*z5^4*conj(z2)*conj(z3) -
  540*conj(z5)^3*z2*z5^4*conj(z2)*conj(z3) + 1080*conj(z4)*conj(z5)^2*z3*z5^4*conj(z2)*conj(z3) - 540*conj(z5)^3*z3*z5^4*conj(z2)*conj(z3) +
  1080*conj(z4)*conj(z5)^2*z4*z5^4*conj(z2)*conj(z3) - 540*conj(z5)^3*z4*z5^4*conj(z2)*conj(z3) - 720*conj(z4)*conj(z5)^2*z5^5*conj(z2)*conj(z3) + 360*conj(z5)^3*z5^5*conj(z2)*conj(z3) -
  1800*conj(z5)^2*conj(z1)*conj(z2)*conj(z3) + 3600*z1*conj(z5)^2*z2*z3*z5^2*conj(z1)*conj(z2)*conj(z3) + 3600*z1*conj(z5)^2*z2*z4*z5^2*conj(z1)*conj(z2)*conj(z3) +
  3600*z1*conj(z5)^2*z3*z4*z5^2*conj(z1)*conj(z2)*conj(z3) + 3600*conj(z5)^2*z2*z3*z4*z5^2*conj(z1)*conj(z2)*conj(z3) - 1800*z1*conj(z5)^2*z2*z5^3*conj(z1)*conj(z2)*conj(z3) -
  1800*z1*conj(z5)^2*z3*z5^3*conj(z1)*conj(z2)*conj(z3) - 1800*conj(z5)^2*z2*z3*z5^3*conj(z1)*conj(z2)*conj(z3) - 1800*z1*conj(z5)^2*z4*z5^3*conj(z1)*conj(z2)*conj(z3) -
  1800*conj(z5)^2*z2*z4*z5^3*conj(z1)*conj(z2)*conj(z3) - 1800*conj(z5)^2*z3*z4*z5^3*conj(z1)*conj(z2)*conj(z3) + 1080*z1*conj(z5)^2*z5^4*conj(z1)*conj(z2)*conj(z3) +
  1080*conj(z5)^2*z2*z5^4*conj(z1)*conj(z2)*conj(z3) + 1080*conj(z5)^2*z3*z5^4*conj(z1)*conj(z2)*conj(z3) + 1080*conj(z5)^2*z4*z5^4*conj(z1)*conj(z2)*conj(z3) -
  720*conj(z5)^2*z5^5*conj(z1)*conj(z2)*conj(z3) - 3600*z6*conj(z6)
h1 = z1*conj(z1) + z2*conj(z2) + z3*conj(z3) + z4*conj(z4) + z5*conj(z5) - n*(1/(n+1))^(2/n)
h2 = 3*z1*z2*z3*z4*z5 + 3*conj(z1)*conj(z2)*conj(z3)*conj(z4)*conj(z5) + 1
h3 = 36*z1*z2*z3*z4*z5*conj(z1)*conj(z2)*conj(z3)*conj(z4)*conj(z5) - 1
cpop = [f; g1; g2; g3; g4; g5; h1; h2; h3]

@time begin
opt,sol,data = complex_tssos_first(cpop, z, 5, numeq=3, TS="block", QUIET=true, solve=false, normality=1)
opt,sol,data = complex_tssos_higher!(data, TS="block", QUIET=true, solve=false)
opt,sol,data = complex_tssos_higher!(data, TS="block", QUIET=true, solve=false)
opt,sol,data = complex_tssos_higher!(data, TS="block", QUIET=true)
end
println(sqrt(-opt))

# Mordell Inequality Conjecture
n = 3
@complex_polyvar z1 z2
z = tuple(z1,z2)
f = (z1-z2)*(conj(z1)-conj(z2))*(2z1+z2)*(2z3+conj(z2))*(2z2+z1)*(2z4+conj(z1))
h = 2*z1*conj(z1) + 2*z2*conj(z2) + z1*conj(z2) + z2*conj(z1) - 3

n = 4
@complex_polyvar z1 z2 z3
z = tuple(z1,z2,z3)
f = (z1-z2)*(conj(z1)-conj(z2))*(z1-z3)*(conj(z1)-conj(z3))*(2z1+z2+z3)*(2z4+conj(z2)+conj(z3))*(z2-z3)*(conj(z2)-conj(z3))*(2z2+z1+z3)*(2z5+conj(z1)+conj(z3))*(2z3+z1+z2)*(2z6+conj(z1)+conj(z2))
h = z1*conj(z1) + z2*conj(z2) + z3*conj(z3) + (z1+z2+z3)*(conj(z1)+conj(z2)+conj(z3)) - 4

n = 5
@complex_polyvar z1 z2 z3 z4
z = tuple(z1,z2,z3,z4)
f = (z1-z2)*(conj(z1)-conj(z2))*(z1-z3)*(conj(z1)-conj(z3))*(z1-z4)*(conj(z1)-conj(z4))*(z2-z3)*(conj(z2)-conj(z3))*(z2-z4)*(conj(z2)-conj(z4))*(z3-z4)*(conj(z3)-conj(z4))*
(2z1+z2+z3+z4)*(2z5+conj(z2)+conj(z3)+conj(z4))*(2z2+z1+z3+z4)*(2z6+conj(z1)+conj(z3)+conj(z4))*(2z3+z1+z2+z4)*(2z7+conj(z1)+conj(z2)+conj(z4))*(2z4+z1+z2+z3)*(2z8+conj(z1)+conj(z2)+conj(z3))
h = z1*conj(z1) + z2*conj(z2) + z3*conj(z3) + z4*conj(z4) + (z1+z2+z3+z4)*(conj(z1)+conj(z2)+conj(z3)+conj(z4)) - 5

cpop = [-f; h]
@time begin
opt,sol,data = complex_tssos_first(cpop, z, 6, numeq=1, TS="block", normality=5, solve=false, QUIET=true)
opt,sol,data = complex_tssos_higher!(data, TS="block", QUIET=true, solve=false)
opt,sol,data = complex_tssos_higher!(data, TS="block", QUIET=true)
end
