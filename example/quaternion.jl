using TSSOS
using TSSOS.QuaternionSOS
using DynamicPolynomials
using COSMO

# z[1:n] are quaternion variables; z[n+1:2n] are their conjugates.
@ncpolyvar z[1:2]
pop = [-z[1]*z[2], 1-z[1]*z[2]]
bound = qtssos(pop, z, 1, 1; TS=false, optimizer=COSMO.Optimizer, QUIET=true)
println("Quaternion bound: ", bound) # -1

# Compare with the real formulation using four coordinates per quaternion.
real_pop, x = quaternion_to_real(pop, z)
real_bound, _, _ = tssos(real_pop, x, 1; TS=false, optimizer=COSMO.Optimizer, QUIET=true)
println("Real bound: ", real_bound)
