using DynamicPolynomials
using TSSOS
using JuMP
using MosekTools
using Random

# Example 4.8
@polyvar x y z
p = [x^2+y^2-y*z, y^2+x^2*z, z^2-x+y]
q = [1+2x^2+y^2+z^2, 1+x^2+2y^2+z^2, 1+x^2+y^2+2z^2]
g = [1-x^2-y^2-z^2]

d = 4
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
h1 = add_poly!(model, [x;y;z], 2d-2, signsymmetry=get_signsymmetry([p[2]; q[2]; g], [x;y;z]))[1]
h2 = add_poly!(model, [x;y;z], 2d-2, signsymmetry=get_signsymmetry([p[3]; q[3]; g], [x;y;z]))[1]
c = @variable(model)
add_psatz!(model, p[1]+(h1+h2-c)*q[1], [x;y;z], g, [], d, QUIET=true, CS=false, TS=false, SO=1, GroebnerBasis=false)
add_psatz!(model, p[2]-h1*q[2], [x;y;z], g, [], d, QUIET=true, CS=false, TS="block", SO=1, GroebnerBasis=false)
add_psatz!(model, p[3]-h2*q[3], [x;y;z], g, [], d, QUIET=true, CS=false, TS="block", SO=1, GroebnerBasis=false)
@objective(model, Max, c)
optimize!(model)
objv = objective_value(model)
@show objv

d = 4
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
h1 = add_poly!(model, [x;y;z], 2d-2, signsymmetry=get_signsymmetry([p[1]; q[1]; g], [x;y;z]))[1]
h2 = add_poly!(model, [x;y;z], 2d-2, signsymmetry=get_signsymmetry([p[3]; q[3]; g], [x;y;z]))[1]
c = @variable(model)
add_psatz!(model, p[1]-h1*q[1], [x;y;z], g, [], d, QUIET=true, CS=false, TS="block", SO=1, GroebnerBasis=false)
add_psatz!(model, p[2]+(h1+h2-c)*q[2], [x;y;z], g, [], d, QUIET=true, CS=false, TS=false, SO=1, GroebnerBasis=false)
add_psatz!(model, p[3]-h2*q[3], [x;y;z], g, [], d, QUIET=true, CS=false, TS="block", SO=1, GroebnerBasis=false)
@objective(model, Max, c)
optimize!(model)
objv = objective_value(model)
@show objv

d = 4
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
h1 = add_poly!(model, [x;y;z], 2d-2, signsymmetry=get_signsymmetry([p[1]; q[1]; g], [x;y;z]))[1]
h2 = add_poly!(model, [x;y;z], 2d-2, signsymmetry=get_signsymmetry([p[2]; q[2]; g], [x;y;z]))[1]
c = @variable(model)
add_psatz!(model, p[1]-h1*q[1], [x;y;z], g, [], d, QUIET=true, CS=false, TS="block", SO=1, GroebnerBasis=false)
add_psatz!(model, p[2]-h2*q[2], [x;y;z], g, [], d, QUIET=true, CS=false, TS="block", SO=1, GroebnerBasis=false)
add_psatz!(model, p[3]+(h1+h2-c)*q[3], [x;y;z], g, [], d, QUIET=true, CS=false, TS=false, SO=1, GroebnerBasis=false)
@objective(model, Max, c)
optimize!(model)
objv = objective_value(model)
@show objv

opt = SumOfRatios(p, q, g, [], [x;y;z], 3, QUIET=true, SignSymmetry=false)

@polyvar w[1:3]
pop = [sum(w), 1-x^2-y^2-z^2, p[1]-w[1]*q[1], p[2]-w[2]*q[2], p[3]-w[3]*q[3]]
opt,sol,data = cs_tssos_first(pop, [x;y;z;w], 3, numeq=3, TS=false, solution=false, QUIET=true)

# Example 5.1
@polyvar x[1:3]
@polyvar a
d = 2
M = 6
p = [a^4*(x[1]^(6d) + x[2]^(6d) + x[3]^(6d)) + (x[1]^(4d)*x[2]^(2d) + x[2]^(4d)*x[3]^(2d) + x[3]^(4d)*x[1]^(2d)) +
a^8*(x[1]^(2d)*x[2]^(4d) + x[2]^(2d)*x[3]^(4d) + x[3]^(2d)*x[1]^(4d)) for a in Vector(1:M-1)/M]
q = [2*a^6*(x[1]^(4d)*x[2]^(2d) + x[2]^(4d)*x[3]^(2d) + x[3]^(4d)*x[1]^(2d)) + 2*a^2*(x[1]^(2d)*x[2]^(4d) + 
x[2]^(2d)*x[3]^(4d) + x[3]^(2d)*x[1]^(4d)) + 3*(1-2*a^2+a^4-2*a^6+a^8)*x[1]^(2d)*x[2]^(2d)*x[3]^(2d) for a in Vector(1:M-1)/M]
h = [3.0 - sum(x.^2)]

@time opt = SumOfRatios(p, q, [], h, x, 3d, QUIET=false, SignSymmetry=true)
@time opt = SumOfRatios(p, q, [], h, x, 3d, QUIET=false, SignSymmetry=false)
@time opt = SumOfRatios(p, q, [], h, x, 3d, QUIET=false, dualize=true, SignSymmetry=false)

@polyvar w[1:M-1]
pop = [sum(w); h; [p[i]-w[i]*q[i] for i=1:M-1]]
@time opt,sol,data = cs_tssos_first(pop, [x;w], 3d+1, numeq=M, TS="block", solution=false, QUIET=true)

# Example 5.2
c = [0.806; 0.517; 0.100; 0.908; 0.965; 0.669; 0.524; 0.902; 0.531; 0.876; 0.462; 0.491; 0.463; 0.714; 0.352; 
0.869; 0.813; 0.811; 0.828; 0.964; 0.789; 0.360; 0.369; 0.992; 0.332; 0.817; 0.632; 0.883; 0.608; 0.326]
A = [9.681 0.667 4.783 9.095 3.517 9.325 6.544 0.211 5.122 2.020;
9.400 2.041 3.788 7.931 2.882 2.672 3.568 1.284 7.033 7.374;
8.025 9.152 5.114 7.621 4.564 4.711 2.996 6.126 0.734 4.982;
2.196 0.415 5.649 6.979 9.510 9.166 6.304 6.054 9.377 1.426;
8.074 8.777 3.467 1.863 6.708 6.349 4.534 0.276 7.633 1.567;
7.650 5.658 0.720 2.764 3.278 5.283 7.474 6.274 1.409 8.208;
1.256 3.605 8.623 6.905 4.584 8.133 6.071 6.888 4.187 5.448;
8.314 2.261 4.224 1.781 4.124 0.932 8.129 8.658 1.208 5.762;
0.226 8.858 1.420 0.945 1.622 4.698 6.228 9.096 0.972 7.637;
7.305 2.228 1.242 5.928 9.133 1.826 4.060 5.204 8.713 8.247;
0.652 7.027 0.508 4.876 8.807 4.632 5.808 6.937 3.291 7.016;
2.699 3.516 5.874 4.119 4.461 7.496 8.817 0.690 6.593 9.789;
8.327 3.897 2.017 9.570 9.825 1.150 1.395 3.885 6.354 0.109;
2.132 7.006 7.136 2.641 1.882 5.943 7.273 7.691 2.880 0.564;
4.707 5.579 4.080 0.581 9.698 8.542 8.077 8.515 9.231 4.670;
8.304 7.559 8.567 0.322 7.128 8.392 1.472 8.524 2.277 7.826;
8.632 4.409 4.832 5.768 7.050 6.715 1.711 4.323 4.405 4.591;
4.887 9.112 0.170 8.967 9.693 9.867 7.508 7.770 8.382 6.740;
2.440 6.686 4.299 1.007 7.008 1.427 9.398 8.480 9.950 1.675;
6.306 8.583 6.084 1.138 4.350 3.134 7.853 6.061 7.457 2.258;
0.652 2.343 1.370 0.821 1.310 1.063 0.689 8.819 8.833 9.070;
5.558 1.272 5.756 9.857 2.279 2.764 1.284 1.677 1.244 1.234;
3.352 7.549 9.817 9.437 8.687 4.167 2.570 6.540 0.228 0.027;
8.798 0.880 2.370 0.168 1.701 3.680 1.231 2.390 2.499 0.064;
1.460 8.057 1.336 7.217 7.914 3.615 9.981 9.198 5.292 1.224;
0.432 8.645 8.774 0.249 8.081 7.461 4.416 0.652 4.002 4.644;
0.679 2.800 5.523 3.049 2.968 7.225 6.730 4.199 9.614 9.229;
4.263 1.074 7.286 5.599 8.291 5.200 9.214 8.272 4.398 4.506;
9.496 4.830 3.150 8.270 5.079 1.231 5.731 9.494 1.883 9.732;
4.138 2.562 2.532 9.661 5.611 5.500 6.886 2.341 9.699 6.500]

N = 30
n = 5
@polyvar x[1:n]
p = -ones(N)
q = [sum((x[j]^2-A[i,j]/10)^2 for j=1:n) + c[i]/100 for i = 1:N]
g = [0.6 - sum((x[i]^2-0.5)^2 for i=1:n)]
k = 3
@time opt = SumOfRatios(p, q, g, [], x, k, QUIET=true, dualize=true, SignSymmetry=false)
@time opt = SumOfRatios(p, q, g, [], x, k, QUIET=true, SignSymmetry=false)
@time opt = SumOfRatios(p, q, g, [], x, k, QUIET=true, SignSymmetry=true)

@polyvar w[1:N]
pop = [sum(w); g; [p[i]-w[i]*q[i] for i=1:N]]
@time opt,sol,data = cs_tssos_first(pop, [x;w], 4, numeq=N, TS="block", solution=false, QUIET=true)

# Example 5.3
using MultivariatePolynomials
Random.seed!(1)
N = 10
n = 6
d = 4
pp = 0.05
s = 3
@polyvar x[1:n]
p = -ones(N)
mons = monomials(x, 0:d)
ind = shuffle(Vector(1:length(mons)))[1:Int(floor(length(mons)*pp))]
q = Vector{Poly{Float64}}(undef, N)
for i = 1:N
	loc = shuffle(Vector(1:length(ind)))[1:s]
	q[i] = 1 + (rand(s)'*mons[ind[loc]])^2
end
g = [1 - sum(x.^2)]
@time opt = SumOfRatios(p, q, g, [], x, d, QUIET=true, dualize=true, SignSymmetry=false)
@time opt = SumOfRatios(p, q, g, [], x, d, QUIET=true, SignSymmetry=false)
@time opt = SumOfRatios(p, q, g, [], x, d, QUIET=true, SignSymmetry=true)

@polyvar w[1:N]
pop = [sum(w); g; [p[i]-w[i]*q[i] for i=1:N]]
@time opt,sol,data = cs_tssos_first(pop, [x;w], d+1, numeq=N, TS="block", solution=false, QUIET=true)


# Example 5.4
N = 5
n = 2*N + 2
@polyvar x[1:n]
p = [sum(x[2*i-1:2*i+1].^2)*prod(x[2*i-1:2*i+1].^2) + x[2*i+2]^8 for i=1:N]
q = [prod(x[2*i-1:2*i+2].^2) for i=1:N]
h = [4 - sum(x[2*i-1:2*i+2].^2) for i=1:N]
d = 5
@time opt = SparseSumOfRatios(p, q, [], h, x, d, QUIET=true, dualize=true, SignSymmetry=false)
@time opt = SparseSumOfRatios(p, q, [], h, x, d, QUIET=true, SignSymmetry=false)
@time opt = SparseSumOfRatios(p, q, [], h, x, d, QUIET=true, SignSymmetry=true)

@polyvar w[1:N]
pop = [sum(w); h; [p[i]-w[i]*q[i] for i=1:N]]
@time opt,sol,data = cs_tssos_first(pop, [x;w], d+1, numeq=2N, TS="block", solution=false, QUIET=true)
@time opt,sol,data = cs_tssos_higher!(data, TS="block", solution=false, QUIET=true)


# Example 5.5
N = 5
d = 2
n = 2*N + 1	
@polyvar x[1:n]
p = [x[2*i-1]^(6d) + x[2i]^(6d) + x[2i+1]^(6d) + 3*x[2i-1]^(2d)*x[2i]^(2d)*x[2i+1]^(2d) for i = 1:N]
q = [x[2i-1]^(4d)*x[2i]^(2d) + x[2i-1]^(2d)*x[2i]^(4d) + x[2i-1]^(4d)*x[2i+1]^(2d) + x[2i-1]^(2d)*x[2i+1]^(4d) + x[2i]^(4d)*x[2i+1]^(2d) + x[2i]^(2d)*x[2i+1]^(4d) for i = 1:N]
h = [3 - sum(x[2*i-1:2*i+1].^2) for i = 1:N]

@time opt = SparseSumOfRatios(p, q, [], h, x, 3d, QUIET=true, SignSymmetry=true)
@time opt = SparseSumOfRatios(p, q, [], h, x, 3d, QUIET=true, SignSymmetry=false)
@time opt = SparseSumOfRatios(p, q, [], h, x, 3d, QUIET=true, dualize=true, SignSymmetry=false)

@polyvar w[1:N]
pop = [sum(w); h; [p[i]-w[i]*q[i] for i=1:N]]
@time opt,sol,data = cs_tssos_first(pop, [x;w], 3d+1, numeq=2N, TS="block", solution=false, QUIET=true)
@time opt,sol,data = cs_tssos_higher!(data, TS="block", solution=false, QUIET=true)


# Example 5.6
n = 51
@polyvar x[1:n]
p = ones(n-1)
q = [100(x[i+1]^2-x[i]^2)^2 + (x[i]^2-1)^2 + 1 for i = 1:n-1]
g = [16-x[i]^2 for i=1:n]
@time opt = SparseSumOfRatios(-p, q, g, [], x, 2, QUIET=true, SignSymmetry=true)
@time opt = SparseSumOfRatios(-p, q, g, [], x, 2, QUIET=true, SignSymmetry=false)
@time opt = SparseSumOfRatios(-p, q, g, [], x, 2, QUIET=true, dualize=true, SignSymmetry=false)
@polyvar w[1:n-1]
pop = [-sum(w); g; [p[i]-w[i]*q[i] for i=1:n-1]]
@time opt,sol,data = cs_tssos_first(pop, [x;w], 4, numeq=n-1, TS="block", solution=false, QUIET=true)

# Example 5.7
N = 20
s = 6
n = N + s
@polyvar x[1:n]
p = [sum(x[i+j-1]*x[i+j] for j=1:s) for i=1:N]
q = [1 + sum(j*x[i+j-1]^2 for j=1:s+1) for i = 1:N]
g = [1 - x[i]^2 for i=1:n]
@time opt = SparseSumOfRatios(p, q, g, [], x, 3, QUIET=true, SignSymmetry=true)
@time opt = SparseSumOfRatios(p, q, g, [], x, 3, QUIET=true, SignSymmetry=false)
@time opt = SparseSumOfRatios(p, q, g, [], x, 3, QUIET=true, dualize=true, SignSymmetry=false)

@polyvar w[1:N]
pop = [sum(w); g; [p[i]-w[i]*q[i] for i=1:N]]
@time begin
opt,sol,data = cs_tssos_first(pop, [x;w], 3, numeq=N, TS="block", solve=false, solution=false, QUIET=true)
opt,sol,data = cs_tssos_higher!(data, TS="block", solution=false, QUIET=true)
end

# Example 6.1
Random.seed!(0)
n = 4
N = 10
@polyvar x[1:2*n]	
p = Poly{Float64}[]
q = Poly{Float64}[]
for i in 1:N
   	A = rand(n, n)
    B = rand(n, n)
    push!(p, x[1:n]'*(A+A')*x[1:n] - 2*x[1:n]'*(B-B')*x[n+1:2*n] + x[n+1:2*n]'*(A+A')*x[n+1:2*n])
	C = rand(n, n)
	D = rand(n, n)
	push!(q, x[1:n]'*(C'*C + D'*D)*x[1:n] - 2*x[1:n]'*(C'*D - D'*C)*x[n+1:2*n] + x[n+1:2*n]'*(C'*C + D'*D)*x[n+1:2*n])
end
h = [1.0 - sum(x.^2)]

d = 2
@time opt = SumOfRatios(-p, q, [], h, x, d, QUIET=true, SignSymmetry=true)
@time opt = SumOfRatios(-p, q, [], h, x, d, QUIET=true, SignSymmetry=false)
@time opt = SumOfRatios(-p, q, [], h, x, d, QUIET=true, dualize=true, SignSymmetry=false)
