using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using TSSOS

n = 3
@polyvar x[1:n]
f = [(x[1]^2+x[2]^2-1/4)*x[1], (x[2]^2+x[3]^2-1/4)*x[2], (x[2]^2+x[3]^2-1/4)*x[3]]
g = [1-x[1]^2, 1-x[2]^2, 1-x[3]^2]
d = 3

model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
v, vc, vb = add_poly!(model, x, 2d-2)
w, wc, wb = add_poly!(model, x, 2d)
Lv = v - sum(f .* differentiate(v, x))
model,info1 = add_psatz!(model, Lv, x, g, [], d, QUIET=true, CS=true, cliques=[[1;2;3]], TS="block", SO=1, Groebnerbasis=false, constrs="con1")
model,info2 = add_psatz!(model, w, x, g, [], d, QUIET=true, CS=true, TS="block", SO=1, Groebnerbasis=false)
model,info3 = add_psatz!(model, w-v-1, x, g, [], d, QUIET=true, CS=true, TS="block", SO=1, Groebnerbasis=false)
supp = get_nbasis(n, 2d, var=Vector(n:-1:1))
moment = get_moment(n, supp, -ones(n), ones(n))
@objective(model, Min, sum(moment.*wc))
optimize!(model)
status = termination_status(model)
if status != MOI.OPTIMAL
    println("termination status: $status")
    status = primal_status(model)
    println("solution status: $status")
end
objv = objective_value(model)

# retrieve Gram matrices
GramMat = Vector{Vector{Vector{Union{Float64,Matrix{Float64}}}}}(undef, info1.cql)
for i = 1:info1.cql
    GramMat[i] = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, 1+length(info1.I[i])+length(info1.J[i]))
    for j = 1:1+length(info1.I[i])+length(info1.J[i])
        GramMat[i][j] = [value.(info1.gram[i][j][l]) for l = 1:info1.cl[i][j]]
    end
end

# retrieve moment matrices
moment = [-dual(constraint_by_name(model, "con1[$i]")) for i=1:size(info1.tsupp, 2)]
MomMat = get_moment_matrix(moment, info1.tsupp, info1.cql, info1.basis)
