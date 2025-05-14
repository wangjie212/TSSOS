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
info1 = add_psatz!(model, Lv, x, g, [], d, TS=false, SO=1, constrs="con1")
info2 = add_psatz!(model, w, x, g, [], d, TS=false, SO=1)
info3 = add_psatz!(model, w-v-1, x, g, [], d, TS=false, SO=1)
moment = get_moment(wb, -ones(n), ones(n))
@objective(model, Min, sum(moment.*wc))
optimize!(model)
objv = objective_value(model)
@show objv
# objv = 3.437648

# retrieve Gram matrices
GramMat = Vector{Vector{Vector{Union{Float64,Matrix{Float64}}}}}(undef, info1.cql)
for i = 1:info1.cql
    GramMat[i] = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, 1+length(info1.I[i]))
    for j = 1:1+length(info1.I[i])
        GramMat[i][j] = [value.(info1.GramMat[i][j][l]) for l = 1:length(info1.GramMat[i][j])]
    end
end

# retrieve moment matrices
MomMat = get_moment_matrix(-dual(constraint_by_name(model, "con1")), info1)
