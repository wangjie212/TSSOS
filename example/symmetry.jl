using PermutationGroups
using DynamicPolynomials
using TSSOS
using JuMP
using MosekTools

@polyvar x[1:4]
f = sum(x) + sum(x.^2)
G = PermGroup([perm"(1,2,3,4)"]) # define the symmetry group
opt,data = tssos_symmetry([f], x, 1, G)
# optimum = -1

@polyvar x[1:4]
f = prod(x) + sum(x.^4)
G = PermGroup([perm"(1,2,3,4)"]) # define the symmetry group
opt,data = tssos_symmetry_first([f], x, 2, G)
opt,data = tssos_symmetry_higher!(data)
# optimum = 0

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry([f], x, 2, G)
# optimum = -1.4174111

@polyvar x[1:3]
f = sum(x) + sum(x.^4) - 4*x[1]*x[2]*x[3]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry([f], x, 2, G)
# optimum = -2.1129138

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry(pop, x, 2, G)
# optimum = -1.3987174

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
opt,sol,data = tssos_first(pop, x, 2, TS="block")
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry_first(pop, x, 2, G, SymmetricConstraint=true)
# optimum = -1.3987174

d = 2
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
lambda = @variable(model)
info = add_psatz_symmetry!(model, f - lambda, x, [1 - sum(x.^2)], [], d, G, SymmetricConstraint=true, TS="block", SO=1)
@objective(model, Max, lambda)
optimize!(model)
objv = objective_value(model)
@show objv
# optimum = -1.3987174

@polyvar x[1:3]
f = sum(x) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry(pop, x, 2, G, numeq=1)
# optimum = -1.3987174

d = 2
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
lambda = @variable(model)
info = add_psatz_symmetry!(model, f - lambda, x, [], [1 - sum(x.^2)], d, G, TS="block", SO=1)
@objective(model, Max, lambda)
optimize!(model)
objv = objective_value(model)
@show objv
# optimum = -1.3987174

@polyvar x[1:3]
f = sum(x.^2) + sum(x.^4)
pop = [f, 1 - sum(x.^2)]
G = PermGroup([perm"(1,2,3)", perm"(1,2)"])
opt,data = tssos_symmetry_first(pop, x, 2, G, numeq=0, TS="block")
# optimum = 0
opt,data = tssos_symmetry_higher!(data, TS="block")
# optimum = 0

using SymbolicWedderburn
import GroupsCore

struct DihedralGroup <: GroupsCore.Group
    n::Int
end

struct DihedralElement <: GroupsCore.GroupElement
    n::Int
    reflection::Bool
    id::Int
end

Base.one(G::DihedralGroup) = DihedralElement(G.n, false, 0)

Base.eltype(::Type{DihedralGroup}) = DihedralElement
Base.IteratorSize(::Type{DihedralGroup}) = Base.HasLength()

function Base.iterate(G::DihedralGroup, prev::DihedralElement=DihedralElement(G.n, false, -1))
    if prev.id + 1 >= G.n
        if prev.reflection
            return nothing
        else
            next = DihedralElement(G.n, true, 0)
        end
    else
        next = DihedralElement(G.n, prev.reflection, prev.id + 1)
    end
    return next, next
end

GroupsCore.order(::Type{T}, G::DihedralGroup) where {T} = convert(T, 2G.n)
GroupsCore.gens(G::DihedralGroup) = [DihedralElement(G.n, false, 1), DihedralElement(G.n, true, 0)]

Base.parent(g::DihedralElement) = DihedralGroup(g.n)
function Base.:(==)(g::DihedralElement, h::DihedralElement)
    return g.n == h.n && g.reflection == h.reflection && g.id == h.id
end

function Base.inv(el::DihedralElement)
    if el.reflection || iszero(el.id)
        return el
    else
        return DihedralElement(el.n, false, el.n - el.id)
    end
end
function Base.:*(a::DihedralElement, b::DihedralElement)
    a.n == b.n || error("Cannot multiply elements from different Dihedral groups")
    id = mod(a.reflection ? a.id - b.id : a.id + b.id, a.n)
    return DihedralElement(a.n, a.reflection != b.reflection, id)
end

Base.copy(a::DihedralElement) = DihedralElement(a.n, a.reflection, a.id)

struct DihedralAction <: OnMonomials end

function SymbolicWedderburn.action(::DihedralAction, el::DihedralElement, mono::AbstractMonomial)
    if iseven(el.reflection + el.id)
        var_x, var_y = x, y
    else
        var_x, var_y = y, x
    end
    sign_x = 1 <= el.id <= 2 ? -1 : 1
    sign_y = 2 <= el.id ? -1 : 1
    return mono([x, y] => [sign_x * var_x, sign_y * var_y])
end
SymbolicWedderburn.coeff_type(::DihedralAction) = Float64

@polyvar x y
f = x^6 + y^6 - x^4 * y^2 - y^4 * x^2 - x^4 - y^4 - x^2 - y^2 + 3x^2 * y^2 + 1
opt,data = tssos_symmetry_first([f], [x;y], 3, DihedralGroup(4), action=DihedralAction(), semisimple=false, TS=false)

G = PermGroup(perm"(1,2)")
opt,data = tssos_symmetry_first([f], [x;y], 3, G)
opt,data = tssos_symmetry_higher!(data)


@complex_polyvar z[1:2]
f = z[1]^2*conj(z[1]^2) + z[2]^2*conj(z[2]^2) + z[1]*conj(z[2]) + z[2]*conj(z[1]) + (1+im)*z[1] + 
(1-im)*conj(z[1]) + (1+im)*z[2] + (1-im)*conj(z[2]) + 1
G = PermGroup(perm"(1,2)")
opt,data = complex_tssos_symmetry_first([f, 1-z[1]*conj(z[1])-z[2]*conj(z[2])], z, 2, G, numeq=1, TS="block", ConjugateBasis=true)
opt,data = complex_tssos_symmetry_higher!(data)
opt,data = complex_tssos_symmetry_first([f, 1-z[1]*conj(z[1]), 1-z[2]*conj(z[2])], z, 2, G, numeq=2, SymmetricConstraint=false, TS="block", ConjugateBasis=true)
