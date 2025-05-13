struct VariablePermutation{V} <: SymbolicWedderburn.ByPermutations
    variables::V
end

function SymbolicWedderburn.action(
    a::VariablePermutation,
    g::AbstractPermutations.AbstractPermutation,
    m::AbstractMonomial,
)
    v = a.variables
    return m(v => SymbolicWedderburn.action(a, g, v))
end

"""
    opt,basis,Gram = tssos_symmetry(pop, x, d, group; numeq=0, QUIET=false)

Compute the symmetry adapted moment-SOS relaxation for polynomial optimization problems.

# Input arguments
- `pop`: polynomial optimization problem
- `x`: POP variables
- `d`: relaxation order
- `group`: permutation group acting on POP variables 
- `numeq`: number of equality constraints
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `opt`: optimum
- `basis`: symmetry adapted basis
- `Gram`: Gram matrix
"""
function tssos_symmetry(pop, x, d, group; numeq=0, QUIET=false)
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    dp = [maxdegree(p) for p in pop] 
    time = @elapsed begin
    action = VariablePermutation(x)
    basis_full = MP.monomials(x, 0:2d)
    basis_half = MP.monomials(x, 0:d)
    wedderburn = WedderburnDecomposition(Float64, group, action, basis_full, basis_half)
    basis = [[wedderburn.Uπs[i].basis[j,:]'*basis_half for j=1:size(wedderburn.Uπs[i].basis, 1)] for i = 1:length(wedderburn.Uπs)] # This is the symmetry adapted basis
    end
    if QUIET == false
        mb = maximum(length.(basis))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    tsupp = [minimum(basis_full[item.nzind]) for item in wedderburn.invariants] # This is the basis for the invariant ring
    sort!(tsupp)
    ltsupp = length(tsupp)
    if QUIET == false
        println("Assembling the SDP...")
        println("There are $ltsupp affine constraints.")
    end
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:length(tsupp)]
    pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, length(wedderburn.Uπs))
    for i = 1:length(wedderburn.Uπs)
        pos[i] = @variable(model, [1:length(basis[i]), 1:length(basis[i])], PSD)
        for j = 1:length(basis[i]), k = j:length(basis[i])
            @inbounds supp,coe = poly_norm(basis[i][j]*basis[i][k], group, action)
            for l = 1:length(supp)
                @inbounds ind = bfind(tsupp, ltsupp, supp[l])
                if j == k
                    @inbounds add_to_expression!(cons[ind], coe[l], pos[i][j,k]) 
                else
                    @inbounds add_to_expression!(cons[ind], 2*coe[l], pos[i][j,k]) 
                end
            end
        end
        for s = 2:length(pop) - numeq
            basis_local = basis[i][maxdegree.(basis[i]) .<= d - Int(ceil(dp[s]/2))]
            if !isempty(basis_local)
                gpos = @variable(model, [1:length(basis_local), 1:length(basis_local)], PSD)
                for j = 1:length(basis_local), k = j:length(basis_local)
                    @inbounds supp,coe = poly_norm(basis_local[j]*basis_local[k]*pop[s], group, action)
                    for l = 1:length(supp)
                        @inbounds ind = bfind(tsupp, ltsupp, supp[l])
                        if j == k
                            @inbounds add_to_expression!(cons[ind], coe[l], gpos[j,k]) 
                        else
                            @inbounds add_to_expression!(cons[ind], 2*coe[l], gpos[j,k]) 
                        end
                    end
                end
            end
        end
    end
    if numeq > 0
        for k = 1:numeq
            # ebasis = basis_full[maxdegree.(basis_full) .<= 2d-dp[length(pop)-numeq+k]]
            ebasis = MP.monomials(x, 0:2d-dp[length(pop)-numeq+k])
            free = @variable(model, [1:length(ebasis)])
            for (i,b) in enumerate(ebasis)
                @inbounds supp,coe = poly_norm(b*pop[length(pop)-numeq+k], group, action)
                for l = 1:length(supp)
                    @inbounds ind = bfind(tsupp, ltsupp, supp[l])
                    @inbounds add_to_expression!(cons[ind], coe[l], free[i])
                end
            end
        end
    end
    supp,coe = poly_norm(pop[1], group, action)
    for l = 1:length(supp)
        @inbounds ind = bfind(tsupp, ltsupp, supp[l])
        @inbounds add_to_expression!(cons[ind], -coe[l])
    end
    lambda = @variable(model)
    @objective(model, Max, lambda)
    cons[1] += lambda
    @constraint(model, cons .== 0)
    if QUIET == false
        println("Solving the SDP...")
    end
    time = @elapsed begin
    optimize!(model)
    end
    if QUIET == false
        println("SDP solving time: $time seconds.")
    end
    SDP_status = termination_status(model)
    if SDP_status != MOI.OPTIMAL
        println("termination status: $SDP_status")
        status = primal_status(model)
        println("solution status: $status")
     end
    optimum = objective_value(model)
    @show optimum
    GramMat = [value.(pos[i]) for i = 1:length(wedderburn.Uπs)]
    return optimum,basis,GramMat
end

"""
    supp,coe = poly_norm(p, group, action)

Compute the normal form of a polynomial under the action of a group.

# Input arguments
- `p`: polynomial
- `group`: group
- `action`: describing how the group acts on monomials

# Output arguments
- `supp`: support of the normal form
- `coe`: coefficients of  the normal form
"""
function poly_norm(p, group, action)
    supp = MP.monomials(p)
    coe = MP.coefficients(p)
    nsupp = [normalform(mon, group, action) for mon in supp]
    supp = copy(nsupp)
    sort!(supp)
    unique!(supp)
    ncoe = zeros(length(supp))
    for i = 1:length(nsupp)
        @inbounds ind = bfind(supp, length(supp), nsupp[i])
        @inbounds ncoe[ind] += coe[i]
    end
    return supp, ncoe
end

"""
    mon = normalform(mon, group, action)

Compute the normal form of a monomial under the action of a group.

# Input arguments
- `mon`: monomial
- `group`: group
- `action`: describing how the group acts on monomials

# Output arguments
- `mon`: normal form of the monomial
"""
function normalform(mon, group, action)
    return minimum([SymbolicWedderburn.action(action, g, mon) for g in group])
end
