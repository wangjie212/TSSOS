# import StarAlgebras as SA

struct VariablePermutation{V} <: SymbolicWedderburn.ByPermutations
    variables::V
end

abstract type OnMonomials <: SymbolicWedderburn.ByLinearTransformation end

function SymbolicWedderburn.action(a::VariablePermutation, g::AbstractPermutations.AbstractPermutation, mon::AbstractMonomial)
    return mon(a.variables => SymbolicWedderburn.action(a, g, a.variables))
end

function action_mon(a::VariablePermutation, g::AbstractPermutations.AbstractPermutation, mon::Vector{UInt16})
    return sort!(SymbolicWedderburn.action(a, g, a.variables)[mon])
end

function action_mon(a::VariablePermutation, g::AbstractPermutations.AbstractPermutation, mon::Tuple{Vector{UInt16},Vector{UInt16}})
    v = SymbolicWedderburn.action(a, g, a.variables)
    return tuple(sort!(v[mon[1]]), sort!(v[mon[2]]))
end

function normalform(mon::Vector{UInt16}, group, action)
    return minimum([action_mon(action, g, mon) for g in group])
end

function normalform(mon::Tuple{Vector{UInt16},Vector{UInt16}}, group, action)
    return minimum([action_mon(action, g, mon) for g in group])
end

function SymbolicWedderburn.decompose(
    k::DP.AbstractPolynomialLike,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    # correct only if basis(hom) == monomials
    I = SymbolicWedderburn._int_type(hom)
    indcs = I[hom[mono] for mono in DP.monomials(k)]
    coeffs = DP.coefficients(k)

    return indcs, coeffs
end

"""
    opt,data = tssos_symmetry(pop, x, d, group; numeq=0, SymmetricConstraint=true, QUIET=false)

Compute the symmetry adapted moment-SOS relaxation for polynomial optimization problems.

# Input arguments
- `pop`: polynomial optimization problem
- `x`: POP variables
- `d`: relaxation order
- `group`: permutation group acting on POP variables 
- `numeq`: number of equality constraints
- `SymmetricConstraint`: whether the constraints are symmetric or not
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `opt`: optimum
- `data`: other auxiliary data 
"""
function tssos_symmetry(pop, x, d, group; numeq=0, action=nothing, semisimple=false, SymmetricConstraint=true, QUIET=false, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    return tssos_symmetry_first(pop, x, d, group, numeq=numeq, action=action, semisimple=semisimple, SymmetricConstraint=SymmetricConstraint, TS=false, QUIET=QUIET, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
end

"""
    opt,data = tssos_symmetry_first(pop, x, d, group; numeq=0, SymmetricConstraint=true, TS="block", QUIET=false)

Compute the symmetry adapted moment-SOS relaxation for polynomial optimization problems.

# Input arguments
- `pop`: polynomial optimization problem
- `x`: POP variables
- `d`: relaxation order
- `group`: permutation group acting on POP variables 
- `numeq`: number of equality constraints
- `SymmetricConstraint`: whether the constraints are symmetric or not
- `TS`: type of term sparsity (`"block"`, `"MD"`, `"MF"`, `false`)
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `opt`: optimum
- `data`: other auxiliary data 
"""
function tssos_symmetry_first(pop, x, d, group; numeq=0, action=nothing, semisimple=false, SymmetricConstraint=true, TS="block", QUIET=false, merge=false, md=3, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    cost = poly(pop[1], x)
    ineq_cons = [poly([UInt16[]], [1]); [poly(p, x) for p in pop[2:length(pop)-numeq]]]
    eq_cons = [poly(p, x) for p in pop[length(pop)-numeq+1:end]]
    if QUIET == false
        println("Starting to compute the symmetry adapted basis...")
    end
    time = @elapsed begin
    if action === nothing
        action = tuple(VariablePermutation(x), VariablePermutation(Vector(1:length(x))))
    end
    supp_d = get_basis(Vector(1:length(x)), d)
    supp_2d = get_basis(Vector(1:length(x)), 2d)
    monos_d = Mono[prod(x[item]) for item in supp_d]
    monos_2d = Mono[prod(x[item]) for item in supp_2d]
    wedderburn = WedderburnDecomposition(Float64, group, action[1], monos_2d, monos_d, semisimple=semisimple)
    basis = Vector{Vector{Vector{poly}}}(undef, length(ineq_cons))
    basis[1] = Vector{Vector{poly}}(undef, length(wedderburn.Uπs))
    for (i, ele) in enumerate(wedderburn.Uπs)
        basis[1][i] = Vector{poly{Float64}}(undef, size(ele.basis,1))
        for j = 1:size(ele.basis,1)
            basis[1][i][j] = poly(supp_d[ele.basis[j,:].nzind], ele.basis[j,:].nzval) # This is the symmetry adapted basis
        end
    end
    for i = 2:length(ineq_cons)
        basis[i] = Vector{poly{Float64}}[]
        if SymmetricConstraint == true
            for bas in basis[1]
                ind = maxdeg.(bas) .<= d - ceil(Int, maxdeg(ineq_cons[i])/2)
                if any(ind)
                    push!(basis[i], bas[ind])
                end
            end
        else
            bas = supp_d[length.(supp_d) .<= d - ceil(Int, maxdeg(ineq_cons[i])/2)]
            push!(basis[i], [poly([ba],[1]) for ba in bas])
        end
    end
    ebasis = Vector{Vector{poly}}(undef, numeq)
    if numeq > 0
        if SymmetricConstraint == true
            ebas = [minimum(supp_2d[item.nzind]) for item in wedderburn.invariants]
        else
            ebas = supp_2d
        end
        for (i, h) in enumerate(eq_cons)
            bas = ebas[length.(ebas) .<= 2d-maxdeg(h)]
            ebasis[i] = [poly([ba],[1]) for ba in bas]
        end
    end
    end
    if QUIET == false
        println("Obtained the symmetry adapted basis in $time seconds.")
    end
    time = @elapsed begin
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    tsupp = nothing
    if TS != false
        tsupp = vcat([[normalform(item, group, action[2]) for item in p.supp] for p in [cost;ineq_cons;eq_cons]]...)
        unique!(tsupp)
        sort!(tsupp)
    end
    blocks,cl,blocksize,eblocks = get_blocks(ineq_cons, eq_cons, tsupp, basis, ebasis, group, action[2], TS=TS, merge=merge, md=md)
    end
    if QUIET == false
        mb = maximum(maximum.(blocksize[1]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        println("Assembling the SDP...")
    end
    optimum,tsupp,GramMat,multiplier,SDP_status = solvesdp(cost, ineq_cons, eq_cons, basis, ebasis, cl, blocksize, blocks, eblocks, group, action[2], QUIET=QUIET, 
    solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
    data = poly_basis(cost, ineq_cons, eq_cons, Float64, group, action[2], basis, ebasis, tsupp, blocksize, blocks, eblocks, GramMat, multiplier, SDP_status)
    return optimum,data
end

"""
    opt,data = complex_tssos_symmetry(pop, x, d, group; numeq=0, SymmetricConstraint=true, ConjugateBasis=false, QUIET=false)

Compute the symmetry adapted moment-HSOS relaxation for complex polynomial optimization problems.

# Input arguments
- `pop`: complex polynomial optimization problem
- `x`: POP variables
- `d`: relaxation order
- `group`: permutation group acting on POP variables 
- `numeq`: number of equality constraints
- `SymmetricConstraint`: whether the constraints are symmetric or not
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `opt`: optimum
- `data`: other auxiliary data 
"""
function complex_tssos_symmetry(pop, x, d, group; numeq=0, action=nothing, semisimple=false, SymmetricConstraint=true, ConjugateBasis=false, QUIET=false, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    return complex_tssos_symmetry_first(pop, x, d, group, numeq=numeq, action=action, semisimple=semisimple, SymmetricConstraint=SymmetricConstraint, ConjugateBasis=ConjugateBasis, TS=false, QUIET=QUIET, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
end

"""
    opt,data = complex_tssos_symmetry_first(pop, x, d, group; numeq=0, SymmetricConstraint=true, ConjugateBasis=false, TS="block", QUIET=false)

Compute the symmetry adapted moment-HSOS relaxation for complex polynomial optimization problems.

# Input arguments
- `pop`: complex polynomial optimization problem
- `x`: POP variables
- `d`: relaxation order
- `group`: permutation group acting on POP variables 
- `numeq`: number of equality constraints
- `SymmetricConstraint`: whether the constraints are symmetric or not
- `TS`: type of term sparsity (`"block"`, `"MD"`, `"MF"`, `false`)
- `QUIET`: run in the quiet mode (`true`, `false`)

# Output arguments
- `opt`: optimum
- `data`: other auxiliary data 
"""
function complex_tssos_symmetry_first(pop::Vector{Poly{T}}, x, d, group; numeq=0, action=nothing, semisimple=false, SymmetricConstraint=true, ConjugateBasis=false, TS="block", QUIET=false, merge=false, md=3, 
    dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para()) where {T<:Number}
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    cost = cpoly(pop[1], x)
    ineq_cons = [cpoly([tuple(UInt16[],UInt16[])], [1]); [cpoly(p, x) for p in pop[2:length(pop)-numeq]]]
    eq_cons = [cpoly(p, x) for p in pop[length(pop)-numeq+1:end]]
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    coe_type = T <: Real ? Float64 : ComplexF64
    time = @elapsed begin
    if action === nothing
        action = tuple(VariablePermutation(x), VariablePermutation(Vector(1:length(x))))
    end
    supp_d = get_basis(Vector(1:length(x)), d)
    if ConjugateBasis == true
        supp_d = filter(item -> sum(length.(item)) <= d, vec([tuple(item1, item2) for item1 in supp_d, item2 in supp_d]))
        supp_2d = unique!(vec([sadd(item1,conj(item2)) for item1 in supp_d, item2 in supp_d]))
        monos_d = [prod(x[item[1]])*MP.conj(prod(x[item[2]])) for item in supp_d]
    else
        monos_d = MonomialVector(Mono[prod(x[item]) for item in supp_d])
        supp_2d = vec([tuple(item1,item2) for item1 in supp_d, item2 in supp_d])
    end
    monos_2d = [prod(x[item[1]])*MP.conj(prod(x[item[2]])) for item in supp_2d]
    wedderburn = WedderburnDecomposition(coe_type, group, action[1], monos_2d, monos_d, semisimple=semisimple)
    basis = Vector{Vector{Vector{cpoly}}}(undef, length(ineq_cons))
    basis[1] = Vector{Vector{cpoly}}(undef, length(wedderburn.Uπs))
    for (i, ele) in enumerate(wedderburn.Uπs)
        basis[1][i] = Vector{cpoly}(undef, size(ele.basis,1))
        for j = 1:size(ele.basis,1)
            if ConjugateBasis == true
                supp = supp_d[ele.basis[j,:].nzind]
            else
                supp = [tuple(item, UInt16[]) for item in supp_d[ele.basis[j,:].nzind]]
            end
            basis[1][i][j] = cpoly(supp, ele.basis[j,:].nzval) # This is the symmetry adapted basis
        end
    end
    for i = 2:length(ineq_cons)
        basis[i] = Vector{cpoly}[]
        if SymmetricConstraint == true
            for bas in basis[1]
                if ConjugateBasis == true
                    ind = maxdeg.(bas) .<= d - ceil(Int, maxdeg(ineq_cons[i])/2)
                else
                    ind = maxcdeg.(bas) .<= d - maxcdeg(ineq_cons[i])
                end
                if any(ind)
                    push!(basis[i], bas[ind])
                end
            end
        else
            if ConjugateBasis == true
                ind = [sum(length.(item)) for item in supp_d] .<= d - ceil(Int, maxdeg(ineq_cons[i])/2)
                push!(basis[i], [cpoly([ba],[1]) for ba in supp_d[ind]])
            else
                ind = length.(supp_d) .<= d - maxcdeg(ineq_cons[i])
                push!(basis[i], [cpoly([tuple(ba,UInt16[])],[1]) for ba in supp_d[ind]])
            end
        end
    end
    ebasis = Vector{Vector{cpoly}}(undef, numeq)
    if numeq > 0
        if SymmetricConstraint == true
            ebas = [minimum(supp_2d[item.nzind]) for item in wedderburn.invariants]
        else
            ebas = supp_2d
        end
        for (i, h) in enumerate(eq_cons)
            if ConjugateBasis == true
                ind = [sum(length.(item)) for item in ebas] .<= 2d-maxdeg(h)
            else
                ind = [max(length(item[1]),length(item[2])) for item in ebas] .<= d-maxcdeg(h)
            end
            ebasis[i] = sort!([cpoly([ba],[1]) for ba in ebas[ind]])
        end
    end
    tsupp = nothing
    if TS != false
        tsupp = vcat([[normalform(item, group, action[2]) for item in p.supp] for p in [cost;ineq_cons;eq_cons]]...)
        unique!(tsupp)
        sort!(tsupp)
    end
    blocks,cl,blocksize,eblocks = get_blocks(ineq_cons, eq_cons, tsupp, basis, ebasis, group, action[2], TS=TS, merge=merge, md=md, field="complex")
    end
    if QUIET == false
        mb = maximum(maximum.(blocksize[1]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        println("Assembling the SDP...")
    end
    optimum,tsupp,GramMat,multiplier,SDP_status = solvesdp(cost, ineq_cons, eq_cons, basis, ebasis, cl, blocksize, blocks, eblocks, group, action[2], QUIET=QUIET, 
    coe_type=coe_type, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
    data = poly_basis(cost, ineq_cons, eq_cons, coe_type, group, action[2], basis, ebasis, tsupp, blocksize, blocks, eblocks, GramMat, multiplier, SDP_status)
    return optimum,data
end

function tssos_symmetry_higher!(data::poly_basis; TS="block", merge=false, md=3, QUIET=false, solver="Mosek", dualize=false, field="real", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    eq_cons = data.eq_cons
    ineq_cons = data.ineq_cons
    group = data.group
    action = data.action
    basis = data.basis
    ebasis = data.ebasis
    ksupp = data.ksupp
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,eblocks = get_blocks(ineq_cons, eq_cons, ksupp, basis, ebasis, group, action, TS=TS, merge=merge, md=md, field=field)
    end
    if blocksize == data.blocksize && eblocks == data.eblocks
        println("No higher TS step of the TSSOS hierarchy!")
        opt = nothing
    else
        if QUIET == false
            mb = maximum(maximum.(blocksize[1]))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,GramMat,multiplier,SDP_status = solvesdp(data.cost, ineq_cons, eq_cons, basis, ebasis, cl, blocksize, blocks, eblocks, group, action, 
        QUIET=QUIET, coe_type=data.coe_type, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
        data.blocks = blocks
        data.eblocks = eblocks
        data.blocksize = blocksize
        data.ksupp = ksupp
        data.GramMat = GramMat
        data.multiplier = multiplier
        data.SDP_status = SDP_status
    end
    return opt,data
end

function complex_tssos_symmetry_higher!(data::poly_basis; TS="block", merge=false, md=3, QUIET=false, solver="Mosek", dualize=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    return tssos_symmetry_higher!(data, TS=TS, merge=merge, md=md, QUIET=QUIET, solver=solver, dualize=dualize, field="complex", cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
end

function solvesdp(cost::poly, ineq_cons, eq_cons, basis, ebasis, cl, blocksize, blocks, eblocks, group, action; QUIET=false, coe_type=Float64, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    if solver == "Mosek"
        if dualize == false
            model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
        else
            model = Model(dual_optimizer(Mosek.Optimizer))
        end
    elseif solver == "COSMO"
        model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter, "time_limit" => cosmo_setting.time_limit))
    elseif solver == "SDPT3"
        model = Model(optimizer_with_attributes(SDPT3.Optimizer))
    elseif solver == "SDPNAL"
        model = Model(optimizer_with_attributes(SDPNAL.Optimizer))
    else
        @error "The solver is currently not supported!"
    end
    time = @elapsed begin
    tsupp = Vector{UInt16}[]
    if all(cl[1] .== 1)
        for (i, bas) in enumerate(basis[1]), j = 1:blocksize[1][i][1], t = j:blocksize[1][i][1]
            append!(tsupp, supp_multi(bas[blocks[1][i][1][j]], bas[blocks[1][i][1][t]], group, action))
        end
    else
        for (k, g) in enumerate(ineq_cons), (i, bas) in enumerate(basis[k]), l = 1:cl[k][i], j = 1:blocksize[k][i][l], t = j:blocksize[k][i][l]
            append!(tsupp, supp_multi(bas[blocks[k][i][l][j]], bas[blocks[k][i][l][t]], group, action, g=g))
        end
        for (k, h) in enumerate(eq_cons), i in eblocks[k]
            append!(tsupp, supp_multi(ebasis[k][i], h, group, action))
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp = length(tsupp)
    if QUIET == false
        println("There are $ltsupp affine constraints.")
    end
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:ltsupp]
    pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, length(ineq_cons))
    for (k, g) in enumerate(ineq_cons)
        pos[k] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[k]))
        for (i, bas) in enumerate(basis[k])
            pos[k][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k][i])
            for (l, bs) in enumerate(blocksize[k][i])
                if bs == 1
                    pos[k][i][l] = @variable(model, lower_bound=0)
                    for (s,item1) in enumerate(bas[blocks[k][i][l][1]].supp), (u,item2) in enumerate(bas[blocks[k][i][l][1]].supp), (v,item3) in enumerate(g.supp)
                        @inbounds ind = bfind(tsupp, ltsupp, normalform(sadd(sadd(item1,item2),item3), group, action))
                        @inbounds add_to_expression!(cons[ind], bas[blocks[k][i][l][1]].coe[s]*bas[blocks[k][i][l][1]].coe[u]*g.coe[v], pos[k][i][l])
                    end
                else
                    pos[k][i][l] = @variable(model, [1:bs, 1:bs], PSD)
                    for j = 1:blocksize[k][i][l], t = j:blocksize[k][i][l], (s,item1) in enumerate(bas[blocks[k][i][l][j]].supp), (u,item2) in enumerate(bas[blocks[k][i][l][t]].supp), (v,item3) in enumerate(g.supp) 
                        @inbounds ind = bfind(tsupp, ltsupp, normalform(sadd(sadd(item1,item2),item3), group, action))
                        if j == t
                            @inbounds add_to_expression!(cons[ind], bas[blocks[k][i][l][j]].coe[s]*bas[blocks[k][i][l][t]].coe[u]*g.coe[v], pos[k][i][l][j,t])
                        else
                            @inbounds add_to_expression!(cons[ind], 2*bas[blocks[k][i][l][j]].coe[s]*bas[blocks[k][i][l][t]].coe[u]*g.coe[v], pos[k][i][l][j,t])
                        end
                    end
                end
            end
        end
    end
    mul = Vector{Vector{VariableRef}}(undef, length(eq_cons))
    for (k, h) in enumerate(eq_cons)
        mul[k] = @variable(model, [1:length(eblocks[k])])
        for (j, i) in enumerate(eblocks[k]), (v, item) in enumerate(h.supp)
            @inbounds ind = bfind(tsupp, ltsupp, normalform(sadd(ebasis[k][i].supp[1],item), group, action))
            @inbounds add_to_expression!(cons[ind], h.coe[v], mul[k][j])
        end
    end
    lambda = @variable(model)
    @objective(model, Max, lambda)
    for (i,item) in enumerate(cost.supp)
        Locb = bfind(tsupp, ltsupp, normalform(item, group, action))
        if Locb === nothing
            @error "The basis is not enough!"
            return nothing,nothing,nothing,nothing,nothing
        else
            cons[Locb] -= cost.coe[i]
        end
    end
    cons[1] += lambda
    @constraint(model, cons==zeros(ltsupp))
    end
    if QUIET == false
        println("SDP assembling time: $time seconds.")
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
    GramMat = Vector{Vector{Vector{Union{Float64,Matrix{Float64}}}}}(undef, length(ineq_cons))
    multiplier = nothing
    for i = 1:length(ineq_cons)
        GramMat[i] = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, length(basis[i]))
        for j = 1:length(basis[i])
            GramMat[i][j] = [value.(pos[i][j][k]) for k = 1:cl[i][j]]
        end
        if length(eq_cons) > 0
            multiplier = [value.(mul[j]) for j = 1:length(eq_cons)]
        end
    end
    return optimum,tsupp,GramMat,multiplier,SDP_status
end

function solvesdp(cost::cpoly, ineq_cons, eq_cons, basis, ebasis, cl, blocksize, blocks, eblocks, group, action; QUIET=false, coe_type=ComplexF64, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    if solver == "Mosek"
        if dualize == false
            model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
        else
            model = Model(dual_optimizer(Mosek.Optimizer))
        end
    elseif solver == "COSMO"
        model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter, "time_limit" => cosmo_setting.time_limit))
    elseif solver == "SDPT3"
        model = Model(optimizer_with_attributes(SDPT3.Optimizer))
    elseif solver == "SDPNAL"
        model = Model(optimizer_with_attributes(SDPNAL.Optimizer))
    else
        @error "The solver is currently not supported!"
    end
    time = @elapsed begin
    tsupp = Tuple{Vector{UInt16},Vector{UInt16}}[]
    if all(cl[1] .== 1)
        for (i, bas) in enumerate(basis[1]), j = 1:blocksize[1][i][1], t = 1:blocksize[1][i][1]
            supp = supp_multi(bas[blocks[1][i][1][j]], conj(bas[blocks[1][i][1][t]]), group, action)
            append!(tsupp, supp[[item[1] <= item[2] for item in supp]])
        end
    else
        for (k, g) in enumerate(ineq_cons), (i, bas) in enumerate(basis[k]), l = 1:cl[k][i], j = 1:blocksize[k][i][l], t = 1:blocksize[k][i][l]
            supp = supp_multi(bas[blocks[k][i][l][j]], conj(bas[blocks[k][i][l][t]]), group, action, g=g)
            append!(tsupp, supp[[item[1] <= item[2] for item in supp]])
        end
        for (k, h) in enumerate(eq_cons), i in eblocks[k]
            supp = supp_multi(ebasis[k][i], h, group, action)
            append!(tsupp, supp[[item[1] <= item[2] for item in supp]])
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp = length(tsupp)
    rcons = [AffExpr(0) for i=1:ltsupp]
    if coe_type == ComplexF64
        icons = [AffExpr(0) for i=1:ltsupp]
    end
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, length(ineq_cons))
    for (k, g) in enumerate(ineq_cons)
        pos[k] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[k]))
        for (i, bas) in enumerate(basis[k])
            pos[k][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k][i])
            for (l, bs) in enumerate(blocksize[k][i])
                if bs == 1
                    pos[k][i][l] = @variable(model, lower_bound=0)
                    for (s,item1) in enumerate(bas[blocks[k][i][l][1]].supp), (u,item2) in enumerate(bas[blocks[k][i][l][1]].supp), (v,item3) in enumerate(g.supp)
                        bi = normalform(sadd(sadd(item1,conj(item2)),item3), group, action)
                        if bi[1] <= bi[2]
                            @inbounds ind = bfind(tsupp, ltsupp, bi)
                            c = bas[blocks[k][i][l][1]].coe[s]*conj(bas[blocks[k][i][l][1]].coe[u])*g.coe[v]
                            @inbounds add_to_expression!(rcons[ind], real(c), pos[k][i][l])
                            if coe_type == ComplexF64
                                @inbounds add_to_expression!(icons[ind], imag(c), pos[k][i][l])
                            end
                        end
                    end
                else
                    if coe_type == ComplexF64
                        pos[k][i][l] = @variable(model, [1:2bs, 1:2bs], PSD)
                    else
                        pos[k][i][l] = @variable(model, [1:bs, 1:bs], PSD)
                    end
                    for j = 1:blocksize[k][i][l], t = 1:blocksize[k][i][l], (s,item1) in enumerate(bas[blocks[k][i][l][j]].supp), (u,item2) in enumerate(bas[blocks[k][i][l][t]].supp), (v,item3) in enumerate(g.supp)
                        bi = normalform(sadd(sadd(item1,conj(item2)),item3), group, action)
                        if bi[1] <= bi[2]
                            @inbounds ind = bfind(tsupp, ltsupp, bi)
                            if coe_type == ComplexF64
                                c = bas[blocks[k][i][l][j]].coe[s]*conj(bas[blocks[k][i][l][t]].coe[u])*g.coe[v]
                                @inbounds add_to_expression!(rcons[ind], real(c), pos[k][i][l][j,t]+pos[k][i][l][j+bs,t+bs])
                                @inbounds add_to_expression!(rcons[ind], -imag(c), pos[k][i][l][j,t+bs]-pos[k][i][l][t,j+bs])
                                if bi[1] != bi[2]
                                    @inbounds add_to_expression!(icons[ind], imag(c), pos[k][i][l][j,t]+pos[k][i][l][j+bs,t+bs])
                                    @inbounds add_to_expression!(icons[ind], real(c), pos[k][i][l][j,t+bs]-pos[k][i][l][t,j+bs])
                                end
                            else
                                @inbounds add_to_expression!(rcons[ind], bas[blocks[k][i][l][j]].coe[s]*bas[blocks[k][i][l][t]].coe[u]*g.coe[v], pos[k][i][l][j,t])
                            end
                        end
                    end
                end
            end
        end
    end
    mul = Vector{Vector{VariableRef}}(undef, length(eq_cons))
    for (k, h) in enumerate(eq_cons)
        temp = ebasis[k][eblocks[k]][[item <= conj(item) for item in ebasis[k][eblocks[k]]]]
        lb = length(temp)
        if coe_type == ComplexF64
            mul[k] = @variable(model, [1:2*lb])
        else
            mul[k] = @variable(model, [1:lb])
        end
        for i in eblocks[k], (v, item) in enumerate(h.supp)
            ba = ebasis[k][i].supp[1]
            bi = normalform(sadd(ba,item), group, action)
            if bi[1] <= bi[2]
                @inbounds ind = bfind(tsupp, ltsupp, bi)
                if ba[1] <= ba[2]
                    loc = bfind(temp, lb, ebasis[k][i])
                    tag = ba[1] == ba[2] ? 0 : 1
                else
                    loc = bfind(temp, lb, conj(ebasis[k][i]))
                    tag = -1
                end
                if coe_type == ComplexF64
                    @inbounds add_to_expression!(rcons[ind], real(h.coe[v])*mul[k][loc] - tag*imag(h.coe[v])*mul[k][loc+lb])
                    if bi[1] != bi[2]
                       @inbounds add_to_expression!(icons[ind], tag*real(h.coe[v])*mul[k][loc+lb] + imag(h.coe[v])*mul[k][loc])
                    end
                else
                    @inbounds add_to_expression!(rcons[ind], h.coe[v], mul[k][loc])
                end
            end
        end
    end
    for (i,item) in enumerate(cost.supp)
        bi = normalform(item, group, action)
        if bi[1] <= bi[2]
            Locb = bfind(tsupp, ltsupp, bi)
            if Locb === nothing
                @error "The basis is not enough!"
                return nothing,nothing,nothing,nothing,nothing
            else
                rcons[Locb] -= real(cost.coe[i])
                if coe_type == ComplexF64 && bi[1] != bi[2]
                    icons[Locb] -= imag(cost.coe[i])
                end
            end
        end
    end
    lambda = @variable(model)
    @objective(model, Max, lambda)
    rcons[1] += lambda
    @constraint(model, rcons .== 0)
    if coe_type == ComplexF64
        icons = icons[[item[1] != item[2] && icons[i] != 0 for (i,item) in enumerate(tsupp)]]
        @constraint(model, icons .== 0)
        ltsupp += length(icons)
    end
    end
    if QUIET == false
        println("There are $ltsupp affine constraints.")
        println("SDP assembling time: $time seconds.")
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
    GramMat = Vector{Vector{Vector{Union{Float64,Matrix{coe_type}}}}}(undef, length(ineq_cons))
    multiplier = nothing
    for i = 1:length(ineq_cons)
        GramMat[i] = Vector{Vector{Union{Float64,Matrix{coe_type}}}}(undef, length(basis[i]))
        for j = 1:length(basis[i])
            if coe_type == Float64
                GramMat[i][j] = [value.(pos[i][j][k]) for k = 1:cl[i][j]]
            else
                GramMat[i][j] = Vector{Union{Float64,Matrix{coe_type}}}(undef, cl[i][j])
                for k = 1:cl[i][j]
                    bs = blocksize[i][j][k]
                    if bs == 1
                        GramMat[i][j][k] = value(pos[i][j][k])
                    else
                        GramMat[i][j][k] = value.(pos[i][j][k][1:bs,1:bs])+value.(pos[i][j][k][bs+1:2*bs,bs+1:2*bs])+value.(pos[i][j][k][bs+1:2*bs,1:bs])*im-value.(pos[i][j][k][1:bs,bs+1:2*bs])*im
                    end
                end
            end
        end
    end
    # if length(eq_cons) > 0
    #     multiplier = Vector{Poly{coe_type}}(undef, length(eq_cons))
    #     for k = 1:length(eq_cons)
    #         temp = ebasis[k][eblocks[k]][[item <= conj(item) for item in ebasis[k][eblocks[k]]]]
    #         lb = length(temp)
    #         mu = [prod(x[item[1]])*conj(prod(x[item[2]])) for item in temp]
    #         if coe_type == Float64
    #             tau = sum(mu .* value.(mul[k]))
    #         else
    #             tau = sum(mu .* (value.(mul[k][1:lb]) .+ value.(mul[k][lb+1:2*lb])*im))
    #         end
    #         multiplier[k] = tau + conj(tau)
    #     end
    # end
    return optimum,tsupp,GramMat,multiplier,SDP_status
end

function get_blocks(ineq_cons, eq_cons, tsupp, basis::Vector{Vector{Vector{T}}}, ebasis::Vector{Vector{T}}, group, action;
     TS="block", merge=false, md=3, field="real") where {T <: Union{poly,cpoly}}
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, length(ineq_cons))
    blocksize = Vector{Vector{Vector{Int}}}(undef, length(ineq_cons))
    cl = Vector{Vector{Int}}(undef, length(ineq_cons))
    if TS == false
        for (k,item) in enumerate(basis)
            blocks[k],blocksize[k],cl[k] = [[Vector(1:length(bas))] for bas in item],[[length(bas)] for bas in item],ones(Int, length(item))
        end
        eblocks = [Vector(1:length(bas)) for bas in ebasis]
    else
        for (k, g) in enumerate(ineq_cons)
            blocks[k] = Vector{Vector{Vector{Int}}}(undef, length(basis[k]))
            blocksize[k] = Vector{Vector{Int}}(undef, length(basis[k]))
            cl[k] = Vector{Int}(undef, length(basis[k]))
            for (i, ba) in enumerate(basis[k])
                G = get_graph(tsupp, ba, group, action, g=g, field=field)
                if TS == "block"
                    blocks[k][i] = connected_components(G)
                    blocksize[k][i] = length.(blocks[k][i])
                    cl[k][i] = length(blocksize[k][i])
                else
                    blocks[k][i],cl[k][i],blocksize[k][i] = chordal_cliques!(G, method=TS)
                    if merge == true
                        blocks[k][i],cl[k][i],blocksize[k][i] = clique_merge!(blocks[k], d=md, QUIET=true)
                    end
                end
            end
        end
        eblocks = [get_eblock(tsupp, h, ebasis[k], group, action) for (k, h) in enumerate(eq_cons)]
    end
    return blocks,cl,blocksize,eblocks
end

function get_graph(tsupp, basis::Vector{T}, group, action; g=poly([UInt16[]], [1]), field="real") where {T <: Union{poly,cpoly}}
    lb = length(basis)
    G = SimpleGraph(lb)
    ltsupp = length(tsupp)
    for i = 1:lb, j = i+1:lb
        if field == "real" 
            supp = supp_multi(basis[i], basis[j], group, action, g=g)
        else
            supp = supp_multi(basis[i], conj(basis[j]), group, action, g=g)
        end
        if findfirst(item -> bfind(tsupp, ltsupp, item) !== nothing, supp) !== nothing
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_eblock(tsupp, h, basis::Vector{T}, group, action) where {T <: Union{poly,cpoly}}
    ltsupp = length(tsupp)
    eblock = Int[]
    for (i, ba) in enumerate(basis)
        supp = supp_multi(ba, h, group, action)
        if findfirst(item -> bfind(tsupp, ltsupp, item) !== nothing, supp) !== nothing
            push!(eblock, i)
        end
    end
    return eblock
end

function add_psatz_symmetry!(model, nonneg::Poly{T}, vars, ineq_cons, eq_cons, order, group; action=nothing, semisimple=false, SymmetricConstraint=true, TS="block", SO=1, merge=false, md=3, QUIET=false) where {T<:Union{Number,AffExpr}}
    obj = poly(nonneg, vars)
    ineq_cons = [poly([UInt16[]], [1]); [poly(p, vars) for p in ineq_cons]]
    eq_cons = [poly(p, vars) for p in eq_cons]
    if QUIET == false
        println("Starting to compute the symmetry adapted basis...")
    end
    time = @elapsed begin
    if action === nothing
        action = tuple(VariablePermutation(vars), VariablePermutation(Vector(1:length(vars))))
    end
    supp_d = get_basis(Vector(1:length(vars)), order)
    supp_2d = get_basis(Vector(1:length(vars)), 2*order)
    monos_d = Mono[prod(vars[item]) for item in supp_d]
    monos_2d = Mono[prod(vars[item]) for item in supp_2d]
    wedderburn = WedderburnDecomposition(Float64, group, action[1], monos_2d, monos_d, semisimple=semisimple)
    basis = Vector{Vector{Vector{poly}}}(undef, length(ineq_cons))
    basis[1] = Vector{Vector{poly}}(undef, length(wedderburn.Uπs))
    for (i, ele) in enumerate(wedderburn.Uπs)
        basis[1][i] = Vector{poly{Float64}}(undef, size(ele.basis,1))
        for j = 1:size(ele.basis,1)
            basis[1][i][j] = poly(supp_d[ele.basis[j,:].nzind], ele.basis[j,:].nzval) # This is the symmetry adapted basis
        end
    end
    for i = 2:length(ineq_cons)
        basis[i] = Vector{poly{Float64}}[]
        if SymmetricConstraint == true
            for bas in basis[1]
                ind = maxdeg.(bas) .<= order - ceil(Int, maxdeg(ineq_cons[i])/2)
                if any(ind)
                    push!(basis[i], bas[ind])
                end
            end
        else
            bas = supp_d[length.(supp_d) .<= order - ceil(Int, maxdeg(ineq_cons[i])/2)]
            push!(basis[i], [poly([ba],[1]) for ba in bas])
        end
    end
    ebasis = Vector{Vector{poly}}(undef,length(eq_cons))
    if length(eq_cons) > 0
        if SymmetricConstraint == true
            ebas = [minimum(supp_2d[item.nzind]) for item in wedderburn.invariants]
        else
            ebas = supp_2d
        end
        for (i, h) in enumerate(eq_cons)
            bas = ebas[length.(ebas) .<= 2*order-maxdeg(h)]
            ebasis[i] = [poly([ba],[1]) for ba in bas]
        end
    end
    end
    if QUIET == false
        println("Obtained the symmetry adapted basis in $time seconds.")
    end
    tsupp = nothing
    if TS != false
        tsupp = vcat([[normalform(item, group, action[2]) for item in p.supp] for p in [obj;ineq_cons;eq_cons]]...)
        unique!(tsupp)
        sort!(tsupp)
    end
    blocks,cl,blocksize,eblocks = get_nblocks(ineq_cons, eq_cons, tsupp, basis, ebasis, group, action[2], TS=TS, SO=SO, merge=merge, md=md)
    tsupp = Vector{UInt16}[]
    if all(cl[1] .== 1)
        for (i, bas) in enumerate(basis[1]), j = 1:blocksize[1][i][1], t = j:blocksize[1][i][1]
            append!(tsupp, supp_multi(bas[blocks[1][i][1][j]], bas[blocks[1][i][1][t]], group, action[2]))
        end
    else
        for (k, g) in enumerate(ineq_cons), (i, bas) in enumerate(basis[k]), l = 1:cl[k][i], j = 1:blocksize[k][i][l], t = j:blocksize[k][i][l]
            append!(tsupp, supp_multi(bas[blocks[k][i][l][j]], bas[blocks[k][i][l][t]], group, action[2], g=g))
        end
        for (k, h) in enumerate(eq_cons), i in eblocks[k]
            append!(tsupp, supp_multi(ebasis[k][i], h, group, action[2]))
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp = length(tsupp)
    cons = [AffExpr(0) for i=1:ltsupp]
    pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, length(ineq_cons))
    for (k, g) in enumerate(ineq_cons)
        pos[k] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[k]))
        for (i, bas) in enumerate(basis[k])
            pos[k][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k][i])
            for (l, bs) in enumerate(blocksize[k][i])
                if bs == 1
                    pos[k][i][l] = @variable(model, lower_bound=0)
                    for (s,item1) in enumerate(bas[blocks[k][i][l][1]].supp), (u,item2) in enumerate(bas[blocks[k][i][l][1]].supp), (v,item3) in enumerate(g.supp)
                        @inbounds ind = bfind(tsupp, ltsupp, normalform(sadd(sadd(item1,item2),item3), group, action[2]))
                        @inbounds add_to_expression!(cons[ind], bas[blocks[k][i][l][1]].coe[s]*bas[blocks[k][i][l][1]].coe[u]*g.coe[v], pos[k][i][l])
                    end
                else
                    pos[k][i][l] = @variable(model, [1:bs, 1:bs], PSD)
                    for j = 1:blocksize[k][i][l], t = j:blocksize[k][i][l], (s,item1) in enumerate(bas[blocks[k][i][l][j]].supp), (u,item2) in enumerate(bas[blocks[k][i][l][t]].supp), (v,item3) in enumerate(g.supp) 
                        @inbounds ind = bfind(tsupp, ltsupp, normalform(sadd(sadd(item1,item2),item3), group, action[2]))
                        if j == t
                            @inbounds add_to_expression!(cons[ind], bas[blocks[k][i][l][j]].coe[s]*bas[blocks[k][i][l][t]].coe[u]*g.coe[v], pos[k][i][l][j,t])
                        else
                            @inbounds add_to_expression!(cons[ind], 2*bas[blocks[k][i][l][j]].coe[s]*bas[blocks[k][i][l][t]].coe[u]*g.coe[v], pos[k][i][l][j,t])
                        end
                    end
                end
            end
        end
    end
    mul = Vector{Vector{VariableRef}}(undef, length(eq_cons))
    for (k, h) in enumerate(eq_cons)
        mul[k] = @variable(model, [1:length(eblocks[k])])
        for (j, i) in enumerate(eblocks[k]), (v, item) in enumerate(h.supp)
            @inbounds ind = bfind(tsupp, ltsupp, normalform(sadd(ebasis[k][i].supp[1],item), group, action[2]))
            @inbounds add_to_expression!(cons[ind], h.coe[v], mul[k][j])
        end
    end
    for (i,item) in enumerate(obj.supp)
        Locb = bfind(tsupp, ltsupp, normalform(item, group, action[2]))
        if Locb === nothing
            @error "The basis is not enough!"
            return nothing
        else
            cons[Locb] -= obj.coe[i]
        end
    end
    @constraint(model, cons .== 0)
    info = poly_basis(obj, ineq_cons, eq_cons, nothing, group, action[2], basis, ebasis, tsupp, blocksize, blocks, eblocks, pos, mul, nothing)
    return info
end

function get_nblocks(ineq_cons, eq_cons, tsupp, basis::Vector{Vector{Vector{poly}}}, ebasis::Vector{Vector{poly}}, group, action; TS="block", SO=1, merge=false, md=3)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, length(ineq_cons))
    blocksize = Vector{Vector{Vector{Int}}}(undef, length(ineq_cons))
    cl = Vector{Vector{Int}}(undef, length(ineq_cons))
    if TS == false
        for (k,item) in enumerate(basis)
            blocks[k],blocksize[k],cl[k] = [[Vector(1:length(bas))] for bas in item],[[length(bas)] for bas in item],ones(Int, length(item))
        end
        eblocks = [Vector(1:length(bas)) for bas in ebasis]
    else
        for i = 1:SO
            if i > 1
                oblocksize = deepcopy(blocksize)
                oeblocks = deepcopy(eblocks)
            end
            for (k, g) in enumerate(ineq_cons)
                blocks[k] = Vector{Vector{Vector{Int}}}(undef, length(basis[k]))
                blocksize[k] = Vector{Vector{Int}}(undef, length(basis[k]))
                cl[k] = Vector{Int}(undef, length(basis[k]))
                for (i, ba) in enumerate(basis[k])
                    G = get_graph(tsupp, ba, group, action, g=g)
                    if TS == "block"
                        blocks[k][i] = connected_components(G)
                        blocksize[k][i] = length.(blocks[k][i])
                        cl[k][i] = length(blocksize[k][i])
                    else
                        blocks[k][i],cl[k][i],blocksize[k][i] = chordal_cliques!(G, method=TS)
                        if merge == true
                            blocks[k][i],cl[k][i],blocksize[k][i] = clique_merge!(blocks[k], d=md, QUIET=true)
                        end
                    end
                end
            end
            eblocks = [get_eblock(tsupp, h, ebasis[k], group, action) for (k, h) in enumerate(eq_cons)]
            if i > 1 && blocksize == oblocksize && eblocks == oeblocks
                println("No higher TS step of the TSSOS hierarchy!")
                break
            end
            if i < SO
                tsupp = Vector{UInt16}[]
                for (i, bas) in enumerate(basis[1]), j = 1:blocksize[1][i][1], t = j:blocksize[1][i][1]
                    append!(tsupp, supp_multi(bas[blocks[1][i][1][j]], bas[blocks[1][i][1][t]], group, action))
                end
                unique!(tsupp)
                sort!(tsupp)
            end
        end
    end
    return blocks,cl,blocksize,eblocks
end
