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
    opt,data = tssos_symmetry(pop, x, d, group; numeq=0, QUIET=false)

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
function tssos_symmetry(pop, x, d, group; numeq=0, SymmetricConstraint=true, QUIET=false, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    return tssos_symmetry_first(pop, x, d, group, numeq=numeq, SymmetricConstraint=SymmetricConstraint, TS=false, QUIET=QUIET, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
end

function tssos_symmetry_first(pop, x, d, group; numeq=0, SymmetricConstraint=true, TS="block", QUIET=false, merge=false, md=3, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    m = length(pop) - 1 - numeq
    time = @elapsed begin
    action = VariablePermutation(x)
    monos_2d = MP.monomials(x, 0:2d)
    monos_d = MP.monomials(x, 0:d)
    wedderburn = WedderburnDecomposition(Float64, group, action, monos_2d, monos_d)
    basis = Vector{Vector{Vector{Poly{Float64}}}}(undef, m+1)
    basis[1] = [[item'*monos_d for item in eachrow(ele.basis)] for ele in wedderburn.Uπs] # This is the symmetry adapted basis
    for i = 2:m+1
        basis[i] = Vector{Poly{Float64}}[]
        if SymmetricConstraint == true
            for bas in basis[1]
                ind = maxdegree.(bas) .<= d - ceil(Int, maxdegree(pop[i])/2)
                if any(ind)
                    push!(basis[i], bas[ind])
                end
            end
        else
            ind = maxdegree.(monos_d) .<= d - ceil(Int, maxdegree(pop[i])/2)
            push!(basis[i], monos_d[ind])
        end
    end
    ebasis = Vector{Vector{Poly{Float64}}}(undef, numeq)
    if numeq > 0
        if SymmetricConstraint == true
            invbas = [minimum(monos_2d[item.nzind]) for item in wedderburn.invariants] # This is the basis for the invariant ring
            for i = 1:numeq
                ebasis[i] = invbas[maxdegree.(invbas) .<= 2d-maxdegree(pop[length(pop)-numeq+i])]
            end
        else
            for i = 1:numeq
                ebasis[i] = monos_2d[maxdegree.(monos_2d) .<= 2d-maxdegree(pop[length(pop)-numeq+i])]
            end
        end
    end
    tsupp = nothing
    if TS != false
        tsupp = vcat([poly_norm(p, group, action)[1] for p in pop]...)
        unique!(tsupp)
        sort!(tsupp)
    end
    blocks,cl,blocksize,eblocks = get_blocks(m, numeq, tsupp, pop, basis, ebasis, group, action, TS=TS, merge=merge, md=md)
    end
    if QUIET == false
        mb = maximum(maximum.(blocksize[1]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        println("Assembling the SDP...")
    end
    optimum,tsupp,GramMat,multiplier,SDP_status = solvesdp(pop, basis, ebasis, cl, blocksize, blocks, eblocks, group, action, numeq=numeq, QUIET=QUIET, 
    solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
    data = poly_basis(pop, numeq, group, action, basis, ebasis, tsupp, blocksize, blocks, eblocks, GramMat, multiplier, SDP_status)
    return optimum,data
end

function tssos_symmetry_higher!(data::poly_basis; TS="block", merge=false, md=3, QUIET=false, solver="Mosek", dualize=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    pop = data.pop
    numeq = data.numeq
    group = data.group
    action = data.action
    basis = data.basis
    ebasis = data.ebasis
    ksupp = data.ksupp
    m = length(pop) - 1 - numeq
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,eblocks = get_blocks(m, numeq, ksupp, pop, basis, ebasis, group, action, TS=TS, merge=merge, md=md)
    end
    if blocksize == data.blocksize && eblocks == data.eblocks
        println("No higher TS step of the TSSOS hierarchy!")
        opt = nothing
    else
        if QUIET == false
            mb = maximum(maximum.(blocksize[1]))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,GramMat,multiplier,SDP_status = solvesdp(pop, basis, ebasis, cl, blocksize, blocks, eblocks, group, action, numeq=numeq, 
        QUIET=QUIET, solver=solver, dualize=dualize, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
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

function solvesdp(pop, basis, ebasis, cl, blocksize, blocks, eblocks, group, action; numeq=0, QUIET=false, dualize=false, solver="Mosek", cosmo_setting=cosmo_para(), mosek_setting=mosek_para())
    m = length(pop) - 1 - numeq
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
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    lambda = @variable(model)
    @objective(model, Max, lambda)
    poly = pop[1] - lambda
    pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, 1+m)
    pos[1] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[1]))
    for (i, bas) in enumerate(basis[1])
        pos[1][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1][i])
        for (l, bs) in enumerate(blocksize[1][i])
            if bs == 1
               pos[1][i][l] = @variable(model, lower_bound=0)
               poly -= bas[blocks[1][i][l][1]]^2 * pos[1][i][l]
            else
               pos[1][i][l] = @variable(model, [1:bs, 1:bs], PSD)
               poly -= bas[blocks[1][i][l]]' * pos[1][i][l] * bas[blocks[1][i][l]]
            end
        end
    end
    for k = 1:m
        pos[k+1] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[k+1]))
        for (i, bas) in enumerate(basis[k+1])
            pos[k+1][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k+1][i])
            for (l, bs) in enumerate(blocksize[k+1][i])
                if bs == 1
                    pos[k+1][i][l] = @variable(model, lower_bound=0)
                    poly -= bas[blocks[k+1][i][l][1]]^2 * pos[k+1][i][l] * pop[k+1]
                else
                    pos[k+1][i][l] = @variable(model, [1:bs, 1:bs], PSD)
                    poly -= bas[blocks[k+1][i][l]]' * pos[k+1][i][l] * bas[blocks[k+1][i][l]] * pop[k+1]
                end
            end
        end
    end
    if numeq > 0
        mul = Vector{Vector{VariableRef}}(undef, numeq)
        for k = 1:numeq
            mul[k] = @variable(model, [1:length(eblocks[k])])
            poly -= ebasis[k][eblocks[k]]' * mul[k] * pop[k+m+1]
        end
    end
    tsupp,coe = poly_norm(poly, group, action)
    @constraint(model, coe .== 0)
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
    GramMat = Vector{Vector{Vector{Union{Float64,Matrix{Float64}}}}}(undef, m+1)
    multiplier = nothing
    for i = 1:m+1
        GramMat[i] = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, length(basis[i]))
        for j = 1:length(basis[i])
            GramMat[i][j] = [value.(pos[i][j][k]) for k = 1:cl[i][j]]
        end
        if numeq > 0
            multiplier = [value.(mul[j]) for j = 1:numeq]
        end
    end
    return optimum,tsupp,GramMat,multiplier,SDP_status
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
    unique!(supp)
    sort!(supp)
    ncoe = zeros(typeof(coe[1]), length(supp))
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

function get_blocks(m::Int, l::Int, tsupp, pop, basis::Vector{Vector{Vector{Poly{Float64}}}}, ebasis::Vector{Vector{Poly{Float64}}}, group, action; TS="block", merge=false, md=3)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, m+1)
    eblocks = Vector{Vector{Int}}(undef, l)
    blocksize = Vector{Vector{Vector{Int}}}(undef, m+1)
    cl = Vector{Vector{Int}}(undef, m+1)
    if TS == false
        for k = 1:m+1
            blocks[k],blocksize[k],cl[k] = [[Vector(1:length(basis[k][i]))] for i = 1:length(basis[k])],[[length(basis[k][i])] for i = 1:length(basis[k])],ones(Int, length(basis[k]))
        end
        for k = 1:l
            eblocks[k] = Vector(1:length(ebasis[k]))
        end
    else
        for k = 1:m+1
            blocks[k] = Vector{Vector{Vector{Int}}}(undef, length(basis[k]))
            blocksize[k] = Vector{Vector{Int}}(undef, length(basis[k]))
            cl[k] = Vector{Int}(undef, length(basis[k]))
            for (i, ba) in enumerate(basis[k])
                if k == 1
                    G = get_graph(tsupp, ba, group, action)
                else
                    G = get_graph(tsupp, ba, group, action, g=pop[k-1])
                end
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
        for k = 1:l
            eblocks[k] = get_eblock(tsupp, pop[k+m+1], ebasis[k], group, action)
        end
    end
    return blocks,cl,blocksize,eblocks
end

function get_graph(tsupp, basis::Vector{Poly{Float64}}, group, action; g=1)
    lb = length(basis)
    G = SimpleGraph(lb)
    ltsupp = length(tsupp)
    for i = 1:lb, j = i+1:lb
        flag = 0
        for item in poly_norm(basis[i] * basis[j] * g, group, action)[1]
            if lbfind(tsupp, ltsupp, item) !== nothing
                flag = 1
                break
            end
        end
        if flag == 1
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_eblock(tsupp, h, basis::Vector{Poly{Float64}}, group, action)
    ltsupp = length(tsupp)
    eblock = Int[]
    for (i, ba) in enumerate(basis)
        flag = 0
        for item in poly_norm(ba * h, group, action)[1]
            if lbfind(tsupp, ltsupp, item) !== nothing
                flag = 1
                break
            end
        end
        if flag == 1
            push!(eblock, i)
        end
    end
    return eblock
end

function add_psatz_symmetry!(model, nonneg::DP.Polynomial{V, M, T}, vars, ineq_cons, eq_cons, order, group; SymmetricConstraint=true, TS="block", SO=1, merge=false, md=3, QUIET=false) where {V, M, T<:Union{Number,AffExpr}}
    m = length(ineq_cons)
    l = length(eq_cons)
    action = VariablePermutation(vars)
    monos_2d = MP.monomials(vars, 0:2*order)
    monos_d = MP.monomials(vars, 0:order)
    wedderburn = WedderburnDecomposition(Float64, group, action, monos_2d, monos_d)
    basis = Vector{Vector{Vector{Poly{Float64}}}}(undef, m+1)
    basis[1] = [[item'*monos_d for item in eachrow(ele.basis)] for ele in wedderburn.Uπs]
    for i = 1:m
        basis[i+1] = Vector{Poly{Float64}}[]
        if SymmetricConstraint == true
            for bas in basis[1]
                ind = maxdegree.(bas) .<= order - ceil(Int, maxdegree(ineq_cons[i])/2)
                if any(ind)
                    push!(basis[i+1], bas[ind])
                end
            end
        else
            ind = maxdegree.(monos_d) .<= order - ceil(Int, maxdegree(ineq_cons[i])/2)
            push!(basis[i+1], monos_d[ind])
        end
    end
    ebasis = Vector{Vector{Poly{Float64}}}(undef, l)
    if l > 0
        if SymmetricConstraint == true
            invbas = [minimum(monos_2d[item.nzind]) for item in wedderburn.invariants] # This is the basis for the invariant ring
            for i = 1:l
                ebasis[i] = invbas[maxdegree.(invbas) .<= 2*order-maxdegree(eq_cons[i])]
            end
        else
            for i = 1:l
                ebasis[i] = monos_2d[maxdegree.(monos_2d) .<= 2*order-maxdegree(eq_cons[i])]
            end
        end
    end
    tsupp = nothing
    if TS != false
        pop = [nonneg]
        if !isempty(ineq_cons)
            pop = [pop; ineq_cons]
        end
        if !isempty(eq_cons)
            pop = [pop; eq_cons]
        end
        tsupp = vcat([poly_norm(p, group, action)[1] for p in pop]...)
        unique!(tsupp)
        sort!(tsupp)
    end
    blocks,cl,blocksize,eblocks = get_blocks(m, l, tsupp, ineq_cons, eq_cons, basis, ebasis, group, action, TS=TS, SO=SO, merge=merge, md=md, QUIET=QUIET)
    poly = nonneg
    pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, 1+m)
    pos[1] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[1]))
    for (i, bas) in enumerate(basis[1])
        pos[1][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1][i])
        for (l, bs) in enumerate(blocksize[1][i])
            if bs == 1
               pos[1][i][l] = @variable(model, lower_bound=0)
               poly -= bas[blocks[1][i][l][1]]^2 * pos[1][i][l]
            else
               pos[1][i][l] = @variable(model, [1:bs, 1:bs], PSD)
               poly -= bas[blocks[1][i][l]]' * pos[1][i][l] * bas[blocks[1][i][l]]
            end
        end
    end
    for k = 1:m
        pos[k+1] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(basis[k+1]))
        for (i, bas) in enumerate(basis[k+1])
            pos[k+1][i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k+1][i])
            for (l, bs) in enumerate(blocksize[k+1][i])
                if bs == 1
                    pos[k+1][i][l] = @variable(model, lower_bound=0)
                    poly -= bas[blocks[k+1][i][l][1]]^2 * pos[k+1][i][l] * ineq_cons[k]
                else
                    pos[k+1][i][l] = @variable(model, [1:bs, 1:bs], PSD)
                    poly -= bas[blocks[k+1][i][l]]' * pos[k+1][i][l] * bas[blocks[k+1][i][l]] * ineq_cons[k]
                end
            end
        end
    end
    mul = nothing
    if l > 0
        mul = Vector{Vector{VariableRef}}(undef, l)
        for k = 1:l
            mul[k] = @variable(model, [1:length(eblocks[k])])
            poly -= ebasis[k][eblocks[k]]' * mul[k] * eq_cons[k]
        end
    end
    tsupp,coe = poly_norm(poly, group, action)
    @constraint(model, coe .== 0)
    info = poly_basis(nothing, l, group, action, basis, ebasis, tsupp, blocksize, blocks, eblocks, pos, mul, nothing)
    return info
end

function get_blocks(m::Int, l::Int, tsupp, ineq_cons, eq_cons, basis, ebasis, group, action; TS="block", SO=1, merge=false, md=3, QUIET=false)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, m+1)
    eblocks = Vector{Vector{Int}}(undef, l)
    blocksize = Vector{Vector{Vector{Int}}}(undef, m+1)
    cl = Vector{Vector{Int}}(undef, m+1)
    if TS == false
        for k = 1:m+1
            blocks[k],blocksize[k],cl[k] = [[Vector(1:length(basis[k][i]))] for i = 1:length(basis[k])],[[length(basis[k][i])] for i = 1:length(basis[k])],ones(Int, length(basis[k]))
        end
        for k = 1:l
            eblocks[k] = Vector(1:length(ebasis[k]))
        end
    else
        for i = 1:SO
            if i > 1
                oblocksize = deepcopy(blocksize)
                oeblocks = deepcopy(eblocks)
            end
            for k = 1:m+1
                blocks[k] = Vector{Vector{Vector{Int}}}(undef, length(basis[k]))
                blocksize[k] = Vector{Vector{Int}}(undef, length(basis[k]))
                cl[k] = Vector{Int}(undef, length(basis[k]))
                for (i, ba) in enumerate(basis[k])
                    if k == 1
                        G = get_graph(tsupp, ba, group, action)
                    else
                        G = get_graph(tsupp, ba, g=ineq_cons[k-1], group, action)
                    end
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
            for k = 1:l
                eblocks[k] = get_eblock(tsupp, eq_cons[k], ebasis[k], group, action)
            end
            if i > 1 && blocksize == oblocksize && eblocks == oeblocks
                println("No higher TS step of the TSSOS hierarchy!")
                break
            end
            if i < SO
                tsupp = DP.Monomial[]
                for t = 1:length(blocks[1]), s = 1:length(blocksize[1][t]), j = 1:blocksize[1][t][s], r = j:blocksize[1][t][s]  
                    append!(tsupp, poly_norm(basis[1][t][blocks[1][t][s][j]]*basis[1][t][blocks[1][t][s][r]], group, action)[1])
                end
                unique!(tsupp)
                sort!(tsupp)
            end
        end
    end
    return blocks,cl,blocksize,eblocks
end
