mutable struct mcpop_data
    n # number of all variables
    nb # number of binary variables
    m # number of all constraints
    numeq # number of equality constraints
    supp # support data
    coe # coefficient data
    basis # monomial bases
    ebasis # monomial bases for equality constraints
    ksupp # extended support at the k-th step
    cql # number of cliques
    cliquesize # sizes of cliques
    cliques # cliques of variables
    I # index sets of inequality constraints
    J # index sets of equality constraints
    ncc # constraints associated to no clique
    blocksize # sizes of blocks
    blocks # block structure
    eblocks # block structrue for equality constraints
    GramMat # Gram matrices
    multiplier # multipliers for equality constraints
    moment # Moment matrix
    solver # SDP solver
    SDP_status
    tol # tolerance to certify global optimality
    flag # 0 if global optimality is certified; 1 otherwise
end

"""
    opt,sol,data = cs_tssos_first(pop, x, d; nb=0, numeq=0, CS="MF", cliques=[], basis=[], ebasis=[], TS="block", merge=false, md=3, solver="Mosek", 
    dualize=false, QUIET=false, solve=true, solution=false, Gram=false, MomentOne=false, tol=1e-4, cosmo_setting=cosmo_para(), mosek_setting=mosek_para())

Compute the first TS step of the CS-TSSOS hierarchy for constrained polynomial optimization.
If `merge=true`, perform the PSD block merging. 
If `solve=false`, then do not solve the SDP.
If `Gram=true`, then output the Gram matrix.
If `MomentOne=true`, add an extra first-order moment PSD constraint to the moment relaxation.

# Input arguments
- `pop`: vector of the objective, inequality constraints, and equality constraints
- `x`: POP variables
- `d`: relaxation order
- `nb`: number of binary variables in `x`
- `numeq`: number of equality constraints
- `CS`: method of chordal extension for correlative sparsity (`"MF"`, `"MD"`, `"NC"`, `false`)
- `cliques`: the set of cliques used in correlative sparsity
- `TS`: type of term sparsity (`"block"`, `"signsymmetry"`, `"MD"`, `"MF"`, `false`)
- `md`: tunable parameter for merging blocks
- `QUIET`: run in the quiet mode (`true`, `false`)
- `tol`: relative tolerance to certify global optimality

# Output arguments
- `opt`: optimum
- `sol`: (near) optimal solution (if `solution=true`)
- `data`: other auxiliary data 
"""
function cs_tssos_first(pop::Vector{P}, x, d; nb=0, numeq=0, CS="MF", cliques=[], basis=[], ebasis=[], TS="block", merge=false, md=3, solver="Mosek", 
    dualize=false, QUIET=false, solve=true, solution=false, Gram=false, MomentOne=false, tol=1e-4, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false) where {P<:AbstractPolynomial}
    n,supp,coe = polys_info(pop, x, nb=nb)
    opt,sol,data = cs_tssos_first(supp, coe, n, d, numeq=numeq, nb=nb, CS=CS, cliques=cliques, basis=basis, ebasis=ebasis, TS=TS,
    merge=merge, md=md, QUIET=QUIET, solver=solver, dualize=dualize, solve=solve, solution=solution, Gram=Gram, MomentOne=MomentOne,
    tol=tol, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, writetofile=writetofile)
    return opt,sol,data
end

"""
    opt,sol,data = cs_tssos_first(supp::Vector{Vector{Vector{UInt16}}}, coe, n, d; nb=0, numeq=0, CS="MF", cliques=[], basis=[], ebasis=[], TS="block", 
    merge=false, md=3, QUIET=false, solver="Mosek", dualize=false, solve=true, solution=false, Gram=false, MomentOne=false, tol=1e-4, 
    cosmo_setting=cosmo_para(), mosek_setting=mosek_para())

Compute the first TS step of the CS-TSSOS hierarchy for constrained polynomial optimization. 
Here the polynomial optimization problem is defined by `supp` and `coe`, corresponding to the supports and coeffients of `pop` respectively.
"""
function cs_tssos_first(supp::Vector{Vector{Vector{UInt16}}}, coe, n, d; numeq=0, nb=0, CS="MF", cliques=[], basis=[], ebasis=[], TS="block", 
    merge=false, md=3, QUIET=false, solver="Mosek", dualize=false, solve=true, solution=false, MomentOne=false, Gram=false, 
    tol=1e-4, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false)
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    m = length(supp) - 1
    supp[1],coe[1] = resort(supp[1], coe[1])
    dc = [maximum(length.(supp[i])) for i=2:m+1]
    if cliques != []
        cql = length(cliques)
        cliquesize = length.(cliques)
    else
        time = @elapsed begin
        CS = CS == true ? "MF" : CS
        cliques,cql,cliquesize = clique_decomp(n, m, numeq, dc, supp, order=d, alg=CS)
        end
        if CS != false && QUIET == false
            mc = maximum(cliquesize)
            println("Obtained the variable cliques in $time seconds. The maximal size of cliques is $mc.")
        end
    end
    I,J,ncc = assign_constraint(m, numeq, supp, cliques, cql)
    if d == "min"
        rlorder = [isempty(I[i]) && isempty(J[i]) ? 1 : ceil(Int, maximum(dc[[I[i]; J[i]]])/2) for i = 1:cql]
    else
        rlorder = d*ones(Int, cql)
    end
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    if isempty(basis)
        basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        ebasis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        for i = 1:cql
            basis[i] = Vector{Vector{Vector{UInt16}}}(undef, length(I[i])+1)
            basis[i][1] = get_basis(cliques[i], rlorder[i], nb=nb)
            ebasis[i] = Vector{Vector{Vector{UInt16}}}(undef, length(J[i]))
            for s = 1:length(I[i])
                basis[i][s+1] = get_basis(cliques[i], rlorder[i]-ceil(Int, dc[I[i][s]]/2), nb=nb)
            end
            for s = 1:length(J[i])
                ebasis[i][s] = get_basis(cliques[i], 2*rlorder[i]-dc[J[i][s]], nb=nb)
            end
        end
    end
    ksupp = nothing
    if TS != false
        ksupp = reduce(vcat, supp)
        for k = 1:cql, item in basis[k][1]
            push!(ksupp, sadd(item, item, nb=nb))
        end
        sort!(ksupp)
        unique!(ksupp)
    end    
    time = @elapsed begin
    ss = nothing
    if TS == "signsymmetry"
        ss = get_signsymmetry(supp, n)
    end
    blocks,eblocks,cl,blocksize = get_blocks(I, J, supp, cliques, cql, ksupp, basis, ebasis, nb=nb, TS=TS, merge=merge, md=md, nvar=n, signsymmetry=ss)
    end
    if QUIET == false
        mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,momone,moment,GramMat,multiplier,SDP_status = solvesdp(m, supp, coe, basis, ebasis, cliques, cql, cliquesize, I, J, ncc, blocks, 
    eblocks, cl, blocksize, numeq=numeq, nb=nb, QUIET=QUIET, TS=TS, solver=solver, dualize=dualize, solve=solve, solution=solution, MomentOne=MomentOne, 
    Gram=Gram, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, writetofile=writetofile)
    data = mcpop_data(n, nb, m, numeq, supp, coe, basis, ebasis, ksupp, cql, cliquesize, cliques, I, J, ncc, blocksize, blocks, eblocks, GramMat, 
    multiplier, moment, solver, SDP_status, tol, 1)
    sol = nothing
    if solution == true
        sol,gap,data.flag = approx_sol(opt, momone, n, cliques, cql, cliquesize, supp, coe, numeq=numeq, tol=tol)
        if data.flag == 1
            sol = gap > 0.5 ? randn(n) : sol
            sol,data.flag = refine_sol(opt, sol, data, QUIET=true, tol=tol)
        end
    end
    return opt,sol,data
end

"""
    opt,sol,data = cs_tssos_higher!(data; TS="block", merge=false, md=3, QUIET=false, solve=true, solution=false, Gram=false, dualize=false, 
    MomentOne=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para())

Compute higher TS steps of the CS-TSSOS hierarchy.
"""
function cs_tssos_higher!(data::mcpop_data; TS="block", merge=false, md=3, QUIET=false, solve=true, solution=false, Gram=false, dualize=false, 
    MomentOne=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false)
    n = data.n
    nb = data.nb
    numeq = data.numeq
    supp = data.supp
    basis = data.basis
    ebasis = data.ebasis
    cql = data.cql
    cliques = data.cliques
    cliquesize = data.cliquesize
    I = data.I
    J = data.J
    tol = data.tol
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,eblocks,cl,blocksize = get_blocks(I, J, supp, cliques, cql, data.ksupp, basis, ebasis, nb=nb, TS=TS, merge=merge, md=md)
    end
    if blocksize == data.blocksize && eblocks == data.eblocks
        println("No higher TS step of the CS-TSSOS hierarchy!")
        opt = sol = nothing
    else
        if TS != false && QUIET == false
            mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,momone,moment,GramMat,multiplier,SDP_status = solvesdp(data.m, supp, data.coe, basis, ebasis, cliques, cql, cliquesize, I, J, data.ncc, blocks, eblocks, cl,
        blocksize, numeq=numeq, nb=nb, QUIET=QUIET, solver=data.solver, solve=solve, solution=solution, dualize=dualize, MomentOne=MomentOne, 
        Gram=Gram, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, writetofile=writetofile)
        sol = nothing
        if solution == true
            sol,gap,data.flag = approx_sol(opt, momone, n, cliques, cql, cliquesize, supp, data.coe, numeq=numeq, tol=tol)
            if data.flag == 1
                sol = gap > 0.5 ? randn(n) : sol
                sol,data.flag = refine_sol(opt, sol, data, QUIET=true, tol=tol)
            end
        end
        data.blocks = blocks
        data.eblocks = eblocks
        data.blocksize = blocksize
        data.ksupp = ksupp
        data.GramMat = GramMat
        data.multiplier = multiplier
        data.moment = moment
        data.SDP_status = SDP_status
    end
    return opt,sol,data
end

function solvesdp(m, supp::Vector{Vector{Vector{UInt16}}}, coe, basis, ebasis, cliques, cql, cliquesize, I, J, ncc, blocks, eblocks, cl, blocksize; 
    numeq=0, nb=0, QUIET=false, TS="block", solver="Mosek", solve=true, solution=false, Gram=false, MomentOne=false, cosmo_setting=cosmo_para(), 
    mosek_setting=mosek_para(), dualize=false, writetofile=false)
    tsupp = Vector{UInt16}[]
    for i = 1:cql, j = 1:cl[i][1], k = 1:blocksize[i][1][j], r = k:blocksize[i][1][j]
        @inbounds bi = sadd(basis[i][1][blocks[i][1][j][k]], basis[i][1][blocks[i][1][j][r]], nb=nb)
        push!(tsupp, bi)
    end
    if TS != false && TS != "signsymmetry"
        for i = 1:cql 
            for (j, w) in enumerate(I[i]), l = 1:cl[i][j+1], t = 1:blocksize[i][j+1][l], r = t:blocksize[i][j+1][l], s = 1:length(supp[w+1])
                ind1 = blocks[i][j+1][l][t]
                ind2 = blocks[i][j+1][l][r]
                @inbounds bi = sadd(sadd(basis[i][j+1][ind1], supp[w+1][s], nb=nb), basis[i][j+1][ind2], nb=nb)
                push!(tsupp, bi)
            end
            for (j, w) in enumerate(J[i]), k in eblocks[i][j], item in supp[w+1]
                @inbounds bi = sadd(ebasis[i][j][k], item, nb=nb)
                push!(tsupp, bi)
            end
        end
        for i ∈ ncc, j = 1:length(supp[i+1])
            push!(tsupp, supp[i+1][j])
        end
    end
    if (MomentOne == true || solution == true) && TS != false
        ksupp = copy(tsupp)
    end
    if (MomentOne == true || solution == true) && TS != false
        for i = 1:cql, j = 1:cliquesize[i]
            push!(tsupp, [cliques[i][j]])
            for k = j+1:cliquesize[i]
                push!(tsupp, [cliques[i][j]; cliques[i][k]])
            end
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    if (MomentOne == true || solution == true) && TS != false
        sort!(ksupp)
        unique!(ksupp)
    else
        ksupp = tsupp
    end
    objv = moment = momone = GramMat = multiplier = SDP_status = nothing
    if solve == true
        ltsupp = length(tsupp)
        if QUIET == false
            println("Assembling the SDP...")
            println("There are $ltsupp affine constraints.")
        end
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
            return nothing,nothing,nothing,nothing,nothing,nothing,nothing
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = [AffExpr(0) for i=1:ltsupp]
        pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, cql)
        for i = 1:cql
            if (MomentOne == true || solution == true) && TS != false
                bs = cliquesize[i]+1
                pos0 = @variable(model, [1:bs, 1:bs], PSD)
                for t = 1:bs, r = t:bs
                    if t == 1 && r == 1
                        bi = UInt16[]
                    elseif t == 1 && r > 1
                        bi = [cliques[i][r-1]]
                    else
                        bi = sadd(cliques[i][t-1], cliques[i][r-1], nb=nb)
                    end
                    Locb = bfind(tsupp, ltsupp, bi)
                    if t == r
                        @inbounds add_to_expression!(cons[Locb], pos0[t,r])
                    else
                        @inbounds add_to_expression!(cons[Locb], 2, pos0[t,r])
                    end
                end
            end
            pos[i] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, 1+length(I[i]))
            pos[i][1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i][1])
            for l = 1:cl[i][1]
                @inbounds bs = blocksize[i][1][l]
                if bs == 1
                    pos[i][1][l] = @variable(model, lower_bound=0)
                    @inbounds bi = sadd(basis[i][1][blocks[i][1][l][1]], basis[i][1][blocks[i][1][l][1]], nb=nb)
                    Locb = bfind(tsupp, ltsupp, bi)
                    @inbounds add_to_expression!(cons[Locb], pos[i][1][l])
                else
                    pos[i][1][l] = @variable(model, [1:bs, 1:bs], PSD)
                    for t = 1:bs, r = t:bs
                        @inbounds ind1 = blocks[i][1][l][t]
                        @inbounds ind2 = blocks[i][1][l][r]
                        @inbounds bi = sadd(basis[i][1][ind1], basis[i][1][ind2], nb=nb)
                        Locb = bfind(tsupp, ltsupp, bi)
                        if t == r
                            @inbounds add_to_expression!(cons[Locb], pos[i][1][l][t,r])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2, pos[i][1][l][t,r])
                        end
                    end
                end
            end
        end
        ## process inequality constraints
        for i = 1:cql, (j, w) in enumerate(I[i])
            pos[i][j+1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i][j+1])
            for l = 1:cl[i][j+1]
                bs = blocksize[i][j+1][l]
                if bs == 1
                    pos[i][j+1][l] = @variable(model, lower_bound=0)
                    ind = blocks[i][j+1][l][1]
                    for s = 1:length(supp[w+1])
                        @inbounds bi = sadd(sadd(basis[i][j+1][ind], supp[w+1][s], nb=nb), basis[i][j+1][ind], nb=nb)
                        Locb = bfind(tsupp, ltsupp, bi)
                        @inbounds add_to_expression!(cons[Locb], coe[w+1][s], pos[i][j+1][l])
                    end
                else
                    pos[i][j+1][l] = @variable(model, [1:bs, 1:bs], PSD)
                    for t = 1:bs, r = t:bs
                        ind1 = blocks[i][j+1][l][t]
                        ind2 = blocks[i][j+1][l][r]
                        for s = 1:length(supp[w+1])
                            @inbounds bi = sadd(sadd(basis[i][j+1][ind1], supp[w+1][s], nb=nb), basis[i][j+1][ind2], nb=nb)
                            Locb = bfind(tsupp, ltsupp, bi)
                            if t == r
                                @inbounds add_to_expression!(cons[Locb], coe[w+1][s], pos[i][j+1][l][t,r])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*coe[w+1][s], pos[i][j+1][l][t,r])
                            end
                        end
                    end
                end
            end
        end
        ## process equality constraints
        if numeq > 0
            free = Vector{Vector{Vector{VariableRef}}}(undef, cql)
            for i = 1:cql
                if !isempty(J[i])
                    free[i] = Vector{Vector{VariableRef}}(undef, length(J[i]))
                    for (j, w) in enumerate(J[i])
                        free[i][j] = @variable(model, [1:length(eblocks[i][j])])
                        for (u,k) in enumerate(eblocks[i][j]), s = 1:length(supp[w+1])
                            @inbounds bi = sadd(ebasis[i][j][k], supp[w+1][s], nb=nb)
                            Locb = bfind(tsupp, ltsupp, bi)
                            @inbounds add_to_expression!(cons[Locb], coe[w+1][s], free[i][j][u])
                        end
                    end
                end
            end
        end
        for i in ncc
            if i <= m-numeq
                pos0 = @variable(model, lower_bound=0)
            else
                pos0 = @variable(model)
            end
            for j = 1:length(supp[i+1])
                Locb = bfind(tsupp, ltsupp, supp[i+1][j])
                @inbounds add_to_expression!(cons[Locb], coe[i+1][j], pos0)
            end
        end
        bc = zeros(ltsupp)
        for i = 1:length(supp[1])
            Locb = bfind(tsupp, ltsupp, supp[1][i])
            if Locb === nothing
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing,nothing,nothing,nothing,nothing
            else
               bc[Locb] = coe[1][i]
            end
        end
        @variable(model, lower)
        cons[1] += lower
        @constraint(model, con, cons==bc)
        @objective(model, Max, lower)
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
        if writetofile != false
            write_to_file(dualize(model), writetofile)
        end
        SDP_status = termination_status(model)
        objv = objective_value(model)
        if SDP_status != MOI.OPTIMAL
           println("termination status: $SDP_status")
           status = primal_status(model)
           println("solution status: $status")
        end
        println("optimum = $objv")
        if Gram == true
            GramMat = Vector{Vector{Vector{Union{Float64,Matrix{Float64}}}}}(undef, cql)
            for i = 1:cql
                GramMat[i] = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, 1+length(I[i]))
                for j = 1:1+length(I[i])
                    GramMat[i][j] = [value.(pos[i][j][l]) for l = 1:cl[i][j]]
                end
            end
            multiplier = Vector{Vector{Vector{Float64}}}(undef, cql)
            for i = 1:cql
                if !isempty(J[i])
                    multiplier[i] = [value.(free[i][j]) for j = 1:length(J[i])]
                end
            end
        end
        measure = -dual(con)
        momone = nothing
        if solution == true  
            momone = get_moment(measure, tsupp, cliques, cql, cliquesize, nb=nb)
        end
        moment = get_moment(measure, tsupp, cliques, cql, cliquesize, basis=basis, nb=nb)
    end
    return objv,ksupp,momone,moment,GramMat,multiplier,SDP_status
end

function get_eblock(tsupp::Vector{Vector{UInt16}}, hsupp::Vector{Vector{UInt16}}, basis::Vector{Vector{UInt16}}; nb=nb, nvar=0, signsymmetry=nothing)
    ltsupp = length(tsupp)
    hlt = length(hsupp)
    eblock = Int[]
    for (i,item) in enumerate(basis)
        if signsymmetry === nothing
            if findfirst(x -> bfind(tsupp, ltsupp, sadd(item, hsupp[x], nb=nb)) !== nothing, 1:hlt) !== nothing
                push!(eblock, i)
            end
        else
            bi = sadd(item, hsupp[1], nb=nb)
            sp = zeros(Int, nvar)
            st = sign_type(bi)
            sp[st] = ones(Int, length(st))
            if all(transpose(signsymmetry)*sp .== 0)
                push!(eblock, i)
            end
        end
    end
    return eblock
end

function get_blocks(I, J, supp::Vector{Vector{Vector{UInt16}}}, cliques, cql, tsupp, basis, ebasis; blocks=[], eblocks=[], cl=[], blocksize=[], TS="block",
    nb=0, merge=false, md=3, nvar=0, signsymmetry=nothing)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    for i = 1:cql
        ksupp = TS == false ? nothing : tsupp[[issubset(tsupp[j], cliques[i]) for j in eachindex(tsupp)]]
        blocks[i],eblocks[i],cl[i],blocksize[i] = get_blocks(length(I[i]), length(J[i]), ksupp, supp[[I[i]; J[i]].+1], basis[i],
        ebasis[i], TS=TS, nb=nb, QUIET=true, merge=merge, md=md, nvar=nvar, signsymmetry=signsymmetry)
    end
    return blocks,eblocks,cl,blocksize
end

function assign_constraint(m, numeq, supp::Vector{Vector{Vector{UInt16}}}, cliques, cql)
    I = [Int[] for i=1:cql]
    J = [Int[] for i=1:cql]
    ncc = Int[]
    for i = 1:m
        ind = findall(k->issubset(unique(reduce(vcat, supp[i+1])), cliques[k]), 1:cql)
        if isempty(ind)
            push!(ncc, i)
        elseif i <= m - numeq
            push!.(I[ind], i)
        else
            push!.(J[ind], i)
        end
    end
    return I,J,ncc
end

function get_graph(tsupp::Vector{Vector{UInt16}}, basis::Vector{Vector{UInt16}}; nb=0, nvar=0, signsymmetry=nothing)
    lb = length(basis)
    G = SimpleGraph(lb)
    ltsupp = length(tsupp)
    for i = 1:lb, j = i+1:lb
        bi = sadd(basis[i], basis[j], nb=nb)
        if signsymmetry === nothing
            if bfind(tsupp, ltsupp, bi) !== nothing
                add_edge!(G, i, j)
            end
        else
            sp = zeros(Int, nvar)
            st = sign_type(bi)
            sp[st] = ones(Int, length(st))
            if all(transpose(signsymmetry)*sp .== 0)
                add_edge!(G, i, j)
            end
        end
    end
    return G
end

function get_graph(tsupp::Vector{Vector{UInt16}}, supp::Vector{Vector{UInt16}}, basis::Vector{Vector{UInt16}}; nb=0, nvar=0, signsymmetry=nothing)
    lb = length(basis)
    ltsupp = length(tsupp)
    G = SimpleGraph(lb)
    for i = 1:lb, j = i+1:lb
        if signsymmetry === nothing
            ind = findfirst(x -> bfind(tsupp, ltsupp, sadd(sadd(basis[i], x, nb=nb), basis[j], nb=nb)) !== nothing, supp)
            if ind !== nothing
                add_edge!(G, i, j)
            end
        else
            bi = sadd(sadd(basis[i], supp[1], nb=nb), basis[j], nb=nb)
            sp = zeros(Int, nvar)
            st = sign_type(bi)
            sp[st] = ones(Int, length(st))
            if all(transpose(signsymmetry)*sp .== 0)
                add_edge!(G, i, j)
            end
        end
    end
    return G
end

function clique_decomp(n, m, numeq, dc, supp::Vector{Vector{Vector{UInt16}}}; order="min", alg="MF")
    if alg == false
        cliques,cql,cliquesize = [UInt16[i for i=1:n]],1,[n]
    else
        G = SimpleGraph(n)
        for i = 1:m+1
            if order == "min" || i == 1 || (order == ceil(Int, dc[i-1]/2) && i <= m-numeq+1) || (2*order == dc[i-1] && i > m-numeq+1)
                foreach(x -> add_clique!(G, unique(x)), supp[i])
            else
                add_clique!(G, unique(reduce(vcat, supp[i])))
            end
        end
        if alg == "NC"
            cliques,cql,cliquesize = max_cliques(G)
        else
            cliques,cql,cliquesize = chordal_cliques!(G, method=alg, minimize=true)
        end
    end
    uc = unique(cliquesize)
    sizes = [sum(cliquesize.== i) for i in uc]
    println("-----------------------------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("-----------------------------------------------------------------------------")
    return cliques,cql,cliquesize
end

function get_moment(measure, tsupp, cliques, cql, cliquesize; basis=[], nb=0)
    moment = Vector{Union{Float64, Symmetric{Float64}, Array{Float64,2}}}(undef, cql)
    ltsupp = length(tsupp)
    for i = 1:cql
        lb = isempty(basis) ? cliquesize[i] + 1 : length(basis[i][1])
        moment[i] = zeros(Float64, lb, lb)
        if basis == []
            for j = 1:lb, k = j:lb
                if j == 1
                    bi = k == 1 ? UInt16[] : [cliques[i][k-1]]
                else
                    bi = sadd(cliques[i][j-1], cliques[i][k-1], nb=nb)
                end
                Locb = bfind(tsupp, ltsupp, bi)
                moment[i][j,k] = measure[Locb]
            end
        else
            for j = 1:lb, k = j:lb
                bi = sadd(basis[i][1][j], basis[i][1][k], nb=nb)
                Locb = bfind(tsupp, ltsupp, bi)
                if Locb !== nothing
                    moment[i][j,k] = measure[Locb]
                end
            end
        end
        moment[i] = Symmetric(moment[i],:U)
    end
    return moment
end
