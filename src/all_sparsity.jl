mutable struct spop_data
    pop # polynomial optimiztion problem
    obj # objective
    ineq_cons # inequality constraints
    eq_cons # equality constraints
    x # variables
    n # number of variables
    nb # number of binary variables
    numeq # number of equality constraints
    basis # monomial bases
    ebasis # monomial bases for equality constraints
    ksupp # extended support at the k-th step
    cliquesize # sizes of cliques
    cliques # cliques of variables
    I # index sets of inequality constraints
    J # index sets of equality constraints
    Iprime # global inequality constraints
    Jprime # global equality constraints
    blocksize # sizes of blocks
    blocks # block structure
    eblocks # block structrue for equality constraints
    GramMat # Gram matrices
    multiplier # multipliers for equality constraints
    moment # Moment matrix
    SDP_status
    rtol # tolerance for rank
    gtol # tolerance for global optimality gap
    ftol # tolerance for feasibility
    flag # 0 if global optimality is certified; 1 otherwise
end

"""
    show_blocks(data)

Display the block structure
"""

function show_blocks(data::spop_data; include_constraints=false)
    for i = 1:data.cql, (j, block) in enumerate(data.blocks[i][1])
        print("clique $i, block $j: ")
        println([prod(data.x[data.basis[i][1][ind]]) for ind in block])
    end
    if include_constraints == true
        for i = 1:data.cql, (l, s) in enumerate(data.I[i]), (j, block) in enumerate(data.blocks[i][l+1])
            print("clique $i, constraint $s, block $j: ")
            println([prod(data.x[data.basis[i][l+1][ind]]) for ind in block])
        end
    end
end

"""
    opt,sol,data = cs_tssos(pop, x, d; nb=0, numeq=0, CS="MF", cliques=[], basis=[], ebasis=[], TS="block", eqTS=TS, merge=false, md=3, 
    dualize=false, QUIET=false, solve=true, solution=false, Gram=false, MomentOne=false, mosek_setting=mosek_para(), model=nothing, 
    rtol=1e-2, gtol=1e-2, ftol=1e-3)

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
- `eqTS`: type of term sparsity for equality constraints (by default the same as `TS`, `false`)
- `md`: tunable parameter for merging blocks
- `QUIET`: run in the quiet mode (`true`, `false`)
- `rtol`: tolerance for rank
- `gtol`: tolerance for global optimality gap
- `ftol`: tolerance for feasibility

# Output arguments
- `opt`: optimum
- `sol`: (near) optimal solution (if `solution=true`)
- `data`: other auxiliary data 
"""
function cs_tssos(pop::Vector{Poly{T}}, x, d; nb=0, numeq=0, CS="MF", cliques=[], basis=[], ebasis=[], TS="block", eqTS=TS, merge=false, md=3,
    dualize=false, QUIET=false, solve=true, solution=false, Gram=false, MomentOne=false, mosek_setting=mosek_para(), model=nothing, 
    writetofile=false, rtol=1e-2, gtol=1e-2, ftol=1e-3) where {T<:Number}
    if nb > 0
        pop = Groebner.normalform(x[1:nb].^2 .- 1, pop)
    end
    npop = [poly(p, x) for p in pop]
    opt,sol,data = cs_tssos(npop, length(x), d, numeq=numeq, nb=nb, CS=CS, cliques=cliques, basis=basis, ebasis=ebasis, TS=TS,
    eqTS=eqTS, merge=merge, md=md, QUIET=QUIET, dualize=dualize, solve=solve, solution=solution, Gram=Gram, MomentOne=MomentOne,
    mosek_setting=mosek_setting, model=model, writetofile=writetofile, rtol=rtol, gtol=gtol, ftol=ftol, pop=pop, x=x)
    return opt,sol,data
end

"""
    opt,sol,data = cs_tssos(npop::Vector{poly{T}}, n, d; nb=0, numeq=0, CS="MF", cliques=[], basis=[], ebasis=[], TS="block", 
    eqTS=TS, merge=false, md=3, QUIET=false, dualize=false, solve=true, solution=false, Gram=false, MomentOne=false, mosek_setting=mosek_para(), 
    model=nothing, rtol=1e-2, gtol=1e-2, ftol=1e-3) where {T<:Number}

Compute the first TS step of the CS-TSSOS hierarchy for constrained polynomial optimization. 
"""
function cs_tssos(npop::Vector{poly{T}}, n, d; numeq=0, nb=0, CS="MF", cliques=[], basis=[], ebasis=[], TS="block", 
    eqTS=TS, merge=false, md=3, QUIET=false, dualize=false, solve=true, solution=false, MomentOne=false, Gram=false, 
    mosek_setting=mosek_para(), model=nothing, writetofile=false, rtol=1e-2, gtol=1e-2, ftol=1e-3, pop=nothing, x=nothing) where {T<:Number}
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    obj = npop[1]
    if pop === nothing
        obj = arrange(obj)
    end
    ineq_cons = [poly{T}([UInt16[]], [1]); npop[2:end-numeq]]
    eq_cons = npop[end-numeq+1:end]
    if !isempty(cliques)
        cql = length(cliques)
        cliquesize = length.(cliques)
    else
        time = @elapsed begin
        CS = CS == true ? "MF" : CS
        cliques,cql,cliquesize = clique_decomp(npop, n, numeq, order=d, alg=CS)
        end
        if CS != false && QUIET == false
            mc = maximum(cliquesize)
            println("Obtained the variable cliques in $time seconds. The maximal size of cliques is $mc.")
        end
    end
    I,J,Iprime,Jprime = assign_constraint(ineq_cons, eq_cons, cliques, cql)
    if d == "min"
        rlorder = [isempty(I[i]) && isempty(J[i]) ? 1 : ceil(Int, maximum([maxdeg.(ineq_cons[I[i]]); maxdeg.(eq_cons[J[i]])])/2) for i = 1:cql]
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
            basis[i] = Vector{Vector{Vector{UInt16}}}(undef, length(I[i]))
            ebasis[i] = Vector{Vector{Vector{UInt16}}}(undef, length(J[i]))
            for (s, k) in enumerate(I[i])
                basis[i][s] = get_basis(cliques[i], rlorder[i]-ceil(Int, maxdeg(ineq_cons[k])/2), nb=nb)
            end
            for (s, k) in enumerate(J[i])
                ebasis[i][s] = get_basis(cliques[i], 2*rlorder[i]-maxdeg(eq_cons[k]), nb=nb)
            end
        end
    end
    ksupp = nothing
    if TS != false && TS != "signsymmetry"
        ksupp = reduce(vcat, [p.supp for p in npop])
        for k = 1:cql, item in basis[k][1]
            push!(ksupp, sadd(item, item, nb=nb))
        end
        sort!(ksupp)
        unique!(ksupp)
    end    
    time = @elapsed begin
    ss = nothing
    if TS == "signsymmetry" || eqTS == "signsymmetry"
        ss = get_signsymmetry(npop, n)
    end
    blocks,cl,blocksize,eblocks = get_blocks(ineq_cons, eq_cons, I, J, cliques, cql, ksupp, basis, ebasis, nb=nb, TS=TS, eqTS=eqTS, merge=merge, md=md, signsymmetry=ss)
    end
    if QUIET == false
        mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,momone,moment,GramMat,multiplier,SDP_status = solvesdp(obj, ineq_cons, eq_cons, basis, ebasis, cliques, cql, cliquesize, I, J, Iprime, Jprime, blocks, 
    eblocks, cl, blocksize, nb=nb, QUIET=QUIET, TS=TS, dualize=dualize, solve=solve, solution=solution, MomentOne=MomentOne, Gram=Gram, mosek_setting=mosek_setting, 
    model=model, writetofile=writetofile)
    data = spop_data(pop, obj, ineq_cons, eq_cons, x, n, nb, numeq, basis, ebasis, ksupp, cliquesize, cliques, I, J, Iprime, Jprime, blocksize, blocks, eblocks, GramMat, 
    multiplier, moment, SDP_status, rtol, gtol, ftol, 1)
    sol = nothing
    if solution == true
        if TS != false
            sol,gap,data.flag = approx_sol(momone, opt, n, cliques, cql, cliquesize, npop, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=true)
            if data.flag == 1
                if gap > 1
                    rsol,status,data.flag = refine_sol(opt, randn(n), data, QUIET=true, gtol=gtol)
                    if status == MOI.LOCALLY_SOLVED
                        sol = rsol
                    end
                else
                    sol,_,data.flag = refine_sol(opt, sol, data, QUIET=true, gtol=gtol)
                end
            end
        else
            sol = extract_solutions_robust(moment, n, d, cliques, cql, cliquesize, pop=pop, x=x, npop=npop, lb=opt, numeq=numeq, check=true, rtol=rtol, gtol=gtol, ftol=ftol, QUIET=QUIET)[1]
            if sol !== nothing
                data.flag = 0
            end
        end
    end
    return opt,sol,data
end

"""
    opt,sol,data = cs_tssos(data; TS="block", eqTS=TS, merge=false, md=3, QUIET=false, solve=true, solution=false, Gram=false, dualize=false, 
    MomentOne=false, mosek_setting=mosek_para(), model=nothing)

Compute higher TS steps of the CS-TSSOS hierarchy.
"""
function cs_tssos(data::spop_data; TS="block", eqTS=TS, merge=false, md=3, QUIET=false, solve=true, solution=false, Gram=false, dualize=false, 
    MomentOne=false, mosek_setting=mosek_para(), model=nothing, writetofile=false)
    obj = data.obj
    ineq_cons = data.ineq_cons
    eq_cons = data.eq_cons
    n = data.n
    nb = data.nb
    numeq = data.numeq
    basis = data.basis
    ebasis = data.ebasis
    cliques = data.cliques
    cliquesize = data.cliquesize
    cql = length(cliques)
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,eblocks = get_blocks(ineq_cons, eq_cons, data.I, data.J, cliques, cql, data.ksupp, basis, ebasis, nb=nb, TS=TS, eqTS=eqTS, merge=merge, md=md)
    end
    if blocksize == data.blocksize && eblocks == data.eblocks
        println("No higher TS step of the CS-TSSOS hierarchy!")
        opt = sol = nothing
    else
        if QUIET == false
            mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,momone,moment,GramMat,multiplier,SDP_status = solvesdp(obj, ineq_cons, eq_cons, basis, ebasis, cliques, cql, cliquesize, data.I, data.J, 
        data.Iprime, data.Jprime, blocks, eblocks, cl, blocksize, nb=nb, QUIET=QUIET, solve=solve, solution=solution, dualize=dualize, MomentOne=MomentOne, 
        Gram=Gram, mosek_setting=mosek_setting, model=model, writetofile=writetofile)
        sol = nothing
        if solution == true
            sol,gap,data.flag = approx_sol(momone, opt, n, cliques, cql, cliquesize, [obj; ineq_cons[2:end]; eq_cons], numeq=numeq, gtol=data.gtol, ftol=data.ftol, QUIET=true)
            if data.flag == 1
                if gap > 1
                    rsol,status,data.flag = refine_sol(opt, randn(n), data, QUIET=true, gtol=data.gtol)
                    if status == MOI.LOCALLY_SOLVED
                        sol = rsol
                    end
                else
                    sol,_,data.flag = refine_sol(opt, sol, data, QUIET=true, gtol=data.gtol)
                end
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

function solvesdp(obj::poly{T}, ineq_cons::Vector{poly{T}}, eq_cons::Vector{poly{T}}, basis, ebasis, cliques, cql, cliquesize, I, J, Iprime, Jprime, 
    blocks, eblocks, cl, blocksize; nb=0, QUIET=false, TS="block", solve=true, solution=false, Gram=false, MomentOne=false, mosek_setting=mosek_para(), 
    model=nothing, dualize=false, writetofile=false) where {T<:Number}
    tsupp = Vector{UInt16}[]
    for i = 1:cql, j = 1:cl[i][1], k = 1:blocksize[i][1][j], r = k:blocksize[i][1][j]
        @inbounds bi = sadd(basis[i][1][blocks[i][1][j][k]], basis[i][1][blocks[i][1][j][r]], nb=nb)
        push!(tsupp, bi)
    end
    if TS != false && TS != "signsymmetry"
        for i = 1:cql
            for (j, w) in enumerate(I[i][2:end]), l = 1:cl[i][j+1], t = 1:blocksize[i][j+1][l], r = t:blocksize[i][j+1][l], item in ineq_cons[w].supp
                @inbounds bi = sadd(basis[i][j+1][blocks[i][j+1][l][t]], item, basis[i][j+1][blocks[i][j+1][l][r]], nb=nb)
                push!(tsupp, bi)
            end
            for (j, w) in enumerate(J[i]), k in eblocks[i][j], item in eq_cons[w].supp
                @inbounds bi = sadd(ebasis[i][j][k], item, nb=nb)
                push!(tsupp, bi)
            end
        end
        for i in Iprime
            append!(tsupp, ineq_cons[i].supp)
        end
        for i in Jprime
            append!(tsupp, eq_cons[i].supp)
        end
    end
    if (MomentOne == true || solution == true) && TS != false
        ksupp = copy(tsupp)
        for i = 1:cql
            append!(tsupp, get_basis(cliques[i], 2, nb=nb))
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
        if QUIET == false
            println("Assembling the SDP...")
            println("There are $(length(tsupp)) affine constraints.")
        end
        if model === nothing
            if dualize == false
                model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
            else
                model = Model(dual_optimizer(Mosek.Optimizer))
            end
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = [AffExpr(0) for i=1:length(tsupp)]
        pos = Vector{Vector{Vector{Symmetric{VariableRef}}}}(undef, cql)
        for i = 1:cql
            if (MomentOne == true || solution == true) && TS != false
                bas = [[UInt16[]]; [UInt16[k] for k in cliques[i]]]
                pos0 = @variable(model, [1:length(bas), 1:length(bas)], PSD)
                for t = 1:length(bas), r = t:length(bas)
                    Locb = bfind(tsupp, sadd(bas[t], bas[r], nb=nb))
                    if t == r
                        @inbounds add_to_expression!(cons[Locb], pos0[t,r])
                    else
                        @inbounds add_to_expression!(cons[Locb], 2, pos0[t,r])
                    end
                end
            end
            pos[i] = Vector{Vector{Symmetric{VariableRef}}}(undef, length(I[i]))
            for (j, w) in enumerate(I[i])
                pos[i][j] = Vector{Symmetric{VariableRef}}(undef, cl[i][j])
                for l = 1:cl[i][j]
                    pos[i][j][l] = @variable(model, [1:blocksize[i][j][l], 1:blocksize[i][j][l]], PSD)
                    for t = 1:blocksize[i][j][l], r = t:blocksize[i][j][l], (s, item) in enumerate(ineq_cons[w].supp)
                        @inbounds bi = sadd(basis[i][j][blocks[i][j][l][t]], item, basis[i][j][blocks[i][j][l][r]], nb=nb)
                        Locb = bfind(tsupp, bi)
                        if t == r
                            @inbounds add_to_expression!(cons[Locb], ineq_cons[w].coe[s], pos[i][j][l][t,r])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2*ineq_cons[w].coe[s], pos[i][j][l][t,r])
                        end
                    end
                end
            end
        end
        free = Vector{Vector{Vector{VariableRef}}}(undef, cql)
        for i = 1:cql
            free[i] = Vector{Vector{VariableRef}}(undef, length(J[i]))
            for (j, w) in enumerate(J[i])
                free[i][j] = @variable(model, [1:length(eblocks[i][j])])
                for (u, k) in enumerate(eblocks[i][j]), (s, item) in enumerate(eq_cons[w].supp)
                    @inbounds bi = sadd(ebasis[i][j][k], item, nb=nb)
                    Locb = bfind(tsupp, bi)
                    @inbounds add_to_expression!(cons[Locb], eq_cons[w].coe[s], free[i][j][u])
                end
            end
        end
        for i in Iprime
            pos0 = @variable(model, lower_bound=0)
            for (j, item) in enumerate(ineq_cons[i].supp)
                Locb = bfind(tsupp, item)
                @inbounds add_to_expression!(cons[Locb], ineq_cons[i].coe[j], pos0)
            end
        end
        for i in Jprime
            pos0 = @variable(model)
            for (j, item) in enumerate(eq_cons[i].supp)
                Locb = bfind(tsupp, item)
                @inbounds add_to_expression!(cons[Locb], eq_cons[i].coe[j], pos0)
            end
        end
        for (i, item) in enumerate(obj.supp)
            Locb = bfind(tsupp, item)
            if Locb === nothing
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing,nothing,nothing,nothing,nothing
            else
               cons[Locb] -= obj.coe[i]
            end
        end
        @variable(model, lower)
        cons[1] += lower
        @constraint(model, con, cons==zeros(length(cons)))
        @objective(model, Max, lower)
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
        end
        if writetofile != false
            write_to_file(dualize(model), writetofile)
        else
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
            objv = objective_value(model)
            if SDP_status != MOI.OPTIMAL
                println("termination status: $SDP_status")
                status = primal_status(model)
                println("solution status: $status")
            end
            println("optimum = $objv")
            if Gram == true
                GramMat = [[[value.(pos[i][j][l]) for l = 1:cl[i][j]] for j = 1:length(I[i])] for i = 1:cql]
                multiplier = [[value.(free[i][j]) for j = 1:length(J[i])] for i = 1:cql]
            end
            dual_var = -dual(con)
            if solution == true && TS != false
                momone = get_moment(dual_var, tsupp, cliques, cql, cliquesize, nb=nb)
            end
            moment = get_moment(dual_var, tsupp, cliques, cql, cliquesize, basis=basis, nb=nb)
        end
    end
    return objv,ksupp,momone,moment,GramMat,multiplier,SDP_status
end

function get_blocks(ineq_cons::Vector{poly{T}}, eq_cons::Vector{poly{T}}, I, J, cliques, cql, tsupp, basis, ebasis; TS="block", eqTS=TS, nb=0, merge=false, md=3, signsymmetry=nothing) where {T<:Number}
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
    for i = 1:cql
        ksupp = (TS == false || TS == "signsymmetry") ? nothing : tsupp[[issubset(item, cliques[i]) for item in tsupp]]
        blocks[i],cl[i],blocksize[i],eblocks[i] = get_blocks(ineq_cons[I[i]], eq_cons[J[i]], ksupp, basis[i],
        ebasis[i], TS=TS, eqTS=eqTS, nb=nb, QUIET=true, merge=merge, md=md, signsymmetry=signsymmetry)
    end
    return blocks,cl,blocksize,eblocks
end

function clique_decomp(npop::Vector{T}, n, numeq; order="min", alg="MF", QUIET=false) where {T<:poly}
    if alg == false
        cliques,cql,cliquesize = [Vector(1:n)],1,[n]
    else
        G = SimpleGraph(n)
        for (i, p) in enumerate(npop)
            if order == "min" || i == 1 || (order == ceil(Int, maxdeg(p)/2) && i <= length(npop)-numeq) || (2*order == maxdeg(p) && i > length(npop)-numeq)
                foreach(item -> add_clique!(G, unique(item)), p.supp)
            else
                add_clique!(G, unique(reduce(vcat, p.supp)))
            end
        end
        if alg == "NC"
            cliques,cql,cliquesize = max_cliques(G)
        else
            cliques,cql,cliquesize = chordal_cliques!(G, method=alg, minimize=true)
        end
    end
    if QUIET == false
        uc = unique(cliquesize)
        sizes = [sum(cliquesize.== i) for i in uc]
        println("-----------------------------------------------------------------------------")
        println("The clique sizes of varibles:\n$uc\n$sizes")
        println("-----------------------------------------------------------------------------")
    end
    return cliques,cql,cliquesize
end

function assign_constraint(ineq_cons::Vector{T1}, eq_cons::Vector{T2}, cliques, cql) where {T1,T2<:poly}
    I = [Int[] for i=1:cql]
    J = [Int[] for i=1:cql]
    Iprime = Int[]
    Jprime = Int[]
    for (i, p) in enumerate(ineq_cons)
        ind = findall(k->issubset(unique(reduce(vcat, p.supp)), cliques[k]), 1:cql)
        if isempty(ind)
            push!(Iprime, i)
        else
            push!.(I[ind], i)
        end
    end
    for (i, p) in enumerate(eq_cons)
        ind = findall(k->issubset(unique(reduce(vcat, p.supp)), cliques[k]), 1:cql)
        if isempty(ind)
            push!(Jprime, i)
        else
            push!.(J[ind], i)
        end
    end
    return I,J,Iprime,Jprime
end

function get_moment(dual_var, tsupp, cliques, cql, cliquesize; basis=[], nb=0)
    moment = Vector{Symmetric{Float64}}(undef, cql)
    for i = 1:cql
        lb = isempty(basis) ? cliquesize[i] + 1 : length(basis[i][1])
        mmat = zeros(Float64, lb, lb)
        if isempty(basis)
            for j = 1:lb, k = j:lb
                if j == 1
                    bi = k == 1 ? UInt16[] : UInt16[cliques[i][k-1]]
                else
                    bi = sadd(UInt16[cliques[i][j-1]], UInt16[cliques[i][k-1]], nb=nb)
                end
                Locb = bfind(tsupp, bi)
                mmat[j,k] = dual_var[Locb]
            end
        else
            for j = 1:lb, k = j:lb
                bi = sadd(basis[i][1][j], basis[i][1][k], nb=nb)
                Locb = bfind(tsupp, bi)
                if Locb !== nothing
                    mmat[j,k] = dual_var[Locb]
                end
            end
        end
        moment[i] = Symmetric(mmat, :U)
    end
    return moment
end
