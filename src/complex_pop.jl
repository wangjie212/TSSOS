mutable struct ccpop_data
    pop # complex polynomial optimiztion problem
    obj # objective
    ineq_cons # inequality constraints
    eq_cons # equality constraints
    z # complex variables
    rlorder # relaxation order
    n # number of complex variables
    nb # number of unit-norm variables
    numeq # number of equality constraints
    ipart # include the imaginary part
    ConjugateBasis # include conjugate variables in monomial bases
    normality # normal order
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
    solver # SDP solver
    SDP_status
    rtol # tolerance for rank
    gtol # tolerance for global optimality gap
    ftol # tolerance for feasibility
    flag # 0 if global optimality is certified; 1 otherwise
end

"""
    opt,sol,data = complex_cs_tssos_first(pop, z, d; nb=0, numeq=0, CS="MF", cliques=[], TS="block", reducebasis=false, 
    merge=false, md=3, solver="Mosek", QUIET=false, solve=true, solution=false, dualize=false, Gram=false, 
    MomentOne=false, ConjugateBasis=false, normality=!ConjugateBasis, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), 
    rtol=1e-2, gtol=1e-2, ftol=1e-3)

Compute the first TS step of the CS-TSSOS hierarchy for constrained complex polynomial optimization. 
If `ConjugateBasis=true`, then include conjugate variables in monomial bases.
If `merge=true`, perform the PSD block merging. 
If `solve=false`, then do not solve the SDP.
If `Gram=true`, then output the Gram matrix.
If `MomentOne=true`, add an extra first-order moment PSD constraint to the moment relaxation.

# Input arguments
- `pop`: vector of the objective, inequality constraints, and equality constraints
- `z`: CPOP variables and their conjugate
- `d`: relaxation order
- `nb`: number of unit-norm variables in `x`
- `numeq`: number of equality constraints
- `CS`: method of chordal extension for correlative sparsity (`"MF"`, `"MD"`, `false`)
- `cliques`: the set of cliques used in correlative sparsity
- `TS`: type of term sparsity (`"block"`, `"MD"`, `"MF"`, `false`)
- `md`: tunable parameter for merging blocks
- `normality`: normal order
- `QUIET`: run in the quiet mode (`true`, `false`)
- `rtol`: tolerance for rank
- `gtol`: tolerance for global optimality gap
- `ftol`: tolerance for feasibility

# Output arguments
- `opt`: optimum
- `sol`: (near) optimal solution (if `solution=true`)
- `data`: other auxiliary data 
"""
function complex_cs_tssos_first(pop::Vector{Poly{T}}, z, d; numeq=0, RemSig=false, nb=0, CS="MF", cliques=[], TS="block", 
    merge=false, md=3, solver="Mosek", reducebasis=false, QUIET=false, solve=true, solution=false, dualize=false, MomentOne=false, 
    ConjugateBasis=false, Gram=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false, normality=!ConjugateBasis, 
    rtol=1e-2, gtol=1e-2, ftol=1e-3) where {T<:Number}
    npop = [cpoly(p, z) for p in pop]
    return complex_cs_tssos_first(npop, length(z), d, numeq=numeq, RemSig=RemSig, nb=nb, CS=CS, cliques=cliques, TS=TS, merge=merge, 
    md=md, solver=solver, reducebasis=reducebasis, QUIET=QUIET, solve=solve, dualize=dualize, solution=solution,
    MomentOne=MomentOne, ConjugateBasis=ConjugateBasis, Gram=Gram, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, 
    writetofile=writetofile, normality=normality, pop=pop, z=z, rtol=rtol, gtol=gtol, ftol=ftol)
end

function complex_tssos_first(pop::Vector{Poly{T}}, z, d; numeq=0, RemSig=false, nb=0, TS="block", 
    merge=false, md=3, solver="Mosek", reducebasis=false, QUIET=false, solve=true, solution=false, dualize=false, MomentOne=false, 
    ConjugateBasis=false, Gram=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false, normality=!ConjugateBasis, 
    rtol=1e-2, gtol=1e-2, ftol=1e-3) where {T<:Number}
    return complex_cs_tssos_first(pop, z, d, numeq=numeq, RemSig=RemSig, nb=nb, CS=false, TS=TS, merge=merge, 
    md=md, solver=solver, reducebasis=reducebasis, QUIET=QUIET, solve=solve, dualize=dualize, solution=solution,
    MomentOne=MomentOne, ConjugateBasis=ConjugateBasis, Gram=Gram, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, 
    writetofile=writetofile, normality=normality, rtol=rtol, gtol=gtol, ftol=ftol)
end

"""
    opt,sol,data = complex_cs_tssos_first(npop::Vector{cpoly{T}}, n, d; nb=0, numeq=0, CS="MF", cliques=[], TS="block", merge=false, md=3, solver="Mosek", solution=false, dualize=false,
    QUIET=false, solve=true, Gram=false, MomentOne=false, normality=!ConjugateBasis, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(),
    rtol=1e-2, gtol=1e-2, ftol=1e-3) where {T<:Number}

Compute the first TS step of the CS-TSSOS hierarchy for constrained complex polynomial optimization. 
"""
function complex_cs_tssos_first(npop::Vector{cpoly{T}}, n::Int, d; numeq=0, RemSig=false, nb=0, CS="MF", cliques=[], 
    TS="block", merge=false, md=3, solver="Mosek", reducebasis=false, QUIET=false, solve=true, solution=false, dualize=false, MomentOne=false, ConjugateBasis=false, 
    Gram=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false, normality=!ConjugateBasis, pop=nothing, z=nothing, 
    rtol=1e-2, gtol=1e-2, ftol=1e-3) where {T<:Number}
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    ipart = T <: Real ? false : true
    obj = npop[1]
    if nb > 0
        obj = arrange(obj, nb)
    end
    ineq_cons = [cpoly{T}([tuple(UInt16[], UInt16[])], [1]); npop[2:end-numeq]]
    eq_cons = npop[end-numeq+1:end]
    if !isempty(cliques)
        cql = length(cliques)
        cliquesize = length.(cliques)
    else
        time = @elapsed begin
        CS = CS == true ? "MF" : CS
        cliques,cql,cliquesize = clique_decomp(npop, n, order=d, alg=CS, QUIET=QUIET)
        end
        if CS != false && QUIET == false
            mc = maximum(cliquesize)
            println("Obtained the variable cliques in $time seconds.\nThe maximal size of cliques is $mc.")
        end
    end
    I,J,Iprime,Jprime = assign_constraint(ineq_cons, eq_cons, cliques, cql)
    if d == "min"
        if ConjugateBasis == false
            rlorder = [isempty(I[i]) && isempty(J[i]) ? 1 : maximum([maxcdeg.(ineq_cons[I[i]]); maxcdeg.(eq_cons[J[i]])]) for i = 1:cql]
        else
            rlorder = [isempty(I[i]) && isempty(J[i]) ? 1 : ceil(Int, maximum([maxdeg.(ineq_cons[I[i]]); maxdeg.(eq_cons[J[i]])])/2) for i = 1:cql]
        end
    else
        rlorder = d*ones(Int, cql)
    end
    ebasis = Vector{Vector{Vector{Tuple{Vector{UInt16},Vector{UInt16}}}}}(undef, cql)
    if ConjugateBasis == false
        if normality == 0
            basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, cql)
        else
            basis = Vector{Vector{Vector{Union{Vector{UInt16}, Tuple{Vector{UInt16},Vector{UInt16}}}}}}(undef, cql)
        end
    else
        basis = Vector{Vector{Vector{Tuple{Vector{UInt16},Vector{UInt16}}}}}(undef, cql)
    end
    for i = 1:cql
        ebasis[i] = Vector{Vector{Tuple{Vector{UInt16},Vector{UInt16}}}}(undef, length(J[i]))
        if ConjugateBasis == false
            if normality == 0
                basis[i] = [get_basis(cliques[i], rlorder[i]-maxcdeg(ineq_cons[j])) for j in I[i]]
            else
                basis[i] = Vector{Vector{Union{Vector{UInt16}, Tuple{Vector{UInt16},Vector{UInt16}}}}}(undef, length(I[i])+cliquesize[i])
                basis[i][1] = get_basis(cliques[i], rlorder[i])
                temp = get_basis(cliques[i], Int(normality))
                for s = 1:cliquesize[i]
                    basis[i][s+1] = [[tuple(item, UInt16[]) for item in temp]; [tuple(item, UInt16[cliques[i][s]]) for item in temp]]
                    if nb > 0
                        basis[i][s+1] = reduce_unitnorm.(basis[i][s+1], nb)
                        unique!(basis[i][s+1])
                    end
                end
                for s = 1:length(I[i])-1
                    basis[i][s+1+cliquesize[i]] = get_basis(cliques[i], rlorder[i]-maxcdeg(ineq_cons[I[i][s+1]]))
                end
                I[i] = [ones(Int, cliquesize[i]); I[i]]
            end
            for s = 1:length(J[i])
                if rlorder[i] < maxcdeg(eq_cons[J[i][s]])
                    @error "The relaxation order is too small!"
                end
                temp = get_basis(cliques[i], rlorder[i]-maxcdeg(eq_cons[J[i][s]]))
                ebasis[i][s] = vec([tuple(item1, item2) for item1 in temp, item2 in temp])
                if nb > 0
                    ebasis[i][s] = reduce_unitnorm.(ebasis[i][s], nb)
                    unique!(ebasis[i][s])
                end
                sort!(ebasis[i][s])
            end
        else
            if normality < d
                basis[i] = Vector{Vector{Tuple{Vector{UInt16},Vector{UInt16}}}}(undef, length(I[i]))
                basis[i][1] = get_conjugate_basis(cliques[i], rlorder[i], nb=nb)
                for s = 1:length(I[i])
                    basis[i][s] = get_conjugate_basis(cliques[i], rlorder[i]-Int(ceil(maxdeg(ineq_cons[I[i][s]])/2)), nb=nb)
                end
            else
                basis[i] = Vector{Tuple{Vector{UInt16},Vector{UInt16}}}(undef, length(I[i])+cliquesize[i])
                basis[i][1] = get_conjugate_basis(cliques[i], rlorder[i], nb=nb)
                temp = get_basis(cliques[i], normality)
                for s = 1:cliquesize[i]
                    basis[i][s+1] = [[tuple(item, UInt16[]) for item in temp]; [tuple(item, UInt16[cliques[i][s]]) for item in temp]]
                    if nb > 0
                        basis[i][s+1] = reduce_unitnorm.(basis[i][s+1], nb)
                        unique!(basis[i][s+1])
                    end
                end
                for s = 1:length(I[i])-1
                    basis[i][s+1+cliquesize[i]] = get_conjugate_basis(cliques[i], rlorder[i]-Int(ceil(maxdeg(ineq_cons[I[i][s+1]])/2)), nb=nb)
                end
                I[i] = [ones(Int, cliquesize[i]); I[i]]
            end
            for s = 1:length(J[i])
                ebasis[i][s] = get_conjugate_basis(cliques[i], 2*rlorder[i]-maxdeg(eq_cons[J[i][s]]), nb=nb)
            end
        end
    end
    tsupp = filter(item -> item[1] <= item[2], vcat([p.supp for p in npop]...))
    sort!(tsupp)
    unique!(tsupp)
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,eblocks = get_blocks(I, J, ineq_cons, eq_cons, cliques, cql, tsupp, basis, ebasis, TS=TS, merge=merge, md=md, nb=nb)
    if RemSig == true
        for i = 1:cql
            basis[i][1] = basis[i][1][union(blocks[i][1][blocksize[i][1] .> 1]...)]
        end
        tsupp = filter(item -> item[1] <= item[2], vcat([p.supp for p in npop]...))
        sort!(tsupp)
        unique!(tsupp)
        blocks,cl,blocksize,eblocks = get_blocks(I, J, ineq_cons, eq_cons, cliques, cql, tsupp, basis, ebasis, TS=TS, merge=merge, md=md, nb=nb)
    end
    if reducebasis == true
        tsupp = get_csupp(rlorder, basis, ebasis, ineq_cons, eq_cons, I, J, Iprime, Jprime, blocks, eblocks, cl, blocksize, cql, cliquesize, norm=true, nb=nb, ConjugateBasis=ConjugateBasis, normality=normality)
        foreach(item -> item[1] == item[2] ? push!(tsupp, item[1]) : nothing, obj.supp)
        sort!(tsupp)
        unique!(tsupp)
        flag = 0
        for i = 1:cql
            ind = [bfind(tsupp, item) !== nothing for item in basis[1][i]]
            if !all(ind)
                basis[1][i] = basis[1][i][ind]
                flag = 1
            end
        end
        if flag == 1
            tsupp = filter(item -> item[1] <= item[2], vcat([p.supp for p in npop]...))
            sort!(tsupp)
            unique!(tsupp)
            blocks,cl,blocksize,eblocks = get_blocks(I, J, ineq_cons, eq_cons, cliques, cql, tsupp, basis, ebasis, TS=TS, nb=nb, merge=merge, md=md)
        end
    end
    end
    if QUIET == false
        mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,sol,GramMat,multiplier,SDP_status = solvesdp(obj, ineq_cons, eq_cons, n, rlorder, basis, ebasis, cliques, cql, cliquesize, I, J, 
    Iprime, Jprime, blocks, eblocks, cl, blocksize, z=z, QUIET=QUIET, TS=TS, ConjugateBasis=ConjugateBasis, solver=solver, solve=solve, MomentOne=MomentOne, 
    ipart=ipart, solution=solution, Gram=Gram, nb=nb, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, dualize=dualize, writetofile=writetofile, normality=normality)
    flag = 1
    if solution == true
        if TS != false
            if pop !== nothing
                asol = check_solution([sol], opt, pop, z, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=true)
                if asol === nothing
                    rsol = [real(sol); imag(sol)]
                    ub = evaluate(obj, sol)
                    gap = abs(opt-ub)/max(1, abs(ub))
                    rsol = gap > 0.5 ? randn(2n) : rsol
                    rsol[abs.(rsol) .< 1e-10] .= 1e-10
                    rpop,x = complex_to_real(pop, z)
                    ub,rsol,status = local_solution([poly(p, x) for p in rpop], 2n, numeq=numeq, startpoint=rsol, QUIET=true)
                    if status == MOI.LOCALLY_SOLVED
                        sol = rsol[1:n] + rsol[n+1:2n]*im
                        ub = evaluate(obj, sol)
                        gap = abs(opt-ub)/max(1, abs(ub))
                        if gap < gtol
                            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
                            println("Successfully extracted one globally optimal solution.")
                        else
                            @printf "Found a locally optimal solution by Ipopt, giving an upper bound: %.8f.\nThe relative optimality gap is: %.6f%%.\n" ub 100*gap
                        end
                    else
                        sol = nothing
                    end
                end
            else
                sol = check_solution([sol], opt, npop, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
            end
        else
            if cql == 1 && pop !== nothing
                sol = extract_solutions_robust(moment[1][1], n, d, type=ComplexF64, pop=pop, x=z, lb=opt, numeq=numeq, check=true, rtol=rtol, gtol=gtol, ftol=ftol, QUIET=QUIET)[1]
            else
                sol = extract_csolutions_robust(moment, n, d, cliques, cql, cliquesize, pop=pop, z=z, npop=npop, lb=opt, numeq=numeq, check=true, rtol=rtol, gtol=gtol, ftol=ftol, QUIET=QUIET)[1]
            end
        end
        if sol !== nothing
            flag = 0
        end
    end
    data = ccpop_data(pop, obj, ineq_cons, eq_cons, z, rlorder, n, nb, numeq, ipart, ConjugateBasis, normality, basis, ebasis, ksupp, cliquesize, cliques, I, J, Iprime, Jprime, blocksize, blocks, eblocks, 
    GramMat, multiplier, moment, solver, SDP_status, rtol, gtol, ftol, flag)
    return opt,sol,data
end

"""
    opt,sol,data = complex_cs_tssos_higher!(data; TS="block", merge=false, md=3, QUIET=false, solve=true, dualize=false,
    solution=false, Gram=false, MomentOne=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para())

Compute higher TS steps of the CS-TSSOS hierarchy.
"""
function complex_cs_tssos_higher!(data::ccpop_data; TS="block", merge=false, md=3, QUIET=false, solve=true, solution=false, Gram=false, dualize=false, 
    MomentOne=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false)
    obj = data.obj
    eq_cons = data.eq_cons
    ineq_cons = data.ineq_cons
    n = data.n
    nb = data.nb
    numeq = data.numeq
    basis = data.basis
    ebasis = data.ebasis
    cliques = data.cliques
    cliquesize = data.cliquesize
    I = data.I
    J = data.J
    cql = length(cliques)
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,eblocks = get_blocks(I, J, ineq_cons, eq_cons, cliques, cql, data.ksupp, basis, ebasis, TS=TS, merge=merge, md=md, nb=nb)
    end
    if blocksize == data.blocksize && eblocks == data.eblocks
        println("No higher TS step of the CS-TSSOS hierarchy!")
        opt = sol = nothing
    else
        if TS != false && QUIET == false
            mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,sol,GramMat,multiplier,SDP_status = solvesdp(obj, ineq_cons, eq_cons, n, data.rlorder, basis, ebasis, cliques, cql, cliquesize, I, J, 
    data.Iprime, data.Jprime, blocks, eblocks, cl, blocksize, z=data.z, QUIET=QUIET, TS=TS, ConjugateBasis=data.ConjugateBasis, solver=data.solver, solve=solve, MomentOne=MomentOne, 
    ipart=data.ipart, solution=solution, Gram=Gram, nb=nb, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, dualize=dualize, writetofile=writetofile, normality=data.normality)
        if solution == true
            if data.pop !== nothing
                asol = check_solution([sol], opt, data.pop, data.z, numeq=numeq, gtol=data.gtol, ftol=data.ftol, QUIET=true)
                if asol === nothing
                    rsol = [real(sol); imag(sol)]
                    ub = evaluate(obj, sol)
                    gap = abs(opt-ub)/max(1, abs(ub))
                    rsol = gap > 0.5 ? randn(2n) : rsol
                    rsol[abs.(rsol) .< 1e-10] .= 1e-10
                    rpop,x = complex_to_real(data.pop, data.z)
                    ub,rsol,status = local_solution([poly(p, x) for p in rpop], 2n, numeq=numeq, startpoint=rsol, QUIET=true)
                    if status == MOI.LOCALLY_SOLVED
                        sol = rsol[1:n] + rsol[n+1:2n]*im
                        ub = evaluate(obj, sol)
                        gap = abs(opt-ub)/max(1, abs(ub))
                        if gap < data.gtol
                            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
                            println("Successfully extracted one globally optimal solution.")
                        else
                            @printf "Found a locally optimal solution by Ipopt, giving an upper bound: %.8f.\nThe relative optimality gap is: %.6f%%.\n" ub 100*gap
                        end
                    else
                        sol = nothing
                    end
                end
            else
                sol = check_solution([sol], opt, [obj; ineq_cons[2:end]; eq_cons], numeq=numeq, gtol=data.gtol, ftol=data.ftol, QUIET=QUIET)
            end
            if sol !== nothing
                data.flag = 0
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

function complex_tssos_higher!(data::ccpop_data; TS="block", merge=false, md=3, QUIET=false, solve=true, solution=false, Gram=false, dualize=false, 
    MomentOne=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false)
    return complex_cs_tssos_higher!(data, TS=TS, merge=merge, md=md, QUIET=QUIET, solve=solve, solution=solution, Gram=Gram, dualize=dualize, 
    MomentOne=MomentOne, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, writetofile=writetofile)
end

function get_csupp(rlorder::Vector{Int}, basis, ebasis, ineq_cons::Vector{T1}, eq_cons::Vector{T2}, I, J, Iprime, Jprime, blocks, eblocks, cl, blocksize, cql, cliquesize; norm=false, nb=0, ConjugateBasis=false, normality=1) where {T1,T2<:cpoly}
    csupp = norm == true ? Vector{UInt16}[] : Tuple{Vector{UInt16},Vector{UInt16}}[]
    for i = 1:cql
        if ConjugateBasis == false
            a = normality > 0 ? 1 + cliquesize[i] : 1
        else
            a = normality >= rlorder[i] ? 1 + cliquesize[i] : 1
        end
        for (j, w) in enumerate(I[i][a+1:end]), l = 1:cl[i][j+a], t = 1:blocksize[i][j+a][l], r = t:blocksize[i][j+a][l], item in ineq_cons[w].supp
            ind1 = blocks[i][j+a][l][t]
            ind2 = blocks[i][j+a][l][r]
            if ConjugateBasis == false
                @inbounds bi = tuple(sadd(basis[i][j+a][ind1], item[1]), sadd(basis[i][j+a][ind2], item[2]))
            else
                @inbounds bi = tuple(sadd(basis[i][j+a][ind1][1], item[1], basis[i][j+a][ind2][2]), sadd(basis[i][j+a][ind1][2], item[2], basis[i][j+a][ind2][1]))
            end
            nb > 0 ? bi = reduce_unitnorm(bi, nb) : nothing
            if norm == true
                bi[1] == bi[2] ? push!(csupp, bi[1]) : nothing
            else
                bi[1] <= bi[2] ? push!(csupp, bi) : push!(csupp, conj(bi))
            end
        end
        for (j, w) in enumerate(J[i]), k in eblocks[i][j], item in eq_cons[w].supp
            @inbounds bi = tuple(sadd(ebasis[i][j][k][1], item[1]), sadd(ebasis[i][j][k][2], item[2]))
            nb > 0 ? bi = reduce_unitnorm(bi, nb) : nothing
            if norm == true
                bi[1] == bi[2] ? push!(csupp, bi[1]) : nothing
            else
                bi[1] <= bi[2] ? push!(csupp, bi) : push!(csupp, conj(bi))
            end
        end
    end
    for i in Iprime, item in ineq_cons[i].supp
        nb > 0 ? item = reduce_unitnorm(item, nb) : nothing
        if norm == true
            item[1] == item[2] ? push!(csupp, item[1]) : nothing
        else
            item[1] <= item[2] ? push!(csupp, item) : nothing
        end
    end
    for i in Jprime, item in eq_cons[i].supp
        nb > 0 ? item = reduce_unitnorm(item, nb) : nothing
        if norm == true
            item[1] == item[2] ? push!(csupp, item[1]) : nothing
        else
            item[1] <= item[2] ? push!(csupp, item) : nothing
        end
    end
    return csupp
end

function solvesdp(obj::T1, ineq_cons::Vector{T2}, eq_cons::Vector{T3}, n, rlorder, basis, ebasis, cliques, cql, cliquesize, I, J, Iprime, Jprime, blocks, eblocks, cl, blocksize; 
    nb=0, z=nothing, QUIET=false, TS="block", ConjugateBasis=false, solver="Mosek", solve=true, dualize=false, Gram=false, MomentOne=false, ipart=true, solution=false, 
    cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false, normality=1) where {T1,T2,T3<:cpoly}
    tsupp = Tuple{Vector{UInt16},Vector{UInt16}}[]
    for i = 1:cql
        if ConjugateBasis == false
            a = normality > 0 ? 1 + cliquesize[i] : 1
        else
            a = normality >= rlorder[i] ? 1 + cliquesize[i] : 1
        end
        for s = 1:a, j = 1:cl[i][s], k = 1:blocksize[i][s][j], r = k:blocksize[i][s][j]
            if ConjugateBasis == false && s == 1
                @inbounds bi = tuple(basis[i][s][blocks[i][s][j][k]], basis[i][s][blocks[i][s][j][r]])
            else
                @inbounds bi = sadd(basis[i][s][blocks[i][s][j][k]], conj(basis[i][s][blocks[i][s][j][r]]))
            end
            nb > 0 ? bi = reduce_unitnorm(bi, nb) : nothing
            bi[1] <= bi[2] ? push!(tsupp, bi) : push!(tsupp, conj(bi))
        end
    end
    if TS != false
        csupp = get_csupp(rlorder, basis, ebasis, ineq_cons, eq_cons, I, J, Iprime, Jprime, blocks, eblocks, cl, blocksize, cql, cliquesize, ConjugateBasis=ConjugateBasis, nb=nb, normality=normality)
        append!(tsupp, csupp)
    end
    if (MomentOne == true || solution == true) && TS != false
        ksupp = copy(tsupp)
        for i = 1:cql, j = 1:cliquesize[i]
            push!(tsupp, tuple(UInt16[], UInt16[cliques[i][j]]))
            for k = j+1:cliquesize[i]
                bi = tuple(UInt16[cliques[i][j]], UInt16[cliques[i][k]])
                nb > 0 ? bi = reduce_unitnorm(bi, nb) : nothing
                push!(tsupp, bi)
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
    objv = moment = sol = GramMat = multiplier = SDP_status= nothing
    if solve == true
        if QUIET == false
            println("Assembling the SDP...")
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
        rcons = [AffExpr(0) for i=1:length(tsupp)]
        if ipart == true
            icons = [AffExpr(0) for i=1:length(tsupp)]
        end
        pos = Vector{Vector{Vector{Symmetric{VariableRef}}}}(undef, cql)
        free = Vector{Vector{Vector{VariableRef}}}(undef, cql)
        for i = 1:cql
            if (MomentOne == true || solution == true) && TS != false
                bs = cliquesize[i] + 1
                if ipart == true
                    pos0 = @variable(model, [1:2*bs, 1:2*bs], PSD)
                else
                    pos0 = @variable(model, [1:bs, 1:bs], PSD)
                end
                for t = 1:bs, r = t:bs
                    if t == 1 && r == 1
                        bi = tuple(UInt16[], UInt16[])
                    elseif t == 1 && r > 1
                        bi = tuple(UInt16[], UInt16[cliques[i][r-1]])
                    else
                        bi = tuple(UInt16[cliques[i][t-1]], UInt16[cliques[i][r-1]])
                        nb > 0 ? bi = reduce_unitnorm(bi, nb) : nothing
                    end
                    Locb = bfind(tsupp, bi)
                    if ipart == true
                        @inbounds add_to_expression!(rcons[Locb], pos0[t,r]+pos0[t+bs,r+bs])
                        @inbounds add_to_expression!(icons[Locb], pos0[t,r+bs]-pos0[r,t+bs])
                    else
                        @inbounds add_to_expression!(rcons[Locb], pos0[t,r])
                    end
                end
            end
            pos[i] = Vector{Vector{Symmetric{VariableRef}}}(undef, length(I[i]))
            for (j, p) in enumerate(ineq_cons[I[i]])
                pos[i][j] = Vector{Symmetric{VariableRef}}(undef, cl[i][j])
                for l = 1:cl[i][j]
                    bs = blocksize[i][j][l]
                    if ipart == true && bs != 1
                        pos[i][j][l] = @variable(model, [1:2bs, 1:2bs], PSD)
                    else
                        pos[i][j][l] = @variable(model, [1:bs, 1:bs], PSD)
                    end
                    for t = 1:bs, r = 1:bs, (s, item) in enumerate(p.supp)
                        if typeof(basis[i][j][1]) == Vector{UInt16}
                            @inbounds bi = tuple(sadd(basis[i][j][blocks[i][j][l][t]], item[1]), sadd(basis[i][j][blocks[i][j][l][r]], item[2]))
                        else
                            @inbounds bi = tuple(sadd(basis[i][j][blocks[i][j][l][t]][1], item[1], basis[i][j][blocks[i][j][l][r]][2]), sadd(basis[i][j][blocks[i][j][l][t]][2], item[2], basis[i][j][blocks[i][j][l][r]][1]))
                        end
                        nb > 0 ? bi = reduce_unitnorm(bi, nb) : nothing
                        if bi[1] <= bi[2]
                            Locb = bfind(tsupp, bi)
                            if ipart == true
                                if bs != 1
                                    @inbounds add_to_expression!(rcons[Locb], real(p.coe[s]), pos[i][j][l][t,r]+pos[i][j][l][t+bs,r+bs])
                                    @inbounds add_to_expression!(rcons[Locb], -imag(p.coe[s]), pos[i][j][l][t,r+bs]-pos[i][j][l][r,t+bs])
                                    if bi[1] != bi[2]
                                        @inbounds add_to_expression!(icons[Locb], imag(p.coe[s]), pos[i][j][l][t,r]+pos[i][j][l][t+bs,r+bs])
                                        @inbounds add_to_expression!(icons[Locb], real(p.coe[s]), pos[i][j][l][t,r+bs]-pos[i][j][l][r,t+bs])
                                    end
                                else
                                    @inbounds add_to_expression!(rcons[Locb], real(p.coe[s]), pos[i][j][l][t,r])
                                    if bi[1] != bi[2]
                                        @inbounds add_to_expression!(icons[Locb], imag(p.coe[s]), pos[i][j][l][t,r])
                                    end
                                end
                            else
                                @inbounds add_to_expression!(rcons[Locb], p.coe[s], pos[i][j][l][t,r])
                            end
                        end
                    end
                end
            end
            free[i] = Vector{Vector{VariableRef}}(undef, length(J[i]))
            for (j, p) in enumerate(eq_cons[J[i]])
                mons = ebasis[i][j][eblocks[i][j]]
                temp = mons[[item[1] <= item[2] for item in mons]]
                lb = length(temp)
                if ipart == true
                    free[i][j] = @variable(model, [1:2*lb])
                else
                    free[i][j] = @variable(model, [1:lb])
                end
                for k in eblocks[i][j], (s, item) in enumerate(p.supp)
                    @inbounds bi = tuple(sadd(ebasis[i][j][k][1], item[1]), sadd(ebasis[i][j][k][2], item[2]))
                    nb > 0 ? bi = reduce_unitnorm(bi, nb) : nothing
                    if bi[1] <= bi[2]
                        Locb = bfind(tsupp, bi)
                        if ebasis[i][j][k][1] <= ebasis[i][j][k][2]
                            loc = bfind(temp, ebasis[i][j][k])
                            tag = 1
                            if ebasis[i][j][k][1] == ebasis[i][j][k][2]
                                tag = 0
                            end
                        else
                            loc = bfind(temp, conj(ebasis[i][j][k]))
                            tag = -1
                        end
                        if ipart == true
                            @inbounds add_to_expression!(rcons[Locb], real(p.coe[s])*free[i][j][loc]-tag*imag(p.coe[s])*free[i][j][loc+lb])
                            if bi[1] != bi[2]
                                @inbounds add_to_expression!(icons[Locb], tag*real(p.coe[s])*free[i][j][loc+lb]+imag(p.coe[s])*free[i][j][loc])
                            end
                        else
                            @inbounds add_to_expression!(rcons[Locb], p.coe[s], free[i][j][loc])
                        end
                    end
                end
            end
        end
        for i in Iprime
            pos0 = @variable(model, lower_bound=0)
            for (j, item) in enumerate(ineq_cons[i].supp)
                if item[1] <= item[2]
                    Locb = bfind(tsupp, item)
                    if ipart == true
                        @inbounds add_to_expression!(icons[Locb], imag(ineq_cons[i].coe[j]), pos0)
                    end
                    @inbounds add_to_expression!(rcons[Locb], real(ineq_cons[i].coe[j]), pos0)
                end
            end
        end
        for i in Jprime
            pos0 = @variable(model)
            for (j, item) in enumerate(eq_cons[i].supp)
                if item[1] <= item[2]
                    Locb = bfind(tsupp, item)
                    if ipart == true
                        @inbounds add_to_expression!(icons[Locb], imag(eq_cons[i].coe[j]), pos0)
                    end
                    @inbounds add_to_expression!(rcons[Locb], real(eq_cons[i].coe[j]), pos0)
                end
            end
        end
        for (i, bi) in enumerate(obj.supp)
            if bi[1] <= bi[2]
                Locb = bfind(tsupp, bi)
                if Locb === nothing
                    @error "The monomial basis is not enough!"
                    return nothing,nothing,nothing,nothing,nothing,nothing,nothing
                else
                    rcons[Locb] -= real(obj.coe[i])
                    if ipart == true && bi[1] != bi[2]
                        icons[Locb] -= imag(obj.coe[i])
                    end
                end
            end
        end
        @variable(model, lower)
        rcons[1] += lower
        @constraint(model, rcon, rcons == zeros(length(tsupp)))
        if ipart == true
            icons = icons[[item[1] != item[2] for item in tsupp]]
            @constraint(model, icon, icons == zeros(length(icons)))
        end
        @objective(model, Max, lower)
        end
        if QUIET == false
            println("There are $(length(tsupp) + length(icons)) affine constraints.")
            println("SDP assembling time: $time seconds.")
        end
        # println(all_constraints(model, include_variable_in_set_constraints=false))
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
                ctype = ipart==true ? ComplexF64 : Float64
                GramMat = Vector{Vector{Vector{Matrix{ctype}}}}(undef, cql)
                multiplier = Vector{Vector{Poly{ctype}}}(undef, cql)
                for i = 1:cql
                    GramMat[i] = Vector{Vector{Matrix{ctype}}}(undef, length(I[i]))
                    for j = 1:length(I[i])
                        GramMat[i][j] = Vector{Matrix{ctype}}(undef, cl[i][j])
                        for l = 1:cl[i][j]
                            if ipart == true && blocksize[i][j][l] > 1
                                bs = blocksize[i][j][l]
                                temp = value.(pos[i][j][l][1:bs,bs+1:2bs])
                                GramMat[i][j][l] = value.(pos[i][j][l][1:bs,1:bs]+pos[i][j][l][bs+1:2bs,bs+1:2bs]) + im*(temp-temp')
                            else
                                GramMat[i][j][l] = value.(pos[i][j][l])
                            end
                        end
                    end
                    if !isempty(J[i])
                        multiplier[i] = Vector{Poly{ctype}}(undef, length(J[i]))
                        for k = 1:length(J[i])
                            temp = ebasis[i][k][eblocks[i][k]][[item[1] <= item[2] for item in ebasis[i][k][eblocks[i][k]]]]
                            lb = length(temp)
                            if ctype == Float64
                                tau = sum(prod(z[temp[j][1]])*conj(prod(z[temp[j][2]]))*value(free[i][k][j]) for j = 1:lb)
                            else
                                tau = sum(prod(z[temp[j][1]])*conj(prod(z[temp[j][2]]))*(value(free[i][k][j])+value(free[i][k][lb+j])*im) for j = 1:lb)
                            end
                            multiplier[i][k] = tau + conj(tau)
                        end
                    end
                end
            end
            rmeasure = -dual(rcon)
            imeasure = nothing
            if ipart == true
                imeasure = -dual(icon)
            end
            if solution == true
                sol = [rmeasure[bfind(tsupp, tuple(UInt16[], UInt16[i]))] for i = 1:n]
                if ipart == true
                    sol += [-imeasure[bfind(tsupp, tuple(UInt16[], UInt16[i]))] for i = 1:n]*im
                end          
            end
            moment = get_cmoment(rmeasure, imeasure, tsupp, cql, blocks, cl, blocksize, basis, ipart=ipart, nb=nb, ConjugateBasis=ConjugateBasis)
        end
    end
    return objv,ksupp,moment,sol,GramMat,multiplier,SDP_status
end

function get_eblock(tsupp::Vector{Tuple{Vector{UInt16},Vector{UInt16}}}, hsupp::Vector{Tuple{Vector{UInt16},Vector{UInt16}}}, basis::Vector{Tuple{Vector{UInt16},Vector{UInt16}}}; nb=0)
    eblock = Int[]
    for (i, item) in enumerate(basis)
        flag = 0
        for temp in hsupp
            bi = tuple(sadd(item[1], temp[1]), sadd(item[2], temp[2]))
            nb > 0 ? bi = reduce_unitnorm(bi, nb) : nothing
            if (bi[1] <= bi[2] && bfind(tsupp, bi) !== nothing) || (bi[1] > bi[2] && bfind(tsupp, conj(bi)) !== nothing)
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

function get_blocks(tsupp, ineq_cons::Vector{T1}, eq_cons::Vector{T2}, basis, ebasis; nb=0, TS="block", merge=false, md=3) where {T1,T2<:cpoly}
    blocks = Vector{Vector{Vector{Int}}}(undef, length(ineq_cons))
    blocksize = Vector{Vector{Int}}(undef, length(ineq_cons))
    cl = Vector{Int}(undef, length(ineq_cons))
    eblocks = Vector{Vector{Int}}(undef, length(eq_cons))
    if TS == false
        for k = 1:length(ineq_cons) 
            blocks[k],blocksize[k],cl[k] = [Vector(1:length(basis[k]))],[length(basis[k])],1
        end
        for k = 1:length(eq_cons)
            eblocks[k] = Vector(1:length(ebasis[k]))
        end
    else
        for (k, p) in enumerate(ineq_cons)
            G = get_graph(tsupp, p.supp, basis[k], nb=nb)
            if TS == "block"
                blocks[k] = connected_components(G)
                blocksize[k] = length.(blocks[k])
                cl[k] = length(blocksize[k])
            else
                blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS)
                if merge == true
                    blocks[k],cl[k],blocksize[k] = clique_merge!(blocks[k], d=md, QUIET=true)
                end
            end
        end
        for (k, p) in enumerate(eq_cons)
            eblocks[k] = get_eblock(tsupp, p.supp, ebasis[k], nb=nb)
        end
    end
    return blocks,cl,blocksize,eblocks
end

function get_blocks(I, J, ineq_cons::Vector{T1}, eq_cons::Vector{T2}, cliques, cql, tsupp, basis, ebasis; TS="block", nb=0, merge=false, md=3) where {T1,T2<:cpoly}
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
    for i = 1:cql
        ksupp = TS == false ? nothing : tsupp[[issubset(union(item[1], item[2]), cliques[i]) for item in tsupp]]
        blocks[i],cl[i],blocksize[i],eblocks[i] = get_blocks(ksupp, ineq_cons[I[i]], eq_cons[J[i]], basis[i], ebasis[i], TS=TS, merge=merge, md=md, nb=nb)
    end
    return blocks,cl,blocksize,eblocks
end

function assign_constraint(ineq_cons::Vector{T1}, eq_cons::Vector{T2}, cliques, cql) where {T1,T2<:cpoly}
    I = [Int[] for i=1:cql]
    J = [Int[] for i=1:cql]
    Iprime = Int[]
    Jprime = Int[]
    for (i, p) in enumerate(ineq_cons)
        ind = findall(k -> issubset(unique(reduce(vcat, [item[1] for item in p.supp])), cliques[k]), 1:cql)
        if isempty(ind)
            push!(Iprime, i)
        else
            push!.(I[ind], i)
        end
    end
    for (i, p) in enumerate(eq_cons)
        ind = findall(k -> issubset(unique(reduce(vcat, [item[1] for item in p.supp])), cliques[k]), 1:cql)
        if isempty(ind)
            push!(Jprime, i)
        else
            push!.(J[ind], i)
        end
    end
    return I,J,Iprime,Jprime
end

function get_graph(tsupp::Vector{Tuple{Vector{UInt16},Vector{UInt16}}}, supp, basis; nb=0)
    lb = length(basis)
    G = SimpleGraph(lb)
    for i = 1:lb, j = i+1:lb
        r = 1
        while r <= length(supp)
            if typeof(basis[i]) == Vector{UInt16}
                bi = tuple(sadd(basis[i], supp[r][1]), sadd(basis[j], supp[r][2]))
            else
                bi = tuple(sadd(basis[i][1], supp[r][1], basis[j][2]), sadd(basis[i][2], supp[r][2], basis[j][1]))
            end
            nb > 0 ? bi = reduce_unitnorm(bi, nb) : nothing
            bi[1] > bi[2] ? bi = conj(bi) : nothing
            if bfind(tsupp, bi) !== nothing
               break
            else
                r += 1
            end
        end
        if r <= length(supp)
            add_edge!(G, i, j)
        end
    end
    return G
end

function clique_decomp(npop::Vector{T}, n; order="min", alg="MF", ReducedCS=true, QUIET=false) where {T<:cpoly}
    if alg == false
        cliques,cql,cliquesize = [Vector(1:n)],1,[n]
    else
        G = SimpleGraph(n)
        for (i, p) in enumerate(npop)
            if i == 1 || order == "min" || (ReducedCS == true && order == maxcdeg(p))
                foreach(a -> add_clique!(G, unique([a[1]; a[2]])), p.supp)
            else
                add_clique!(G, unique(vcat([item[1] for item in p.supp]...)))
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

function get_cmoment(rmeasure, imeasure, tsupp, cql, blocks, cl, blocksize, basis; ipart=true, nb=0, ConjugateBasis=false)
    ctype = ipart == true ? ComplexF64 : Float64
    moment = Vector{Vector{Matrix{ctype}}}(undef, cql)
    if ipart == true
        itsupp = tsupp[[item[1] != item[2] for item in tsupp]]
    end
    for i = 1:cql
        moment[i] = Vector{Matrix{ctype}}(undef, cl[i][1])
        for l = 1:cl[i][1]
            bs = blocksize[i][1][l]
            rtemp = zeros(Float64, bs, bs)
            if ipart == true
                itemp = zeros(Float64, bs, bs)
            end
            for t = 1:bs, r = t:bs
                if ConjugateBasis == false
                    bi = tuple(basis[i][1][blocks[i][1][l][t]], basis[i][1][blocks[i][1][l][r]])
                else
                    bi = sadd(basis[i][1][blocks[i][1][l][t]], conj(basis[i][1][blocks[i][1][l][r]]))
                end
                nb > 0 ? bi = reduce_unitnorm(bi, nb) : nothing
                if ipart == true
                    if bi[1] < bi[2]
                        Locb = bfind(itsupp, bi)
                        itemp[t,r] = imeasure[Locb]
                    elseif bi[1] > bi[2]
                        Locb = bfind(itsupp, conj(bi))
                        itemp[t,r] = -imeasure[Locb]
                    end
                end
                bi[1] > bi[2] ? bi = conj(bi) : nothing
                Locb = bfind(tsupp, bi)
                rtemp[t,r] = rmeasure[Locb]
            end
            rtemp = (rtemp + rtemp')/2
            if ipart == true
                itemp = (itemp - itemp')/2
                moment[i][l] = rtemp - itemp*im
            else
                moment[i][l] = rtemp
            end
        end
    end
    return moment
end
