mutable struct pop_data
    pop # polynomial optimiztion problem
    obj # objective
    ineq_cons # inequality constraints
    eq_cons # equality constraints
    x # variables
    n # number of variables
    nb # number of binary variables
    numeq # number of equality constraints
    gb # Grobner basis
    leadsupp # leader terms of the Grobner basis
    basis # monomial bases
    ebasis # monomial bases for equality constraints
    ksupp # extended support at the k-th step
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
function show_blocks(data::pop_data; include_constraints=false)
    for (j, block) in enumerate(data.blocks[1])
        print("block $j: ")
        println([prod(data.x[data.basis[1][item]]) for item in block])
    end
    if include_constraints == true
        for (i, blocks) in enumerate(data.blocks[2:end]), (j, block) in enumerate(blocks)
            print("constraint $i, block $j: ")
            println([prod(data.x[data.basis[i+1][item]]) for item in block])
        end
    end
end

"""
    opt,sol,data = tssos(pop, x, d; nb=0, numeq=0, GroebnerBasis=true, basis=[], reducebasis=false, TS="block", 
    merge=false, md=3, QUIET=false, solve=true, MomentOne=false, Gram=false, solution=false, model=nothing,
    mosek_setting=mosek_para(), normality=false, rtol=1e-2, gtol=1e-2, ftol=1e-3)

Compute the first TS step of the TSSOS hierarchy for constrained polynomial optimization.
If `reducebasis=true`, then remove monomials from the monomial basis by diagonal inconsistency.
If `GroebnerBasis=true`, then exploit the quotient ring structure defined by the equality constraints.
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
- `TS`: type of term sparsity (`"block"`, `"signsymmetry"`, `"MD"`, `"MF"`, `false`)
- `md`: tunable parameter for merging blocks
- `normality`: impose the normality condtions (`true`, `false`)
- `QUIET`: run in the quiet mode (`true`, `false`)
- `rtol`: tolerance for rank
- `gtol`: tolerance for global optimality gap
- `ftol`: tolerance for feasibility

# Output arguments
- `opt`: optimum
- `sol`: (near) optimal solution (if `solution=true`)
- `data`: other auxiliary data 
"""
function tssos(pop::Vector{Poly{T}}, x, d; nb=0, numeq=0, newton=false, feasibility=false, GroebnerBasis=false, basis=[], reducebasis=false, TS="block", 
    merge=false, md=3, QUIET=false, solve=true, dualize=false, MomentOne=false, Gram=false, solution=false, normality=false, mosek_setting=mosek_para(), 
    writetofile=false, model=nothing, rtol=1e-2, gtol=1e-2, ftol=1e-3) where {T<:Number}
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    n = length(x)
    if nb > 0
        pop = Groebner.normalform(x[1:nb].^2 .- 1, pop)
    end
    if numeq > 0 && GroebnerBasis == true
        println("Computing the Gröbner basis...")
        println("This might be slow. You can set GroebnerBasis=false to close it.")
        eq_cons = poly{T}[]
        gb = groebner(convert.(Poly{Float64}, pop[end-numeq+1:end]), ordering=DegRevLex())
        obj = poly(Groebner.normalform(gb, pop[1], ordering=DegRevLex()), x)
        lead = leading_ideal(gb, ordering=DegRevLex())
        leadsupp = [UInt16[] for i=1:length(lead)]
        for (i, mon) in enumerate(lead)
            mon = convert(DP.Monomial, mon)
            ind = mon.z .> 0
            vars = mon.vars[ind]
            exp = mon.z[ind]
            for j in eachindex(vars)
                append!(leadsupp[i], bfind_rev(x, vars[j])*ones(UInt16, exp[j]))
            end
        end
    else
        obj = poly(pop[1], x)
        eq_cons = [poly(p, x) for p in pop[end-numeq+1:end]]
        gb = leadsupp = []
    end
    ineq_cons = [poly{T}([UInt16[]], [1]); [poly(p, x) for p in pop[2:end-numeq]]]
    ss = nothing
    if normality == true || TS == "signsymmetry"
        ss = get_signsymmetry([obj; ineq_cons; eq_cons], n)
    end
    if isempty(basis)
        if newton == true
            supp = obj.supp
            if !isempty(obj.supp[1]) && feasibility == false
                supp = [supp; [UInt16[]]]
            end
            basis = [newton_basis(n, d, supp)]
        else
            basis = [get_basis(n, d-Int(ceil(maxdeg(p)/2)), nb=nb, lead=leadsupp) for p in ineq_cons]
        end
        ebasis = [get_basis(n, 2*d-maxdeg(p), nb=nb, lead=leadsupp) for p in eq_cons]
    end
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    ksupp = nothing
    if TS != false && TS != "signsymmetry"
        ksupp = reduce(vcat, [p.supp for p in [obj; ineq_cons; eq_cons]])
        for item in basis[1]
            push!(ksupp, sadd(item, item, nb=nb))
        end
        sort!(ksupp)
        unique!(ksupp)
    end
    time = @elapsed begin
    blocks,cl,blocksize,eblocks = get_blocks(ineq_cons, eq_cons, ksupp, basis, ebasis, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md, signsymmetry=ss)
    if reducebasis == true && GroebnerBasis == false
        csupp = get_csupp(ineq_cons, eq_cons, basis, ebasis, blocks, eblocks, cl, blocksize, nb=nb)
        psupp = [obj.supp; [UInt16[]]; csupp]
        basis[1],flag = reducebasis!(psupp, basis[1], blocks[1], cl[1], blocksize[1], nb=nb)
        if flag == 1
            ksupp = reduce(vcat, [p.supp for p in [obj; ineq_cons; eq_cons]])
            for item in basis[1]
                push!(ksupp, sadd(item, item, nb=nb))
            end
            sort!(ksupp)
            unique!(ksupp)
            blocks,cl,blocksize,eblocks = get_blocks(ineq_cons, eq_cons, ksupp, basis, ebasis, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md, signsymmetry=ss)
        end
    end
    end
    if QUIET == false
        mb = maximum(maximum.(blocksize))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,momone,GramMat,multiplier,SDP_status = solvesdp(obj, ineq_cons, eq_cons, n, basis, ebasis, blocks, eblocks, cl, blocksize, nb=nb, gb=gb, x=x, dualize=dualize, TS=TS,
    lead=leadsupp, QUIET=QUIET, solve=solve, solution=solution, MomentOne=MomentOne, Gram=Gram, mosek_setting=mosek_setting, writetofile=writetofile, signsymmetry=ss, normality=normality, model=model)
    data = pop_data(pop, obj, ineq_cons, eq_cons, x, n, nb, numeq, gb, leadsupp, basis, ebasis, ksupp, blocksize, blocks, eblocks, GramMat, multiplier, moment, SDP_status, rtol, gtol, ftol, 1)
    sol = nothing
    if solution == true
        if TS != false || !isempty(gb)
            sol,gap,data.flag = extract_solution(momone, opt, pop, x, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=true)
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
            sol = extract_solutions_robust(moment[1], n, d, pop=pop, x=x, lb=opt, numeq=numeq, basis=basis[1], check=true, rtol=rtol, gtol=gtol, ftol=ftol, QUIET=QUIET)[1]
            if sol !== nothing
                data.flag = 0
            end
        end
    end
    return opt,sol,data
end

"""
    opt,sol,data = tssos(f, x; newton=true, reducebasis=false, TS="block", merge=false, md=3, feasibility=false, QUIET=false, solve=true, 
    dualize=false, MomentOne=false, Gram=false, solution=false, normality=false, model=nothing, mosek_setting=mosek_para(), rtol=1e-2, gtol=1e-2, ftol=1e-3)

Compute the first TS step of the TSSOS hierarchy for unconstrained polynomial optimization.
If `newton=true`, then compute a monomial basis by the Newton polytope method.
If `reducebasis=true`, then remove monomials from the monomial basis by diagonal inconsistency.
If `TS="block"`, use maximal chordal extensions; if `TS="MD"`, use approximately smallest chordal extensions. 
If `merge=true`, perform the PSD block merging. 
If `feasibility=true`, then solve the feasibility problem.
If `solve=false`, then do not solve the SDP.
If `Gram=true`, then output the Gram matrix.
If `MomentOne=true`, add an extra first-order moment PSD constraint to the moment relaxation.

# Input arguments
- `f`: objective
- `x`: POP variables
- `TS`: type of term sparsity (`"block"`, `"signsymmetry"`, `"MD"`, `"MF"`, `false`)
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
function tssos(f::Poly{T}, x; newton=true, reducebasis=false, TS="block", merge=false, md=3, feasibility=false, QUIET=false, solve=true, 
    dualize=false, MomentOne=false, Gram=false, solution=false, mosek_setting=mosek_para(), model=nothing, rtol=1e-2, gtol=1e-2, ftol=1e-3) where {T<:Number}
    return tssos([f], x, Int(ceil(MP.maxdegree(f)/2)), newton=newton, feasibility=feasibility, GroebnerBasis=false, reducebasis=reducebasis, TS=TS, merge=merge, md=md, 
    QUIET=QUIET, solve=solve, dualize=dualize, MomentOne=MomentOne, Gram=Gram, solution=solution, rtol=rtol, gtol=gtol, ftol=ftol, mosek_setting=mosek_setting, model=model)
end

"""
    opt,sol,data = tssos(data; TS="block", merge=false, md=3, QUIET=false, solve=true, feasibility=false, dualize=false, Gram=false,
    MomentOne=false, normality=false, solution=false, mosek_setting=mosek_para(), model=nothing)

Compute higher TS steps of the TSSOS hierarchy.
"""
function tssos(data::pop_data; TS="block", merge=false, md=3, QUIET=false, solve=true, feasibility=false, dualize=false, MomentOne=false, Gram=false,
    solution=false, normality=false, mosek_setting=mosek_para(), model=nothing, writetofile=false)
    obj = data.obj
    ineq_cons = data.ineq_cons
    eq_cons = data.eq_cons
    n = data.n
    nb = data.nb
    numeq = data.numeq
    x = data.x
    gb = data.gb
    leadsupp = data.leadsupp
    basis = data.basis
    ebasis = data.ebasis
    ksupp = data.ksupp
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,eblocks = get_blocks(ineq_cons, eq_cons, ksupp, basis, ebasis, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md)
    end
    if blocksize == data.blocksize && eblocks == data.eblocks
        println("No higher TS step of the TSSOS hierarchy!")
        opt = sol = nothing
    else
        if QUIET == false
            mb = maximum(maximum.(blocksize))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,momone,GramMat,multiplier,SDP_status = solvesdp(obj, ineq_cons, eq_cons, n, basis, ebasis, blocks, eblocks, cl, blocksize, nb=nb, gb=gb, x=x, 
        lead=leadsupp, TS=TS, feasibility=feasibility, QUIET=QUIET, solve=solve, dualize=dualize, solution=solution, MomentOne=MomentOne, Gram=Gram, mosek_setting=mosek_setting, 
        model=model, normality=normality, writetofile=writetofile)
        sol = nothing
        if solution == true
            sol,gap,data.flag = extract_solution(momone, opt, data.pop, x, numeq=numeq, QUIET=true, gtol=data.gtol, ftol=data.ftol)
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

function get_csupp(ineq_cons::Vector{poly{T}}, eq_cons::Vector{poly{T}}, basis, ebasis, blocks, eblocks, cl, blocksize; nb=0) where {T<:Number}
    csupp = Vector{UInt16}[]
    for (k, p) in enumerate(ineq_cons[2:end]), i = 1:cl[k+1], j = 1:blocksize[k+1][i], r = j:blocksize[k+1][i], item in p.supp
        @inbounds bi = sadd(basis[k+1][blocks[k+1][i][j]], item, basis[k+1][blocks[k+1][i][r]], nb=nb)
        push!(csupp, bi)
    end
    for (k, p) in enumerate(eq_cons), i in eblocks[k], item in p.supp
        @inbounds push!(csupp, sadd(ebasis[k][i], item, nb=nb))
    end
    return csupp
end

function reducebasis!(supp, basis, blocks, cl, blocksize; nb=0)
    esupp = supp[map(item -> all(i -> iseven(count(==(i), item)), unique(item)), supp)]
    init,flag,check = 0,0,0
    while init == 0 || check > 0
        init,check = 1,0
        tsupp = esupp
        for i = 1:cl
            if blocksize[i] > 1
                for j = 1:blocksize[i], r = j+1:blocksize[i]
                    @inbounds bi = sadd(basis[blocks[i][j]], basis[blocks[i][r]], nb=nb)
                    push!(tsupp, bi)
                end
            end
        end
        unique!(tsupp)
        sort!(tsupp)
        for i = 1:cl
            lo = blocksize[i]
            indexb = Vector(1:lo)
            j = 1
            while lo >= j
                bi = sadd(basis[blocks[i][indexb[j]]], basis[blocks[i][indexb[j]]], nb=nb)
                Locb = bfind(tsupp, bi)
                if Locb === nothing
                   check,flag = 1,1
                   deleteat!(indexb, j)
                   lo -= 1
                else
                   j += 1
                end
            end
            blocks[i] = blocks[i][indexb]
            blocksize[i] = lo
        end
    end
    if flag == 1
       indexb = blocks[1]
       for i = 2:cl
           indexb = append!(indexb, blocks[i])
       end
       sort!(indexb)
       unique!(indexb)
       return basis[indexb],flag
    else
       return basis,flag
    end
end

function get_graph(tsupp, supp::Vector{Vector{UInt16}}, basis::Vector{Vector{UInt16}}; nb=0, signsymmetry=nothing)
    lb = length(basis)
    G = SimpleGraph(lb)
    for i = 1:lb, j = i+1:lb
        if signsymmetry === nothing
            ind = findfirst(item -> bfind(tsupp, sadd(basis[i], item, basis[j], nb=nb)) !== nothing, supp)
            if ind !== nothing
                add_edge!(G, i, j)
            end
        else
            bi = sadd(basis[i], supp[1], basis[j], nb=nb)
            sp = zeros(Int, size(signsymmetry, 1))
            sp[sign_type(bi)] .= 1
            if all(transpose(signsymmetry)*sp .== 0)
                add_edge!(G, i, j)
            end
        end
    end
    return G
end

function get_eblock(tsupp, hsupp, basis::Vector{Vector{UInt16}}; nb=0, signsymmetry=nothing)
    eblock = Int[]
    for (i, item) in enumerate(basis)
        if signsymmetry === nothing
            if findfirst(nitem -> bfind(tsupp, sadd(item, nitem, nb=nb)) !== nothing, hsupp) !== nothing
                push!(eblock, i)
            end
        else
            bi = sadd(item, hsupp[1], nb=nb)
            sp = zeros(Int, size(signsymmetry, 1))
            sp[sign_type(bi)] .= 1
            if all(transpose(signsymmetry)*sp .== 0)
                push!(eblock, i)
            end
        end
    end
    return eblock
end

function get_blocks(ineq_cons::Vector{poly{T}}, eq_cons::Vector{poly{T}}, tsupp, basis, ebasis; nb=0, TS="block", QUIET=true, merge=false, md=3, signsymmetry=nothing) where {T<:Number}
    blocks = Vector{Vector{Vector{Int}}}(undef, length(ineq_cons))
    blocksize = Vector{Vector{Int}}(undef, length(ineq_cons))
    cl = Vector{Int}(undef, length(ineq_cons))
    if TS == false
        for k = 1:length(ineq_cons)
            blocks[k],blocksize[k],cl[k] = [Vector(1:length(basis[k]))],[length(basis[k])],1
        end
        eblocks = [Vector(1:length(ebasis[k])) for k = 1:length(eq_cons)]
    else
        for k = 1:length(ineq_cons)
            G = get_graph(tsupp, ineq_cons[k].supp, basis[k], nb=nb, signsymmetry=signsymmetry)
            if TS == true || TS == "block" || TS == "signsymmetry"
                blocks[k] = connected_components(G)
                blocksize[k] = length.(blocks[k])
                cl[k] = length(blocksize[k])
            else
                blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS)
                if merge == true
                    blocks[k],cl[k],blocksize[k] = clique_merge!(blocks[k], d=md, QUIET=true)
                end
            end
            if QUIET == false
                sb = sort(unique(blocksize[k]), rev=true)
                numb = [sum(blocksize[k].== i) for i in sb]
                k0 = k - 1
                println("-----------------------------------------------------------------------------")
                println("The sizes of PSD blocks for the $k0-th SOS multiplier:\n$sb\n$numb")
                println("-----------------------------------------------------------------------------")
            end
        end
        eblocks = Vector{Vector{Int}}(undef, length(eq_cons))
        for k = 1:length(eq_cons)
            eblocks[k] = get_eblock(tsupp, eq_cons[k].supp, ebasis[k], nb=nb, signsymmetry=signsymmetry)
        end
    end
    return blocks,cl,blocksize,eblocks
end

function solvesdp(obj, ineq_cons::Vector{poly{T}}, eq_cons::Vector{poly{T}}, n, basis, ebasis, blocks, eblocks, cl, blocksize; nb=0, gb=[], x=[], lead=[], TS="block",
    QUIET=true, solve=true, feasibility=false, dualize=false, solution=false, MomentOne=false, Gram=false, mosek_setting=mosek_para(), signsymmetry=false, writetofile=false, 
    normality=false, model=nothing) where {T<:Number}
    ksupp = Vector{UInt16}[]
    for i = 1:cl[1], j = 1:blocksize[1][i], r = j:blocksize[1][i]
        @inbounds bi = sadd(basis[1][blocks[1][i][j]], basis[1][blocks[1][i][r]], nb=nb)
        push!(ksupp, bi)
    end
    if TS != false && TS != "signsymmetry"
        csupp = get_csupp(ineq_cons, eq_cons, basis, ebasis, blocks, eblocks, cl, blocksize, nb=nb)
        append!(ksupp, csupp)
    end
    sort!(ksupp)
    unique!(ksupp)
    objv = moment = momone = GramMat = multiplier = SDP_status = nothing
    if solve == true
        tsupp = copy(ksupp)
        if normality == true
            wbasis = basis[1]
            bs = length(wbasis)  
            hyblocks = Vector{Vector{Vector{Int}}}(undef, n)
            for i = 1:n
                G = SimpleGraph(2bs)
                for j = 1:bs, k = j:bs
                    bi = sadd(wbasis[j], wbasis[k], nb=nb)
                    sp = zeros(Int, size(signsymmetry, 1))
                    sp[sign_type(bi)] .= 1
                    if all(transpose(signsymmetry)*sp .== 0)
                        add_edge!(G, j, k)
                    end
                    bi = sadd(wbasis[j], wbasis[k], UInt16[i;i], nb=nb)
                    sp = zeros(Int, size(signsymmetry, 1))
                    sp[sign_type(bi)] .= 1
                    if all(transpose(signsymmetry)*sp .== 0)
                        add_edge!(G, j+bs, k+bs)
                    end
                    bi = sadd(wbasis[j], wbasis[k], UInt16[i], nb=nb)
                    sp = zeros(Int, size(signsymmetry, 1))
                    sp[sign_type(bi)] .= 1
                    if all(transpose(signsymmetry)*sp .== 0)
                        add_edge!(G, j, k+bs)
                    end
                end
                hyblocks[i] = connected_components(G)
                for block in hyblocks[i], j = 1:length(block), k = j:length(block)
                    if block[j] <= bs && block[k] > bs
                        bi = sadd(wbasis[block[j]], wbasis[block[k]-bs], UInt16[i], nb=nb)
                        push!(tsupp, bi)
                    elseif block[j] > bs
                        bi = sadd(wbasis[block[j]-bs], wbasis[block[k]-bs], UInt16[i;i], nb=nb)
                        push!(tsupp, bi)
                    end
                end
            end
        end
        if (MomentOne == true || solution == true) && TS != false
            append!(tsupp, get_basis(n, 2, nb=nb, lead=lead))
        end
        if !isempty(gb)
            unique!(tsupp)
            nsupp = Vector{UInt16}[]
            for item in tsupp
                if divide(item, lead)
                    append!(nsupp, reminder(item, x, gb).supp)
                else
                    push!(nsupp, item)
                end
            end
            tsupp = nsupp
        end
        sort!(tsupp)
        unique!(tsupp)
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
        if normality == true
            for i = 1:n, block in hyblocks[i]
                hnom = @variable(model, [1:length(block), 1:length(block)], PSD)
                for j = 1:length(block), k = j:length(block)
                    if block[k] <= bs
                        bi = sadd(wbasis[block[j]], wbasis[block[k]], nb=nb)
                    elseif block[j] <= bs && block[k] > bs
                        bi = sadd(wbasis[block[j]], wbasis[block[k]-bs], UInt16[i], nb=nb)
                    else
                        bi = sadd(wbasis[block[j]-bs], wbasis[block[k]-bs], UInt16[i;i], nb=nb)
                    end
                    if !isempty(gb) && divide(bi, lead)
                        rem = reminder(bi, x, gb)
                        for (l, item) in enumerate(rem.supp)
                            Locb = bfind(tsupp, item)
                            if j == k
                                @inbounds add_to_expression!(cons[Locb], rem.coe[l], hnom[j,k])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*rem.coe[l], hnom[j,k])
                            end
                        end
                    else
                        Locb = bfind(tsupp, bi)
                        if j == k
                            @inbounds add_to_expression!(cons[Locb], hnom[j,k])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2, hnom[j,k])
                        end
                    end                
                end
            end
        end
        if (MomentOne == true || solution == true) && TS != false
            pos0 = @variable(model, [1:n+1, 1:n+1], PSD)
            bas = [[UInt16[]]; [UInt16[i] for i = 1:n]]
            for j = 1:n+1, k = j:n+1
                @inbounds bi = sadd(bas[j], bas[k], nb=nb)
                if !isempty(gb) && divide(bi, lead)
                    rem = reminder(bi, x, gb)
                    for (l, item) in enumerate(rem.supp)
                        Locb = bfind(tsupp, item)
                        if j == k
                            @inbounds add_to_expression!(cons[Locb], rem.coe[l], pos0[j,k])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2*rem.coe[l], pos0[j,k])
                        end
                    end
                else
                    Locb = bfind(tsupp, bi)
                    if j == k
                        @inbounds add_to_expression!(cons[Locb], pos0[j,k])
                    else
                        @inbounds add_to_expression!(cons[Locb], 2, pos0[j,k])
                    end
                end
            end
        end
        pos = Vector{Vector{Symmetric{VariableRef}}}(undef, length(ineq_cons))
        for (k, p) in enumerate(ineq_cons)
            pos[k] = Vector{Symmetric{VariableRef}}(undef, cl[k])
            for i = 1:cl[k]
                pos[k][i] = @variable(model, [1:blocksize[k][i], 1:blocksize[k][i]], PSD)
                for j = 1:blocksize[k][i], r = j:blocksize[k][i], (s, it) in enumerate(p.supp)
                    @inbounds bi = sadd(basis[k][blocks[k][i][j]], it, basis[k][blocks[k][i][r]], nb=nb)
                    if !isempty(gb) && divide(bi, lead)
                        rem = reminder(bi, x, gb)
                        for (l, item) in enumerate(rem.supp)
                            Locb = bfind(tsupp, item)
                            if j == r
                                @inbounds add_to_expression!(cons[Locb], p.coe[s]*rem.coe[l], pos[k][i][j,r])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*p.coe[s]*rem.coe[l], pos[k][i][j,r])
                            end
                        end
                    else
                        Locb = bfind(tsupp, bi)
                        if j == r
                            @inbounds add_to_expression!(cons[Locb], p.coe[s], pos[k][i][j,r])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2*p.coe[s], pos[k][i][j,r])
                        end
                    end
                end
            end
        end
        free = Vector{Vector{VariableRef}}(undef, length(eq_cons))
        for (k, p) in enumerate(eq_cons)
            free[k] = @variable(model, [1:length(eblocks[k])])
            for (i, j) in enumerate(eblocks[k]), (s, item) in enumerate(p.supp)
                @inbounds bi = sadd(ebasis[k][j], item, nb=nb)
                Locb = bfind(tsupp, bi)
                @inbounds add_to_expression!(cons[Locb], p.coe[s], free[k][i])
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
        if feasibility == false
            @variable(model, lower)
            cons[1] += lower
            @objective(model, Max, lower)
        else
            tnorm = AffExpr(0)
            for i = 1:cl[1]
                if blocksize[1][i] == 1
                    tnorm += pos[1][i]
                else
                    tnorm += tr(pos[1][i])
                end
            end
            @objective(model, Min, tnorm)
        end
        @constraint(model, con, cons==zeros(length(cons)))
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
                GramMat = [[value.(pos[k][i]) for i = 1:cl[k]] for k = 1:length(ineq_cons)]
                multiplier = [value.(free[j]) for j = 1:length(eq_cons)]
            end
            dual_var = -dual(con)
            moment = Vector{Symmetric{Float64}}(undef, cl[1])
            for i = 1:cl[1]
                mmat = zeros(blocksize[1][i], blocksize[1][i])
                for j = 1:blocksize[1][i], k = j:blocksize[1][i]
                    bi = sadd(basis[1][blocks[1][i][j]], basis[1][blocks[1][i][k]], nb=nb)
                    if !isempty(gb) && divide(bi, lead)
                        rem = reminder(bi, x, gb)
                        for (l, item) in enumerate(rem.supp)
                            Locb = bfind(tsupp, item)
                            mmat[j,k] += rem.coe[l]*dual_var[Locb]
                        end
                    else
                        Locb = bfind(tsupp, bi)
                        mmat[j,k] = dual_var[Locb]
                    end
                end
                moment[i] = Symmetric(mmat, :U)
            end
            if solution == true && (TS != false || !isempty(gb))
                momone = zeros(Float64, n+1, n+1)
                bas = [[UInt16[]]; [UInt16[i] for i = 1:n]]
                for j = 1:n+1, k = j:n+1
                    bi = sadd(bas[j], bas[k], nb=nb)
                    if !isempty(gb) && divide(bi, lead)
                        rem = reminder(bi, x, gb)
                        for (l, item) in enumerate(rem.supp)
                            Locb = bfind(tsupp, item)
                            momone[j,k] += rem.coe[l]*dual_var[Locb]
                        end
                    else
                        Locb = bfind(tsupp, bi)
                        momone[j,k] = dual_var[Locb]
                    end
                end
                momone = Symmetric(momone, :U)
            end
        end
    end
    return objv,ksupp,moment,momone,GramMat,multiplier,SDP_status
end
