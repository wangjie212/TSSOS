mutable struct cpop_data
    n # number of variables
    nb # number of binary variables
    m # number of constraints
    numeq # number of equality constraints
    x # set of variables
    pop # polynomial optimization problem
    gb # Grobner basis
    leadsupp # leader terms of the Grobner basis
    supp # support data
    coe # coefficient data
    basis # monomial bases
    ksupp # extended support at the k-th step
    sb # sizes of different blocks
    numb # numbers of different blocks
    blocks # block structure
    cl # numbers of blocks
    blocksize # sizes of blocks
    GramMat # Gram matrix
    moment # Moment matrix
    solver # SDP solver
    SDP_status
    tol # tolerance to certify global optimality
    flag # 0 if global optimality is certified; 1 otherwise
end

"""
    opt,sol,data = tssos_first(pop, x, d; nb=0, numeq=0, quotient=true, basis=[],
    reducebasis=false, TS="block", merge=false, md=3, solver="Mosek", QUIET=false, solve=true,
    MomentOne=false, Gram=false, solution=false, tol=1e-4)

Compute the first TS step of the TSSOS hierarchy for constrained polynomial optimization.
If `quotient=true`, then exploit the quotient ring structure defined by the equality constraints.
If `merge=true`, perform the PSD block merging. 
If `solve=false`, then do not solve the SDP.
If `Gram=true`, then output the Gram matrix.
If `MomentOne=true`, add an extra first order moment matrix to the moment relaxation.

# Input arguments
- `pop`: vector of the objective, inequality constraints, and equality constraints
- `x`: POP variables
- `d`: relaxation order
- `nb`: number of binary variables in `x`
- `numeq`: number of equality constraints
- `TS`: type of term sparsity (`"block"`, `"MD"`, `"MF"`, `false`)
- `md`: tunable parameter for merging blocks
- `QUIET`: run in the quiet mode or not (`true`, `false`)
- `tol`: relative tolerance to certify global optimality

# Output arguments
- `opt`: optimum
- `sol`: (near) optimal solution (if `solution=true`)
- `data`: other auxiliary data 
"""
function tssos_first(pop, x, d; nb=0, numeq=0, quotient=true, basis=[], reducebasis=false, TS="block", merge=false, md=3, solver="Mosek", 
    QUIET=false, solve=true, dualize=false, MomentOne=false, Gram=false, solution=false, tol=1e-4, cosmo_setting=cosmo_para())
    println("*********************************** TSSOS ***********************************")
    println("Version 1.0.0, developed by Jie Wang, 2020--2023")
    println("TSSOS is launching...")
    n = length(x)
    if nb > 0
        gb = x[1:nb].^2 .- 1
        for i in eachindex(pop)
            pop[i] = rem(pop[i], gb)
        end
    end
    if numeq > 0 && quotient == true
        cpop = copy(pop)
        gb = convert.(Polynomial{true,Float64}, cpop[end-numeq+1:end])
        cpop = cpop[1:end-numeq]
        if QUIET == false
            println("Starting to compute the Gröbner basis...")
            println("This might take much time. You can set quotient=false to close it.")
        end
        SemialgebraicSets.gröbnerbasis!(gb)
        cpop[1] = rem(cpop[1], gb)
        lead = leadingmonomial.(gb)
        llead = length(lead)
        leadsupp = zeros(UInt8, n, llead)
        for i = 1:llead, j = 1:n
            @inbounds leadsupp[j,i] = MultivariatePolynomials.degree(lead[i], x[j])
        end
    else
        cpop = pop
        gb = []
        leadsupp = []
    end
    m = length(cpop)-1
    dg = zeros(Int, m)
    coe = Vector{Vector{Float64}}(undef, m+1)
    supp = Vector{Array{UInt8,2}}(undef, m+1)
    for k = 1:m+1
        mons = monomials(cpop[k])
        coe[k] = coefficients(cpop[k])
        supp[k] = zeros(UInt8,n,length(mons))
        for i in eachindex(mons), j = 1:n
            @inbounds supp[k][j,i] = MultivariatePolynomials.degree(mons[i],x[j])
        end
    end
    isupp = supp[1]
    for i = 2:m+1
        dg[i-1] = maxdegree(pop[i])
        isupp = [isupp supp[i]]
    end
    if basis == []
        basis = Vector{Array{UInt8,2}}(undef, m+1)
        basis[1] = get_basis(n, d, nb=nb, lead=leadsupp)
        for k = 1:m
            basis[k+1] = get_basis(n, d-Int(ceil(dg[k]/2)), nb=nb, lead=leadsupp)
        end
    end
    tsupp = [isupp bin_add(basis[1],basis[1],nb)]
    tsupp = sortslices(tsupp, dims=2)
    tsupp = unique(tsupp, dims=2)
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,sb,numb,_ = get_cblocks(m, tsupp, supp[2:end], basis, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md)
    if reducebasis == true && quotient == false
        gsupp = get_gsupp(n, m, supp, basis[2:end], blocks[2:end], cl[2:end], blocksize[2:end], nb=nb)
        psupp = [supp[1] zeros(UInt8,n)]
        psupp = [psupp gsupp]
        basis[1],flag = reducebasis!(psupp, basis[1], blocks[1], cl[1], blocksize[1], nb=nb)
        if flag == 1
            tsupp = [isupp bin_add(basis[1],basis[1],nb)]
            tsupp = sortslices(tsupp, dims=2)
            tsupp = unique(tsupp, dims=2)
            blocks,cl,blocksize,sb,numb,_ = get_cblocks(m, tsupp, supp[2:end], basis, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md)
        end
    end
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,momone,GramMat,SDP_status = blockcpop(n, m, supp, coe, basis, blocks, cl, blocksize, nb=nb, numeq=numeq, gb=gb, x=x, dualize=dualize,
    lead=leadsupp, solver=solver, QUIET=QUIET, solve=solve, solution=solution, MomentOne=MomentOne, Gram=Gram, cosmo_setting=cosmo_setting)
    data = cpop_data(n, nb, m, numeq, x, pop, gb, leadsupp, supp, coe, basis, ksupp, sb, numb, blocks, cl, blocksize, GramMat, moment, solver, SDP_status, tol, 1)
    sol = nothing
    if solution == true
        sol,gap,data.flag = extract_solution(momone, opt, pop, x, numeq=numeq, tol=tol)
        if data.flag == 1
            sol = gap > 0.5 ? randn(n) : sol
            sol,data.flag = refine_sol(opt, sol, data, QUIET=true, tol=tol)
        end
    end
    return opt,sol,data
end

function tssos_higher!(data::cpop_data; TS="block", merge=false, md=3, QUIET=false, solve=true, dualize=false, MomentOne=false, Gram=false,
    solution=false, cosmo_setting=cosmo_para())
    n = data.n
    nb = data.nb
    m = data.m
    numeq = data.numeq
    x = data.x
    pop = data.pop
    gb = data.gb
    leadsupp = data.leadsupp
    supp = data.supp
    coe = data.coe
    basis = data.basis
    ksupp = data.ksupp
    sb = data.sb
    numb = data.numb
    blocks = data.blocks
    cl = data.cl
    blocksize = data.blocksize
    solver = data.solver
    tol = data.tol
    ksupp = sortslices(ksupp, dims=2)
    ksupp = unique(ksupp, dims=2)
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,sb,numb,status = get_cblocks(m, ksupp, supp[2:end], basis, blocks=blocks, cl=cl, blocksize=blocksize, sb=sb, numb=numb, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md)
    end
    opt = nothing
    sol = nothing
    if status == 1
        if QUIET == false
            mb = maximum(maximum.(sb))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,momone,GramMat,SDP_status = blockcpop(n, m, supp, coe, basis, blocks, cl, blocksize, nb=nb, numeq=numeq, gb=gb, x=x, lead=leadsupp,
        solver=solver, QUIET=QUIET, solve=solve, dualize=dualize, solution=solution, MomentOne=MomentOne, Gram=Gram, cosmo_setting=cosmo_setting)
        if solution == true
            sol,gap,data.flag = extract_solution(momone, opt, pop, x, numeq=numeq, tol=tol)
            if data.flag == 1
                sol = gap > 0.5 ? randn(n) : sol
                sol,data.flag = refine_sol(opt, sol, data, QUIET=true, tol=tol)
            end
        end
        data.ksupp = ksupp
        data.sb = sb
        data.numb = numb
        data.blocks = blocks
        data.cl = cl
        data.blocksize = blocksize
        data.GramMat = GramMat
        data.moment = moment
        data.SDP_status = SDP_status
    end
    return opt,sol,data
end

function get_gsupp(n, m, supp, gbasis, blocks, cl, blocksize; nb=0)
    gsupp = zeros(UInt8, n, sum(size(supp[k+1],2)*Int(sum(Int.(blocksize[k]).^2+blocksize[k])/2) for k=1:m))
    l = 1
    for k = 1:m, i = 1:cl[k], j = 1:blocksize[k][i], r = j:blocksize[k][i], s = 1:size(supp[k+1],2)
        @inbounds bi = bin_add(gbasis[k][:,blocks[k][i][j]], gbasis[k][:,blocks[k][i][r]], nb)
        @inbounds bi = bin_add(bi, supp[k+1][:,s], nb)
        @inbounds gsupp[:,l] = bi
        l += 1
    end
    return gsupp
end

function reducebasis!(supp, basis, blocks, cl, blocksize; nb=0)
    esupp = supp[:, all.(iseven, eachcol(supp))]
    init,flag,check = 0,0,0
    while init==0 || check>0
        init,check = 1,0
        tsupp = esupp
        for i = 1:cl
            if blocksize[i] > 1
                for j = 1:blocksize[i], r = j+1:blocksize[i]
                    @inbounds bi = bin_add(basis[:,blocks[i][j]], basis[:,blocks[i][r]], nb)
                    tsupp = [tsupp bi]
                end
            end
        end
        tsupp = unique(tsupp, dims=2)
        tsupp = sortslices(tsupp, dims=2)
        ltsupp = size(tsupp, 2)
        for i = 1:cl
            lo = blocksize[i]
            indexb = [k for k=1:lo]
            j = 1
            while lo >= j
                bi = bin_add(basis[:,blocks[i][indexb[j]]], basis[:,blocks[i][indexb[j]]], nb)
                Locb = bfind(tsupp, ltsupp, bi)
                if Locb == 0
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
       return basis[:,indexb],flag
    else
       return basis,flag
    end
end

function get_cgraph(tsupp::Array{UInt8, 2}, supp::Array{UInt8, 2}, basis::Array{UInt8, 2}; nb=0, balanced=false)
    lb = size(basis, 2)
    G = SimpleGraph(lb)
    ltsupp = size(tsupp, 2)
    for i = 1:lb, j = i+1:lb
        r = 1
        while r <= size(supp, 2)
            bi = bin_add(basis[:,i], basis[:,j], nb)
            bi = bin_add(bi, supp[:,r], nb)
            if bfind(tsupp, ltsupp, bi) != 0
               break
            else
                r += 1
            end
        end
        if r <= size(supp, 2)
           add_edge!(G, i, j)
        end
    end
    return G
end

function get_cblocks(m, tsupp, supp, basis; blocks=[], cl=[], blocksize=[], sb=[], numb=[], nb=0,
    TS="block", balanced=false, QUIET=true, merge=false, md=3)
    if isempty(blocks)
        blocks = Vector{Vector{Vector{UInt16}}}(undef, m+1)
        blocksize = Vector{Vector{UInt16}}(undef, m+1)
        cl = Vector{UInt16}(undef, m+1)
    end
    if TS == false
        for k = 1:m+1
            lb = ndims(basis[k])==1 ? length(basis[k]) : size(basis[k],2)
            blocks[k] = [[i for i=1:lb]]
            blocksize[k] = [lb]
            cl[k] = 1
        end
        status = 1
        nsb = Int.(blocksize[1])
        nnumb = [1]
    else
        G = get_graph(tsupp, basis[1], nb=nb, balanced=balanced)
        if TS == "block"
            blocks[1] = connected_components(G)
            blocksize[1] = length.(blocks[1])
            cl[1] = length(blocksize[1])
        else
            blocks[1],cl[1],blocksize[1] = chordal_cliques!(G, method=TS, minimize=false)
            if merge == true
                blocks[1],cl[1],blocksize[1] = clique_merge!(blocks[1], d=md, QUIET=true)
            end
        end
        nsb = sort(Int.(unique(blocksize[1])), rev=true)
        nnumb = [sum(blocksize[1].== i) for i in nsb]
        if isempty(sb) || nsb!=sb || nnumb!=numb
            status = 1
            if QUIET == false
                println("-----------------------------------------------------------------------------")
                println("The sizes of PSD blocks:\n$nsb\n$nnumb")
                println("-----------------------------------------------------------------------------")
            end
            for k = 1:m
                G = get_cgraph(tsupp, supp[k], basis[k+1], nb=nb, balanced=balanced)
                if TS == "block"
                    blocks[k+1] = connected_components(G)
                    blocksize[k+1] = length.(blocks[k+1])
                    cl[k+1] = length(blocksize[k+1])
                else
                    blocks[k+1],cl[k+1],blocksize[k+1] = chordal_cliques!(G, method=TS, minimize=false)
                    if merge == true
                        blocks[k+1],cl[k+1],blocksize[k+1] = clique_merge!(blocks[k+1], d=md, QUIET=true)
                    end
                end
            end
        else
            status = 0
            if QUIET == false
                println("No higher TS step of the TSSOS hierarchy!")
            end
        end
    end
    return blocks,cl,blocksize,nsb,nnumb,status
end

function blockcpop(n, m, supp, coe, basis, blocks, cl, blocksize; nb=0, numeq=0, gb=[], x=[], lead=[], solver="Mosek", 
    QUIET=true, solve=true, dualize=false, solution=false, MomentOne=false, Gram=false, cosmo_setting=cosmo_para())
    ksupp = zeros(UInt8, n, Int(sum(Int.(blocksize[1]).^2+blocksize[1])/2))
    k = 1
    for i = 1:cl[1], j = 1:blocksize[1][i], r = j:blocksize[1][i]
        @inbounds bi = bin_add(basis[1][:,blocks[1][i][j]], basis[1][:,blocks[1][i][r]], nb)
        @inbounds ksupp[:,k] = bi
        k += 1
    end
    objv = moment = momone = GramMat = SDP_status = nothing
    if solve == true
        tsupp = ksupp
        if m > 0
            gsupp = get_gsupp(n, m, supp, basis[2:end], blocks[2:end], cl[2:end], blocksize[2:end], nb=nb)
            tsupp = [tsupp gsupp]
        end
        if MomentOne == true || solution == true
            tsupp = [tsupp get_basis(n, 2, nb=nb)]
        end
        if !isempty(gb)
            tsupp = unique(tsupp, dims=2)
            nsupp = zeros(UInt8, n)
            llead = size(lead, 2)
            for col in eachcol(tsupp)
                if divide(col, lead, n, llead)
                    temp = reminder(col, x, gb, n)[2]
                    nsupp = [nsupp temp]
                else
                    nsupp = [nsupp col]
                end
            end
            tsupp = nsupp
        end
        tsupp = sortslices(tsupp, dims=2)
        tsupp = unique(tsupp, dims=2)
        ltsupp = size(tsupp, 2)
        if QUIET == false
            println("Assembling the SDP...")
            println("There are $ltsupp affine constraints.")
        end
        if solver == "Mosek"
            if dualize == false
                model = Model(optimizer_with_attributes(Mosek.Optimizer))
            else
                model = Model(dual_optimizer(Mosek.Optimizer))
            end
        elseif solver == "COSMO"
            model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter))
        elseif solver == "SDPT3"
            model = Model(optimizer_with_attributes(SDPT3.Optimizer))
        elseif solver == "SDPNAL"
            model = Model(optimizer_with_attributes(SDPNAL.Optimizer))
        else
            @error "The solver is currently not supported!"
            return nothing,nothing,nothing,nothing,nothing,nothing
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = [AffExpr(0) for i=1:ltsupp]
        pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1])
        for i = 1:cl[1]
            if MomentOne == true || solution == true
                pos0 = @variable(model, [1:n+1, 1:n+1], PSD)
                for j = 1:n+1, k = j:n+1
                    @inbounds bi = bin_add(basis[1][:,j], basis[1][:,k], nb)
                    if !isempty(gb) && divide(bi, lead, n, llead)
                        bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                        for l = 1:bi_lm
                            Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                            if j == k
                               @inbounds add_to_expression!(cons[Locb], bi_coe[l], pos0[j,k])
                            else
                               @inbounds add_to_expression!(cons[Locb], 2*bi_coe[l], pos0[j,k])
                            end
                        end
                    else
                        Locb = bfind(tsupp,ltsupp,bi)
                        if j == k
                           @inbounds add_to_expression!(cons[Locb], pos0[j,k])
                        else
                           @inbounds add_to_expression!(cons[Locb], 2, pos0[j,k])
                        end
                    end
                end
            end
            bs = blocksize[1][i]
            if bs == 1
               @inbounds pos[i] = @variable(model, lower_bound=0)
               @inbounds bi = bin_add(basis[1][:,blocks[1][i][1]], basis[1][:,blocks[1][i][1]], nb)
               if !isempty(gb) && divide(bi, lead, n, llead)
                   bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                   for l = 1:bi_lm
                       Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                       @inbounds add_to_expression!(cons[Locb], bi_coe[l], pos[i])
                   end
               else
                   Locb = bfind(tsupp,ltsupp,bi)
                   @inbounds add_to_expression!(cons[Locb], pos[i])
               end
            else
               @inbounds pos[i] = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:bs, r = j:bs
                   @inbounds bi = bin_add(basis[1][:,blocks[1][i][j]], basis[1][:,blocks[1][i][r]], nb)
                   if !isempty(gb) && divide(bi, lead, n, llead)
                       bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                       for l = 1:bi_lm
                           Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                           if j == r
                              @inbounds add_to_expression!(cons[Locb], bi_coe[l], pos[i][j,r])
                           else
                              @inbounds add_to_expression!(cons[Locb], 2*bi_coe[l], pos[i][j,r])
                           end
                       end
                   else
                       Locb = bfind(tsupp, ltsupp, bi)
                       if j == r
                          @inbounds add_to_expression!(cons[Locb], pos[i][j,r])
                       else
                          @inbounds add_to_expression!(cons[Locb], 2, pos[i][j,r])
                       end
                   end
               end
            end
        end
        gpos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m)
        for k = 1:m
            gpos[k] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k+1])
            for i = 1:cl[k+1]
                bs = blocksize[k+1][i]
                if bs == 1
                    if !isempty(gb) || k <= m-numeq
                        gpos[k][i] = @variable(model, lower_bound=0)
                    else
                        gpos[k][i] = @variable(model)
                    end
                    for s = 1:size(supp[k+1],2)
                        @inbounds bi = bin_add(basis[k+1][:,blocks[k+1][i][1]], basis[k+1][:,blocks[k+1][i][1]], nb)
                        @inbounds bi = bin_add(bi, supp[k+1][:,s], nb)
                        if !isempty(gb) && divide(bi, lead, n, llead)
                            bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                            for l = 1:bi_lm
                                Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                                @inbounds add_to_expression!(cons[Locb], coe[k+1][s]*bi_coe[l], gpos[k][i])
                            end
                        else
                            Locb = bfind(tsupp, ltsupp, bi)
                            @inbounds add_to_expression!(cons[Locb], coe[k+1][s], gpos[k][i])
                        end
                    end
                else
                    if !isempty(gb) || k <= m-numeq
                       gpos[k][i] = @variable(model, [1:bs, 1:bs], PSD)
                    else
                       gpos[k][i] = @variable(model, [1:bs, 1:bs], Symmetric)
                    end
                    for j = 1:bs, r = j:bs, s = 1:size(supp[k+1],2)
                        @inbounds bi = bin_add(basis[k+1][:,blocks[k+1][i][j]], basis[k+1][:,blocks[k+1][i][r]], nb)
                        @inbounds bi = bin_add(bi,supp[k+1][:,s], nb)
                        if !isempty(gb) && divide(bi, lead, n, llead)
                            bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                            for l = 1:bi_lm
                                Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                                if j == r
                                   @inbounds add_to_expression!(cons[Locb], coe[k+1][s]*bi_coe[l], gpos[k][i][j,r])
                                else
                                   @inbounds add_to_expression!(cons[Locb], 2*coe[k+1][s]*bi_coe[l], gpos[k][i][j,r])
                                end
                            end
                        else
                            Locb = bfind(tsupp, ltsupp, bi)
                            if j == r
                               @inbounds add_to_expression!(cons[Locb], coe[k+1][s], gpos[k][i][j,r])
                            else
                               @inbounds add_to_expression!(cons[Locb], 2*coe[k+1][s], gpos[k][i][j,r])
                            end
                        end
                    end
                end
            end
        end
        bc = zeros(ltsupp)
        for i = 1:size(supp[1],2)
            Locb = bfind(tsupp, ltsupp, supp[1][:,i])
            if Locb == 0
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing,nothing,nothing,nothing
            else
               bc[Locb] = coe[1][i]
           end
        end
        @variable(model, lower)
        cons[1] += lower
        @constraint(model, con[i=1:ltsupp], cons[i]==bc[i])
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
        SDP_status = termination_status(model)
        objv = objective_value(model)
        if SDP_status != MOI.OPTIMAL
           println("termination status: $SDP_status")
           status = primal_status(model)
           println("solution status: $status")
        end
        println("optimum = $objv")
        if Gram == true
            GramMat = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, m+1)
            GramMat[1] = [value.(pos[i]) for i = 1:cl[1]]
            for k = 1:m
                GramMat[k+1] = [value.(gpos[k][i]) for i = 1:cl[k+1]]
            end
        end
        dual_var = -dual.(con)
        moment = Vector{Matrix{Float64}}(undef, cl[1])
        for i = 1:cl[1]
            moment[i] = zeros(blocksize[1][i],blocksize[1][i])
            for j = 1:blocksize[1][i], k = j:blocksize[1][i]
                bi = bin_add(basis[1][:,blocks[1][i][j]], basis[1][:,blocks[1][i][k]], nb)
                if !isempty(gb) && divide(bi, lead, n, llead)
                    bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                    moment[i][j,k] = 0
                    for l = 1:bi_lm
                        Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                        moment[i][j,k] += bi_coe[l]*dual_var[Locb]
                    end
                else
                    Locb = bfind(tsupp, ltsupp, bi)
                    moment[i][j,k] = dual_var[Locb]
                end
            end
            moment[i] = Symmetric(moment[i],:U)
        end
        if solution == true
            momone = zeros(Float64, n+1, n+1)
            for j = 1:n+1, k = j:n+1
                bi = bin_add(basis[1][:,j], basis[1][:,k], nb)
                if !isempty(gb) && divide(bi, lead, n, llead)
                    bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                    momone[j,k] = 0
                    for l = 1:bi_lm
                        Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                        momone[j,k] += bi_coe[l]*dual_var[Locb]
                    end
                else
                    Locb = bfind(tsupp, ltsupp, bi)
                    momone[j,k] = dual_var[Locb]
                end
            end
            momone = Symmetric(momone,:U)
        end
    end
    return objv,ksupp,moment,momone,GramMat,SDP_status
end
