mutable struct upop_data
    n # number of variables
    nb # number of binary variables
    x # set of variables
    f # objective function
    supp # support data
    coe # coefficient data
    basis # monomial basis
    ksupp # extending support at the k-th step
    blocks # the block structrue
    sb # sizes of different blocks
    numb # numbers of different blocks
    GramMat # Gram matrix
    moment # Moment matrix
    solver # SDP solver
    tol # tolerance to certify global optimality
    flag # 0 if global optimality is certified; 1 otherwise
end

mutable struct cosmo_para
    eps_abs::Float64
    eps_rel::Float64
    max_iter::Int64
end

cosmo_para() = cosmo_para(1e-5, 1e-5, 1e4)

"""
    opt,sol,data = tssos_first(f, x; nb=0, newton=true, reducebasis=false, TS="block", merge=false,
    md=3, solver="Mosek", QUIET=false, solve=true, MomentOne=false, Gram=false, solution=false, tol=1e-4)

Compute the first step of the TSSOS hierarchy for unconstrained polynomial optimization.
If `newton=true`, then compute a monomial basis by the Newton polytope method.
If `reducebasis=true`, then remove monomials from the monomial basis by diagonal inconsistency.
If `TS="block"`, use maximal chordal extensions; if `TS="MD"`, use approximately smallest chordal
extensions. If `merge=true`, perform the PSD block merging. If `MomentOne=true`, add an extra first
order moment matrix to the moment relaxation.
Return the optimum, the (near) optimal solution (if `solution=true`) and other auxiliary data.

# Arguments
- `f`: the objective function for unconstrained polynomial optimization.
- `x`: the set of variables.
- `nb`: the number of binary variables in `x`.
- `md`: the tunable parameter for merging blocks.
- `tol`: the relative tolerance to certify global optimality.
"""
function tssos_first(f, x; nb=0, order=0, newton=true, reducebasis=false, TS="block", merge=false, feasible=false, md=3, solver="Mosek", 
    QUIET=false, solve=true, dualize=false, MomentOne=false, Gram=false, solution=false, tol=1e-4, cosmo_setting=cosmo_para())
    println("*********************************** TSSOS ***********************************")
    println("Version 1.0.0, developed by Jie Wang, 2020--2023")
    println("TSSOS is launching...")
    n = length(x)
    if nb > 0
        f = rem(f, x[1:nb].^2 .- 1)
    end
    mon = monomials(f)
    coe = coefficients(f)
    lm = length(mon)
    supp = zeros(UInt8, n, lm)
    for i = 1:lm, j = 1:n
        @inbounds supp[j,i] = MultivariatePolynomials.degree(mon[i], x[j])
    end
    if order == 0
        d = Int(maxdegree(f)/2)
    else
        d = order
    end
    if order == 0 && newton == true
       if sum(supp[:,end]) != 0 && feasible == false
          supp = [supp zeros(UInt8, n)]
          coe = [coe; 0]
       end
       basis = newton_basis(n, d, supp, solver=solver)
    else
       basis = get_basis(n, d, nb=nb)
    end
    tsupp = [supp bin_add(basis, basis, nb)]
    tsupp = sortslices(tsupp, dims=2)
    tsupp = unique(tsupp, dims=2)
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    blocks,cl,blocksize,sb,numb,_ = get_blocks(tsupp, basis, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md)
    if reducebasis == true
        psupp = [supp zeros(UInt8, n)]
        basis,flag = reducebasis!(psupp, basis, blocks, cl, blocksize, nb=nb)
        if flag == 1
            tsupp = [supp bin_add(basis, basis, nb)]
            tsupp = sortslices(tsupp, dims=2)
            tsupp = unique(tsupp, dims=2)
            blocks,cl,blocksize,sb,numb,_ = get_blocks(tsupp, basis, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md)
        end
    end
    if TS != false && QUIET == false
        mb = maximum(maximum.(sb))
        println("Obtained the block structure. The maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,momone,GramMat = blockupop(n, supp, coe, basis, blocks, cl, blocksize, nb=nb, solver=solver, feasible=feasible,
    QUIET=QUIET, solve=solve, dualize=dualize, solution=solution, MomentOne=MomentOne, Gram=Gram, cosmo_setting=cosmo_setting)
    data = upop_data(n, nb, x, f, supp, coe, basis, ksupp, blocks, sb, numb, GramMat, moment, solver, tol, 1)
    sol = nothing
    if solution == true
        sol,gap,data.flag = extract_solution(momone, opt, [f], x, tol=tol)
        if data.flag == 1
            sol = gap > 0.5 ? randn(n) : sol
            sol,data.flag = refine_sol(opt, sol, data, QUIET=true, tol=tol)
        end
    end
    return opt,sol,data
end

"""
    opt,sol,data = tssos_higher!(data; TS="block", merge=false, md=3, QUIET=false, solve=true,
    MomentOne=false, solution=false, tol=1e-4)

Compute higher steps of the TSSOS hierarchy.
Return the optimum, the (near) optimal solution (if `solution=true`) and other auxiliary data.
"""
function tssos_higher!(data::upop_data; TS="block", merge=false, md=3, QUIET=false, solve=true, feasible=false, MomentOne=false, Gram=false, 
    solution=false, cosmo_setting=cosmo_para(), dualize=false)
    n = data.n
    nb = data.nb
    x = data.x
    f = data.f
    supp = data.supp
    coe = data.coe
    basis = data.basis
    ksupp = data.ksupp
    sb = data.sb
    numb = data.numb
    solver = data.solver
    tol = data.tol
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    blocks,cl,blocksize,sb,numb,status = get_blocks(ksupp, basis, sb=sb, numb=numb, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md)
    opt = nothing
    sol = nothing
    if status == 1
        if QUIET == false
            mb = maximum(maximum.(sb))
            println("Obtained the block structure. The maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,momone,GramMat = blockupop(n, supp, coe, basis, blocks, cl, blocksize, nb=nb, solver=solver, feasible=feasible, 
        QUIET=QUIET, solve=solve, dualize=dualize, solution=solution, MomentOne=MomentOne, Gram=Gram, cosmo_setting=cosmo_setting)
        if solution == true
            sol,gap,data.flag = extract_solution(momone, opt, [f], x, tol=tol)
            if data.flag == 1
                sol = gap > 0.5 ? randn(n) : sol
                sol,data.flag = refine_sol(opt, sol, data, QUIET=true, tol=tol)
            end
        end
        data.ksupp = ksupp
        data.GramMat = GramMat
        data.sb = sb
        data.numb = numb
        data.blocks = blocks
        data.moment = moment
    end
    return opt,sol,data
end

function bin_add(bi, bj, nb)
    bs = bi + bj
    if nb > 0
        bs[1:nb,:] = isodd.(bs[1:nb,:])
    end
    return bs
end

function get_basis(n, d; nb=0, lead=[])
    lb = binomial(n+d, d)
    basis = zeros(UInt8, n, lb)
    i = 0
    t = 1
    while i < d+1
        t += 1
        if basis[n,t-1] == i
           if i < d
              basis[1,t] = i+1
           end
           i += 1
        else
            j = findfirst(x->basis[x,t-1]!=0, 1:n)
            basis[:,t] = basis[:,t-1]
            if j == 1
               basis[1,t] -= 1
               basis[2,t] += 1
            else
               basis[1,t] = basis[j,t] - 1
               basis[j,t] = 0
               basis[j+1,t] += 1
            end
        end
    end
    if nb > 0
        basis_bin = basis[1:nb,:]
        basis_valid = all.(x->x<=1, eachcol(basis_bin))
        basis = basis[:, basis_valid]
    end
    if !isempty(lead)
        basis_valid = map(a->!divide(a, lead, n, size(lead,2)), eachcol(basis))
        basis = basis[:, basis_valid]
    end
    return basis
end

function divide(a, lead, n, llead)
    return any(j->all(i->lead[i,j]<=a[i], 1:n), 1:llead)
end

function reminder(a, x, gb, n)
    remind = rem(prod(x.^a), gb)
    mon = monomials(remind)
    coe = coefficients(remind)
    lm = length(mon)
    supp = zeros(UInt8,n,lm)
    for i = 1:lm, j = 1:n
        @inbounds supp[j,i]=MultivariatePolynomials.degree(mon[i],x[j])
    end
    return lm,supp,coe
end

function newton_basis(n, d, supp; e=1e-5, solver="Mosek")
    lsupp = size(supp,2)
    basis = get_basis(n, d)
    lb = size(basis,2)
    A0 = [-1/2*supp' ones(lsupp,1)]
    t = 1
    indexb = [i for i=1:lb]
    temp = sortslices(supp, dims=2)
    while t <= lb
          i = indexb[t]
          if bfind(temp, lsupp, UInt8(2)*basis[:,i]) != 0
             t += 1
          else
             if solver == "Mosek"
                model = Model(optimizer_with_attributes(Mosek.Optimizer))
             elseif solver == "SDPT3"
                model = Model(optimizer_with_attributes(SDPT3.Optimizer))
             elseif solver == "SDPNAL"
                model = Model(optimizer_with_attributes(SDPNAL.Optimizer))
             else
                @error "The solver is currently not supported!"
                return nothing
             end
             set_optimizer_attribute(model, MOI.Silent(), true)
             @variable(model, x[1:n+1], lower_bound=-10, upper_bound=10)
             @constraint(model, [A0; [basis[:,i]' -1]]*x .<= zeros(lsupp+1))
             @objective(model, Min, [basis[:,i]' -1]*x)
             optimize!(model)
             vx = value.(x)
             if abs(objective_value(model)) <= e && sum(abs.(vx)) <= e
                t += 1
             else
                if abs(objective_value(model)) <= e && sum(abs.(vx)) > e
                   t += 1
                else
                   lb -= 1
                   indexb = deleteat!(indexb, t)
                end
                r = t
                while lb >= r
                      j = indexb[r]
                      if [basis[:,j]' -1]*vx <= -e
                         lb -= 1
                         indexb = deleteat!(indexb, r)
                      else
                         r += 1
                      end
                end
             end
          end
    end
    return basis[:,indexb]
end

function generate_basis!(supp, basis)
    supp = sortslices(supp, dims=2)
    supp = unique(supp, dims=2)
    lsupp = size(supp, 2)
    lb = size(basis, 2)
    indexb = UInt32[]
    for i = 1:lb, j = i:lb
        bi = basis[:,i] + basis[:,j]
        if bfind(supp, lsupp, bi) != 0
             push!(indexb, i, j)
        end
    end
    sort!(indexb)
    unique!(indexb)
    return basis[:,indexb]
end

function bfind(A, l, a)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if ndims(A) == 2
            temp = A[:, mid]
        else
            temp = A[mid]
        end
        if temp == a
           return mid
        elseif temp < a
           low = mid+1
        else
           high = mid-1
        end
    end
    return 0
end

function get_graph(tsupp::Array{UInt8, 2}, basis::Array{UInt8, 2}; nb=0, balanced=false)
    lb = size(basis,2)
    G = SimpleGraph(lb)
    ltsupp = size(tsupp,2)
    for i = 1:lb, j = i+1:lb
        bi = bin_add(basis[:,i], basis[:,j], nb)
        if bfind(tsupp, ltsupp, bi) != 0
           add_edge!(G, i, j)
        end
    end
    return G
end

function get_blocks(tsupp, basis; sb=[], numb=[], nb=0, TS="block", minimize=false,
    QUIET=true, merge=false, md=3)
    if TS == false
        blocksize = [size(basis,2)]
        blocks = [[i for i=1:size(basis,2)]]
        cl = 1
    else
        G = get_graph(tsupp, basis, nb=nb)
        if TS == "block"
            blocks = connected_components(G)
            blocksize = length.(blocks)
            cl = length(blocksize)
        else
            blocks,cl,blocksize = chordal_cliques!(G, method=TS, minimize=minimize)
            if merge == true
                blocks,cl,blocksize = clique_merge!(blocks, d=md, QUIET=true)
            end
        end
    end
    nsb = sort(unique(blocksize), rev=true)
    nnumb = [sum(blocksize.== i) for i in nsb]
    if isempty(sb) || nsb!=sb || nnumb!=numb
        status = 1
        if QUIET == false
            println("-----------------------------------------------------------------------------")
            println("The sizes of PSD blocks:\n$nsb\n$nnumb")
            println("-----------------------------------------------------------------------------")
        end
    else
        status = 0
        println("No higher TS step of the TSSOS hierarchy!")
    end
    return blocks,cl,blocksize,nsb,nnumb,status
end

function blockupop(n, supp, coe, basis, blocks, cl, blocksize; nb=0, solver="Mosek", feasible=false, QUIET=true, solve=true, 
    solution=false, MomentOne=false, Gram=false, cosmo_setting=cosmo_para(), dualize=false)
    tsupp = zeros(UInt8, n, Int(sum(Int.(blocksize).^2+blocksize)/2))
    k = 1
    for i = 1:cl, j = 1:blocksize[i], r = j:blocksize[i]
        @inbounds bi = bin_add(basis[:,blocks[i][j]], basis[:,blocks[i][r]], nb)
        @inbounds tsupp[:,k] = bi
        k += 1
    end
    if MomentOne == true || solution == true
        ksupp = copy(tsupp)
        ksupp = unique(ksupp, dims=2)
        ksupp = sortslices(ksupp, dims=2)
        tsupp = [tsupp get_basis(n,2,nb=nb)]
        tsupp = unique(tsupp, dims=2)
        tsupp = sortslices(tsupp, dims=2)
    else
        tsupp = unique(tsupp, dims=2)
        tsupp = sortslices(tsupp, dims=2)
        ksupp = tsupp
    end
    objv,moment,momone,GramMat = nothing,nothing,nothing,nothing
    if solve == true
        ltsupp = size(tsupp,2)
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
            return nothing,nothing,nothing,nothing
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons = [AffExpr(0) for i=1:ltsupp]
        if MomentOne == true || solution == true
            pos0 = @variable(model, [1:n+1, 1:n+1], PSD)
            for j = 1:n+1, k = j:n+1
                @inbounds bi = bin_add(basis[:,j],basis[:,k],nb)
                Locb = bfind(tsupp, ltsupp, bi)
                if j == k
                   @inbounds add_to_expression!(cons[Locb], pos0[j,k])
                else
                   @inbounds add_to_expression!(cons[Locb], 2, pos0[j,k])
                end
            end
        end
        pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl)
        for i = 1:cl
            bs = blocksize[i]
            if bs == 1
               @inbounds pos[i] = @variable(model, lower_bound=0)
               @inbounds bi = bin_add(basis[:,blocks[i][1]], basis[:,blocks[i][1]], nb)
               Locb = bfind(tsupp, ltsupp, bi)
               @inbounds add_to_expression!(cons[Locb], pos[i])
            else
               @inbounds pos[i] = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:blocksize[i], r = j:blocksize[i]
                   @inbounds bi = bin_add(basis[:,blocks[i][j]], basis[:,blocks[i][r]], nb)
                   Locb = bfind(tsupp, ltsupp, bi)
                   if j == r
                       @inbounds add_to_expression!(cons[Locb], pos[i][j,r])
                   else
                       @inbounds add_to_expression!(cons[Locb], 2, pos[i][j,r])
                   end
               end
            end
        end
        bc = zeros(ltsupp)
        for i in axes(supp,2)
            Locb = bfind(tsupp, ltsupp, supp[:,i])
            if Locb == 0
               @error "The monomial basis is not enough!"
               return nothing,nothing,nothing,nothing
            else
               bc[Locb] = coe[i]
            end
        end
        if feasible == false
            @variable(model, lower)
            @objective(model, Max, lower)
            cons[1] += lower
        else
            tnorm = AffExpr(0)
            for i = 1:cl
                if blocksize[i] == 1
                    tnorm += pos[i]
                else
                    tnorm += tr(pos[i])
                end
            end
            @objective(model, Min, tnorm)
        end
        @constraint(model, con[i=1:ltsupp], cons[i]==bc[i])
        if QUIET == false
            println("Solving the SDP...")
        end
        time = @elapsed begin
        optimize!(model)
        end
        if QUIET == false
            println("SDP solving time: $time seconds.")
        end
        status = termination_status(model)
        objv = objective_value(model)
        if status != MOI.OPTIMAL
           println("termination status: $status")
           status = primal_status(model)
           println("solution status: $status")
        end
        println("optimum = $objv")
        if Gram == true
            GramMat = [value.(pos[i]) for i = 1:cl]
        end
        dual_var = -dual.(con)
        moment = Vector{Matrix{Float64}}(undef, cl)
        for i = 1:cl
            moment[i] = zeros(blocksize[i],blocksize[i])
            for j = 1:blocksize[i], k = j:blocksize[i]
                bi = bin_add(basis[:,blocks[i][j]], basis[:,blocks[i][k]], nb)
                Locb = bfind(tsupp, ltsupp, bi)
                moment[i][j,k] = dual_var[Locb]
            end
            moment[i] = Symmetric(moment[i],:U)
        end
        if solution == true
            momone = zeros(Float64, n+1, n+1)
            for j = 1:n+1, k = j:n+1
                bi = bin_add(basis[:,j], basis[:,k], nb)
                Locb = bfind(tsupp, ltsupp, bi)
                momone[j,k] = dual_var[Locb]
            end
            momone = Symmetric(momone,:U)
        end
    end
    return objv,ksupp,moment,momone,GramMat
end

function extract_solution(moment, opt, pop, x; numeq=0, tol=1e-4)
    n = length(x)
    m = length(pop) - 1
    F = eigen(moment, n+1:n+1)
    sol = sqrt(F.values[1])*F.vectors[:,1]
    if abs(sol[1]) < 1e-8
        return nothing,1,1
    else
        sol = sol[2:end]/sol[1]
        ub = polynomial(pop[1])(x => sol)
        gap = abs(opt-ub)/max(1, abs(ub))
        flag = gap >= tol ? 1 : 0
        for i = 1:m-numeq
            if polynomial(pop[i+1])(x => sol) <= -tol
                flag = 1
            end
        end
        for i = m-numeq+1:m
            if abs(polynomial(pop[i+1])(x => sol)) >= tol
                flag = 1
            end
        end
        if flag == 0
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
        end
        return sol,gap,flag
    end
end
