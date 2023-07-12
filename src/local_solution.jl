function local_solution(n, m, supp::Vector{Vector{Vector{UInt16}}}, coe; nb=0, numeq=0,
    startpoint=[], QUIET=false)
    model = Model(optimizer_with_attributes(Ipopt.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    if QUIET == true
        set_optimizer_attribute(model, "print_level", 0)
    end
    if isempty(startpoint)
        @variable(model, x[1:n], start = 0)
    else
        @variable(model, x[i=1:n], start = startpoint[i])
    end
    @NLobjective(model, Min, sum(coe[1][j]*prod(x[supp[1][j][k]] for k=1:length(supp[1][j])) for j=1:length(supp[1])))
    for i = 1:nb
        @NLconstraint(model, x[i]^2-1==0)
    end
    for i in 1:m-numeq
        @NLconstraint(model, sum(coe[i+1][j]*prod(x[supp[i+1][j][k]] for k=1:length(supp[i+1][j])) for j=1:length(supp[i+1]))>=0)
    end
    for i in m-numeq+1:m
        @NLconstraint(model, sum(coe[i+1][j]*prod(x[supp[i+1][j][k]] for k=1:length(supp[i+1][j])) for j=1:length(supp[i+1]))==0)
    end
    optimize!(model)
    status = termination_status(model)
    objv = objective_value(model)
    if QUIET == false
        println("optimum = $objv")
    end
    return objv,value.(x),status
end

function local_solution(n, m, supp::Vector{Array{UInt8, 2}}, coe; nb=0, numeq=0,
    startpoint=[], QUIET=false)
    model = Model(optimizer_with_attributes(Ipopt.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    if QUIET == true
        set_optimizer_attribute(model, "print_level", 0)
    end
    if isempty(startpoint)
        @variable(model, x[1:n], start = 0)
    else
        @variable(model, x[i=1:n], start = startpoint[i])
    end
    @NLobjective(model, Min, sum(coe[1][j]*prod(x[k]^supp[1][k,j] for k=1:n) for j=1:size(supp[1], 2)))
    for i = 1:nb
        @NLconstraint(model, x[i]^2-1==0)
    end
    for i in 1:m-numeq
        @NLconstraint(model, sum(coe[i+1][j]*prod(x[k]^supp[i+1][k,j] for k=1:n) for j=1:size(supp[i+1], 2))>=0)
    end
    for i in m-numeq+1:m
        @NLconstraint(model, sum(coe[i+1][j]*prod(x[k]^supp[i+1][k,j] for k=1:n) for j=1:size(supp[i+1], 2))==0)
    end
    optimize!(model)
    status = termination_status(model)
    objv = objective_value(model)
    if QUIET == false
        println("optimum = $objv")
    end
    return objv,value.(x),status
end

"""
    ref_sol,ub,rel_gap = refine_sol(opt, sol, data, QUIET=false, tol=1e-4)

Refine the obtained solution by a local solver.
Return the refined solution, the upper bound given by the local solver and the
relative optimality gap.
"""
function refine_sol(opt, sol, data::upop_data; QUIET=false, tol=1e-4)
    n = data.n
    nb = data.nb
    supp = data.supp
    coe = data.coe
    for i = 1:n
        if abs(sol[i]) < 1e-10
            sol[i] = 1e-10
        end
    end
    ub,rsol,status = local_solution(n, 0, [supp], [coe], nb=nb, numeq=0, startpoint=sol, QUIET=QUIET)
    if status == MOI.LOCALLY_SOLVED
        gap = abs(opt-ub)/max(1, abs(ub))
        if gap < tol
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
            return rsol,0
        else
            @printf "Found a locally optimal solution by Ipopt, giving an upper bound: %.8f.\nThe relative optimality gap is: %.6f%%.\n" ub 100*gap
            return rsol,1
        end
    else
        println("The local solver failed refining the solution!")
        return sol,1
    end
end

function refine_sol(opt, sol, data::Union{cpop_data,mcpop_data}; QUIET=false, tol=1e-4)
    n = data.n
    nb = data.nb
    numeq = data.numeq
    if typeof(data) == cpop_data && !isempty(data.gb)
        m = length(data.pop)-1
        supp = Vector{Array{UInt8,2}}(undef, m+1)
        coe = Vector{Vector{Float64}}(undef, m+1)
        supp[2:m+1-numeq] = data.supp[2:end]
        coe[2:m+1-numeq] = data.coe[2:end]
        for k in [1; [k for k=m+2-numeq:m+1]]
            mons = monomials(data.pop[k])
            coe[k] = coefficients(data.pop[k])
            supp[k] = zeros(UInt8, n, length(mons))
            for i in eachindex(mons), j = 1:n
                @inbounds supp[k][j,i] = MultivariatePolynomials.degree(mons[i], data.x[j])
            end
        end
    else
        m = data.m
        supp = data.supp
        coe = data.coe
    end
    for i = 1:n
        if abs(sol[i]) < 1e-10
            sol[i] = 1e-10
        end
    end
    ub,rsol,status = local_solution(n, m, supp, coe, nb=nb, numeq=numeq, startpoint=sol, QUIET=QUIET)
    if status == MOI.LOCALLY_SOLVED
        gap = abs(opt-ub)/max(1, abs(ub))
        if gap < tol
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
            return rsol,0
        else
            @printf "Found a locally optimal solution by Ipopt, giving an upper bound: %.8f.\nThe relative optimality gap is: %.6f%%.\n" ub 100*gap
            return rsol,1
        end
    else
        println("The local solver failed refining the solution!")
        return sol,1
    end
end
