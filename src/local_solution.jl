function local_solution(n, m, supp::Vector{Vector{Vector{UInt16}}}, coe; numeq=0,
    startpoint=[], QUIET=false)
    model=Model(optimizer_with_attributes(Ipopt.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    if QUIET==true
        set_optimizer_attribute(model, "print_level", 0)
    end
    if isempty(startpoint)
        @variable(model, x[1:n], start = 0)
    else
        @variable(model, x[i=1:n], start = startpoint[i])
    end
    @NLobjective(model, Min, sum(coe[1][j]*prod(x[supp[1][j][k]] for k=1:length(supp[1][j])) for j=1:length(supp[1])))
    for i in 1:m-numeq
        @NLconstraint(model, sum(coe[i+1][j]*prod(x[supp[i+1][j][k]] for k=1:length(supp[i+1][j])) for j=1:length(supp[i+1]))>=0)
    end
    for i in m-numeq+1:m
        @NLconstraint(model, sum(coe[i+1][j]*prod(x[supp[i+1][j][k]] for k=1:length(supp[i+1][j])) for j=1:length(supp[i+1]))==0)
    end
    optimize!(model)
    status=termination_status(model)
    objv = objective_value(model)
    if QUIET==false
        println("optimum = $objv")
    end
    return objv,value.(x),status
end

function local_solution(n, m, supp::Vector{Array{UInt8, 2}}, coe; numeq=0,
    startpoint=[], QUIET=false)
    model=Model(optimizer_with_attributes(Ipopt.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    if QUIET==true
        set_optimizer_attribute(model, "print_level", 0)
    end
    if isempty(startpoint)
        @variable(model, x[1:n], start = 0)
    else
        @variable(model, x[i=1:n], start = startpoint[i])
    end
    @NLobjective(model, Min, sum(coe[1][j]*prod(x[k]^supp[1][k,j] for k=1:n) for j=1:size(supp[1], 2)))
    for i in 1:m-numeq
        @NLconstraint(model, sum(coe[i+1][j]*prod(x[k]^supp[i+1][k,j] for k=1:n) for j=1:size(supp[i+1], 2))>=0)
    end
    for i in m-numeq+1:m
        @NLconstraint(model, sum(coe[i+1][j]*prod(x[k]^supp[i+1][k,j] for k=1:n) for j=1:size(supp[i+1], 2))==0)
    end
    optimize!(model)
    status=termination_status(model)
    objv = objective_value(model)
    if QUIET==false
        println("optimum = $objv")
    end
    return objv,value.(x),status
end

"""
    ref_sol,upper_bound,rel_gap = refine_sol(opt, sol, data, QUIET=false, tol=1e-4)

Refine the obtained solution by a local solver.
Return the refined solution, the upper bound given by the local solver and the
relative optimality gap.
"""
function refine_sol(opt, sol, data::upop_data; QUIET=false, tol=1e-4)
    n=data.n
    supp=data.supp
    coe=data.coe
    upper_bound,rsol,status=local_solution(n, 0, [supp], [coe], numeq=0, startpoint=sol, QUIET=QUIET)
    if status==MOI.LOCALLY_SOLVED
        gap=abs(upper_bound)>1 ? abs((opt-upper_bound)/upper_bound) : abs(opt-upper_bound)
        if gap < tol
            println("Global optimality certified!")
        end
    else
        rsol,upper_bound,gap=sol,nothing,nothing
        println("The local solver failed!")
    end
    return rsol,upper_bound,gap
end

function refine_sol(opt, sol, data::Union{cpop_data,mcpop_data}; QUIET=false, tol=1e-4)
    n=data.n
    m=data.m
    numeq=data.numeq
    supp=data.supp
    coe=data.coe
    upper_bound,rsol,status=local_solution(n, m, supp, coe, numeq=numeq, startpoint=sol, QUIET=QUIET)
    if status==MOI.LOCALLY_SOLVED
        gap=abs(upper_bound)>1 ? abs((opt-upper_bound)/upper_bound) : abs(opt-upper_bound)
        if gap < tol
            println("Global optimality certified!")
        end
    else
        rsol,upper_bound,gap=sol,nothing,nothing
        println("The local solver failed!")
    end
    return rsol,upper_bound,gap
end
