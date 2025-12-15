"""
    obj,sol,status = local_solution(npop, n; nb=0, numeq=0,
    startpoint=[], QUIET=false)

Compute a local solution by a local solver.
# Input arguments
- `npop`: polynomial optimiztion problem
- `n`: number of POP variables
- `nb`: number of binary variables
- `numeq`: number of equality constraints
- `startpoint`: provide a start point
- `QUIET`: run in the quiet mode or not (`true`, `false`)

# Output arguments
- `obj`: local optimum
- `sol`: local solution
- `status`: solver termination status
"""
function local_solution(npop, n; nb=0, numeq=0, startpoint=[], QUIET=false)
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
    @NLobjective(model, Min, sum(npop[1].coe[j]*prod(x[k] for k in item) for (j, item) in enumerate(npop[1].supp)))
    for i = 1:nb
        @NLconstraint(model, x[i]^2-1==0)
    end
    for (i, p) in enumerate(npop[2:end])
        if i < length(npop) - numeq
            @NLconstraint(model, sum(p.coe[j]*prod(x[k] for k in item) for (j, item) in enumerate(p.supp)) >= 0)
        else
            @NLconstraint(model, sum(p.coe[j]*prod(x[k] for k in item) for (j, item) in enumerate(p.supp)) == 0)
        end
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
    ref_sol,status,flag = refine_sol(opt, sol, data, QUIET=false, tol=1e-2)

Refine the obtained solution by a local solver.
Return the refined solution, and `flag=0` if global optimality is certified, `flag=1` otherwise.
"""
function refine_sol(opt, sol, data::Union{cpop_data,mcpop_data}; QUIET=false, gtol=1e-2)
    numeq = data.numeq
    if typeof(data) == cpop_data && !isempty(data.gb)
        npop = [data.obj; data.ineq_cons[2:end]; [poly(p, data.x) for p in data.pop[end-numeq+1:end]]]
    else
        npop = [data.obj; data.ineq_cons[2:end]; data.eq_cons]
    end
    sol[abs.(sol) .< 1e-10] .= 1e-10
    ub,rsol,status = local_solution(npop, data.n, nb=data.nb, numeq=numeq, startpoint=sol, QUIET=QUIET)
    if status == MOI.LOCALLY_SOLVED
        gap = abs(opt-ub)/max(1, abs(ub))
        if gap < gtol
            println("------------------------------------------------")
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
            println("Successfully extracted one globally optimal solution.")
            println("------------------------------------------------")
            return rsol,status,0
        else
            @printf "Found a locally optimal solution by Ipopt, giving an upper bound: %.8f.\nThe relative optimality gap is: %.6f%%.\n" ub 100*gap
            return rsol,status,1
        end
    else
        println("The local solver failed refining the solution!")
        return sol,status,1
    end
end
