"""
    obj,sol,status = local_solution(n, m, supp::Vector{Union{Vector{Vector{UInt16}}, Array{UInt8, 2}}}, coe; nb=0, numeq=0,
    startpoint=[], QUIET=false)

Compute a local solution by a local solver.
# Input arguments
- `n`: number of POP variables
- `m`: number of POP constraints
- `supp`: supports of the POP
- `coe`: coefficients of the POP
- `nb`: number of binary variables
- `numeq`: number of equality constraints
- `startpoint`: provide a start point
- `QUIET`: run in the quiet mode or not (`true`, `false`)

# Output arguments
- `obj`: local optimum
- `sol`: local solution
- `status`: solver termination status
"""
function local_solution(n, m, supp::Vector{Vector{Vector{UInt16}}}, coe; nb=0, numeq=0, startpoint=[], QUIET=false)
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

function local_solution(n, m, supp::Vector{Array{UInt8, 2}}, coe; nb=0, numeq=0, startpoint=[], QUIET=false)
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
    ref_sol,flag = refine_sol(opt, sol, data, QUIET=false, tol=1e-2)

Refine the obtained solution by a local solver.
Return the refined solution, and `flag=0` if global optimality is certified, `flag=1` otherwise.
"""
function refine_sol(opt, sol, data::Union{cpop_data,mcpop_data}; QUIET=false, gtol=1e-2)
    n = data.n
    nb = data.nb
    numeq = data.numeq
    if typeof(data) == cpop_data && !isempty(data.gb)
        m = length(data.pop) - 1
        supp = Vector{Array{UInt8,2}}(undef, m+1)
        coe = Vector{Vector{Float64}}(undef, m+1)
        supp[2:m+1-numeq] = data.supp[2:end]
        coe[2:m+1-numeq] = data.coe[2:end]
        for k in [1; [k for k=m+2-numeq:m+1]]
            mons = MP.monomials(data.pop[k])
            coe[k] = MP.coefficients(data.pop[k])
            supp[k] = zeros(UInt8, n, length(mons))
            for i in eachindex(mons), j = 1:n
                @inbounds supp[k][j,i] = MP.degree(mons[i], data.x[j])
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
        if gap < gtol
            println("------------------------------------------------")
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
            println("Successfully extracted one globally optimal solution.")
            println("------------------------------------------------")
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
