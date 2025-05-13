# Homogenize the polynomial f with the homogenization variable z
function homogenize(f, z)
    d = maxdegree(f)
    ts = [term*z^(d-maxdegree(term)) for term in MP.terms(f)]
    return sum(ts)
end

function solve_hpop(cost, vars, ineq_cons, eq_cons, order; QUIET=false, CS="MF", type=2, ε=0, TS="block", SO=1, nnhomovar=false, GroebnerBasis=false, Mommat=false)
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    if CS != false
        fsupp = poly_info(cost, vars)[1]
        if ineq_cons != []
            gsupp = npolys_info(ineq_cons, vars)[1]
        else
            gsupp = Matrix{UInt8}[]
        end
        if eq_cons != []
            hsupp = npolys_info(eq_cons, vars)[1]
        else
            hsupp = Matrix{UInt8}[]
        end
        cliques,cql,_ = clique_decomp(length(vars), length(ineq_cons), length(eq_cons), fsupp, gsupp, hsupp, alg="NC", QUIET=QUIET)  
    end
    @polyvar z
    if ineq_cons != []
        ineq_cons = homogenize.(ineq_cons, z)
        dg = maxdegree.(ineq_cons)
    else
        dg = [0]
    end
    if eq_cons != []
        eq_cons = homogenize.(eq_cons, z)
    end
    cost = homogenize(cost, z)
    d = maxdegree(cost)
    if nnhomovar == true || isodd(d) || any(isodd.(dg))
        push!(ineq_cons, z)
    end
    if CS != false && cql > 1
        if type == 1
            @polyvar y[1:cql]
            for i = 1:cql
                push!(eq_cons, sum(vars[cliques[i]].^2) + z^2 + y[i]^2 - 1)
            end
        else
            freq = zeros(length(vars))
            for clique in cliques
                freq[clique] .+= 1
            end
            if type == 2
                @polyvar y[1:cql-1]
                push!(ineq_cons, 2 - sum(vars[cliques[1]].^2) - z^2 - y[1]^2)
                push!(eq_cons, sum(1 ./ freq[cliques[1]] .* vars[cliques[1]].^2) + 1/cql*z^2 - y[1]^2)
                for i = 2:cql-1
                    push!(ineq_cons, 3 - sum(vars[cliques[i]].^2) - z^2 - y[i-1]^2 - y[i]^2)
                    push!(eq_cons, sum(1 ./ freq[cliques[i]] .* vars[cliques[i]].^2) + 1/cql*z^2 + y[i-1]^2 - y[i]^2)
                end
                push!(ineq_cons, 2 - sum(vars[cliques[end]].^2) - z^2 - y[end]^2)
                push!(eq_cons, sum(1 ./ freq[cliques[end]] .* vars[cliques[end]].^2) + 1/cql*z^2 + y[end]^2 - 1)
            else
                @polyvar y[1:cql]
                for i = 1:cql
                    push!(eq_cons, sum(1 ./ freq[cliques[i]] .* vars[cliques[i]].^2) + 1/cql*z^2 + y[i]^2 - 1)
                end
                push!(eq_cons, sum(y.^2) - cql + 1)
            end
            # append!(ineq_cons, 1 .- vars.^2)
            # push!(ineq_cons, 1 - z^2)
            # append!(ineq_cons, 1 .- y.^2)
        end
        nvars = [vars; z; y]
    else
        push!(eq_cons, sum(vars.^2) + z^2 - 1)
        nvars = [vars; z]
    end
    if ε > 0
        d0 = iseven(d) ? d : d+1
        cost += ε*sum(nvars.^d0)
    end
    time = @elapsed begin
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    @variable(model, λ)
    info = add_psatz!(model, cost-λ*z^d, nvars, ineq_cons, eq_cons, order, QUIET=QUIET, CS=CS, TS=TS, SO=SO, GroebnerBasis=GroebnerBasis, constrs="con")
    @objective(model, Max, λ)
    optimize!(model)
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("termination status: $status")
        status = primal_status(model)
        println("solution status: $status")
    end
    optimum = objective_value(model)
    @show optimum
    end
    println("POP solving time: $time seconds.")
    moment = -dual(constraint_by_name(model, "con"))
    MomMat = get_moment_matrix(moment, info)
    return optimum,MomMat,info
end