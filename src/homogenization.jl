# Homogenize the polynomial f with the homogenization variable z
function homogenize(f, z)
    d = MP.maxdegree(f)
    ts = [term*z^(d-MP.maxdegree(term)) for term in MP.terms(f)]
    return sum(ts)
end

function solve_hpop(cost, x, ineq_cons, eq_cons, order; QUIET=false, CS="MF", type=2, ε=0, TS="block", eqTS=TS, SO=1, nnhomovar=false, GroebnerBasis=false)
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    if CS != false
        f = poly(cost, x)
        g = isempty(ineq_cons) ? poly[] : [poly(p, x) for p in ineq_cons]
        h = isempty(eq_cons) ? poly[] : [poly(p, x) for p in eq_cons]
        cliques,cql,_ = clique_decomp([f; g; h], length(x), length(h), alg="NC", QUIET=QUIET)  
    end
    @polyvar z
    if !isempty(ineq_cons)
        ineq_cons = homogenize.(ineq_cons, z)
        dg = MP.maxdegree.(ineq_cons)
    else
        dg = [0]
    end
    if !isempty(eq_cons)
        eq_cons = homogenize.(eq_cons, z)
    end
    cost = homogenize(cost, z)
    d = MP.maxdegree(cost)
    if nnhomovar == true || isodd(d) || any(isodd.(dg))
        push!(ineq_cons, z)
    end
    if CS != false && cql > 1
        if type == 1
            @polyvar y[1:cql]
            for i = 1:cql
                push!(eq_cons, sum(x[cliques[i]].^2) + z^2 + y[i]^2 - 1)
            end
        else
            freq = zeros(length(x))
            for clique in cliques
                freq[clique] .+= 1
            end
            if type == 2
                @polyvar y[1:cql-1]
                push!(ineq_cons, 2 - sum(x[cliques[1]].^2) - z^2 - y[1]^2)
                push!(eq_cons, sum(1 ./ freq[cliques[1]] .* x[cliques[1]].^2) + 1/cql*z^2 - y[1]^2)
                for i = 2:cql-1
                    push!(ineq_cons, 3 - sum(x[cliques[i]].^2) - z^2 - y[i-1]^2 - y[i]^2)
                    push!(eq_cons, sum(1 ./ freq[cliques[i]] .* x[cliques[i]].^2) + 1/cql*z^2 + y[i-1]^2 - y[i]^2)
                end
                push!(ineq_cons, 2 - sum(x[cliques[end]].^2) - z^2 - y[end]^2)
                push!(eq_cons, sum(1 ./ freq[cliques[end]] .* x[cliques[end]].^2) + 1/cql*z^2 + y[end]^2 - 1)
            else
                @polyvar y[1:cql]
                for i = 1:cql
                    push!(eq_cons, sum(1 ./ freq[cliques[i]] .* x[cliques[i]].^2) + 1/cql*z^2 + y[i]^2 - 1)
                end
                push!(eq_cons, sum(y.^2) - cql + 1)
            end
            # append!(ineq_cons, 1 .- x.^2)
            # push!(ineq_cons, 1 - z^2)
            # append!(ineq_cons, 1 .- y.^2)
        end
        nvars = [x; z; y]
    else
        push!(eq_cons, sum(x.^2) + z^2 - 1)
        nvars = [x; z]
    end
    if ε > 0
        d0 = iseven(d) ? d : d+1
        cost += ε*sum(nvars.^d0)
    end
    time = @elapsed begin
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    @variable(model, λ)
    info = add_psatz!(model, cost-λ*z^d, nvars, ineq_cons, eq_cons, order, QUIET=QUIET, CS=CS, TS=TS, eqTS=eqTS, SO=SO, GroebnerBasis=GroebnerBasis, constrs="con")
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
