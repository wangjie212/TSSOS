"""
    optimum = SumOfRatios(p, q, g, h, x, d; QUIET=false, SignSymmetry=true, GroebnerBasis=false)

Minimizing the sum of ratios p1/q1 + ... + pN/qN on the set defined by g >= 0 and h == 0.

# Input arguments
- `p`: vector of denominators
- `q`: vector of numerator
- `g`: inequality constraints
- `h`: equality constraints
- `x`: vector of variables
- `d`: relaxation order

# Output arguments
- `SignSymmetry`: exploit sign symmetries or not (`true`, `false`)
- `GroebnerBasis`: exploit the quotient ring structure or not (`true`, `false`)
"""
function SumOfRatios(p, q, g, h, x, d; QUIET=false, dualize=false, SignSymmetry=true, mosek_setting=mosek_para(), GroebnerBasis=false)
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    N = length(p)
    if SignSymmetry == true
        polys = [p; q]
        if !isempty(g)
            polys = [polys; g]
        end
        if !isempty(h)
            polys = [polys; h]
        end
        ss = get_signsymmetry(polys, x)
        TS = "block"
    else
        ss = TS = false
    end
    dq = maxdegree.(q)
    time = @elapsed begin
    if dualize == false
        model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
    else
        model = Model(dual_optimizer(Mosek.Optimizer))
    end
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    hh = Vector{Polynomial{true, AffExpr}}(undef, N-1)
    for i = 1:N-1
        hh[i] = add_poly!(model, x, 2d-max(dq[i], dq[N]), signsymmetry=ss)[1]
        add_psatz!(model, p[i]-hh[i]*q[i], x, g, h, d, QUIET=QUIET, CS=false, TS=TS, GroebnerBasis=GroebnerBasis)
    end
    c = @variable(model)
    add_psatz!(model, p[N]+(sum(hh)-c)*q[N], x, g, h, d, QUIET=QUIET, CS=false, TS=TS, GroebnerBasis=GroebnerBasis)
    @objective(model, Max, c)
    optimize!(model)
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("termination status: $status")
        status = primal_status(model)
        println("solution status: $status")
    end
    optimum = objective_value(model)
    end
    @show optimum
    println("solving time: $time seconds.")
    return optimum
end

"""
    optimum = SparseSumOfRatios(p, q, g, h, x, d; QUIET=false, SignSymmetry=true, GroebnerBasis=false)

Minimizing the sum of sparse ratios p1/q1 + ... + pN/qN on the set defined by g >= 0 and h == 0.

# Input arguments
- `p`: vector of denominators
- `q`: vector of numerator
- `g`: inequality constraints
- `h`: equality constraints
- `x`: vector of variables
- `d`: relaxation order

# Output arguments
- `SignSymmetry`: exploit sign symmetries or not (`true`, `false`)
- `GroebnerBasis`: exploit the quotient ring structure or not (`true`, `false`)
"""
function SparseSumOfRatios(p, q, g, h, x, d; QUIET=false, dualize=false, SignSymmetry=true, mosek_setting=mosek_para(), GroebnerBasis=false)
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    N = length(p)
    if SignSymmetry == true
        polys = [p; q]
        if !isempty(g)
            polys = [polys; g]
        end
        if !isempty(h)
            polys = [polys; h]
        end
        ss = get_signsymmetry(polys, x)
    end
    TS = SignSymmetry == true ? "block" : false
    dq = maxdegree.(q)
    vp = Vector{Vector{PolyVar{true}}}(undef, N)
    for i = 1:N
        if typeof(p[i]) == Polynomial
            vp[i] = variables(p[i])
        else
            vp[i] = PolyVar{true}[]
        end
    end
    I = [unique([vp[i]; variables(q[i])]) for i=1:N]
    vg = variables.(g)
    J = [findall(j->issubset(vg[j], I[i]), 1:length(g)) for i=1:N]
    vh = variables.(h)
    K = [findall(j->issubset(vh[j], I[i]), 1:length(h)) for i=1:N]
    U = [findall(j->!isempty(intersect(I[i], I[j])), i+1:N) .+ i for i=1:N-1]
    V = [findall(j->!isempty(intersect(I[i], I[j])), 1:i-1) for i=2:N]
    time = @elapsed begin
    if dualize == false
        model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
    else
        model = Model(dual_optimizer(Mosek.Optimizer))
    end
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    hh = Vector{Vector{Polynomial{true, AffExpr}}}(undef, N)
    for i = 1:N-1
        hh[i] = Vector{Polynomial{true, AffExpr}}(undef, length(U[i]))
        for (ind, j) in enumerate(U[i])
            Iij = intersect(I[i], I[j])
            if SignSymmetry == true
                temp = [ncbfind(x, length(x), item) for item in Iij]
                hh[i][ind] = add_poly!(model, Iij, 2d-max(dq[i], dq[j]), signsymmetry=ss[temp, :])[1]
            else
                hh[i][ind] = add_poly!(model, Iij, 2d-max(dq[i], dq[j]))[1]
            end
        end
    end
    hh[N] = Polynomial{true, AffExpr}[]
    c = @variable(model, [1:N])
    for i = 1:N
        if i == 1
            temp = 0
        else
            ind = [bfind(U[j], length(U[j]), i) for j in V[i-1]]
            temp = sum(hh[V[i-1][j]][ind[j]] for j=1:length(V[i-1]))
        end
        add_psatz!(model, p[i]-(c[i]+sum(hh[i])-temp)*q[i], I[i], g[J[i]], h[K[i]], d, QUIET=QUIET, CS=false, TS=TS, GroebnerBasis=GroebnerBasis)
    end 
    @objective(model, Max, sum(c))
    end
    println("SDP assembling time: $time seconds.")
    println("Solving the SDP...")
    time = @elapsed begin
    optimize!(model)
    end
    println("SDP solving time: $time seconds.")
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("termination status: $status")
        status = primal_status(model)
        println("solution status: $status")
    end
    optimum = objective_value(model)
    @show optimum
    return optimum
end