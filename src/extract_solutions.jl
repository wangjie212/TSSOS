using RowEchelon

# extract a solution from the eigenvector associated with the maximal eigenvalue of the moment matrix
function extract_solution(moment, lb, pop, x; numeq=0, gtol=1e-2, ftol=1e-3, QUIET=true)
    F = eigen(moment, length(x)+1:length(x)+1)
    sol = sqrt(F.values[1])*F.vectors[:,1]
    if abs(sol[1]) < 1e-6
        return nothing,1,1
    else
        sol = sol[2:end]/sol[1]
        ub = MP.polynomial(pop[1])(x => sol)
        gap = abs(lb-ub)/max(1, abs(ub))
        flag = 1
        if check_optimality(sol, lb, pop, x, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
            flag = 0
            println("------------------------------------------------")
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
            println("Successfully extracted one globally optimal solution.")
            println("------------------------------------------------")
        end
        return sol,gap,flag
    end
end

"""
    sol = extract_solutions(moment, d; pop=nothing, x=nothing, lb=nothing, numeq=0, nb=0, rtol=1e-2, gtol=1e-2, ftol=1e-3)

Extract a tuple of solutions from the moment matrix using Henrion-Lasserre's algorithm.

# Input arguments
- `moment`: moment matrix
- `d`: relaxation order
- `pop`: polynomial optimization problem
- `x`: set of variables
- `lb`: SDP lower bound
- `numeq`: number of equality constraints
- `nb`: number of binary variables
- `rtol`: tolerance for rank
- `gtol`: tolerance for global optimality gap
- `ftol`: tolerance for feasibility

# Output arguments
- `sol`: a tuple of solutions
"""
function extract_solutions(moment, d; pop=nothing, x=nothing, lb=nothing, numeq=0, nb=0, basis=[], check=false, rtol=1e-2, gtol=1e-2, ftol=1e-3, QUIET=true)
    n = length(x)
    if rank(moment, rtol) == 1
        sol = [moment[2:n+1, 1]]
        if check == true && check_optimality(sol[1], lb, pop, x, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
            ub = MP.polynomial(pop[1])(x => sol[1])
            gap = abs(lb-ub)/max(1, abs(ub))
            println("------------------------------------------------")
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
            println("Successfully extracted one globally optimal solution.")
            println("------------------------------------------------")
        end
    else
        U,pivots = rref_with_pivots!(Matrix(moment), rtol)
        U = U[1:length(pivots), :]'
        if isempty(basis)
            basis = get_basis(n, d)
        end
        w = basis[:, pivots]
        lw = size(w, 2)
        if any([sum(w[:,j]) == d for j = 1:lw])
            return nothing
        end
        println("Rank of the moment matrix = ", lw)
        N = Vector{Matrix{Float64}}(undef, n)
        for i = 1:n
            kk = UInt16[]
            temp = zeros(UInt8, n)
            temp[i] = 1
            for j = 1:lw
                if w[i,j] == 1 && i <= nb
                    xwj = w[:,j] - temp
                else
                    xwj = w[:,j] + temp
                end
                locb = bfind_to(basis, size(basis, 2), xwj, n)
                kk = push!(kk, locb)
            end
            N[i] = U[kk,:]
        end
        rands = rand(n)
        rands = rands/sum(rands)
        L = schur(sum(rands[i]*N[i] for i in 1:n)).Z
        sol = [[L[:,i]'*N[j]*L[:,i] for j = 1:n] for i = 1:lw]
        if check == true
            sol = check_solution(sol, lb, pop, x, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
        end
    end
    return sol
end

"""
    sol = extract_solutions_robust(moment, n, d; pop=nothing, x=nothing, lb=nothing, numeq=0, nb=0, rtol=1e-2, gtol=1e-2, ftol=1e-3)

Extract a tuple of solutions from the moment matrix using the GNS robust algorithm.

# Input arguments
- `moment`: moment matrix
- `n`: number of variables
- `d`: relaxation order
- `pop`: polynomial optimization problem
- `x`: set of variables
- `lb`: SDP lower bound
- `numeq`: number of equality constraints
- `nb`: number of binary variables
- `rtol`: tolerance for rank
- `gtol`: tolerance for global optimality gap
- `ftol`: tolerance for feasibility

# Output arguments
- `sol`: a tuple of solutions
- `w`: weight vector
"""
function extract_solutions_robust(moment, n, d; type=Float64, pop=nothing, x=nothing, lb=nothing, numeq=0, basis=[], check=false, rtol=1e-2, gtol=1e-2, ftol=1e-3, QUIET=true)
    if isempty(basis)
        basis = get_basis(n, d)
    end
    ls = binomial(n+d-1, n)
    N = Vector{Matrix{type}}(undef, n)
    for i = 1:n
        N[i] = zeros(type, ls, ls)
        temp = zeros(UInt8, n)
        temp[i] = 1
        for j = 1:ls, k = 1:ls
            loc = bfind_to(basis, size(basis, 2), basis[:,k] + temp, n)
            N[i][k,j] = moment[loc, j]
        end
    end
    F = svd(moment[1:ls, 1:ls])
    S = sqrt.(F.S[F.S/F.S[1] .> rtol])
    for i = 1:n
        N[i] = Diagonal(S.^(-1))*F.Vt[1:length(S),:]*N[i]*F.U[:,1:length(S)]*Diagonal(S.^(-1))
    end
    rands = rand(n)
    rands = rands/sum(rands)
    L = schur(sum(rands[i]*N[i] for i in 1:n)).Z
    sol = [[L[:,i]'*N[j]*L[:,i] for j = 1:n] for i = 1:length(S)]
    W = L'*Diagonal(S)*F.Vt[1:length(S),:]
    w = abs.(W[:,1]).^2
    if check == true
        sol = check_solution(sol, lb, pop, x, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
    end
    return sol,w
end

function extract_solutions_robust(moment, n, d, cliques, cql, cliquesize; pop=nothing, x=nothing, supp=[], coe=[], lb=nothing, numeq=0, check=false, rtol=1e-2, gtol=1e-2, ftol=1e-3, QUIET=true)
    ssol = Vector{Vector{Vector{Float64}}}(undef, cql)
    sol = zeros(n)
    freq = zeros(n)
    for i = 1:cql
        ssol[i],w = extract_solutions_robust(moment[i], cliquesize[i], d, check=false, rtol=rtol, gtol=gtol, ftol=ftol)
        sol[cliques[i]] += ssol[i][argmax(w)]
        freq[cliques[i]] .+= 1
    end
    sol ./= freq
    if check == true
        if pop !== nothing
            sol = check_solution([sol], lb, pop, x, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
        else
            sol = check_solution([sol], lb, supp, coe, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
        end
    end
    return sol,ssol
end

function extract_csolutions_robust(moment, n, d, cliques, cql, cliquesize; pop=nothing, z=nothing, supp=[], coe=[], lb=nothing, numeq=0, check=false, rtol=1e-2, gtol=1e-2, ftol=1e-3, QUIET=true)
    ssol = Vector{Vector{Vector{ComplexF64}}}(undef, cql)
    sol = zeros(ComplexF64, n)
    freq = zeros(n)
    for i = 1:cql
        ssol[i],w = extract_solutions_robust(moment[i][1], cliquesize[i], d, type=ComplexF64, check=false, rtol=rtol, gtol=gtol, ftol=ftol)
        sol[cliques[i]] += ssol[i][argmax(w)]
        freq[cliques[i]] .+= 1
    end
    sol ./= freq
    if check == true
        if pop !== nothing
            sol = check_solution([sol], lb, pop, z, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
        else
            sol = check_solution([sol], lb, supp, coe, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
        end
    end
    return sol,ssol
end

"""
    sol = extract_solutions_pmo(n, d, p, moment; rtol=1e-2)

Extract a tuple of solutions for a polynomial matrix optimization problem using Henrion-Lasserre's algorithm.

# Input arguments
- `n`: number of variables
- `d`: relaxation order
- `p`: size of the objective matrix
- `moment`: moment matrix
- `rtol`: tolerance for rank

# Output arguments
- `sol`: a tuple of solutions
"""
function extract_solutions_pmo(n, d, p, moment; basis=[], rtol=1e-2)
    U,pivots = rref_with_pivots!(Matrix(moment), rtol)
    U = U[1:length(pivots), :]'
    if isempty(basis)
        basis = get_basis(n, d)
    end
    w = basis[:, ceil.(Int, pivots./p)]
    ind = cmod.(pivots, p)
    lw = size(w, 2)
    if any([sum(w[:,j]) == d for j = 1:lw])
        return nothing
    end
    N = Vector{Matrix{Float64}}(undef, n)
    for i = 1:n
        kk = UInt16[]
        temp = zeros(UInt8, n)
        temp[i] = 1
        for j = 1:lw
            xwj = w[:,j] + temp
            loc = bfind_to(basis, size(basis, 2), xwj, n)
            kk = push!(kk, (loc-1)*p + ind[j])
        end
        N[i] = U[kk,:]
    end
    rands = rand(n)
    rands = rands/sum(rands)
    L = schur(sum(rands[i]*N[i] for i in 1:n)).Z
    sol = [[L[:,i]'*N[j]*L[:,i] for j = 1:n] for i = 1:lw]
    return sol
end

"""
    sol = extract_solutions_pmo_robust(n, d, p, moment; rtol=1e-2)

Extract a tuple of solutions for a polynomial matrix optimization problem using the GNS robust algorithm.

# Input arguments
- `n`: number of variables
- `d`: relaxation order
- `p`: size of the objective matrix
- `moment`: moment matrix
- `rtol`: tolerance for rank

# Output arguments
- `sol`: a tuple of solutions
"""
function extract_solutions_pmo_robust(n, d, p, moment; basis=[], rtol=1e-2)
    if isempty(basis)
        basis = get_basis(n, d)
    end
    ls = p*binomial(n+d-1, n)
    N = Vector{Matrix{Float64}}(undef, n)
    for i = 1:n
        N[i] = zeros(Float64, ls, ls)
        temp = zeros(UInt8, n)
        temp[i] = 1
        for j = 1:ls, k = j:ls
            loc = bfind_to(basis, size(basis, 2), basis[:,ceil(Int, k/p)] + temp, n)
            N[i][j,k] = moment[j, (loc-1)*p + cmod(k, p)]
        end
        N[i] = Symmetric(N[i],:U)
    end
    F = svd(moment[1:ls, 1:ls])
    S = F.S[F.S/F.S[1] .> rtol]
    S = sqrt.(S).^(-1)
    for i = 1:n
        N[i] = Diagonal(S)*F.Vt[1:length(S),:]*N[i]*F.U[:,1:length(S)]*Diagonal(S)
    end
    rands = rand(n)
    rands = rands/sum(rands)
    L = schur(sum(rands[i]*N[i] for i in 1:n)).Z
    sol = [[L[:,i]'*N[j]*L[:,i] for j = 1:n] for i = 1:length(S)]
    return sol
end

function extract_weight_matrix(n, d, q, sol, moment; tol=1e-2)
    ind = [1]
    for i = 2:length(sol)
        if abs(norm(sol[i]-sol[i-1])) > tol
            push!(ind, i)
        end
    end
    sol = sol[ind]
    basis = get_basis(n, d)
    A = kron([prod(v.^item) for item in eachcol(basis), v in sol], Diagonal(ones(q)))
    ind = rref_with_pivots!(Matrix(A'), tol)[2]
    W = A[ind, :]^(-1)*moment[ind, 1:q]
    return [W[(i-1)*q+1:i*q, :] for i = 1:length(sol)]
end

function bfind_to(A, l, a, n)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        order = comp_to(A[:,mid], a, n)
        if order == 0
           return mid
        elseif order < 0
           low = mid + 1
        else
           high = mid - 1
        end
    end
    return nothing
end

function comp_to(a, b, n)
    if sum(a) < sum(b)
        return -1
    elseif sum(a) > sum(b)
        return 1
    else
        i = 1
        while i <= n
            if a[end+1-i] < b[end+1-i]
                return -1
            elseif a[end+1-i] > b[end+1-i]
                return 1
            else
                i += 1
            end
        end
    end
    return 0
end

# extract an approximate solution from the moment matrix
function approx_sol(moment, lb, n, cliques, cql, cliquesize, supp, coe; numeq=0, gtol=1e-2, ftol=1e-3, QUIET=true)
    qsol = Float64[]
    A = zeros(sum(cliquesize), n)
    q = 1
    for k = 1:cql
        cqs = cliquesize[k]
        F = eigen(moment[k], cqs+1:cqs+1)
        temp = sqrt(F.values[1])*F.vectors[:,1]
        if temp[1] == 0
            temp = zeros(cqs)
        else
            temp = temp[2:cqs+1]./temp[1]
        end
        append!(qsol, temp)
        for j = 1:cqs
            A[q,cliques[k][j]] = 1
            q += 1
        end
    end
    sol = (A'*A)\(A'*qsol)
    flag = 1
    ub = eval(supp[1], coe[1], sol)
    gap = abs(lb-ub)/max(1, abs(ub))
    if check_optimality(sol, lb, supp, coe, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
        flag = 0
        println("------------------------------------------------")
        @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
        println("Successfully extracted one globally optimal solution.")
        println("------------------------------------------------")
    end
    return sol,gap,flag
end

function check_optimality(sol, lb, pop::Vector{DP.Polynomial{V, M, T}}, x; numeq=0, gtol=1e-2, ftol=1e-3, QUIET=true) where {V, M, T<:Number}
    flag = 1
    ub = pop[1](x => sol)
    gap = abs(lb-ub)/max(1, abs(ub))
    if QUIET == false
        @printf "Global optimality gap = %.6f%%!\n" 100*gap
    end
    if gap > gtol
        flag = 0
    end
    m = length(pop) - 1 - numeq
    if m > 0
        vio = [pop[j](x => sol) for j = 2:m+1]
        mvio = max(-minimum(real.(vio)), 0)
        if mvio > ftol
            flag = 0
        end
        if QUIET == false
            println("Maximal inequality violation = ", mvio)
        end
    end
    if numeq > 0
        vio = [pop[j](x => sol) for j = m+2:length(pop)]
        mvio = maximum(abs.(vio))
        if mvio > ftol
            flag = 0
        end
        if QUIET == false
            println("Maximal equality violation = ", mvio)
        end
    end
    return flag == 1
end

function check_optimality(sol, lb, supp::Union{Vector{Vector{Vector{UInt16}}}, Vector{Vector{Vector{Vector{UInt16}}}}}, coe; numeq=0, gtol=1e-2, ftol=1e-3, QUIET=true)
    flag = 1
    ub = eval(supp[1], coe[1], sol)
    gap = abs(lb-ub)/max(1, abs(ub))
    if QUIET == false
        @printf "Global optimality gap = %.6f%%!\n" 100*gap
    end
    if gap > gtol
        flag = 0
    end
    m = length(supp) - 1 - numeq
    if m > 0
        vio = [eval(supp[j], coe[j], sol) for j = 2:m+1]
        mvio = max(-minimum(vio), 0)
        if mvio > ftol
            flag = 0
        end
        if QUIET == false
            println("Maximal inequality violation = ", mvio)
        end
    end
    if numeq > 0
        vio = [eval(supp[j], coe[j], sol) for j = m+2:length(supp)]
        mvio = maximum(abs.(vio))
        if mvio > ftol
            flag = 0
        end
        if QUIET == false
            println("Maximal equality violation = ", mvio)
        end
    end
    return flag == 1
end

function check_solution(candidate_sol, lb, pop::Vector{DP.Polynomial{V, M, T}}, x; numeq=0, gtol=1e-2, ftol=1e-3, QUIET=true) where {V, M, T<:Number}
    sol = typeof(candidate_sol[1])[]
    for (i, cand) in enumerate(candidate_sol)
        if QUIET == false
            println("------------------------------------------------")
            println("Check the $i-th candidate solution:")
        end
        if check_optimality(cand, lb, pop, x, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
            push!(sol, cand)
        end
    end
    nsol = length(sol)
    if nsol > 0
        println("------------------------------------------------")
        if nsol == 1
            ub = MP.polynomial(pop[1])(x => sol[1])
            gap = abs(lb-ub)/max(1, abs(ub))
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
            println("Successfully extracted one globally optimal solution.")
        else
            ub = minimum([real(MP.polynomial(pop[1])(x => s)) for s in sol])
            gap = abs(lb-ub)/max(1, abs(ub))
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
            println("Successfully extracted ", nsol, " globally optimal solutions.")
        end
        println("------------------------------------------------")
    else
        sol = nothing
    end
    return sol
end

function check_solution(candidate_sol, lb, supp::Union{Vector{Vector{Vector{UInt16}}}, Vector{Vector{Vector{Vector{UInt16}}}}}, coe; numeq=0, gtol=1e-2, ftol=1e-3, QUIET=true)
    sol = typeof(candidate_sol[1])[]
    for (i, cand) in enumerate(candidate_sol)
        if QUIET == false
            println("------------------------------------------------")
            println("Check the $i-th candidate solution:")
        end
        if check_optimality(cand, lb, supp, coe, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=QUIET)
            push!(sol, cand)
        end
    end
    nsol = length(sol)
    if nsol > 0
        println("------------------------------------------------")
        if nsol == 1
            ub = eval(supp[1], coe[1], sol[1])
            gap = abs(lb-ub)/max(1, abs(ub))
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
            println("Successfully extracted one globally optimal solution.")
        else
            ub = minimum([eval(supp[1], coe[1], s) for s in sol])
            gap = abs(lb-ub)/max(1, abs(ub))
            @printf "Global optimality certified with relative optimality gap %.6f%%!\n" 100*gap
            println("Successfully extracted ", nsol, " globally optimal solutions.")
        end
        println("------------------------------------------------")
    else
        sol = nothing
    end
    return sol
end
