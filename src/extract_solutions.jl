using RowEchelon

"""
    sol = extract_solutions(pop, x, d, opt, moment; numeq=0, nb=0, tol=1e-4)

Extract a set of solutions for polynomial optimization.

# Input arguments
- `pop`: polynomial optimization problem
- `x`: set of variables
- `d`: relaxation order
- `opt`: optimum
- `moment`: moment matrix
- `numeq`: number of equality constraints
- `nb`: number of binary variables
- `tol`: tolerance to obtain the column echelon form

# Output arguments
- `sol`: a set of solutions
"""
function extract_solutions(pop, x, d, opt, moment; numeq=0, nb=0, tol=1e-2)
    n = length(x)
    if rank(moment, tol) == 1
        sol = moment[2:n+1, 1]
        println("------------------------------------------------")
        println("Global optimality certified!")
        println("Successfully extracted ", 1 ," globally optimal solution.")
        println("------------------------------------------------")
    else
        U,pivots = rref_with_pivots!(Matrix(moment), tol)
        U = U[1:length(pivots), :]'
        basis = get_basis(n, d)
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
        sol = Vector{Float64}[]
        m = length(pop) - 1
        for i = 1:lw
            atom = [L[:,i]'*N[j]*L[:,i] for j = 1:n]
            flag = 1
            println("------------------------------------------------")
            println("Check atom ", i)
            gap = pop[1](x => atom) - opt
            println("Global optimality gap = ", gap)
            if abs(gap) > 1e-4
                flag = 0
            end
            if m - numeq > 0
                vio = [pop[j](x => atom) for j = 2:m+1-numeq]
                mvio = max(-minimum(vio), 0)
                if mvio > 1e-3
                    flag = 0
                end
                println("Maximal inequality violation = ", mvio)
            end
            if numeq > 0
                vio = [pop[j](x => atom) for j = m-numeq+2:m+1]
                mvio = maximum(abs.(vio))
                if mvio > 1e-3
                    flag = 0
                end
                println("Maximal equality violation = ", mvio)
            end
            if flag == 1
                push!(sol, atom)
            end
        end
        nsol = length(sol)
        if nsol > 0
            println("------------------------------------------------")
            println("Global optimality certified!")
            println("Successfully extracted ", nsol, " globally optimal solutions.")
            println("------------------------------------------------")
        else
            sol = nothing
        end
    end
    return sol
end

function extract_solutions_robust(n, d, moment; tol=1e-2)
    basis = get_basis(n, d)
    ls = binomial(n+d-1, n)
    N = Vector{Matrix{Float64}}(undef, n)
    for i = 1:n
        N[i] = zeros(Float64, ls, ls)
        temp = zeros(UInt8, n)
        temp[i] = 1
        for j = 1:ls, k = j:ls
            loc = bfind_to(basis, size(basis, 2), basis[:,k] + temp, n)
            N[i][j,k] = moment[j, loc]
        end
        N[i] = Symmetric(N[i],:U)
    end
    F = svd(moment[1:ls, 1:ls])
    S = F.S[F.S .> tol]
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

"""
    sol = extract_solutions_pmo(n, d, p, moment; tol=1e-4)

Extract a set of solutions for polynomial matrix optimization.

# Input arguments
- `n`: number of variables
- `d`: relaxation order
- `p`: size of the objective matrix
- `moment`: moment matrix
- `tol`: tolerance to obtain the column echelon form

# Output arguments
- `sol`: a set of solutions
"""
function extract_solutions_pmo(n, d, p, moment; tol=1e-2)
    U,pivots = rref_with_pivots!(Matrix(moment), tol)
    U = U[1:length(pivots), :]'
    basis = get_basis(n, d)
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

function extract_solutions_robust_pmo(n, d, p, moment; tol=1e-2)
    basis = get_basis(n, d)
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
    S = F.S[F.S .> tol]
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
