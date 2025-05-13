using Random

"""
    get_basis_indices(n::Int64, d::Int64) -> Matrix{Int64}

Generate the multi-index representation of the monomial basis for `n` variables 
of degree `d`.
    
# Arguments
- `n::Int64`: The number of variables.
- `d::Int64`: The maximum degree of the monomials.
    
# Returns
- `Matrix{Int64}`: A matrix where each column represents a multi-index for a 
monomial. The rows correspond to the variables, and the columns correspond 
to the exponents in each monomial. The total degree of each column sums to is always smaller or equal to d. 
""" 
function get_basis_indices(n::Int64, d::Int64) 
    sd = binomial(n + d, d)
    basis = zeros(Int64, n, sd)  
    
    i = Int64(0)  # Tracks the degree being processed
    t = Int64(1)  # Index of the current basis element
    
    while i <= d
        if basis[n, t] == i
            if i < d
                @inbounds begin
                    t += 1
                    basis[1, t] = i + 1  
                    i += 1
                end
            else
                i += 1
            end
        else
            j = Int64(1)
            while basis[j, t] == 0
                @inbounds j += 1
            end
            
            @inbounds begin
                t += 1
                basis[:, t] = basis[:, t - 1]  
                
                if j == 1
                    basis[1, t] -= 1
                    basis[2, t] += 1
                else
                    basis[1, t] = basis[j, t] - 1
                    basis[j, t] = 0
                    basis[j + 1, t] += 1
                end
            end
        end
    end
    
    return basis
end

"""
   Generate the monomial basis of degree `d` for a given set of variables.
    
# Arguments
- `vars::Vector`: A vector of symbolic variables.
- `d::Int64`: The total degree of the monomials to generate.
    
# Returns
- `Vector`: A vector containing the monomial basis of degree `d` of size binomial(length(vars)+d,d). 
Each monomial is represented as a product of the variables raised to their corresponding powers.
"""
function get_basis(vars::Vector{DP.Variable}, d::Int64)
    b_indices = get_basis_indices(length(vars),d)

    monomial_basis = [] 
    for col in 1:size(b_indices, 2)
        monomials = vars .^ b_indices[:, col]
        push!(monomial_basis, monomials)
    end
    
    monob = prod.(monomial_basis)
    return monob
end

"""
    get_monomial_exponents(monomial, vars) -> Vector{Int}

Compute the exponent vector for a given monomial with respect to a specified set of variables.
    
# Arguments
- `monomial`: The input monomial for which the exponents are to be computed. This should be a symbolic term.
- `vars::Vector`: A vector containing the variables (e.g., `[x1, x2, ..., xn]`) to be used for indexing the exponents.
    
# Returns
- `Vector{Int}`: A vector of integers representing the exponents of the monomial for each variable in `vars`.    
"""
function get_monomial_exponents(monomial, vars)
    exponents = zeros(Int, length(vars))
    
    for (i, var) in enumerate(vars)
        exponents[i] = degree(monomial, var)
    end
    
    return exponents
end

"""
    find_moment_matrix_entry(vars, d, target_exponents) -> Tuple{Int, Int}

Find the position of a monomial (i.e., its associated moment) in the moment matrix of order `d`.
    
# Arguments
- `vars::Vector`: A vector of variables (e.g., `[x1, x2, ..., xn]`).
- `d::Int`: The order of the moment matrix.
- `target_exponents::Vector{Int}`: A vector representing the exponents of the target monomial. The length of this vector should match the number of variables in `vars`.
    
# Returns
- `Tuple{Int, Int}`: A pair `(i, j)` representing the row and column indices in the moment matrix where the target monomial (i.e., its associated moment) appears.       
"""
function find_moment_matrix_entry(vars, d, target_exponents)
    basis_d = get_basis_indices(length(vars), d)     # Basis for indexing the matrix
    basis_2d = get_basis_indices(length(vars), 2*d)   # Basis for the moment matrix entries

    # Step 1: Find the index of the target monomial in the 2d basis
    target_index = nothing
    for col in 1:size(basis_2d, 2)
        if basis_2d[:, col] == target_exponents
            target_index = col
            break
        end
    end
    if target_index === nothing
        error("Target monomial not found in basis for degree 2d.")
    end

    # Step 2: Map the target monomial to pairs of basis monomials
    for i in 1:size(basis_d, 2)
        for j in i:size(basis_d, 2) 
            if basis_d[:, i] + basis_d[:, j] == target_exponents
                return i, j 
            end
        end
    end

    error("Target monomial cannot be represented as a product of basis monomials.")
end

"""
    extract_moment_vector(vars, d, monomials, Mm) -> Vector

Extract the moment vector of order 2d corresponding to a given set of monomials from a moment matrix.
    
# Arguments
- `vars::Vector`: A vector of variables (e.g., `[x1, x2, ..., xn]`).
- `d::Int`: The order of the moment matrix.
- `monomials::Vector`: A vector of monomials for which the corresponding moments are to be extracted. Each monomial should be expressed in terms of the variables in `vars`.
- `Mm::Matrix`: The moment matrix of order `d`, containing moments for the given variables.

# Returns
- `Vector`: A vector containing the extracted moments, where the i-th entry corresponds to the moment associated to the i-th monomial in `monomials`.
"""
function extract_moment_vector(vars, d, monomials, Mm)    
    momemnt_vec_2d=[]

    for i=1:length(monomials)
        temp_exponent = get_monomial_exponents(monomials[i], vars)
        temp_row,temp_column=find_moment_matrix_entry(vars, d, temp_exponent)
        push!(momemnt_vec_2d, Mm[temp_row,temp_column])
    end
    
    return momemnt_vec_2d
end

"""
    Constructs Christoffel polynomial of order d (CDK) using Singular Value Decomposition (SVD) for the input moment matrix `Mm` of order smaller or equal than d. 
    
# Arguments
- `vars`: vector of variables, e.g., @polyvar x[1:5]
- `dc`:  degree of the CDK (has to be smaller than or equal to the order of the moment matrix)
- `Mm::Matrix`: The input moment matrix from which which the CDK is constructed.
- `threshold::Float64` (optional): A threshold value for filtering eigenvalues (i.e., deciding the numerical rank of `Mm`) . Default is `0.001`.
    
# Returns
- `p_alpha_squared[1:negativeEV]`: Polynomials in the kernel of the moment matrix.
- `cdk`: The constructed CDK (positive part of the moment matrix).
- `positiveEV`: The number of positive eigenvalues of the moment matrix.
- `minimum(Qval)`: The minimum eigenvalue of the moment matrix.
- `negativeEV`: The number of negative eigenvalues of the moment matrix.    
"""
function construct_CDK(vars, dc, Mm; threshold=0.001)
    sd=binomial(length(vars)+dc,dc)

    eigen_result = eigen(Mm[1:sd,1:sd])

    Qvec = eigen_result.vectors
    Qval = eigen_result.values

    p_alpha_squared=[]   #Vector of orhonormal polynomials 
    
    for i=1:length(Qval)
        push!(p_alpha_squared, (get_basis(vars, dc)'*Qvec[:,i])^2)
    end

    threshold=0.001

    positiveEV = count(x -> abs(x) > threshold, Qval)
    negativeEV=count(x -> abs(x) <= threshold, Qval)

    if positiveEV==length(Qval)
        cdk=sum(p_alpha_squared[i]/(Qval[i]) for i=1:length(Qval))
    else
        cdk=sum(p_alpha_squared[i]/(Qval[i]) for i=(length(Qval)-positiveEV+1):length(Qval))
    end 

    return p_alpha_squared[1:negativeEV], cdk, positiveEV,minimum(Qval),negativeEV
        
end

"""
    construct_CDK_cs(x, dc, Mm, cliques, threshold=0.001) -> Tuple
    Construct the `Cdk`, kernel, and eigenvalue-related statistics for a given set of cliques and moment matrices.
    
# Arguments
- `x::Vector`: A vector of polynomial variables.
- `dc::Int`: The degree of the Christoffel polynomial constraints. It must be  <= than the order of relaxation d.
- `Mm::Vector`: A vector of moment matrices, one for each clique.
- `cliques::Vector{Vector{Int}}`: A list of cliques, where each clique is a vector of indices corresponding to variables in `x`.
- `threshold::Float64` (optional): A small positive value used to distinguish between eigenvalues considered positive and negative. Defaults to `0.001`.
    
# Returns
- `Tuple`: A tuple containing:
- `Kernel::Vector`: A vector of squared orthonormal polynomial values corresponding to the polynomials in the kernel of clique-based moment matrices.
- `Cdk::Vector`: A vector of Christoffel poolynomials of degree 2*dc, one for each clique.
- `posEV::Vector`: A vector of counts of eigenvalues above the threshold for each clique.
- `minEV::Vector`: A vector of the smallest eigenvalue for each clique.
- `negEV::Vector`: A vector of counts of eigenvalues below or equal to the threshold for each clique. 
"""
function construct_CDK_cs(x, dc, Mm, cliques, threshold=0.001)  
    Cdk = []
    Kernel = []
    posEV = []
    negEV = [] 
    minEV = []

    cliques_basis = [x[cliques[i]] for i=1:length(cliques)] 
    local_basis = []
    for i in 1:length(cliques)
        push!(local_basis, get_basis(cliques_basis[i], dc))
    end

    for i = 1:length(cliques)
        sd=binomial(length(cliques_basis[i])+dc,dc)

        eigen_result = eigen(Mm[i][1:sd,1:sd])
        
        Qvec = eigen_result.vectors
        Qval = eigen_result.values

        p_alpha_squared = []  # Vector of orthonormal polynomials 
    
        for j = 1:length(Qval)
            push!(p_alpha_squared, (local_basis[i]' * Qvec[:, j])^2)
        end

        positiveEV = count(x -> abs(x) > threshold, Qval)
        negativeEV = count(x -> abs(x) <= threshold, Qval)

        if positiveEV == length(Qval)
            cdk = sum(p_alpha_squared[j] / Qval[j] for j in 1:length(Qval))
        else
            cdk = sum(p_alpha_squared[j] / Qval[j] for j in (length(Qval) - positiveEV + 1):length(Qval))
        end

        push!(Cdk, cdk)
        push!(Kernel, p_alpha_squared[1:negativeEV])
        push!(posEV, positiveEV)
        push!(negEV, negativeEV)
        push!(minEV, minimum(Qval))
    end

    return Kernel, Cdk, posEV, minEV, negEV
end

"""
    Performs iterative bound strengthening (i.e, implements H1) for a given POP.

# Arguments
- `x::Vector`: Initial vector of decision variables.
- `n::Int`: Dimension or size of the problem.
- `d::Int`: Degree of the relaxation.
- `dc::Int`: Degree of the Christoffel polynomials relaxation. dc must be smaller than or equal to d.
- `N::Int`: Maximum number of iterations to perform.
- `eps::Float64`: Perturbation factor deciding the volume of CDK sublevel sets.
- `gap_tol::Float64` (optional): Threshold for optimality gap tolerance. Default is `0.5`.
- `beta::Float64` (optional): Threshold beta controlling the strength of the Tikhonov regularization of CDK. Default is `0.001`.

# Returns
A list with the following elements:
1. `OPT`: Relaxation bound from each iteration.
2. `gap_history`: A history of the optimality gaps at each iteration.
3. `Moment_matrices`: Moment matrices computed during each iterations used to construct CDK sublevel constraints.
4. `K`: Parameter to keep track of the rank of the moment matrices for each iteration.
5. `Certain_overrestriction`: A measure or flag for over-restriction (when H1 is sure that further feasible set reductions would yield an invalid bound).
6. `ubRef`: Upper bounds with respect to which the optimality gap is measured.
"""
function run_H1(pop,x,d,dc,N,eps,gap_tol=0.5,beta=0.001)
    Moment_matrices=[]
    Certain_overrestriction=false
    
    f=pop[1]
    normalization_factor=maximum(abs.(coefficients(f)))


    println("_________________________ SOLVING THE INITIAL RELAXATION: _________________________")
    println()

    opt, sol, data = cs_tssos_first(pop, x, d, TS=false,CS=false, solution=true)
    
    println()
    println("________________________________________________________________________________________________ ")

    push!(Moment_matrices,data.moment[1])


    ubRef=f(sol)

    initial_gap=(abs((ubRef-opt))/abs(ubRef))*100

    current_gap=initial_gap
    current_opt=opt

    println()
    println("Value of the initial relaxation = ", opt, ".")
    println("Global upper bound = ", f(sol), ".")
    println("Initial optimality gap = ", round(initial_gap,digits=5), "%.")
    println()

    
    OPT=[round(opt,digits=6)]
    gap_history=[round(current_gap,digits=4)]

    i=0

    k=1                                         #To keep track of the rank reduction of the moment matrix
    K=Int.(ones(size(data.moment[1], 1)))



    CDK_positive=Any[]
    CDK_zero=Any[]

    while i < N 
        
        if current_gap < gap_tol
            println("RELAXATION IS GOOD ENOUGH.")
            println()
            break
        end
    
        println()
        printstyled(">>>>>>>>>>>> ___________________ Iteration ", i+1, " ______________________ <<<<<<<<<<<<", color=:red, bold=true)
        println()
        temp = construct_CDK(x, dc, data.moment[1])
        println()


        if temp[3] == size(data.moment[1], 1)       #For full-rank matrices, there are no kernel constraints

            K[temp[3]]=K[temp[3]]+1

            cdk = temp[2]

            coeffs = coefficients(cdk)
            moment_vector = extract_moment_vector(x, dc, monomials(cdk), data.moment[1])
            
            L_ystar = coeffs'*moment_vector
            
            println("L_ystar(CDK) = ",  L_ystar, " and (1-eps)*L_ystar(CDK) = ", (1-eps)*L_ystar, ".")

            popcdk = [pop[1]/maximum(abs.(coefficients(pop[1])))]

            push!(CDK_positive, (1-eps)*L_ystar - cdk)

            for i=1:length(CDK_positive)
                push!(popcdk, CDK_positive[i]/maximum(abs.(coefficients(CDK_positive[i]))))    #Normalization
            end
            

            for i=2:length(pop)
                push!(popcdk, pop[i]/maximum(abs.(coefficients(pop[i]))))    #Normalization - initial constraints are added at the end
                                                                             #In case there are some equality constraints in the original POP, t
                                                                             #his preserves the structure needed to use the parameter "numeq" in TSSOS
            end

            println("Total number of constraints at iteration ",i+1, " is ",  length(popcdk)-1, ".")
            println()
            
            rng = Random.default_rng() 
            Random.seed!(rng, nothing)
            optcdk,solcdk,datacdk=cs_tssos_first(popcdk,x,d,solution=true,TS=false,CS=false)

            if f(solcdk)< ubRef
                ubRef=f(solcdk)
                println("\e[1;32mGlobal upper bound updated to $ubRef.\e[0m")
            end

            current_opt = optcdk*normalization_factor

            if current_opt > ubRef+1e-6                    #to avoid some innacuracies
                println()
                println(" Overrestriction !!! ")
                Certain_overrestriction=true
                break
            end

            data=datacdk
            push!(Moment_matrices,data.moment[1])


            push!(OPT,round(current_opt, digits=6))
            println()
            
            current_gap= (abs((ubRef-current_opt))/abs(ubRef))*100
            println("Current gap = ", round(current_gap,digits=4),"%.")
            
            push!(gap_history,round(current_gap,digits=4))
            println("Gap history = ", gap_history ,"%")
            println("History of lower bounds = ", OPT)
            




        else

            K[temp[3]]=K[temp[3]]+1
                                    
            cdk = temp[2]            
            coeffs = coefficients(cdk)
            
            moment_vector = extract_moment_vector(x, dc, monomials(cdk), data.moment[1])
            
            L_ystar = coeffs'*moment_vector
            
            println("L_ystar(CDK) = ", L_ystar, " and the chosen level set, namely (1-eps)*L_ystar(CDK) = ", (1-eps)*L_ystar, ".")

            Gammas=Float64[]

            for j=1:temp[5]
                push!(Gammas, coefficients(temp[1][j])'* extract_moment_vector(x, dc, monomials(temp[1][j]), data.moment[1]))
            end
            
            popcdk = [pop[1]/maximum(abs.(coefficients(pop[1])))]
            
            push!(CDK_positive, (1-eps)*L_ystar - cdk)

            for i=1:length(CDK_positive)
                push!(popcdk, CDK_positive[i]/maximum(abs.(coefficients(CDK_positive[i]))))    #Normalization
            end
                
            thresholds_kernel=Float64[]

            sum_kernel=sum([beta-temp[1][j] for j=1:temp[5]])

            push!(CDK_zero, sum_kernel)

            for i=1:length(CDK_zero)
                 push!(popcdk, CDK_zero[i]/maximum(abs.(coefficients(CDK_zero[i]))))         #Normalization
            end

            for i=2:length(pop)
                push!(popcdk, pop[i]/maximum(abs.(coefficients(pop[i]))))    #Normalization - initial constraints are added at the end
                                                                             #In case there are some equality constraints in the original POP, t
                                                                             #his preserves the structure needed to use the parameter "numeq" in TSSOS
            end
            
            println("Total number of constraints at iteration $(i+1) = ", length(popcdk)-1, ".")
            println()
                
            rng = Random.default_rng() 
            Random.seed!(rng, nothing)
            optcdk,solcdk,datacdk=cs_tssos_first(popcdk,x,d,solution=true,TS=false,CS=false)

            if f(solcdk) < ubRef
                ubRef=f(solcdk)
                println("\e[1;32mGlobal upper bound updated to $ubRef.\e[0m")
            end

            current_opt = optcdk*normalization_factor
            
            println()
            println("Current opt is ", current_opt, ".")
            println("Current best upper bound is ", ubRef, ".")


            if current_opt > ubRef+1e-6                    #to avoid some innacuracies
                println()
                println(" Overrestriction !!! ")
                Certain_overrestriction=true
                break
            end


            data=datacdk
            push!(Moment_matrices,data.moment[1])
  
            push!(OPT, round(current_opt, digits=6))
            println()
            
            current_gap= (abs((ubRef-current_opt))/abs(ubRef))*100
            println("Current gap = ", round(current_gap,digits=4),"%.")
            
            push!(gap_history,round(current_gap,digits=4))
            println("Gap history = ", gap_history ,"%")
            println("History of lower bounds = ", OPT)
            
            println(" ________________________________________________________________________________________________ ")


            
                
        end
        i=i+1
        
    end
                

    return OPT,gap_history, Moment_matrices, K, Certain_overrestriction,ubRef

            
            
end 

function run_H1CS(pop,x,d,dc,N,eps,gap_tol=0.5,beta=0.0001)

    UBH=[]
    Moment_matrices=[]
    Certain_overrestriction=false

    f=pop[1]
    normalization_factor=maximum(abs.(coefficients(f)))

   
    println("_________________________ SOLVING THE INITIAL RELAXATION: _________________________")
    println()

    opt, sol, data = cs_tssos_first(pop, x, d, TS=false, solution=true)

    
    println(" ________________________________________________________________________________________________ ")

    push!(Moment_matrices,data.moment)
    push!(UBH,f(sol))

    ubRef=f(sol)

    initial_gap=(abs((ubRef-opt))/abs(ubRef))*100

    current_gap=initial_gap
    current_opt=opt

    println()
    println("Value of the initial relaxation = ", opt, ".")
    println("Global upper bound = ", f(sol), ".")
    println("Initial optimality gap = ", round(initial_gap,digits=5), "%.")
    println()

    
    OPT=[round(opt,digits=6)]
    gap_history=[round(current_gap,digits=4)]

    i=0
                                            
    K=[Int.(ones(size(data.moment[i], 1))) for i=1:length(data.cliques)]  #To keep track of the rank reduction of the moment matrix



    CDK_positive=Any[]
    CDK_zero=Any[]

    while i < N 
        
        if current_gap < gap_tol
            println("RELAXATION IS GOOD ENOUGH.")
            println()
            break
        end
    
        println()
        printstyled(">>>>>>>>>>>> ___________________ Iteration ", i+1, " ______________________ <<<<<<<<<<<<", color=:red, bold=true)
        println()
        println()
        

        temp = construct_CDK_cs(x, dc, data.moment, data.cliques)
        

        for k = 1:length(data.cliques)    #Update the CDK constraints for each clique

            if temp[3][k] == size(data.moment[k], 1)      #HIf the moment matrix of the clique k is full rank, there are no kernel constraints

                K[k][temp[3][k]]=K[k][temp[3][k]]+1
                                        
                cdk = temp[2][k]           
                coeffs = coefficients(cdk)
                moment_vector = extract_moment_vector(x[data.cliques[k]], dc, monomials(cdk),data.moment[k])
                
                L_ystar = coeffs'*moment_vector
        
                println("Clique $k: L_ystar(CDK_k) = ", L_ystar, " and the chosen level set, namely (1-eps)*L_ystar(CDK_k) = ", (1-eps)*L_ystar, ".")
        
                push!(CDK_positive, (1-eps)*L_ystar - cdk)

            else 
                    
                K[k][temp[3][k]]=K[k][temp[3][k]]+1
                                        
                cdk = temp[2][k]           
                coeffs = coefficients(cdk)
                moment_vector = extract_moment_vector(x[data.cliques[k]], dc, monomials(cdk),data.moment[k])
                
                L_ystar = coeffs'*moment_vector
        
                println("Clique $k: L_ystar(CDK_$k) = ", L_ystar, " and the chosen level set, namely (1-eps)*L_ystar(CDK_$k) = ", (1-eps)*L_ystar, ".")
        
                push!(CDK_positive, (1-eps)*L_ystar - cdk)
            
                sum_kernel=sum([beta-temp[1][k][j] for j=1:temp[5][k]])
        
                push!(CDK_zero, sum_kernel)
                    
            
                        
            end 
                
        end
            

        popcdk = [pop[1]/maximum(abs.(coefficients(pop[1])))]     #Construct the updated POP
                
        for i=1:length(CDK_positive)
            push!(popcdk, CDK_positive[i]/maximum(abs.(coefficients(CDK_positive[i]))))    #Normalization
        end

        if length(CDK_zero)!=0
            for i=1:length(CDK_zero)
                 push!(popcdk, CDK_zero[i]/maximum(abs.(coefficients(CDK_zero[i]))))         #Normalization of kernel constraints, in case they exist
            end
        end

        for i=2:length(pop)
            push!(popcdk,pop[i]/maximum(abs.(coefficients(pop[i]))))    #Normalization; Again, put all the initial constraints in the end
        end
        
                
    
        
        println("Total number of constraints at iteration $(i+1) = ", length(popcdk)-1, ".")
        println()

    
        rng = Random.default_rng() 
        Random.seed!(rng, nothing)
        optcdk,solcdk,datacdk=cs_tssos_first(popcdk,x,d,solution=true,TS=false)

        push!(UBH,f(solcdk))
    
        if f(solcdk) < ubRef
            ubRef = f(solcdk)
            println("\e[1;32mGlobal upper bound updated to $ubRef.\e[0m")
        end

        current_opt = optcdk*normalization_factor
        
        println()
        println("Current opt is ", current_opt, ".")
        println("Current best upper bound is ", ubRef, ".")
    
    
        if current_opt > ubRef+1e-6                    #to avoid some innacuracies
            println()
            println(" Overrestriction !!! ")
            Certain_overrestriction=true
            break
        end
    
    
        data=datacdk
        push!(Moment_matrices,data.moment)

        push!(OPT, round(current_opt, digits=6))
        println()
        
        current_gap= (abs((ubRef-current_opt))/abs(ubRef))*100
        println("Current gap = ", round(current_gap,digits=4),"%.")
        
        push!(gap_history,round(current_gap,digits=4))
        println("Gap history = ", gap_history ,"%")
        println("History of lower bounds = ", OPT)
        
        println(" ________________________________________________________________________________________________ ")
    
                
                    
        i=i+1
            
        
    end
                

    return OPT, gap_history, Moment_matrices, K, Certain_overrestriction, UBH

            
end

"""
    extract_marginal_indices(k::Int, indices::Matrix{Int}) -> Vector{Int}
    
Extract the column indices corresponding to the marginal monomials for a given variable  x_k.
    
# Arguments
- `k::Int`: The index of the variable x_k  for which marginal monomial indices are to be extracted.
- `indices::Matrix{Int}`: A matrix of monomial basis indices, where each column represents a monomial, and each row represents the exponents of a variable in that monomial.

# Returns
- `Vector{Int}`: A vector containing the indices of columns in `indices` that correspond to marginal monomials for x_k.
"""
function extract_marginal_indices(k::Int,indices::Matrix{Int})
    n, sd = size(indices)  
    
    marginal_indices = []

    for col in 1:sd
        if indices[k, col] > 0 && all(indices[i, col] == 0 for i in 1:n if i != k)
            push!(marginal_indices, col)  
        end
    end

    return  marginal_indices  
end

"""
    find_clique_and_local_index(k, cliques) -> Vector{Int}

Find the clique and local index of an element in a list of cliques.
    
# Arguments
- `k::Int`: The element/variable to search for.
- `cliques::Vector{Vector{Int}}`: A list of cliques, where each clique is a vector of integers representing decision variable indices.
    
# Returns
- `Vector{Int}`: A two-element vector `[global_clique_index, local_index]`, where:
- `global_clique_index` is the index of the clique (in `cliques`) containing the element `k`.
- `local_index` is the position of `k` within the identified clique.    
"""  
function find_clique_and_local_index(k, cliques)  
    for i in 1:length(cliques)
        clique = cliques[i]

        for local_idx in 1:length(clique)
            if clique[local_idx] == k
                return [i, local_idx]
            end
        end
    end
    error("The index k=$k does not belong to any clique!")
end

"""
    Constructs marginal Christoffel polynomials (CDK) using Singular Value Decomposition (SVD) for the input matrix `Mm`. 
    
# Arguments
- `vars`: vector of variables, e.g., @polyvar x[1:5].
- `k::vector` - Index of the decision variable for which CDK needs to be constructed.
- `dc::Int`: - Degree of the marginal Christoffel polynomials (has to be smaller than or equal to d - order of the moment matrix).
- `Mm::Matrix`: The input moment matrix from which the marginal CDK is constructed.
- `threshold::Float64` (optional): A threshold value for filtering eigenvalue (i.e., deciding the numerical rank of `Mm`) . Default is `0.001`.
    
# Returns
- `cdk`: The constructed marginal CDK.
- `p_alpha_squared[1:negativeEV]`: Polynomials in the kernel of the marginal moment matrix.
- `positiveEV`: The set of positive eigenvalues of the marginal moment matrix.
- `negativeEV`: The set of negative eigenvalues of the marginal moment matrix.
- `minimum(Qval)`: The minimum eigenvalue of the marginal moment matrix.    
"""
function construct_marginal_CDK(vars, k, dc, Mm, threshold=0.001)
    local_indices=extract_marginal_indices(k,get_basis_indices(length(vars),dc))
    eigen_result = eigen(Mm[[1,local_indices...],[1,local_indices...]])

    Qvec = eigen_result.vectors
    Qval = eigen_result.values

    p_alpha_squared=[]   #Vector of orhonormal polynomials 
    
    for j=1:length(Qval)
        push!(p_alpha_squared, (get_basis([vars[k]], dc)'*Qvec[:,j])^2)
    end

    positiveEV = count(x -> abs(x) > threshold, Qval)
    negativeEV=count(x -> abs(x) <= threshold, Qval)

    if positiveEV==length(Qval)
        cdk=sum(p_alpha_squared[i]/(Qval[i]) for i=1:length(Qval))
    else
        cdk=sum(p_alpha_squared[i]/(Qval[i]) for i=(length(Qval)-positiveEV+1):length(Qval))
    end 

    return p_alpha_squared[1:negativeEV], cdk, positiveEV, negativeEV, minimum(Qval)
        
end

"""
    construct_marginal_CDK_cs(vars, k, dc, Mm, cliques; threshold=0.001) -> Tuple

Constructs the marginal Christoffel polynomials (CDK) for a given variable x_k using eigenvalue decomposition, using clique-based moment matrices.
    
# Arguments
- `vars::Vector`: A vector of polynomial variables (e.g., created with `@polyvar x[1:n]`).
- `k::Int`: The global index of the variable x_k  for which the marginal CDK is to be constructed.
- `dc::Int`: The degree of the marginal Christoffel polynomials. Must be less than or equal to the order of the associated moment matrix.
- `Mm::Vector{Matrix}`: A vector of moment matrices, where each matrix corresponds to a clique.
- `cliques::Vector{Vector{Int}}`: A list of cliques, where each clique is a subset of the indices of `vars`.
- `threshold::Float64` (optional): A threshold for filtering eigenvalues, used to determine the numerical rank of the moment matrix. Default is `0.001`.
    
# Returns
A tuple containing:
1. `p_alpha_squared[1:negativeEV]::Vector`: Polynomials in the kernel of the marginal moment matrix.
2. `cdk::Float64`: The constructed marginal Christoffel polynomial.
3. `positiveEV::Int`: The count of positive eigenvalues of the marginal moment matrix.
4. `negativeEV::Int`: The count of negative eigenvalues of the marginal moment matrix.
5. `minimum(Qval)::Float64`: The smallest eigenvalue of the marginal moment matrix.
"""
function construct_marginal_CDK_cs(vars, k, dc, Mm, cliques, threshold=0.001)
    clique_index, local_index = find_clique_and_local_index(k, cliques)
    
    local_indices=extract_marginal_indices(local_index, get_basis_indices(length(cliques[clique_index]),dc))
    eigen_result = eigen(Mm[clique_index][[1,local_indices...],[1,local_indices...]])
    
    Qvec = eigen_result.vectors
    Qval = eigen_result.values

    p_alpha_squared=[]   #Vector of orhonormal polynomials 
    
    for j=1:length(Qval)
        push!(p_alpha_squared, (get_basis([vars[k]], dc)'*Qvec[:,j])^2)
    end

    positiveEV = count(x -> abs(x) > threshold, Qval)
    negativeEV=count(x -> abs(x) <= threshold, Qval)

    if positiveEV==length(Qval)
        cdk=sum(p_alpha_squared[i]/(Qval[i]) for i=1:length(Qval))
    else
        cdk=sum(p_alpha_squared[i]/(Qval[i]) for i=(length(Qval)-positiveEV+1):length(Qval))
    end 

    return p_alpha_squared[1:negativeEV], cdk, positiveEV, negativeEV, minimum(Qval)
        
end

"""
    Applies H2 algorithm to a given POP instance
    
# Arguments
- `pop::Vector`: Vector defining the problem, the first element being the objective function, the others being constraints.
- `x::Vector`: Initial vector of decision variables 
- `d::Int`: - Degree of the relaxation
- `dc::Int`: - Degree of the marginal Christoffel polynomials (has to be smaller than or equal to d - order of the moment matrix)
- `local_sol::Vector`: Available solution (here obtained from an external local optimizer IPOPT).
- `tau::Float64`: Filtering parameter 
- `beta::Float64` (optional): Threshold beta controling the strength of the Tikhonov regularization of marginal Christoffel polynomials. Default is `0.001`.
    
# Returns
A tuple of the following:
1. `[opt, optcdk]`: A pair containing the optimal solutions of moment relaxation before and after CDK strengthening via H2.
2. `[initial_gap, gapcdk]`: Relative optimality gaps before and after CDK strengthening via H2.
3. `Gammas`: Thresholds used to construct sublevel sets of marginal Christoffel polynomials (CDK).
4. `[ubRef, f(solcdk)]`: The available upper bound before and after CDK strengthening via H2.
5. `data.moment[1]`: The moment matrix from which marginal CDK polynomials were constructed    
"""
function run_H2(pop,x,d,dc,local_sol,tau,beta=0.001)          
    f=pop[1]
    
    println("_________________________ SOLVING THE INITIAL RELAXATION: _________________________")
    println()

    
    opt, sol, data = cs_tssos_first(pop, x, d, TS=false,CS=false, solution=true)

    ubRef=minimum([f(sol), f(local_sol)])
    initial_gap=(abs((ubRef-opt))/abs(ubRef))*100

    current_gap=initial_gap
    current_opt=opt

    println()
    println("Value of the initial relaxation = ", opt, ".")
    println("Global upper bound = ", ubRef, ".")
    println("Initial optimality gap = ", round(initial_gap,digits=5), "%.")
    println()

    println("_________________________________________________________________________________________________ ")

    CDK_marg_SVD=[]
    for i=1:n
        push!(CDK_marg_SVD,construct_marginal_CDK(x, i, dc, data.moment[1]))
    end
    
    popcdk = [f]
    Gammas=[]
    gk=[]

    for i=1:n
        temp=CDK_marg_SVD[i]
        gamma_i=temp[2](local_sol[i])+1e-6 # to avoid numerical issues
        push!(gk,[temp[1][l](local_sol[i]) for l=1:length(temp[1])])
        push!(Gammas,gamma_i)   

        if gamma_i < tau
            if temp[4] > 0
                push!(popcdk,sum([beta-temp[1][j] for j=1:length(temp[1])]), gamma_i-temp[2])
            else
                push!(popcdk, gamma_i-temp[2])
            end
        end
        
                            
    end
    for i=2:length(pop)
             push!(popcdk,pop[i])
    end
    for i=1:length(popcdk)
        popcdk[i]=popcdk[i]/maximum(abs.(coefficients(popcdk[i])))
    end 
        
    println()
    printstyled("_______________________________ APPLYING H2 _______________________________", color=:red, bold=true)
    println()
    println()
    println("New total number of constraints: ", length(popcdk)-1)
    println()

    rng = Random.default_rng() 
    Random.seed!(rng, nothing)
    optcdk,solcdk,datacdk=cs_tssos_first(popcdk,x,d,solution=true,TS=false,CS=false)
        
    println()
    println("Output: _______________________________________________________________")
    
    optcdk=optcdk*maximum(abs.(coefficients(pop[1])))
    gapcdk= abs((ubRef-optcdk)/ubRef)*100
    println("Thresholds 'gamma_i' = ", Gammas)
    println("\e[1mNEW BOUND\e[0m = ", optcdk)
    println()

    
    return [opt,optcdk],[initial_gap,gapcdk],Gammas,[ubRef,f(solcdk)],data.moment[1]
end 

"""
   Applies H2CS algorithm to a given POP instance

# Arguments
- `pop::Vector`: Vector defining the problem, the first element being the objective function, the others being constraints.
- `x::Vector`: Initial vector of decision variables 
- `d::Int`: - Degree of the relaxation
- `dc::Int`: - Degree of the marginal Christoffel polynomials (has to be smaller than or equal to d - order of the moment matrix)
- `local_sol::Vector`: Available solution (here obtained from an external local optimizer IPOPT).
- `tau::Float64`: Filtering parameter 
- `beta::Float64` (optional): Threshold beta controling the strength of the Tikhonov regularization of marginal Christoffel polynomials. Default is `0.001`.

# Returns
A tuple of the following:
1. `[opt, optcdk]`: A pair containing the optimal solutions of the sparse moment relaxation before and after CDK strengthening via H2CS.
2. `[initial_gap, gapcdk]`: Relative optimality gaps before and after CDK strengthening via H2CS.
3. `Gammas`: Thresholds used to construct sublevel sets of marginal Christoffel polynomials (CDK).
4. `[ubRef, f(solcdk)]`: The available upper bound before and after CDK strengthening via H2CS.
5. `data.moment[1]`: The moment matrix from which marginal CDK polynomials were constructed.
"""
function run_H2CS(pop,x,d,dc,local_sol,tau,beta=0.001)
    f=pop[1]
    
    println("_________________________ SOLVING THE INITIAL RELAXATION: _________________________")
    println()

    opt, sol, data = cs_tssos_first(pop, x, d, TS=false, solution=true)

    ubRef=minimum([f(sol), f(local_sol)])
    initial_gap=(abs((ubRef-opt))/abs(ubRef))*100

    current_gap=initial_gap
    current_opt=opt

    println()
    println("Value of the initial relaxation = ", opt, ".")
    println("Global upper bound = ", ubRef, ".")
    println("Initial optimality gap = ", round(initial_gap,digits=5), "%.")
    println()

    println("_________________________________________________________________________________________________ ")

    CDK_marg_SVD=[]
    
    for i=1:length(x)
        push!(CDK_marg_SVD,construct_marginal_CDK_cs(x, i, dc, data.moment, data.cliques))
    end
    
    popcdk = [f]
    Gammas=[]
    gk=[]

    for i=1:length(x)
        temp=CDK_marg_SVD[i]
        gamma_i=temp[2](local_sol[i])+1e-6 # to avoid numerical issues
        push!(gk,[temp[1][l](local_sol[i]) for l=1:length(temp[1])])
        push!(Gammas,gamma_i)   

        if gamma_i < tau
            if temp[4] > 0
                push!(popcdk,sum([beta-temp[1][j] for j=1:length(temp[1])]), gamma_i-temp[2])
            else
                push!(popcdk, gamma_i-temp[2])
            end
        end
        
                            
    end
    for i=2:length(pop)
             push!(popcdk,pop[i])
    end
    for i=1:length(popcdk)
        popcdk[i]=popcdk[i]/maximum(abs.(coefficients(popcdk[i])))
    end 
    
    println()
    printstyled("_______________________________ APPLYING H2 _______________________________", color=:red, bold=true)
    println()
    println()
    println("New total number of constraints: ", length(popcdk)-1)
    println()

    rng = Random.default_rng() 
    Random.seed!(rng, nothing)
    optcdk,solcdk,datacdk=cs_tssos_first(popcdk,x,d,solution=true,TS=false)
        
    println()
    println("Output: _______________________________________________________________")
    
    optcdk=optcdk*maximum(abs.(coefficients(pop[1])))
    gapcdk= abs((ubRef-optcdk)/ubRef)*100
    println("Thresholds 'gamma_i' = ", Gammas)
    println("\e[1mNEW BOUND\e[0m = ", optcdk)
    println()
    
    return [opt,optcdk],[initial_gap,gapcdk],Gammas, [ubRef,f(solcdk)],[data.moment,datacdk.moment]     
            
end 
