# find the location of an entry a in the sorted sequence A
function bfind(A, a)
    low = 1
    high = ndims(A) == 2 ? size(A, 2) : length(A)
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        temp = ndims(A) == 2 ? A[:, mid] : A[mid]
        if temp == a
           return mid
        elseif temp < a
           low = mid + 1
        else
           high = mid - 1
        end
    end
    return nothing
end

function bfind_rev(A, a)
    low = 1
    high = length(A)
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if A[mid] == a
           return mid
        elseif A[mid] < a
            high = mid - 1
        else
            low = mid + 1
        end
    end
    return nothing
end

function nbfind(A, a)
    low = 1
    high = length(A)
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        order = comp(A[mid], a)
        if order == 0
           return mid
        elseif order == -1
            low = mid + 1
        else
            high = mid - 1
        end
    end
    return nothing
end

function comp(a, b)
    if length(a) < length(b)
        return -1
    elseif length(a) > length(b)
        return 1
    elseif a < b
        return -1
    elseif a > b
        return 1
    else
        return 0
    end
end

function sadd(a::Vector{UInt16}...; nb=0)
    w = sort!(vcat(a...))
    if nb > 0
        i = 1
        while i < length(w)
            if w[i] <= nb
                if w[i] == w[i+1]
                    deleteat!(w, i:i+1)
                else
                    i += 1
                end
            else
                break
            end
        end
    end
    return w
end

function sadd(a::Tuple{Vector{UInt16},Vector{UInt16}}, b::Tuple{Vector{UInt16},Vector{UInt16}})
    return tuple(sort!([a[1]; b[1]]), sort!([a[2]; b[2]]))
end

function sadd(a::Tuple{Vector{UInt16},Vector{UInt16}}, b::Tuple{Vector{UInt16},Vector{UInt16}}, c::Tuple{Vector{UInt16},Vector{UInt16}})
    return tuple(sort!([a[1]; b[1]; c[1]]), sort!([a[2]; b[2]; c[2]]))
end

function get_signsymmetry(npop::Vector{poly{T}}, n) where {T<:Union{Number,AffExpr,GenericAffExpr}}
    supp = zeros(UInt8, sum(length(p.supp) for p in npop), n)
    l = 1
    for p in npop, bi in p.supp
        for k in bi
            supp[l, k] += 1
        end
        l += 1
    end
    return nullspace(matrix(GF(2), unique(supp, dims=1)))[2]
end

function get_signsymmetry(pop::Vector{Poly{T}}, x) where {T<:Union{Number,AffExpr,GenericAffExpr}}
    return get_signsymmetry([poly(p, x) for p in pop], length(x))
end

# generate the standard monomial basis
function get_basis(var::Vector{Int}, d::Int; nb=0, lead=[])
    lb = binomial(length(var)+d, d)
    basis = Vector{Vector{UInt16}}(undef, lb)
    basis[1] = UInt16[]
    i = 0
    t = 1
    while i < d + 1
        t += 1
        if sum(basis[t-1]) == var[end]*i
           if i < d
               basis[t] = var[1]*ones(UInt16, i+1)
           end
           i += 1
        else
            j = bfind(var, basis[t-1][1])
            basis[t] = copy(basis[t-1])
            ind = findfirst(x->basis[t][x]!=var[j], 1:length(basis[t]))
            if ind === nothing
                ind = length(basis[t]) + 1
            end
            if j != 1
                basis[t][1:ind-2] = var[1]*ones(UInt16, ind-2)
            end
            basis[t][ind-1] = var[j+1]
        end
    end
    if nb > 0
        ind = [!any([item[j] == item[j+1] && item[j] <= nb for j = 1:length(item)-1]) item in basis]
        basis = basis[ind]
    end
    if !isempty(lead)
        ind = map(item -> !divide(item, lead), basis)
        basis = basis[ind]
    end
    return basis
end

function get_basis(n::Int, d::Int; nb=0, lead=[])
    return get_basis(Vector(1:n), d, nb=nb, lead=lead)
end

function get_conjugate_basis(var::Vector{Int}, d::Int; nb=0)
    temp = get_basis(var, d)
    basis = vec([tuple(item1, item2) for item1 in temp, item2 in temp])
    basis = basis[[length(item[1]) + length(item[2]) <= d for item in basis]]
    if nb > 0
        basis = reduce_unitnorm.(basis, nb)
        unique!(basis)
    end
    sort!(basis)
    return basis
end

function newton_basis(n, d, supp; tol=1e-5, solver="Mosek")
    nsupp = zeros(Int, n, length(supp))
    for i = 1:length(supp), j in supp[i]
        nsupp[j, i] += 1
    end
    basis = get_basis(n, d)
    lb = length(basis)
    A0 = [-1/2*nsupp' ones(length(supp),1)]
    t = 1
    ind = Vector(1:lb)
    nsupp = sortslices(nsupp, dims=2)
    while t <= lb
        i = ind[t]
        bi = zeros(Int, n)
        for j in basis[i]
            bi[j] += 1
        end
        if bfind(nsupp, 2*bi) !== nothing
            t += 1
        else
            if solver == "Mosek"
                model = Model(optimizer_with_attributes(Mosek.Optimizer))
            elseif solver == "SDPT3"
                model = Model(optimizer_with_attributes(SDPT3.Optimizer))
            elseif solver == "SDPNAL"
                model = Model(optimizer_with_attributes(SDPNAL.Optimizer))
            elseif solver == "COSMO"
                model = Model(optimizer_with_attributes(COSMO.Optimizer))
            else
                @error "The solver is currently not supported!"
                return nothing
            end
            set_optimizer_attribute(model, MOI.Silent(), true)
            @variable(model, x[1:n+1], lower_bound=-10, upper_bound=10)
            @constraint(model, [A0; [bi' -1]]*x .<= zeros(length(supp)+1))
            @objective(model, Min, [bi' -1]*x)
            optimize!(model)
            vx = value.(x)
            if abs(objective_value(model)) <= tol && sum(abs.(vx)) <= tol
                t += 1
            else
                if abs(objective_value(model)) <= tol && sum(abs.(vx)) > tol
                   t += 1
                else
                   lb -= 1
                   deleteat!(ind, t)
                end
                r = t
                while lb >= r
                    j = ind[r]
                    bj = zeros(Int, n)
                    for k in basis[j]
                        bj[k] += 1
                    end
                    if [bj' -1]*vx <= -tol
                        lb -= 1
                        deleteat!(ind, r)
                    else
                        r += 1
                    end
                end
            end
        end
    end
    return basis[ind]
end

function generate_basis!(supp, basis)
    sort!(supp)
    unique!(supp)
    ind = Int[]
    for i = 1:length(basis), j = i:length(basis)
        if bfind(supp, sadd(basis[i], basis[j])) !== nothing
             push!(ind, i, j)
        end
    end
    sort!(ind)
    unique!(ind)
    return basis[ind]
end

function arrange(p)
    nsupp = copy(p.supp)
    sort!(nsupp)
    unique!(nsupp)
    ncoe = zeros(typeof(p.coe[1]), length(nsupp))
    for i in eachindex(p.supp)
        locb = bfind(nsupp, p.supp[i])
        ncoe[locb] += p.coe[i]
    end
    return poly(nsupp, ncoe)
end

function arrange(p, nb)
    nsupp = reduce_unitnorm.(p.supp, nb)
    sort!(nsupp)
    unique!(nsupp)
    ncoe = zeros(typeof(p.coe[1]), length(nsupp))
    for i in eachindex(p.supp)
        locb = bfind(nsupp, p.supp[i])
        ncoe[locb] += p.coe[i]
    end
    return cpoly(nsupp, ncoe)
end

function divide(mon, lead)
    return any(item -> all(i -> count(==(i), item) <= count(==(i), mon), unique(item)), lead)
end

function reminder(a, x, gb)
    exp = zeros(Int, length(x))
    for i in a
        exp[i] += 1
    end
    rem = Groebner.normalform(gb, prod(x.^exp), ordering=DegRevLex())
    return poly(rem, x)
end

function sign_type(a::Vector{UInt16})
    st = UInt16[]
    if length(a) == 1
        push!(st, a[1])
    elseif length(a) > 1
        r = 1
        for i = 2:length(a)
            if a[i] == a[i-1]
                r += 1
            else
                if isodd(r)
                    push!(st, a[i-1])
                end
                r = 1
            end
        end
        if isodd(r)
            push!(st, a[end])
        end
    end
    return st
end

"""
    p,coe,mon = add_poly!(model, x, degree)

Generate an unknown polynomial of given degree whose coefficients are from the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `x`: set of variables
- `degree`: degree of the polynomial

# Output arguments
- `p`: the polynomial 
- `coe`: coefficients of the polynomial 
- `mon`: monomials of the polynomial 
"""
function add_poly!(model, x, degree::Int; signsymmetry=false)
    mon = vcat([MP.monomials(x, i) for i = 0:degree]...)
    if signsymmetry != false
        ind = [all(transpose(signsymmetry)*exponents(item) .== 0) for item in mon]
        mon = mon[ind]
    end
    coe = @variable(model, [1:length(mon)])
    p = coe'*mon
    return p,coe,mon
end

"""
    p,coe,basis = add_poly_cheby!(model, x, degree)

Generate an unknown polynomial of given degree in the Chebyshev basis whose coefficients are from the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `x`: set of variables
- `degree`: degree of the polynomial

# Output arguments
- `p`: the polynomial 
- `coe`: coefficients of the polynomial 
- `basis`: Chebyshev basis 
"""
function add_poly_cheby!(model, x, degree::Int)
    basis = basis_covering_monomials(ChebyshevBasis, MP.monomials(x, 0:degree))
    coe = @variable(model, [1:length(basis)])
    p = coe'*basis
    return p,coe,basis
end

function add_poly!(model, x, supp::Vector{Vector{UInt16}})
    mon = Mono[prod(x[item]) for item in supp]
    coe = @variable(model, [1:length(mon)])
    p = coe'*mon
    return p,coe,mon
end

function complex_to_real(cpop, z)
    n = length(z)
    @polyvar x[1:2n]
    pop = Vector{Poly{Float64}}(undef, length(cpop))
    for (i, cp) in enumerate(cpop)
        temp = cp(z => x[1:n]+im*x[n+1:2n])
        pop[i] = real.(MP.coefficients(temp))'*MP.monomials(temp)
    end
    return pop,x
end

function cmod(a, m)
    s = mod(a, m)
    return s == 0 ? m : s
end

# generate an SOS polynomial with variables x and degree 2d
"""
    sos = add_SOS!(model, x, d)

Generate an SOS polynomial of degree 2d whose coefficients are from the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `x`: set of variables
- `d`: half degree of the SOS polynomial

# Output arguments
- `sos`: the sos polynomial 
"""
function add_SOS!(model, x, d)
    basis = vcat([MP.monomials(x, i) for i = 0:d]...)
    sos = 0
    pos = @variable(model, [1:length(basis), 1:length(basis)], PSD)
    for j = 1:length(basis), k = j:length(basis)
        if j == k
            @inbounds sos += pos[j,k]*basis[j]*basis[k]
        else
            @inbounds sos += 2*pos[j,k]*basis[j]*basis[k]
        end
    end
    return sos
end

# generate an SOS polynomial with monomial `basis`
"""
    sos = add_SOS!(model, basis)

Generate an SOS polynomial with monomial `basis` whose coefficients are from the JuMP `model`.

# Input arguments
- `model`: a JuMP optimization model
- `basis`: monomial basis

# Output arguments
- `sos`: the sos polynomial 
"""
function add_SOS!(model, basis)
    sos = 0
    pos = @variable(model, [1:length(basis), 1:length(basis)], PSD)
    for j = 1:length(basis), k = j:length(basis)
        if j == k
            @inbounds sos += pos[j,k]*basis[j]*basis[k]
        else
            @inbounds sos += 2*pos[j,k]*basis[j]*basis[k]
        end
    end
    return sos
end
