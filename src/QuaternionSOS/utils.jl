# A term consists of paired norm factors and an ordered noncommutative word.
function standardterm(b::Vector{Vector{UInt16}}, n)
    a = deepcopy(b)
    i = 1
    while i < length(a[3])
        x, y = a[3][i:i+1]
        if abs(Int(x)-Int(y)) == n
            push!(a[1], min(x, y))
            push!(a[2], max(x, y))
            deleteat!(a[3], i:i+1)
            i = max(1, i-1)
        else
            i += 1
        end
    end
    sort!(a[1])
    sort!(a[2])
    return a
end

function star(b::Vector{Vector{UInt16}}, n)
    word = UInt16[v <= n ? v+n : v-n for v in reverse(b[3])]
    return standardterm([copy(b[1]), copy(b[2]), word], n)
end

# Products never mutate their input terms.
function qproduct(terms, n)
    return standardterm([vcat((term[i] for term in terms)...) for i in 1:3], n)
end
qtermadd4(a, b, c, d, n) = qproduct((a, b, c, d), n)
qtermadd(a, b, n) = qproduct((a, star(b, n)), n)
qtermaddleft(a, b, n) = qproduct((star(a, n), b), n)
qtermadd3(a, b, c, n) = qproduct((a, b, star(c, n)), n)
qtermadd3left(a, b, c, n) = qproduct((star(a, n), b, c), n)

function mono_to_term(m, q, n)
    varmap = Dict(v => UInt16(i) for (i, v) in enumerate(q))
    word = UInt16[]
    for (v, exponent) in zip(m.vars, m.z)
        iszero(exponent) && continue
        haskey(varmap, v) || throw(ArgumentError("Polynomial contains a variable outside q."))
        append!(word, fill(varmap[v], exponent))
    end
    return standardterm([UInt16[], UInt16[], word], n)
end
mono_to_term(m::Number, q, n) = [UInt16[], UInt16[], UInt16[]]

get_qncbasis(n::Int, d; kwargs...) = get_qncbasis(collect(1:n), n, d; kwargs...)

function get_qncbasis(vars::AbstractVector{<:Integer}, n, d; conjubasis=false)
    d >= 0 || throw(ArgumentError("Basis degree must be nonnegative."))
    alphabet = conjubasis ? [vars; vars .+ n] : vars
    basis = [[UInt16[], UInt16[], UInt16[]]]
    degree = copy(basis)
    for _ in 1:d
        degree = unique([standardterm([copy(term[1]), copy(term[2]), [term[3]; UInt16(v)]], n)
                         for term in degree for v in alphabet])
        append!(basis, degree)
    end
    return basis
end

function qreduce_unitnorm(b, n; nb=0)
    a = standardterm(b, n)
    # Only the first nb variables are constrained to unit norm.
    keep = findall(>(nb), a[1])
    return [a[1][keep], a[2][keep], a[3]]
end

function qresort(supp, coe, n; nb=0)
    reduced = qreduce_unitnorm.(supp, n; nb=nb)
    nsupp = sort!(unique(reduced))
    ncoe = zeros(eltype(coe), length(nsupp))
    for (term, coefficient) in zip(reduced, coe)
        ncoe[bfind(nsupp, term)] += coefficient
    end
    return nsupp, ncoe
end

function reduced_terms(term, n, nb, ipart)
    term = qreduce_unitnorm(term, n; nb)
    cyclic = qreduce_unitnorm(canonical_qmono(term, n), n; nb)
    return ipart ? term : cyclic, cyclic
end

# Re(a* H (c*w) b) = Re((H*c) w b a*). Its adjoint has the
# reversed/conjugated word and the conjugate of the ENTIRE coefficient H*c.
function quaternion_right_product(h, c)
    cr = real(c)
    ci, cj, ck = imag_part(c)
    hr, hi, hj, hk = h
    return (hr*cr-hi*ci-hj*cj-hk*ck,
            hr*ci+hi*cr+hj*ck-hk*cj,
            hr*cj-hi*ck+hj*cr+hk*ci,
            hr*ck+hi*cj-hj*ci+hk*cr)
end

function add_quaternion_equality_pair!(rcons, icons, jcons, kcons, tsupp, ttsupp,
                                      word, coefficient, n, nb)
    for (term, components) in ((word, coefficient),
                               (star(word, n), (coefficient[1], -coefficient[2],
                                                -coefficient[3], -coefficient[4])))
        term, cyclic = reduced_terms(term, n, nb, true)
        loc = bfind(tsupp, term)
        real_loc = bfind(ttsupp, cyclic)
        add_to_expression!(rcons[real_loc], components[1])
        add_to_expression!(icons[loc], components[2])
        add_to_expression!(jcons[loc], components[3])
        add_to_expression!(kcons[loc], components[4])
    end
end

function conjugate_representatives(supp, n; cyclic=false)
    keep = trues(length(supp))
    visited = Set{Int}()
    for i in eachindex(supp)
        i in visited && continue
        adjoint = star(supp[i], n)
        cyclic && (adjoint = canonical_qmono(adjoint, n))
        j = bfind(supp, adjoint)
        push!(visited, i)
        if j !== nothing && j != i
            push!(visited, j)
            keep[max(i, j)] = false
        end
    end
    return findall(keep)
end

"""
    canonical_qmono(m, n)

Choose a representative under cyclic invariance of the real part of a quaternion word.
"""
function canonical_qmono(m::Vector{Vector{UInt16}}, n::Int)

    r  = sort(deepcopy(m[1]))
    im = sort(deepcopy(m[2]))
    v  = deepcopy(m[3])

    nn = UInt16(n)

    while true
        old_r  = deepcopy(r)
        old_im = deepcopy(im)
        old_v  = deepcopy(v)

        temp = standardterm([r, im, v], n)

        r  = temp[1]
        im = temp[2]
        v  = temp[3]

        if length(v) >= 2

            if v[1] == v[end] + nn
                push!(r, v[end])
                push!(im, v[1])

                v = v[2:end-1]

            elseif v[1] == v[end] - nn
                push!(im, v[end])
                push!(r, v[1])

                v = v[2:end-1]
            end
        end

        sort!(r)
        sort!(im)

        if r == old_r && im == old_im && v == old_v
            break
        end
    end

    if length(v) >= 2

        rotations = Vector{Vector{UInt16}}()

        d = length(v)

        for k in 0:d-1
            rot = vcat(v[k+1:end], v[1:k])
            push!(rotations, rot)
        end

        v = minimum(rotations)
    end

    return standardterm([r, im, v], n)
end

function generate_rcomm_constraints(q, n)
    comm_constraints = typeof(1.0*q[1]+2.0*q[2])[]

    for i in 1:n
        for j in 1:n
            if i != j
                expr = q[i]*q[j] + q[i+n]*q[j] -
                    q[j]*q[i] - q[j]*q[i+n]
                push!(comm_constraints, expr)
            end
        end
    end

    return comm_constraints
end

function generate_icomm_constraints(n)

    Qi = QuaternionF64(0,1,0,0)
    Qj = QuaternionF64(0,0,1,0)
    Qk = QuaternionF64(0,0,0,1)

    pop_num = n*(n-1)

    supp = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef,pop_num)
    coe  = Vector{Vector{QuaternionF64}}(undef,3)

    cnt = 1

    for i in 1:n
        for j in 1:n

            i == j && continue

            supp[cnt] = [

                [[UInt16[], UInt16[], [UInt16(i)]],[UInt16[], UInt16[], [UInt16(j)]]],

                [[UInt16[], UInt16[], UInt16[]],[UInt16[], UInt16[], [UInt16(i+n), UInt16(j)]]],

                [[UInt16[], UInt16[], [UInt16(j), UInt16(i)]],[UInt16[], UInt16[], UInt16[]]],

                [[UInt16[], UInt16[], [UInt16(j)]],[UInt16[], UInt16[], [UInt16(i+n)]]]
            ]
            cnt += 1

        end
    end
    coe[1] = QuaternionF64[-Qi,Qi,Qi,-Qi]
    coe[2] = QuaternionF64[-Qj,Qj,Qj,-Qj]
    coe[3] = QuaternionF64[-Qk,Qk,Qk,-Qk]
    return supp, coe

end
