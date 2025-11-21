# include("symmetric_tssos_benchmark_tools.jl");

using PermutationGroups, DynamicPolynomials, TSSOS, JuMP, MosekTools;


polydeg(p) = (ts = terms(p); isempty(ts) ? 0 : maximum(degree, ts))


# 1D ring Ising quartic

function dihedral_group(n::Int)
    cyc = "(" * join(1:n, ",") * ")"
    r = eval(Meta.parse("perm\"$cyc\""))

    pairs = ["($(i),$(n+1-i))" for i in 1:div(n,2)]
    sstr = isempty(pairs) ? "()" : join(pairs, "")
    s = eval(Meta.parse("perm\"$sstr\""))

    return PermGroup([r, s])  # |D_n| = 2n
end

function build_ising_quartic(n::Int; a=1.0, b=0.5, c=-1.0, d=0.2)
n = n
@polyvar X[1:n]
nxt(i) = i == n ? 1 : i + 1
f = sum(a*X[i]^4 + b*X[i]^2*X[nxt(i)]^2 + c*X[i]*X[nxt(i)] + d*X[i]^2 for i in 1:n)
g = n - sum((X[i] - X[i==n ? 1 : i+1])^2 for i in 1:n)
G = dihedral_group(n)
d = ceil(Int, maximum([polydeg(f); map(polydeg, [g])...])/2)
return X, f, g, G, d
end


# 2D torus grid quartic

function cyclic_group(n::Int)
    @assert n >= 1
    cycle_str = "(" * join(1:n, ",") * ")"
    gen = parse(Perm, cycle_str)
    return PermGroup([gen])
end


function Cpq_group(p::Int, q::Int)
    idx(i,j) = (i-1)*q + j
    r = if q == 1
        perm"()"
    else
        cyc_r = join(["(" * join((idx(i,j) for j in 1:q), ",") * ")" for i in 1:p], "")
        eval(Meta.parse("perm\"$cyc_r\""))
    end
    s = if p == 1
        perm"()"
    else
        cyc_s = join(["(" * join((idx(i,j) for i in 1:p), ",") * ")" for j in 1:q], "")
        eval(Meta.parse("perm\"$cyc_s\""))
    end
    return PermGroup([r, s])  # ≅ C_p × C_q
end

function build_grid_quartic(p::Int, q::Int; a=1.0, b=0.5, c=-1.0, d=0.2, opt::Union{Symbol,String}=:smooth)
    n = p*q
    @polyvar X[1:n]
    idx(i,j) = (i-1)*q + j
    nxti(i) = i == p ? 1 : i+1
    nxtj(j) = j == q ? 1 : j+1

    f = zero(X[1])
    for i in 1:p, j in 1:q
        Xi = X[idx(i,j)]
        Xr = X[idx(i,nxtj(j))]
        Xu = X[idx(nxti(i),j)]
        f += a*Xi^4 + b*Xi^2*Xr^2 + b*Xi^2*Xu^2 + c*Xi*Xr + c*Xi*Xu + d*Xi^2
    end

    osym = Symbol(opt)
    if osym === :smooth
        S2 = 2*n
        g = S2 - sum( (X[idx(i,j)] - X[idx(i,nxtj(j))])^2 +
                      (X[idx(i,j)] - X[idx(nxti(i),j)])^2 for i in 1:p, j in 1:q )
    elseif osym === :magnet
        g = n^2 - sum(X)^2
    elseif osym === :spin
        g = sum( (x^2 - 1)^2 for x in X )
    else
        error("Unknown opt=$opt. Use :smooth, :magnet, or :spin.")
    end

    G = Cpq_group(p,q)
    d = ceil(Int, maximum([polydeg(f); map(polydeg, [g])...])/2)

    return X, f, g, G, d
end


# symmetric chordal extension

function symmetric_group(n::Int)
    @assert n >= 2
    gens = [parse(Perm, "($i,$(i+1))") for i in 1:(n-1)]
    return PermGroup(gens)
end

function build_symmetric_quartic(n::Int)

    @polyvar X[1:n]

    t1 = (1/n) * sum(X[i]^4 for i in 1:n)

    t2 = (1/binomial(n,4)) * sum(sum(sum(sum(
        X[i] * X[j] * X[k] * X[l]
        for l in (k+1):n)
        for k in (j+1):(n-1))
        for j in (i+1):(n-2))
        for i in 1:(n-3)
        )

    t3 = (1/binomial(n,3)) * sum(sum(sum(
        X[i] * X[j] * X[k]
        for k in (j+1):n)
        for j in (i+1):(n-1))
        for i in 1:(n-2)
        )

    t4 = (1/n) * sum(X[i] for i in 1:n)

    f = t1 + t2 + t3 + t4

    G=symmetric_group(n)

    r=2

    return X, f, G, r
end