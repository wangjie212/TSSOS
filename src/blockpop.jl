mutable struct cpop_data
    pop # polynomial optimization problem
    x # variables
    n # number of variables
    nb # number of binary variables
    m # number of constraints
    numeq # number of equality constraints
    gb # Grobner basis
    leadsupp # leader terms of the Grobner basis
    supp # support data
    coe # coefficient data
    basis # monomial bases
    ebasis # monomial bases for equality constraints
    ksupp # extended support at the k-th step
    blocksize # sizes of blocks
    blocks # block structure
    eblocks # block structrue for equality constraints
    GramMat # Gram matrices
    multiplier # multipliers for equality constraints
    moment # Moment matrix
    solver # SDP solver
    SDP_status
    rtol # tolerance for rank
    gtol # tolerance for global optimality gap
    ftol # tolerance for feasibility
    flag # 0 if global optimality is certified; 1 otherwise
end

"""
    show_blocks(data)

Display the block structure
"""
function show_blocks(data::cpop_data)
    for j = 1:length(data.blocks[1])
        print("block $j: ")
        println([prod(data.x.^data.basis[1][:, data.blocks[1][j][k]]) for k = 1:data.blocksize[1][j]])
    end
end

"""
    opt,sol,data = tssos_first(pop, x, d; nb=0, numeq=0, GroebnerBasis=true, basis=[], reducebasis=false, TS="block", 
    merge=false, md=3, solver="Mosek", QUIET=false, solve=true, MomentOne=false, Gram=false, solution=false, 
    cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), normality=false, rtol=1e-2, gtol=1e-2, ftol=1e-3)

Compute the first TS step of the TSSOS hierarchy for constrained polynomial optimization.
If `reducebasis=true`, then remove monomials from the monomial basis by diagonal inconsistency.
If `GroebnerBasis=true`, then exploit the quotient ring structure defined by the equality constraints.
If `merge=true`, perform the PSD block merging. 
If `solve=false`, then do not solve the SDP.
If `Gram=true`, then output the Gram matrix.
If `MomentOne=true`, add an extra first-order moment PSD constraint to the moment relaxation.

# Input arguments
- `pop`: vector of the objective, inequality constraints, and equality constraints
- `x`: POP variables
- `d`: relaxation order
- `nb`: number of binary variables in `x`
- `numeq`: number of equality constraints
- `TS`: type of term sparsity (`"block"`, `"signsymmetry"`, `"MD"`, `"MF"`, `false`)
- `md`: tunable parameter for merging blocks
- `normality`: impose the normality condtions (`true`, `false`)
- `QUIET`: run in the quiet mode (`true`, `false`)
- `rtol`: tolerance for rank
- `gtol`: tolerance for global optimality gap
- `ftol`: tolerance for feasibility

# Output arguments
- `opt`: optimum
- `sol`: (near) optimal solution (if `solution=true`)
- `data`: other auxiliary data 
"""
function tssos_first(pop::Vector{Poly{T}}, x, d; nb=0, numeq=0, newton=false, feasibility=false, GroebnerBasis=false, basis=[], reducebasis=false, TS="block", merge=false, md=3, solver="Mosek", 
    QUIET=false, solve=true, dualize=false, MomentOne=false, Gram=false, solution=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false, normality=false, 
    rtol=1e-2, gtol=1e-2, ftol=1e-3) where {T<:Number}
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    n = length(x)
    if nb > 0
        pop = Groebner.normalform(x[1:nb].^2 .- 1, pop)
    end
    if numeq > 0 && GroebnerBasis == true
        cpop = copy(pop)
        gb = convert.(Poly{Float64}, cpop[end-numeq+1:end])
        cpop = cpop[1:end-numeq]
        println("Computing the Gröbner basis...")
        println("This might be slow. You can set GroebnerBasis=false to close it.")
        gb = groebner(gb, ordering=DegRevLex())
        cpop[1] = Groebner.normalform(gb, cpop[1], ordering=DegRevLex())
        lead = leading_ideal(gb, ordering=DegRevLex())
        leadsupp = zeros(UInt8, n, length(lead))
        for i = 1:length(lead), j = 1:n
            @inbounds leadsupp[j,i] = MP.degree(convert(DP.Monomial, lead[i]), x[j])
        end
    else
        cpop = pop
        gb = leadsupp = []
    end
    ss = nothing
    if normality == true || TS == "signsymmetry"
        ss = get_signsymmetry(pop, x)
    end
    m = length(cpop) - 1
    supp,coe = npolys_info(cpop, x)
    isupp = reduce(hcat, supp)
    neq = isempty(gb) ? numeq : 0
    if isempty(basis)
        basis = Vector{Array{UInt8,2}}(undef, m-neq+1)
        if newton == true
            if sum(supp[1][:,1]) != 0 && feasibility == false
               supp[1] = [zeros(UInt8, n) supp[1]]
               coe[1] = [0; coe[1]]
            end
            basis[1] = newton_basis(n, d, supp[1], solver=solver)
        else
            basis[1] = get_basis(n, d, nb=nb, lead=leadsupp)
        end
        for k = 1:m-neq
            basis[k+1] = get_basis(n, d-Int(ceil(MP.maxdegree(pop[k+1])/2)), nb=nb, lead=leadsupp)
        end
        ebasis = nothing
        if isempty(gb) && numeq > 0
            ebasis = Vector{Array{UInt8,2}}(undef, numeq)
            for k = 1:numeq
                ebasis[k] = get_basis(n, 2*d-MP.maxdegree(pop[k+1+m-numeq]), nb=nb, lead=leadsupp)
            end
        end
    end
    tsupp = [isupp bin_add(basis[1], basis[1], nb)]
    tsupp = sortslices(tsupp, dims=2)
    tsupp = unique(tsupp, dims=2)
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,cl,blocksize,eblocks = get_blocks(m-neq, neq, tsupp, supp[2:end], basis, ebasis, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md, signsymmetry=ss)
    if reducebasis == true && GroebnerBasis == false
        gsupp = get_csupp(n, m, numeq, supp, basis[2:end], ebasis, blocks[2:end], eblocks, cl[2:end], blocksize[2:end], nb=nb)
        psupp = [supp[1] zeros(UInt8,n)]
        psupp = [psupp gsupp]
        basis[1],flag = reducebasis!(psupp, basis[1], blocks[1], cl[1], blocksize[1], nb=nb)
        if flag == 1
            tsupp = [isupp bin_add(basis[1], basis[1], nb)]
            tsupp = sortslices(tsupp, dims=2)
            tsupp = unique(tsupp, dims=2)
            blocks,cl,blocksize,eblocks = get_blocks(m-numeq, numeq, tsupp, supp[2:end], basis, ebasis, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md, signsymmetry=ss)
        end
    end
    end
    if QUIET == false
        mb = maximum(maximum.(blocksize))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt,ksupp,moment,momone,GramMat,multiplier,SDP_status = solvesdp(n, m, supp, coe, basis, ebasis, blocks, eblocks, cl, blocksize, nb=nb, numeq=numeq, gb=gb, x=x, dualize=dualize, TS=TS,
    lead=leadsupp, solver=solver, QUIET=QUIET, solve=solve, solution=solution, MomentOne=MomentOne, Gram=Gram, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, writetofile=writetofile, 
    signsymmetry=ss, normality=normality)
    data = cpop_data(pop, x, n, nb, m, numeq, gb, leadsupp, supp, coe, basis, ebasis, ksupp, blocksize, blocks, eblocks, GramMat, multiplier, moment, solver, SDP_status, rtol, gtol, ftol, 1)
    sol = nothing
    if solution == true
        if TS != false || (numeq > 0 && GroebnerBasis == true)
            sol,gap,data.flag = extract_solution(momone, opt, pop, x, numeq=numeq, gtol=gtol, ftol=ftol, QUIET=true)
            if data.flag == 1
                sol = gap > 0.5 ? randn(n) : sol
                sol,data.flag = refine_sol(opt, sol, data, QUIET=true, gtol=gtol)
            end
        else
            sol = extract_solutions_robust(moment[1], n, d, pop=pop, x=x, lb=opt, numeq=numeq, basis=basis[1], check=true, rtol=rtol, gtol=gtol, ftol=ftol, QUIET=QUIET)[1]
            if sol !== nothing
                data.flag = 0
            end
        end
    end
    return opt,sol,data
end

"""
    opt,sol,data = tssos_first(f, x; newton=true, reducebasis=false, TS="block", merge=false, md=3, feasibility=false, solver="Mosek", QUIET=false, solve=true, 
    dualize=false, MomentOne=false, Gram=false, solution=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), rtol=1e-2, gtol=1e-2, ftol=1e-3)

Compute the first TS step of the TSSOS hierarchy for unconstrained polynomial optimization.
If `newton=true`, then compute a monomial basis by the Newton polytope method.
If `reducebasis=true`, then remove monomials from the monomial basis by diagonal inconsistency.
If `TS="block"`, use maximal chordal extensions; if `TS="MD"`, use approximately smallest chordal extensions. 
If `merge=true`, perform the PSD block merging. 
If `feasibility=true`, then solve the feasibility problem.
If `solve=false`, then do not solve the SDP.
If `Gram=true`, then output the Gram matrix.
If `MomentOne=true`, add an extra first-order moment PSD constraint to the moment relaxation.

# Input arguments
- `f`: objective
- `x`: POP variables
- `TS`: type of term sparsity (`"block"`, `"signsymmetry"`, `"MD"`, `"MF"`, `false`)
- `md`: tunable parameter for merging blocks
- `QUIET`: run in the quiet mode (`true`, `false`)
- `rtol`: tolerance for rank
- `gtol`: tolerance for global optimality gap
- `ftol`: tolerance for feasibility

# Output arguments
- `opt`: optimum
- `sol`: (near) optimal solution (if `solution=true`)
- `data`: other auxiliary data 
"""
function tssos_first(f::Poly{T}, x; newton=true, reducebasis=false, TS="block", merge=false, md=3, feasibility=false, solver="Mosek", QUIET=false, solve=true, 
    dualize=false, MomentOne=false, Gram=false, solution=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), rtol=1e-2, gtol=1e-2, ftol=1e-3) where {T<:Number}
    return tssos_first([f], x, Int(ceil(MP.maxdegree(f)/2)), newton=newton, feasibility=feasibility, GroebnerBasis=false, reducebasis=reducebasis, TS=TS, merge=merge, md=md, solver=solver, 
    QUIET=QUIET, solve=solve, dualize=dualize, MomentOne=MomentOne, Gram=Gram, solution=solution, rtol=rtol, gtol=gtol, ftol=ftol, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting)
end

"""
    opt,sol,data = tssos_higher!(data; TS="block", merge=false, md=3, QUIET=false, solve=true, feasibility=false, dualize=false, Gram=false,
    MomentOne=false, solution=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), normality=false)

Compute higher TS steps of the TSSOS hierarchy.
"""
function tssos_higher!(data::cpop_data; TS="block", merge=false, md=3, QUIET=false, solve=true, feasibility=false, dualize=false, MomentOne=false, Gram=false,
    solution=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false, normality=false)
    n = data.n
    nb = data.nb
    m = data.m
    numeq = data.numeq
    x = data.x
    gb = data.gb
    leadsupp = data.leadsupp
    supp = data.supp
    coe = data.coe
    basis = data.basis
    ebasis = data.ebasis
    ksupp = data.ksupp
    solver = data.solver
    ksupp = sortslices(ksupp, dims=2)
    ksupp = unique(ksupp, dims=2)
    if QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    neq = isempty(gb) ? numeq : 0
    blocks,cl,blocksize,eblocks = get_blocks(m-neq, neq, ksupp, supp[2:end], basis, ebasis, nb=nb, TS=TS, QUIET=QUIET, merge=merge, md=md)
    end
    if blocksize == data.blocksize && eblocks == data.eblocks
        println("No higher TS step of the TSSOS hierarchy!")
        opt = sol = nothing
    else
        if QUIET == false
            mb = maximum(maximum.(blocksize))
            println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
        end
        opt,ksupp,moment,momone,GramMat,multiplier,SDP_status = solvesdp(n, m, supp, coe, basis, ebasis, blocks, eblocks, cl, blocksize, nb=nb, numeq=numeq, gb=gb, x=x, lead=leadsupp, TS=TS,
        solver=solver, feasibility=feasibility, QUIET=QUIET, solve=solve, dualize=dualize, solution=solution, MomentOne=MomentOne, Gram=Gram, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, 
        normality=normality, writetofile=writetofile)
        sol = nothing
        if solution == true
            sol,gap,data.flag = extract_solution(momone, opt, data.pop, x, numeq=numeq, QUIET=true, gtol=data.gtol, ftol=data.ftol)
            if data.flag == 1
                sol = gap > 0.5 ? randn(n) : sol
                sol,data.flag = refine_sol(opt, sol, data, QUIET=true, gtol=data.gtol)
            end
        end
        data.blocks = blocks
        data.eblocks = eblocks
        data.blocksize = blocksize
        data.ksupp = ksupp
        data.GramMat = GramMat
        data.multiplier = multiplier
        data.moment = moment
        data.SDP_status = SDP_status
    end
    return opt,sol,data
end

function get_csupp(n, m, numeq, supp, gbasis, ebasis, blocks, eblocks, cl, blocksize; nb=0)
    s = 0
    if m > numeq
        s += sum(size(supp[k+1],2)*Int(sum(Int.(blocksize[k]).^2+blocksize[k])/2) for k=1:m-numeq)
    end
    if numeq > 0
        s += sum(size(supp[k+m-numeq+1],2)*length(eblocks[k]) for k=1:numeq)
    end
    gsupp = zeros(UInt8, n, s)
    l = 1
    for k = 1:m-numeq, i = 1:cl[k], j = 1:blocksize[k][i], r = j:blocksize[k][i], col in eachcol(supp[k+1])
        @inbounds bi = bin_add(gbasis[k][:,blocks[k][i][j]], gbasis[k][:,blocks[k][i][r]], nb)
        @inbounds gsupp[:,l] = bin_add(bi, col, nb)
        l += 1
    end
    for k = 1:numeq, i in eblocks[k], col in eachcol(supp[k+m-numeq+1])
        @inbounds gsupp[:,l] = bin_add(ebasis[k][:,i], col, nb)
        l += 1
    end
    return gsupp
end

function reducebasis!(supp, basis, blocks, cl, blocksize; nb=0)
    esupp = supp[:, all.(iseven, eachcol(supp))]
    init,flag,check = 0,0,0
    while init==0 || check>0
        init,check = 1,0
        tsupp = esupp
        for i = 1:cl
            if blocksize[i] > 1
                for j = 1:blocksize[i], r = j+1:blocksize[i]
                    @inbounds bi = bin_add(basis[:,blocks[i][j]], basis[:,blocks[i][r]], nb)
                    tsupp = [tsupp bi]
                end
            end
        end
        tsupp = unique(tsupp, dims=2)
        tsupp = sortslices(tsupp, dims=2)
        ltsupp = size(tsupp, 2)
        for i = 1:cl
            lo = blocksize[i]
            indexb = [k for k=1:lo]
            j = 1
            while lo >= j
                bi = bin_add(basis[:,blocks[i][indexb[j]]], basis[:,blocks[i][indexb[j]]], nb)
                Locb = bfind(tsupp, ltsupp, bi)
                if Locb === nothing
                   check,flag = 1,1
                   deleteat!(indexb, j)
                   lo -= 1
                else
                   j += 1
                end
            end
            blocks[i] = blocks[i][indexb]
            blocksize[i] = lo
        end
    end
    if flag == 1
       indexb = blocks[1]
       for i = 2:cl
           indexb = append!(indexb, blocks[i])
       end
       sort!(indexb)
       unique!(indexb)
       return basis[:,indexb],flag
    else
       return basis,flag
    end
end

function get_graph(tsupp::Array{UInt8, 2}, basis::Array{UInt8, 2}; nb=0, nvar=0, signsymmetry=nothing)
    lb = size(basis,2)
    G = SimpleGraph(lb)
    ltsupp = size(tsupp,2)
    for i = 1:lb, j = i+1:lb
        bi = bin_add(basis[:,i], basis[:,j], nb)
        if signsymmetry === nothing
            if bfind(tsupp, ltsupp, bi) !== nothing
                add_edge!(G, i, j)
            end
        else
            if all(transpose(signsymmetry)*bi .== 0)
                add_edge!(G, i, j)
            end
        end
    end
    return G
end

function get_graph(tsupp::Array{UInt8, 2}, supp::Array{UInt8, 2}, basis::Array{UInt8, 2}; nb=0, nvar=0, signsymmetry=nothing)
    lb = size(basis, 2)
    G = SimpleGraph(lb)
    ltsupp = size(tsupp, 2)
    for i = 1:lb, j = i+1:lb
        if signsymmetry === nothing
            ind = findfirst(x -> bfind(tsupp, ltsupp, bin_add(bin_add(basis[:,i], basis[:,j], nb), supp[:,x], nb)) !== nothing, 1:size(supp, 2))
            if ind !== nothing
                add_edge!(G, i, j)
            end
        else
            bi = bin_add(bin_add(basis[:,i], basis[:,j], nb), supp[:,1], nb)
            if all(transpose(signsymmetry)*bi .== 0)
                add_edge!(G, i, j)
            end
        end
    end
    return G
end

function get_blocks(tsupp, basis; nb=0, TS="block", QUIET=true, merge=false, md=3, signsymmetry=nothing)
    if TS == false
        blocksize = [size(basis,2)]
        blocks = [[i for i=1:size(basis,2)]]
        cl = 1
    else
        G = get_graph(tsupp, basis, nb=nb, signsymmetry=signsymmetry)
        if TS == "block" || TS == "signsymmetry"
            blocks = connected_components(G)
            blocksize = length.(blocks)
            cl = length(blocksize)
        else
            blocks,cl,blocksize = chordal_cliques!(G, method=TS)
            if merge == true
                blocks,cl,blocksize = clique_merge!(blocks, d=md, QUIET=true)
            end
        end
    end
    if QUIET == false
        sb = sort(unique(blocksize), rev=true)
        numb = [sum(blocksize.== i) for i in sb]
        println("-----------------------------------------------------------------------------")
        println("The sizes of PSD blocks:\n$sb\n$numb")
        println("-----------------------------------------------------------------------------")
    end
    return blocks,cl,blocksize
end

function get_eblock(tsupp::Array{UInt8, 2}, hsupp::Array{UInt8, 2}, basis::Array{UInt8, 2}; nb=0, nvar=0, signsymmetry=nothing)
    ltsupp = size(tsupp, 2)
    hlt = size(hsupp, 2)
    eblock = Int[]
    for i = 1:size(basis, 2)
        if signsymmetry === nothing
            if findfirst(x -> bfind(tsupp, ltsupp, bin_add(basis[:,i], hsupp[:,x], nb)) !== nothing, 1:hlt) !== nothing
                push!(eblock, i)
            end
        else
            bi = bin_add(basis[:,i], hsupp[:,1], nb)
            if all(transpose(signsymmetry)*bi .== 0)
                push!(eblock, i)
            end
        end
    end
    return eblock
end

function get_blocks(m, l, tsupp, supp, basis, ebasis; nb=0, TS="block", QUIET=true, merge=false, md=3, nvar=0, signsymmetry=nothing)
    blocks = Vector{Vector{Vector{Int}}}(undef, m+1)
    eblocks = Vector{Vector{Int}}(undef, l)
    blocksize = Vector{Vector{Int}}(undef, m+1)
    cl = Vector{Int}(undef, m+1)
    if TS == false
        for k = 1:m+1
            lb = ndims(basis[k])==1 ? length(basis[k]) : size(basis[k], 2)
            blocks[k],blocksize[k],cl[k] = [Vector(1:lb)],[lb],1
        end
        for k = 1:l
            lb = ndims(ebasis[k])==1 ? length(ebasis[k]) : size(ebasis[k], 2)
            eblocks[k] = Vector(1:lb)
        end
    else
        for k = 1:m+1
            if k == 1
                G = get_graph(tsupp, basis[1], nb=nb, nvar=nvar, signsymmetry=signsymmetry)
            else
                G = get_graph(tsupp, supp[k-1], basis[k], nb=nb, nvar=nvar, signsymmetry=signsymmetry)
            end
            if TS == "block" || TS == "signsymmetry"
                blocks[k] = connected_components(G)
                blocksize[k] = length.(blocks[k])
                cl[k] = length(blocksize[k])
            else
                blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS)
                if merge == true
                    blocks[k],cl[k],blocksize[k] = clique_merge!(blocks[k], d=md, QUIET=true)
                end
            end
            if QUIET == false
                sb = sort(Int.(unique(blocksize[k])), rev=true)
                numb = [sum(blocksize[k].== i) for i in sb]
                k0 = k - 1
                println("-----------------------------------------------------------------------------")
                println("The sizes of PSD blocks for the $k0-th SOS multiplier:\n$sb\n$numb")
                println("-----------------------------------------------------------------------------")
            end
        end
        for k = 1:l
            eblocks[k] = get_eblock(tsupp, supp[k+m], ebasis[k], nb=nb, nvar=nvar, signsymmetry=signsymmetry)
        end
    end
    return blocks,cl,blocksize,eblocks
end

function solvesdp(n, m, supp, coe, basis, ebasis, blocks, eblocks, cl, blocksize; nb=0, numeq=0, gb=[], x=[], lead=[], solver="Mosek", TS="block",
    QUIET=true, solve=true, feasibility=false, dualize=false, solution=false, MomentOne=false, Gram=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), 
    signsymmetry=false, writetofile=false, normality=false)
    ksupp = zeros(UInt8, n, numele(blocksize[1]))
    k = 1
    for i = 1:cl[1], j = 1:blocksize[1][i], r = j:blocksize[1][i]
        @inbounds bi = bin_add(basis[1][:,blocks[1][i][j]], basis[1][:,blocks[1][i][r]], nb)
        @inbounds ksupp[:,k] = bi
        k += 1
    end
    neq = isempty(gb) ? numeq : 0
    if TS != false && TS != "signsymmetry"
        gsupp = get_csupp(n, m, neq, supp, basis[2:end], ebasis, blocks[2:end], eblocks, cl[2:end], blocksize[2:end], nb=nb)
        ksupp = [ksupp gsupp]
    end
    objv = moment = momone = GramMat = multiplier = SDP_status = nothing
    if solve == true
        tsupp = ksupp
        if normality == true
            wbasis = basis[1]
            bs = size(wbasis, 2)  
            hyblocks = Vector{Vector{Vector{Int}}}(undef, n)
            for i = 1:n
                G = SimpleGraph(2bs)
                for j = 1:bs, k = j:bs
                    bi = bin_add(wbasis[:, j], wbasis[:, k], nb)
                    if all(transpose(signsymmetry)*bi .== 0)
                        add_edge!(G, j, k)
                    end
                    temp = zeros(UInt8, n)
                    temp[i] = 2
                    bi = bin_add(bin_add(wbasis[:, j], wbasis[:, k], nb), temp, nb)
                    if all(transpose(signsymmetry)*bi .== 0)
                        add_edge!(G, j+bs, k+bs)
                    end
                    temp[i] = 1
                    bi = bin_add(bin_add(wbasis[:, j], wbasis[:, k], nb), temp, nb)
                    if all(transpose(signsymmetry)*bi .== 0)
                        add_edge!(G, j, k+bs)
                    end
                end
                hyblocks[i] = connected_components(G)
                for l = 1:length(hyblocks[i])
                    for j = 1:length(hyblocks[i][l]), k = j:length(hyblocks[i][l])
                        if hyblocks[i][l][j] <= bs && hyblocks[i][l][k] > bs
                            temp = zeros(UInt8, n)
                            temp[i] = 1
                            bi = bin_add(bin_add(wbasis[:, hyblocks[i][l][j]], wbasis[:, hyblocks[i][l][k]-bs], nb), temp, nb)
                            tsupp = [tsupp bi]
                        elseif hyblocks[i][l][j] > bs
                            temp = zeros(UInt8, n)
                            temp[i] = 2
                            bi = bin_add(bin_add(wbasis[:, hyblocks[i][l][j]-bs], wbasis[:, hyblocks[i][l][k]-bs], nb), temp, nb)
                            tsupp = [tsupp bi]
                        end
                    end
                end
            end
        end
        if (MomentOne == true || solution == true) && TS != false
            tsupp = [tsupp get_basis(n, 2, nb=nb, lead=lead)]
        end
        if !isempty(gb)
            tsupp = unique(tsupp, dims=2)
            nsupp = zeros(UInt8, n)
            llead = size(lead, 2)
            for col in eachcol(tsupp)
                if divide(col, lead, n, llead)
                    temp = reminder(col, x, gb, n)[2]
                    nsupp = [nsupp temp]
                else
                    nsupp = [nsupp col]
                end
            end
            tsupp = nsupp
        end
        tsupp = sortslices(tsupp, dims=2)
        tsupp = unique(tsupp, dims=2)
        ltsupp = size(tsupp, 2)
        if QUIET == false
            println("Assembling the SDP...")
            println("There are $ltsupp affine constraints.")
        end
        if solver == "Mosek"
            if dualize == false
                model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
            else
                model = Model(dual_optimizer(Mosek.Optimizer))
            end
        elseif solver == "COSMO"
            model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter, "time_limit" => cosmo_setting.time_limit))
        elseif solver == "SDPT3"
            model = Model(optimizer_with_attributes(SDPT3.Optimizer))
        elseif solver == "SDPNAL"
            model = Model(optimizer_with_attributes(SDPNAL.Optimizer))
        else
            @error "The solver is currently not supported!"
            return nothing,nothing,nothing,nothing,nothing,nothing,nothing
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        cons = [AffExpr(0) for i=1:ltsupp]
        if normality == true
            for i = 1:n, l = 1:length(hyblocks[i])
                hbs = length(hyblocks[i][l])
                hnom = @variable(model, [1:hbs, 1:hbs], PSD)
                for j = 1:hbs, k = j:hbs
                    temp = zeros(UInt8, n)
                    if hyblocks[i][l][k] <= bs
                        bi = bin_add(wbasis[:, hyblocks[i][l][j]], wbasis[:, hyblocks[i][l][k]], nb)
                    elseif hyblocks[i][l][j] <= bs && hyblocks[i][l][k] > bs
                        temp[i] = 1
                        bi = bin_add(bin_add(wbasis[:, hyblocks[i][l][j]], wbasis[:, hyblocks[i][l][k]-bs], nb), temp, nb)
                    else
                        temp[i] = 2
                        bi = bin_add(bin_add(wbasis[:, hyblocks[i][l][j]-bs], wbasis[:, hyblocks[i][l][k]-bs], nb), temp, nb)
                    end
                    if !isempty(gb) && divide(bi, lead, n, llead)
                        bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                        for l = 1:bi_lm
                            Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                            if j == k
                                @inbounds add_to_expression!(cons[Locb], bi_coe[l], hnom[j,k])
                            else
                                @inbounds add_to_expression!(cons[Locb], 2*bi_coe[l], hnom[j,k])
                            end
                        end
                    else
                        Locb = bfind(tsupp, ltsupp, bi)
                        if j == k
                            @inbounds add_to_expression!(cons[Locb], hnom[j,k])
                        else
                            @inbounds add_to_expression!(cons[Locb], 2, hnom[j,k])
                        end
                    end                
                end
            end
        end
        pos = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[1])
        for i = 1:cl[1]
            if MomentOne == true || solution == true
                pos0 = @variable(model, [1:n+1, 1:n+1], PSD)
                for j = 1:n+1, k = j:n+1
                    @inbounds bi = bin_add(basis[1][:,j], basis[1][:,k], nb)
                    if !isempty(gb) && divide(bi, lead, n, llead)
                        bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                        for l = 1:bi_lm
                            Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                            if j == k
                               @inbounds add_to_expression!(cons[Locb], bi_coe[l], pos0[j,k])
                            else
                               @inbounds add_to_expression!(cons[Locb], 2*bi_coe[l], pos0[j,k])
                            end
                        end
                    else
                        Locb = bfind(tsupp,ltsupp,bi)
                        if j == k
                           @inbounds add_to_expression!(cons[Locb], pos0[j,k])
                        else
                           @inbounds add_to_expression!(cons[Locb], 2, pos0[j,k])
                        end
                    end
                end
            end
            bs = blocksize[1][i]
            if bs == 1
               @inbounds pos[i] = @variable(model, lower_bound=0)
               @inbounds bi = bin_add(basis[1][:,blocks[1][i][1]], basis[1][:,blocks[1][i][1]], nb)
               if !isempty(gb) && divide(bi, lead, n, llead)
                   bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                   for l = 1:bi_lm
                       Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                       @inbounds add_to_expression!(cons[Locb], bi_coe[l], pos[i])
                   end
               else
                   Locb = bfind(tsupp, ltsupp, bi)
                   @inbounds add_to_expression!(cons[Locb], pos[i])
               end
            else
               @inbounds pos[i] = @variable(model, [1:bs, 1:bs], PSD)
               for j = 1:bs, r = j:bs
                   @inbounds bi = bin_add(basis[1][:,blocks[1][i][j]], basis[1][:,blocks[1][i][r]], nb)
                   if !isempty(gb) && divide(bi, lead, n, llead)
                       bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                       for l = 1:bi_lm
                           Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                           if j == r
                              @inbounds add_to_expression!(cons[Locb], bi_coe[l], pos[i][j,r])
                           else
                              @inbounds add_to_expression!(cons[Locb], 2*bi_coe[l], pos[i][j,r])
                           end
                       end
                   else
                       Locb = bfind(tsupp, ltsupp, bi)
                       if j == r
                          @inbounds add_to_expression!(cons[Locb], pos[i][j,r])
                       else
                          @inbounds add_to_expression!(cons[Locb], 2, pos[i][j,r])
                       end
                   end
               end
            end
        end
        if m > neq
            gpos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m-neq)
        end
        for k = 1:m-neq
            gpos[k] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[k+1])
            for i = 1:cl[k+1]
                bs = blocksize[k+1][i]
                if bs == 1
                    gpos[k][i] = @variable(model, lower_bound=0)
                    for s = 1:size(supp[k+1],2)
                        @inbounds bi = bin_add(basis[k+1][:,blocks[k+1][i][1]], basis[k+1][:,blocks[k+1][i][1]], nb)
                        @inbounds bi = bin_add(bi, supp[k+1][:,s], nb)
                        if !isempty(gb) && divide(bi, lead, n, llead)
                            bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                            for l = 1:bi_lm
                                Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                                @inbounds add_to_expression!(cons[Locb], coe[k+1][s]*bi_coe[l], gpos[k][i])
                            end
                        else
                            Locb = bfind(tsupp, ltsupp, bi)
                            @inbounds add_to_expression!(cons[Locb], coe[k+1][s], gpos[k][i])
                        end
                    end
                else
                    gpos[k][i] = @variable(model, [1:bs, 1:bs], PSD)
                    for j = 1:bs, r = j:bs, s = 1:size(supp[k+1],2)
                        @inbounds bi = bin_add(basis[k+1][:,blocks[k+1][i][j]], basis[k+1][:,blocks[k+1][i][r]], nb)
                        @inbounds bi = bin_add(bi, supp[k+1][:,s], nb)
                        if !isempty(gb) && divide(bi, lead, n, llead)
                            bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                            for l = 1:bi_lm
                                Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                                if j == r
                                   @inbounds add_to_expression!(cons[Locb], coe[k+1][s]*bi_coe[l], gpos[k][i][j,r])
                                else
                                   @inbounds add_to_expression!(cons[Locb], 2*coe[k+1][s]*bi_coe[l], gpos[k][i][j,r])
                                end
                            end
                        else
                            Locb = bfind(tsupp, ltsupp, bi)
                            if j == r
                               @inbounds add_to_expression!(cons[Locb], coe[k+1][s], gpos[k][i][j,r])
                            else
                               @inbounds add_to_expression!(cons[Locb], 2*coe[k+1][s], gpos[k][i][j,r])
                            end
                        end
                    end
                end
            end
        end
        if isempty(gb) && numeq > 0
            free = Vector{Vector{VariableRef}}(undef, numeq)
            for k = 1:numeq
                free[k] = @variable(model, [1:length(eblocks[k])])
                for (i,j) in enumerate(eblocks[k]), s = 1:size(supp[k+m-numeq+1], 2)
                    @inbounds bi = bin_add(ebasis[k][:,j], supp[k+m-numeq+1][:,s], nb)
                    Locb = bfind(tsupp, ltsupp, bi)
                    @inbounds add_to_expression!(cons[Locb], coe[k+m-numeq+1][s], free[k][i])
                end
            end
        end
        for i = 1:size(supp[1], 2)
            Locb = bfind(tsupp, ltsupp, supp[1][:,i])
            if Locb === nothing
                @error "The monomial basis is not enough!"
                return nothing,nothing,nothing,nothing,nothing,nothing,nothing
            else
               cons[Locb] -= coe[1][i]
            end
        end
        if feasibility == false
            @variable(model, lower)
            cons[1] += lower
            @objective(model, Max, lower)
        else
            tnorm = AffExpr(0)
            for i = 1:cl[1]
                if blocksize[1][i] == 1
                    tnorm += pos[1][i]
                else
                    tnorm += tr(pos[1][i])
                end
            end
            @objective(model, Min, tnorm)
        end
        @constraint(model, con, cons==zeros(ltsupp))
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
            println("Solving the SDP...")
        end
        time = @elapsed begin
        optimize!(model)
        end
        if QUIET == false
            println("SDP solving time: $time seconds.")
        end
        if writetofile != false
            write_to_file(dualize(model), writetofile)
        end
        SDP_status = termination_status(model)
        objv = objective_value(model)
        if SDP_status != MOI.OPTIMAL
           println("termination status: $SDP_status")
           status = primal_status(model)
           println("solution status: $status")
        end
        println("optimum = $objv")
        if Gram == true
            GramMat = Vector{Vector{Union{Float64,Matrix{Float64}}}}(undef, m-neq+1)
            GramMat[1] = [value.(pos[i]) for i = 1:cl[1]]
            for k = 1:m-neq
                GramMat[k+1] = [value.(gpos[k][i]) for i = 1:cl[k+1]]
            end
            if neq > 0
                multiplier = [value.(free[j]) for j = 1:neq]
            end
        end
        dual_var = -dual(con)
        moment = Vector{Matrix{Float64}}(undef, cl[1])
        for i = 1:cl[1]
            moment[i] = zeros(blocksize[1][i], blocksize[1][i])
            for j = 1:blocksize[1][i], k = j:blocksize[1][i]
                bi = bin_add(basis[1][:,blocks[1][i][j]], basis[1][:,blocks[1][i][k]], nb)
                if !isempty(gb) && divide(bi, lead, n, llead)
                    bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                    moment[i][j,k] = 0
                    for l = 1:bi_lm
                        Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                        moment[i][j,k] += bi_coe[l]*dual_var[Locb]
                    end
                else
                    Locb = bfind(tsupp, ltsupp, bi)
                    moment[i][j,k] = dual_var[Locb]
                end
            end
            moment[i] = Symmetric(moment[i],:U)
        end
        if solution == true
            momone = zeros(Float64, n+1, n+1)
            for j = 1:n+1, k = j:n+1
                bi = bin_add(basis[1][:,j], basis[1][:,k], nb)
                if !isempty(gb) && divide(bi, lead, n, llead)
                    bi_lm,bi_supp,bi_coe = reminder(bi, x, gb, n)
                    momone[j,k] = 0
                    for l = 1:bi_lm
                        Locb = bfind(tsupp, ltsupp, bi_supp[:,l])
                        momone[j,k] += bi_coe[l]*dual_var[Locb]
                    end
                else
                    Locb = bfind(tsupp, ltsupp, bi)
                    momone[j,k] = dual_var[Locb]
                end
            end
            momone = Symmetric(momone,:U)
        end
    end
    return objv,ksupp,moment,momone,GramMat,multiplier,SDP_status
end
