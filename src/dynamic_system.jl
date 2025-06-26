function get_dynamic_sparsity(f, g, x, d; TS=["block","block"], SO=[1,1], merge=false, md=3, QUIET=false)
    n = length(x)
    m = length(g)
    fsupp = npolys_info(f, x)[1]
    flt = size.(fsupp, 2)
    df = MP.maxdegree.(f)
    gsupp = npolys_info(g, x)[1]
    glt = size.(gsupp, 2)
    dg = MP.maxdegree.(g)
    basis = Vector{Array{UInt8,2}}(undef, m+1)
    basis[1] = get_basis(n, d)
    for i = 1:m
        basis[i+1] = get_basis(n, d-Int(ceil(dg[i]/2)))
    end
    dv = 2d + 1 - maximum(df)
    if TS[1] != false
        tsupp = hcat(gsupp...)
        tsupp = [tsupp get_Lsupp(n, tsupp, fsupp, flt)]
        tsupp = sortslices(tsupp, dims=2)
        tsupp = unique(tsupp, dims=2)
        vsupp = tsupp[:, [sum(item) <= dv for item in eachcol(tsupp)]]
        tsupp1 = [vsupp get_Lsupp(n, vsupp, fsupp, flt)]
        tsupp1 = sortslices(tsupp1, dims=2)
        tsupp1 = unique(tsupp1, dims=2)
    else
        vsupp = get_basis(n, dv)
        tsupp1 = tsupp = nothing
    end
    vblocks,vsupp,tsupp,status = get_vblocks(m, dv, tsupp1, tsupp, vsupp, fsupp, flt, gsupp, glt, basis, TS=TS[1], SO=SO[1], merge=merge, md=md, QUIET=QUIET)
    wsupp = wblocks = nothing
    if status == 1
        if TS[1] != false
            tsupp = sortslices(tsupp, dims=2)
            tsupp = unique(tsupp, dims=2)
        end
        wblocks,status = get_blocks(m, tsupp, gsupp, glt, basis, TS=TS[2], SO=SO[2], merge=merge, md=md, QUIET=QUIET)
        if status == 1
            wsupp = get_tsupp(n, m, gsupp, glt, basis, wblocks)
        end
    end
    return vsupp,vblocks,wsupp,wblocks,status
end

function get_Lsupp(n, vsupp, fsupp, flt)
    Lsupp = zeros(UInt8, n, 1)
    for i = 1:length(flt)
        temp = zeros(UInt8, n)
        temp[i] = 1
        for j = 1:size(vsupp, 2)
            if vsupp[i, j] > 0
                for k = 1:flt[i]
                    Lsupp = [Lsupp vsupp[:,j]-temp+fsupp[i][:,k]]
                end
            end
        end
    end
    Lsupp = sortslices(Lsupp, dims=2)
    Lsupp = unique(Lsupp, dims=2)
    return Lsupp
end

function get_tsupp(n, m, gsupp, glt, basis, blocks)
    blocksize = [length.(blocks[i]) for i = 1:m+1]
    cl = length.(blocksize)
    supp1 = zeros(UInt8, n, Int(sum(Int.(blocksize[1]).^2+blocksize[1])/2))
    k = 1
    for i = 1:cl[1], j = 1:blocksize[1][i], r = j:blocksize[1][i]
        supp1[:,k] = basis[1][:,blocks[1][i][j]] + basis[1][:,blocks[1][i][r]]
        k += 1
    end
    supp2 = zeros(UInt8, n, sum(glt[i]*Int(sum(Int.(blocksize[i+1]).^2+blocksize[i+1])/2) for i=1:m))
    l = 1
    for k = 1:m, i = 1:cl[k+1], j = 1:blocksize[k+1][i], r = j:blocksize[k+1][i], s = 1:glt[k]
        supp2[:,l] = basis[k+1][:,blocks[k+1][i][j]] + basis[k+1][:,blocks[k+1][i][r]] + gsupp[k][:,s]
        l += 1
    end
    tsupp = [supp1 supp2]
    tsupp = sortslices(tsupp,dims=2)
    tsupp = unique(tsupp,dims=2)
    return tsupp
end

function get_blocks(m::Int, tsupp, gsupp::Vector{Array{UInt8, 2}}, glt, basis::Vector{Array{UInt8, 2}}; TS="block", SO=1, merge=false, md=3, QUIET=false)
    blocks = Vector{Vector{Vector{Int}}}(undef, m+1)
    if TS == false
        blocks = [[Vector(1:size(basis[k],2))] for k = 1:m+1]
        status = 1
    else
        status = 1
        for i = 1:SO
            if i > 1
                oblocks = deepcopy(blocks)
            end
            for k = 1:m+1
                if k == 1
                    G = get_graph(tsupp, basis[1])
                else
                    G = get_graph(tsupp, gsupp[k-1], basis[k])
                end
                if TS == "block"
                    blocks[k] = connected_components(G)
                else
                    blocks[k] = chordal_cliques!(G, method=TS)[1]
                    if merge == true
                        blocks[k] = clique_merge!(blocks[k], QUIET=true, d=md)[1]
                    end
                end
            end
            if i == 1 || oblocks != blocks
                if i < SO
                    blocksize = length.(blocks[1])
                    tsupp = zeros(UInt8, size(gsupp[1], 1), numele(blocksize))
                    s = 1
                    for t = 1:length(blocks[1]), j = 1:blocksize[t], r = j:blocksize[t]
                        tsupp[:,s] = basis[1][:,blocks[1][t][j]] + basis[1][:,blocks[1][t][r]]
                        s += 1
                    end
                    tsupp = sortslices(tsupp, dims=2)
                    tsupp = unique(tsupp, dims=2)
                end
            else
                println("No higher TS step of the TSSOS hierarchy!")
                status = 0
                break
            end
        end
    end
    return blocks,status
end

function get_vblocks(m::Int, dv, tsupp1, tsupp, vsupp::Array{UInt8, 2}, fsupp, flt, gsupp::Vector{Array{UInt8, 2}}, glt, basis::Vector{Array{UInt8, 2}}; TS="block", SO=1, merge=false, md=3, QUIET=false)
    blocks = Vector{Vector{Vector{Int}}}(undef, m+1)
    if TS == false
        blocks = [[Vector(1:size(basis[k],2))] for k = 1:m+1]
        status = 1
    else
        status = 1
        qvsupp = vsupp
        for i = 1:SO
            if i > 1
                oblocks = deepcopy(blocks)
            end
            for k = 1:m+1
                if k == 1
                    G = get_graph(tsupp1, basis[1])
                else
                    G = get_graph(tsupp1, gsupp[k-1], basis[k])
                end
                if TS == "block"
                    blocks[k] = connected_components(G)
                else
                    blocks[k] = chordal_cliques!(G, method=TS)[1]
                    if merge == true
                        blocks[k] = clique_merge!(blocks[k], QUIET=true, d=md)[1]
                    end
                end
            end
            if i == 1 || oblocks != blocks || size(vsupp, 2) != size(qvsupp, 2)
                vsupp = qvsupp
                if i < SO
                    blocksize = length.(blocks[1])
                    tsupp = zeros(UInt8, size(gsupp[1], 1), Int(sum(Int.(blocksize).^2+blocksize)/2))
                    s = 1
                    for t = 1:length(blocks[1]), j = 1:blocksize[t], r = j:blocksize[t]
                        tsupp[:,s] = basis[1][:,blocks[1][t][j]] + basis[1][:,blocks[1][t][r]]
                        s += 1
                    end
                    tsupp = sortslices(tsupp, dims=2)
                    tsupp = unique(tsupp, dims=2)
                    qvsupp = tsupp[:, [sum(item) <= dv for item in eachcol(tsupp)]]
                    tsupp1 = [tsupp get_Lsupp(size(gsupp[1], 1), qvsupp, fsupp, flt)]
                    tsupp1 = sortslices(tsupp1, dims=2)
                    tsupp1 = unique(tsupp1, dims=2)
                end
            else
                println("No higher TS step of the TSSOS hierarchy!")
                status = 0
                break
            end
        end
    end
    return blocks,vsupp,tsupp,status
end
