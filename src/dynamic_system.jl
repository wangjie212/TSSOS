function get_dynamic_sparsity(f, g, x, d; TS=["block","block"], SO=[1,1], merge=false, md=3)
    f = [poly(p, x) for p in f]
    g = [poly([UInt16[]], [1]); [poly(p, x) for p in g]]
    basis = [get_basis(length(x), d-Int(ceil(maxdeg(p)/2))) for p in g]
    dv = 2d + 1 - maximum(maxdeg.(f))
    if TS[1] != false
        supp = vcat([p.supp for p in g]...)
        append!(supp, get_Lsupp(supp, f))
        unique!(supp)
        vsupp = supp[[length(item) <= dv for item in supp]]
        supp = [vsupp; get_Lsupp(vsupp, f)]
        sort!(supp)
        unique!(supp)
    else
        vsupp = get_basis(length(x), dv)
        supp = nothing
    end
    vblocks,vsupp,tsupp,status = get_vblocks(dv, supp, vsupp, f, g, basis, TS=TS[1], SO=SO[1], merge=merge, md=md)
    wsupp = wblocks = nothing
    if status == 1
        if TS[1] != false
            sort!(tsupp)
            unique!(tsupp)
        end
        wblocks,status = get_blocks(tsupp, g, basis, TS=TS[2], SO=SO[2], merge=merge, md=md)
        if status == 1
            wsupp = get_tsupp(g, basis, wblocks)
        end
    end
    return vsupp,vblocks,wsupp,wblocks,status
end

function get_Lsupp(vsupp, f)
    Lsupp = Vector{UInt16}[]
    for (i, p) in enumerate(f), a in vsupp
        loc = bfind(a, i)
        if loc !== nothing
            ca = copy(a)
            deleteat!(ca, loc)
            append!(Lsupp, [sadd(ca, item) for item in p.supp])
        end
    end
    sort!(Lsupp)
    unique!(Lsupp)
    return Lsupp
end

function get_tsupp(g, basis, blocks)
    tsupp = Vector{UInt16}[]
    for (k, p) in enumerate(g), block in blocks[k], i = 1:length(block), j = i:length(block), item in p.supp
        push!(tsupp, sadd(basis[k][block[i]], item, basis[k][block[j]]))
    end
    sort!(tsupp)
    unique!(tsupp)
    return tsupp
end

function get_blocks(tsupp, g, basis; TS="block", SO=1, merge=false, md=3)
    blocks = Vector{Vector{Vector{Int}}}(undef, length(g))
    if TS == false
        blocks = [[Vector(1:length(basis[k]))] for k = 1:length(g)]
        status = 1
    else
        status = 1
        for i = 1:SO
            if i > 1
                oblocks = deepcopy(blocks)
            end
            for (k, p) in enumerate(g)
                G = get_graph(tsupp, p.supp, basis[k])
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
                    tsupp = Vector{UInt16}[]
                    for block in blocks[1], j = 1:length(block), k = j:length(block)
                        push!(tsupp, sadd(basis[1][block[j]], basis[1][block[k]]))
                    end
                    sort!(tsupp)
                    unique!(tsupp)
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

function get_vblocks(dv, supp, vsupp, f, g, basis; TS="block", SO=1, merge=false, md=3)
    tsupp = Vector{UInt16}[]
    if TS == false
        blocks = [[Vector(1:length(basis[k]))] for k = 1:length(g)]
        status = 1
    else
        blocks = Vector{Vector{Vector{Int}}}(undef, length(g))
        status = 1
        qvsupp = vsupp
        for i = 1:SO
            if i > 1
                oblocks = deepcopy(blocks)
            end
            for (k, p) in enumerate(g)
                G = get_graph(supp, p.supp, basis[k])
                if TS == "block"
                    blocks[k] = connected_components(G)
                else
                    blocks[k] = chordal_cliques!(G, method=TS)[1]
                    if merge == true
                        blocks[k] = clique_merge!(blocks[k], QUIET=true, d=md)[1]
                    end
                end
            end
            if i == 1 || oblocks != blocks || length(vsupp) != length(qvsupp)
                vsupp = qvsupp
                if i < SO
                    tsupp = Vector{UInt16}[]
                    for block in blocks[1], j = 1:length(block), k = j:length(block)
                        push!(tsupp, sadd(basis[1][block[j]], basis[1][block[k]]))
                    end
                    sort!(tsupp)
                    unique!(tsupp)
                    qvsupp = tsupp[[length(item) <= dv for item in tsupp]]
                    supp = [tsupp; get_Lsupp(qvsupp, f)]
                    sort!(supp)
                    unique!(supp)
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
