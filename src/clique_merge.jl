function fweight(a, b; d=3)
    return length(a)^d + length(b)^d - length(union(a, b))^d
end

function mergeclique!(cliques, cliqueG, i, j; d=3)
    cliques[i] = unique(sort([cliques[i]; cliques[j]]))
    cliques[j] = cliques[end]
    cliques = cliques[1:end-1]
    for neighbor in neighbors(cliqueG, j)
        if neighbor != i
            add_edge!(cliqueG, i, neighbor)
        end
    end
    rem_vertex!(cliqueG, j)
    for neighbor in neighbors(cliqueG, i)
        weight = fweight(cliques[i], cliques[neighbor], d=d)
        set_prop!(cliqueG, i, neighbor, :weight, weight)
    end
    return cliques,cliqueG
end

function clique_merge!(cliques; d=3, QUIET=true)
    cql = length(cliques)
    cliqueG = MetaGraph(cql)
    for i = 1:cql, j = i+1:cql
        if intersect(cliques[i], cliques[j]) != []
            add_edge!(cliqueG, i, j)
            weight = fweight(cliques[i], cliques[j], d=d)
            set_prop!(cliqueG, i, j, :weight, weight)
        end
    end
    merge = true
    while merge == true
        mweight = 0
        medge = 0
        for edge in edges(cliqueG)
            weight = get_prop(cliqueG, edge, :weight)
            if weight > mweight
                mweight = weight
                medge = edge
            end
        end
        if mweight > 0
            cliques,cliqueG = mergeclique!(cliques, cliqueG, src(medge), dst(medge), d=d)
        else
            break
        end
    end
    cliquesize = length.(cliques)
    cql = length(cliquesize)
    if QUIET == false
        uc = sort(unique(cliquesize), rev=true)
        sizes = [sum(cliquesize.== i) for i in uc]
        println("-------------------------------------------")
        println("The size of blocks after merging:\n$uc\n$sizes")
        println("-------------------------------------------")
    end
    return cliques,cql,cliquesize
end
