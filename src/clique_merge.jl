function fweight(a,b;d=3)
    return length(a)^d+length(b)^d-length(union(a,b))^d
end

function mergeclique!(cliques,cliqueG,i,j;d=3)
    sort!(union!(cliques[i], cliques[j]))
    deleteat!(cliques,j)
    for neighbor in neighbors(cliqueG, j)
        if neighbor!=i
            add_edge!(cliqueG, i, neighbor)
        end
    end
    rem_vertex!(cliqueG, j)
    for neighbor in neighbors(cliqueG, i)
        weight=fweight(cliques[i],cliques[neighbor],d=d)
        set_prop!(cliqueG, i, neighbor, :weight, weight)
    end
end

function clique_merge!(cliques,cql;QUIET=true,d=3)
    cliqueG=MetaGraph(cql)
    for i=1:cql
        for j=i+1:cql
            if intersect(cliques[i], cliques[j])!=[]
                add_edge!(cliqueG, i, j)
                weight=fweight(cliques[i],cliques[j],d=d)
                set_prop!(cliqueG, i, j, :weight, weight)
            end
        end
    end
    merge=true
    while merge==true
        mweight=0
        medge=0
        for edge in edges(cliqueG)
            weight=get_prop(cliqueG, edge, :weight)
            if weight>mweight
                mweight=weight
                medge=edge
            end
        end
        if mweight>0
            if src(medge)<dst(medge)
                mergeclique!(cliques,cliqueG,src(medge),dst(medge),d=d)
            else
                mergeclique!(cliques,cliqueG,dst(medge),src(medge),d=d)
            end
        else
            break
        end
    end
    cliquesize=length.(cliques)
    cql=UInt16(length(cliquesize))
    if QUIET==false
        uc=unique(cliquesize)
        sizes=[sum(cliquesize.== i) for i in uc]
        println("-------------------------------------------")
        println("The size of blocks after merging:\n$uc\n$sizes")
    end
    return cliques,cql,cliquesize
end
