function chordal_cliques!(G; method="MF", minimize=false)
    # choose algorithm
    alg = method == "MF" && minimize ? CliqueTrees.MinimalChordal(CliqueTrees.MF())  :
          method == "MD" && minimize ? CliqueTrees.MinimalChordal(CliqueTrees.MMD()) :
          method == "MF" && !minimize ? CliqueTrees.MF()                              :
          method == "MD" && !minimize ? CliqueTrees.MMD()                             :
          error()

    # compute maximal cliques
    label, tree = CliqueTrees.cliquetree(G; alg)
    
    # triangulate graph
    F = CliqueTrees.FilledGraph(tree)
    
    for edge in edges(F)
        add_edge!(G, label[src(edge)], label[dst(edge)])
    end
    
    # return maximal cliques
    maximal_cliques = Vector{Vector{UInt16}}(undef, length(tree))
    
    for (i, clique) in enumerate(tree)
        maximal_cliques[i] = sort!(label[clique])
    end
    
    cliquesize = length.(maximal_cliques)
    cql = length(cliquesize)
    return maximal_cliques, cql, cliquesize
end

function add_clique!(G, nodes)
    for i in 1:length(nodes)-1, j in i+1:length(nodes)
        add_edge!(G, nodes[i], nodes[j])
    end
end

function max_cliques(G)
    cliques = convert(Vector{Vector{UInt16}}, maximal_cliques(G))
    sort!.(cliques)
    cliquesize = length.(cliques)
    cql = length(cliquesize)
    return cliques,cql,cliquesize
end
