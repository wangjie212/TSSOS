struct CGraph
    neighbors
    disabled
end

CGraph() = CGraph(Set{UInt16}[], BitSet())

function cadd_node!(G::CGraph)
    push!(G.neighbors, Set{UInt16}())
    return length(G.neighbors)
end

num_nodes(G::CGraph) = length(G.neighbors)

function cadd_edge!(G::CGraph, i, j)
    if !(j in G.neighbors[i])
        push!(G.neighbors[i], j)
        push!(G.neighbors[j], i)
    end
    return (i, j)
end

function cadd_clique!(G::CGraph, nodes)
    for i in 1:(length(nodes) - 1)
        for j in (i + 1):length(nodes)
            cadd_edge!(G, nodes[i], nodes[j])
        end
    end
end

disable_node!(G::CGraph, node) = push!(G.disabled, node)
is_enabled(G::CGraph, node) = !(node in G.disabled)
enable_all_nodes!(G::CGraph) = empty!(G.disabled)

function Cneighbors(G::CGraph, node)
    return G.neighbors[node]
end

function _num_edges_subgraph(G::CGraph, nodes, node)
    neighs = Cneighbors(G, node)
    return count(nodes) do node
        is_enabled(G, node) && node in neighs
    end
end
function num_edges_subgraph(G::CGraph, nodes)
    return mapreduce(+, nodes; init = 0) do node
        is_enabled(G, node) ? _num_edges_subgraph(G, nodes, node) : 0
    end
end

function num_missing_edges_subgraph(G::CGraph, nodes)
    n = count(node -> is_enabled(G, node), nodes)
    return div(n * (n - 1) - num_edges_subgraph(G, nodes), 2)
end

function fill_in(G::CGraph, node)
    return num_missing_edges_subgraph(G, Cneighbors(G, node))
end

function is_clique(G::CGraph, nodes)
    return iszero(num_missing_edges_subgraph(G, nodes))
end

struct FillInCache
    graph::CGraph
    fill_in
end

FillInCache(graph::CGraph) = FillInCache(graph, [fill_in(graph, i) for i in 1:num_nodes(graph)])
num_nodes(G::FillInCache) = num_nodes(G.graph)
Cneighbors(G::FillInCache, node) = Cneighbors(G.graph, node)

function cadd_edge!(G::FillInCache, i, j)
    ni = Cneighbors(G, i)
    nj = Cneighbors(G, j)
    for node in ni
        if node in nj
            G.fill_in[node] -= 1
        end
    end
    G.fill_in[i] += count(k -> is_enabled(G.graph, k), ni) - _num_edges_subgraph(G.graph, ni, j)
    G.fill_in[j] += count(k -> is_enabled(G.graph, k), nj) - _num_edges_subgraph(G.graph, nj, i)
    cadd_edge!(G.graph, i, j)
end

fill_in(G::FillInCache, node) = G.fill_in[node]

function disable_node!(G::FillInCache, node)
    for neighbor in Cneighbors(G, node)
        nodes = Cneighbors(G, neighbor)
        G.fill_in[neighbor] -= (length(nodes) - 1) - _num_edges_subgraph(G.graph, nodes, node)
    end
    disable_node!(G.graph, node)
end

is_enabled(G::FillInCache, node) = is_enabled(G.graph, node)

abstract type AbstractGreedyAlgorithm end

struct GreedyFillIn <: AbstractGreedyAlgorithm end
cache(G::CGraph, ::GreedyFillIn) = FillInCache(G)
heuristic_value(G::FillInCache, node, ::GreedyFillIn) = fill_in(G, node)

function _greedy_triangulation!(G, algo::AbstractGreedyAlgorithm)
    candidate_cliques = Vector{UInt16}[]
    for i in 1:num_nodes(G)
        node = argmin(map(1:num_nodes(G)) do node
            if is_enabled(G, node)
                heuristic_value(G, node, algo)
            else
                typemax(Int)
            end
        end)
	neighbor_nodes = [neighbor for neighbor in Cneighbors(G, node) if is_enabled(G, neighbor)]
	push!(candidate_cliques, [node, neighbor_nodes...])
	for i in eachindex(neighbor_nodes)
		for j in (i+1):length(neighbor_nodes)
			cadd_edge!(G, neighbor_nodes[i], neighbor_nodes[j])
		end
	end
    disable_node!(G, node)
    end
    return candidate_cliques
end

function chordal_extension(G::CGraph, algo::AbstractGreedyAlgorithm)
    candidate_cliques= _greedy_triangulation!(cache(G, algo), algo)
    enable_all_nodes!(G)
    sort!.(candidate_cliques)
    unique!(candidate_cliques)
    candidate_cliques = candidate_cliques[[is_clique(G, clique) for clique in candidate_cliques]]
    sort!(candidate_cliques, by = x -> length(x))
    reverse!(candidate_cliques)
    maximal_cliques = [first(candidate_cliques)]
    for clique in Iterators.drop(candidate_cliques, 1)
        if all(other_clique -> !(clique ⊆ other_clique), maximal_cliques)
            push!(maximal_cliques, clique)
        end
    end
    cliquesize=length.(maximal_cliques)
    cql=UInt16(length(cliquesize))
    return maximal_cliques,cql,cliquesize
end

function max_cliques(A)
    cliques=maximal_cliques(Graph(A))
    sort!.(cliques)
    cliquesize=length.(cliques)
    cql=UInt16(length(cliquesize))
    return cliques,cql,cliquesize
end
