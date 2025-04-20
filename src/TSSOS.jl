module TSSOS

using MosekTools
using JuMP
using Graphs
using DynamicPolynomials
using MultivariatePolynomials
using Ipopt
using LinearAlgebra
using MetaGraphs
using SemialgebraicSets
using COSMO
using Dualization
using Printf
using AbstractAlgebra
using Random
using SymbolicWedderburn
# using AbstractPermutations
# using Hypatia
# using SDPT3
# using SDPNAL

import CliqueTrees

export tssos_first, tssos_higher!, cs_tssos_first, cs_tssos_higher!, local_solution, refine_sol,
cosmo_para, mosek_para, add_psatz!, add_poly!, get_nbasis, get_moment, get_moment_matrix, get_cmoment, homogenize, 
solve_hpop, get_signsymmetry, SumOfRatios, SparseSumOfRatios, LinearPMI_first, LinearPMI_higher!, 
show_blocks, complex_to_real, add_SOSMatrix!, sparseobj, get_mmoment, extract_solutions, extract_solutions_robust, 
extract_solutions_pmo, extract_solutions_robust_pmo, extract_weight_matrix, add_SOS!, tssos_symmetry,
run_H1, run_H1CS, run_H2, run_H2CS, construct_CDK, construct_marginal_CDK, construct_CDK_cs, construct_marginal_CDK_cs

mutable struct cosmo_para
    eps_abs::Float64
    eps_rel::Float64
    max_iter::Int64
    time_limit::Float64
end

cosmo_para() = cosmo_para(1e-5, 1e-5, 1e4, 0)

mutable struct mosek_para
    tol_pfeas::Float64
    tol_dfeas::Float64
    tol_relgap::Float64
    time_limit::Int64
    num_threads::Int64
end

mosek_para() = mosek_para(1e-8, 1e-8, 1e-8, -1, 0)

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

include("clique_merge.jl")
include("blockpop.jl")
include("nblockmix.jl")
include("complex.jl")
include("utils.jl")
include("local_solution.jl")
include("extract_solutions.jl")
include("add_psatz.jl")
include("homogenize.jl")
include("matrixsos.jl")
include("sum_of_ratios.jl")
include("CDK.jl")
include("symmetry.jl")

end
