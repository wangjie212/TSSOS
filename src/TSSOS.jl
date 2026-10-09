module TSSOS

using Base.Threads
using JuMP
using Graphs
using DynamicPolynomials
using MultivariatePolynomials
using Ipopt
using LinearAlgebra
using MetaGraphs
using SemialgebraicSets
using Groebner
using Dualization
using Printf
using AbstractAlgebra
using Random
using SymbolicWedderburn
using AbstractPermutations
import Reexport
Reexport.@reexport using MultivariateBases
import DynamicPolynomials as DP
import MultivariatePolynomials as MP
import CliqueTrees

export tssos, cs_tssos, complex_tssos, complex_cs_tssos, LinearPMI, sparseobj
export arrange, bfind, MosekParameters
export local_solution, refine_sol, extract_solutions, extract_solutions_robust, extract_solutions_pmo, extract_solutions_pmo_robust, extract_weight_matrix
export add_SOS!, add_SOSMatrix!, add_poly!, add_psatz!, add_complex_psatz!, add_psatz_cheby!, add_poly_cheby!
export OnMonomials, tssos_symmetry, complex_tssos_symmetry, get_signsymmetry, add_psatz_symmetry!
export homogenize, solve_hpop, SumOfRatios, SparseSumOfRatios, get_dynamic_sparsity
export show_blocks, complex_to_real, get_mmoment, get_basis, get_moment, get_moment_matrix, get_cmoment
export run_H1, run_H1CS, run_H2, run_H2CS, construct_CDK, construct_marginal_CDK, construct_CDK_cs, construct_marginal_CDK_cs

mutable struct MosekParameters
    tol_pfeas::Float64
    tol_dfeas::Float64
    tol_relgap::Float64
    time_limit::Int64
    num_threads::Int64
end

MosekParameters() = MosekParameters(1e-8, 1e-8, 1e-8, -1, 0)

#temporarily receive mosek_setting parameter so as to not break the interface.
#in future work, remove mosek_setting entirely in favor of MosekExt's mosek_optimizer_from_settings
function default_optimizer(mosek_setting=nothing)
    mosek_extension = Base.get_extension(TSSOS, :MosekExt)
    cosmo_extension = Base.get_extension(TSSOS, :COSMOExt)
    if !isnothing(mosek_extension)
        @info "Using Mosek as default solver" maxlog=1 
        return mosek_extension.mosek_optimizer(mosek_setting)
    elseif !isnothing(cosmo_extension)
        @info "Using COSMO as default solver since MosekTools was not loaded" maxlog=1
        return cosmo_extension.cosmo_optimizer(mosek_setting)
    else
        error("No supported SDP solver found. Currently supported default solvers are Mosek and COSMO: import at least one before TSSOS. This is necessary because there are currently multiple internal SDP solver calls not available for configuring (change soon)")
    end
end

include("polynomial.jl")
include("utils.jl")
include("chordal_extension.jl")
include("clique_merge.jl")
include("term_sparsity.jl")
include("all_sparsity.jl")
include("local_solution.jl")
include("extract_solutions.jl")
include("add_psatz.jl")
include("homogenization.jl")
include("sum_of_ratios.jl")
include("matrixsos.jl")
include("dynamic_system.jl")
include("CDK.jl")
include("Chebyshev_basis.jl")
include("complex_pop.jl")
include("symmetry.jl")

end
