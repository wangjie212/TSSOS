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
using Groebner
using COSMO
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

export tssos_first, tssos_higher!, cs_tssos_first, cs_tssos_higher!, complex_tssos_first, complex_tssos_higher!, complex_cs_tssos_first, 
complex_cs_tssos_higher!, LinearPMI_first, LinearPMI_higher!, sparseobj
export cosmo_para, mosek_para
export local_solution, refine_sol, extract_solutions, extract_solutions_robust, extract_solutions_pmo, extract_solutions_pmo_robust, extract_weight_matrix
export add_SOS!, add_SOSMatrix!, add_poly!, add_psatz!, add_complex_psatz!, add_psatz_cheby!, add_poly_cheby!
export OnMonomials, tssos_symmetry, complex_tssos_symmetry, get_signsymmetry, tssos_symmetry_first, complex_tssos_symmetry_first, tssos_symmetry_higher!, complex_tssos_symmetry_higher!, add_psatz_symmetry!
export homogenize, solve_hpop, SumOfRatios, SparseSumOfRatios, get_dynamic_sparsity
export show_blocks, complex_to_real, get_mmoment, get_basis, get_moment, get_moment_matrix, get_cmoment
export run_H1, run_H1CS, run_H2, run_H2CS, construct_CDK, construct_marginal_CDK, construct_CDK_cs, construct_marginal_CDK_cs

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

const Mono = DP.Monomial{DP.Commutative{DP.CreationOrder}, Graded{LexOrder}}
const Poly{T} = DP.Polynomial{DP.Commutative{DP.CreationOrder}, Graded{LexOrder}, T}
const PolyLike = Union{Mono,Term,Poly}
export Mono,Poly

include("chordal_extension.jl")
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
include("dynamic_system.jl")
include("Chebyshev_basis.jl")
include("symmetry.jl")

end
