module QuaternionSOS

using DynamicPolynomials
using MultivariatePolynomials
using Graphs
using JuMP
using LinearAlgebra
using Quaternions
import DynamicPolynomials as DP

using ..TSSOS: bfind, add_clique!, max_cliques, chordal_cliques!, clique_merge!, _sdp_model, has_sdp_point

export qtssos, quaternion_to_real

const Poly = DP.Polynomial{DP.Commutative{DP.CreationOrder}, Graded{LexOrder}, Float64}
include("utils.jl")
include("optimization.jl")

end
