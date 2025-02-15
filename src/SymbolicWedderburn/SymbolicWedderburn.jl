module SymbolicWedderburn

using LinearAlgebra
using SparseArrays
using Primes

using Cyclotomics
using GroupsCore
import AbstractPermutations as AP
import AbstractPermutations: degree
import PermutationGroups as PG
using StarAlgebras
import StarAlgebras as SA

export symmetry_adapted_basis, WedderburnDecomposition, VariablePermutation
export basis, direct_summands, invariant_vectors, issimple, multiplicity

include("Characters/Characters.jl")
using .Characters
import .Characters: row_echelon_form!
import .Characters.FiniteFields

include("ext_homomorphisms.jl")
include("ext_hom_schreier.jl")
include("actions.jl")
include("group_action_error.jl")
include("action_characters.jl")
include("matrix_projections.jl")
include("image_basis.jl")
include("minimal_projections.jl")
include("direct_summands.jl")
include("sa_basis.jl")
include("wedderburn_decomposition.jl")
include("action_polynomials.jl")

end # module
