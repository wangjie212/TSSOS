module Characters
using LinearAlgebra
using SparseArrays
using Primes
using GroupsCore
using Cyclotomics

import AbstractPermutations as AP
import AbstractPermutations: degree
import PermutationGroups as PG

export AbstractClassFunction, Character, CharacterTable

export conjugacy_classes,
    multiplicities,
    degree,
    irreducible_characters,
    isirreducible,
    table

include("gf.jl")

include("echelon_form.jl")
include("eigenspacedecomposition.jl")
include("cmmatrix.jl")

include("powermap.jl")
include("character_tables.jl")
include("class_functions.jl")
include("chars.jl")
include("io.jl")

include("dixon.jl")
end
